# Shared data-loading functions, included by all plot scripts.
# No top-level execution — only function definitions.

using DataManipulation
using AccessorsExtra
using StructArrays
using FlexiJoins
using VirtualObservatory
using UncertainSkyCoords
using SkyCoords: separation, position_angle
import MonteCarloMeasurements as MCM
using Distributions: Poisson, ccdf
using DirectionalStatistics
using Unitful, UnitfulAngles, UnitfulAstro
using StaticArrays
using LinearAlgebra
using Statistics
using DateFormats
using IntervalSets
using Random: seed!
using VOTables, DictArrays
import CodecZstd
import FixedWidthTables as FWT
import QuackIO
using AstroAngles: hms2rad, dms2rad
using Unitful: NoUnits
using Dates: DateTime

# ── Flare model ───────────────────────────────────────────────────────

const GAIA_DR3_EPOCH = 2016.0

#  PM-only: Z is 2N x 4, y = [VG_alpha, VG_delta, mu_alpha*, mu_delta]
function build_Z_pm(Nt, dt)
    vcat(
        hcat(ones(Nt), zeros(Nt), dt,        zeros(Nt)),
        hcat(zeros(Nt), ones(Nt), zeros(Nt),  dt),
    )
end

function build_Sigma_inv_pm(Sigma_VG, Sigma_mu)
    S = zeros(4, 4)
    S[1:2, 1:2] .= inv(Symmetric(Sigma_VG))
    S[3:4, 3:4] .= inv(Symmetric(Sigma_mu))
    S
end

#  Physical model:  position(t) = x_flare + A / S(t)
#
#  In 2D stacked form:
#    D = [I_2, q*I_2]  (2N x 4)  where q = 1/S per epoch
#    theta = [x_flare_alpha, x_flare_delta, A_alpha, A_delta]
#
#  Gaia basis Z compresses the epoch data into K astrometric coefficients y:
#    G = (Z'Z)^{-1}Z'    (K x 2N, deterministic)
#    M = G * D            (K x 4)
#
#  GLS solution of y = M theta:
#    (M'Sigma^{-1}M) theta = M'Sigma^{-1}y    (always 4x4)
#
function estimate_flare_2d(y, q, Z, Sigma_inv)
    Nt = length(q)
    z = zero(eltype(q))
    zv = fill(z, Nt)

    D = vcat(
        hcat(ones(typeof(z), Nt), zv, q,  zv),
        hcat(zv, ones(typeof(z), Nt), zv,  q),
    )

    G = (Z' * Z) \ Z'
    M = G * D

    MtS = M' * Sigma_inv
    MtSM = MtS * M
    MtSy = MtS * y

    theta = MtSM \ MtSy
    return (; x_flare=SVector(theta[1], theta[2]), A=SVector(theta[3], theta[4]))
end

function flare_major_std(x_flare)
    px = x_flare[1].particles
    py = x_flare[2].particles
    C = cov(hcat(px, py))
    sqrt(maximum(eigvals(Symmetric(C))))
end

# ── load_vgmatches: RFC×Gaia cross-match ────

function load_vgmatches(root)
    raw_dir = joinpath(root, "data/raw")

    # RFC catalog
    rfc_tbl = let
        raw = FWT.read(joinpath(raw_dir, "rfc_2025d_cat.txt"), (
            J2000_name = (5:14, String),
            IVS_name   = (17:24, String),
            ra         = (27:41, String),
            dec        = (43:57, String),
            ra_err     = (60:65, Float64),
            dec_err    = (67:72, Float64),
            err_corr   = (76:81, Float64),
        ); skiprows_startwith=["#"])
        map(raw) do r
            dec = dms2rad(r.dec)
            coords = U.Value(ICRSCoords(hms2rad(r.ra), dec),
                             U.CovMat(σx=r.ra_err * cos(dec) * u"mas", σy=r.dec_err * u"mas", ρ=r.err_corr))
            (; names=(IVS=r.IVS_name, J2000=r.J2000_name), coords)
        end |> StructArray
    end

    # Source density counts
    counts_tbl = QuackIO.read_csv(StructArray, joinpath(raw_dir, "gaia_rfc_counts.csv.zst"))

    # Gaia matches VOTable
    raw_gaia = open(joinpath(raw_dir, "gaia_rfc_matches.vot.zst")) do io
        VOTables.read(DictArray, CodecZstd.ZstdDecompressorStream(io); unitful=true, quiet=true)
    end

    # Build structured Gaia DR3 rows (only fields used by paper scripts)
    gaia_rows = map(raw_gaia) do r
        coord = U.Value(ICRSCoords(r.ra, r.dec), U.CovMat(; σx=r.ra_error, σy=r.dec_error, ρ=r.ra_dec_corr))
        # PM columns are either all present or all missing
        pm = ismissing(r.pmra) ?
            U.Value(SVector(NaN, NaN) .* u"mas/yr", U.CovMat(; σx=NaN32 * u"mas/yr", σy=NaN32 * u"mas/yr", ρ=NaN32)) :
            U.Value(SVector(r.pmra, r.pmdec), U.CovMat(; σx=r.pmra_error, σy=r.pmdec_error, ρ=Float32(r.pmra_pmdec_corr)))
        (; r.rfc_name, G = (; r.source_id, coord, pm))
    end

    # Join Gaia matches with RFC catalog, grouped by RFC source
    MR = innerjoin(
        (; M=gaia_rows, R=rfc_tbl),
        by_key(:rfc_name, @o(_.names.IVS));
        groupby=:R, cardinality=(M=*, R=1),
    )

    # Join with source density counts
    MRC = innerjoin(
        (; __=MR, C=counts_tbl),
        by_key(@o(_.R.names.IVS), :rfc_name);
        cardinality=(__=0:1, C=1),
    )

    # Per-group: keep closest Gaia match, compute Poisson p_chance, apply Bonferroni
    p_max = 1.0 / length(rfc_tbl)  # Bonferroni correction (n_chance=1)

    _COUNTS_RADIUS = 300u"arcsecond"
    @p MRC |> filtermap() do (; R, M, C)
        @assert C.n_gaia > 0
        closest = argmin(m -> norm(U.value(separation(SphericalOffsetFlat, R.coords, m.G.coord))), M)
        offset = U.value(separation(SphericalOffsetFlat, R.coords, closest.G.coord))
        λ = NoUnits(C.n_gaia / (π * _COUNTS_RADIUS^2) * π * norm(offset)^2)
        p_chance = ccdf(Poisson(λ), 0)  # P(X ≥ 1)
        p_chance ≥ p_max && return nothing
        (; R, closest.G, p_chance,
           VG=(vec=separation(SphericalOffsetFlat, R.coords, closest.G.coord),
               pa=position_angle(R.coords, closest.G.coord)))
    end |> StructArray
end

# ── load_glightcurves: Gaia epoch photometry ──────────────────────────

function load_glightcurves(root)
    raw = open(joinpath(root, "data/raw/gaia_epoch_photometry.vot.zst")) do io
        VOTables.read(StructArray, CodecZstd.ZstdDecompressorStream(io); unitful=true, quiet=true)
    end

    @p raw |> group_vg(_.source_id) |> map() do rows
        measurements = @p let
            StructArray((
                time   = rows.g_transit_time,
                flux   = rows.g_transit_flux .±ᵤ rows.g_transit_flux_error,
                reject = rows.variability_flag_g_reject,
            ))
            filter(!_.reject && !isnan(_.flux))
            map((; _.time, _.flux))
            collect
        end
        (; source_id=key(rows), measurements)
    end |> StructArray
end

# ── load_jetdirs: VLBI jet directions from VizieR ────────────────────

load_jetdirs() = @p let
    execute(TAPService(:vizier), """select * from "J/ApJS/260/4/table2" """; cache=true)
    map((J2000=_.Name, pa=_.PA))
end

# ── load_allsrc: VGmatches + jetdirs ─────────────────────────────────

function load_allsrc(root)
    VGmatches = load_vgmatches(root)
    jetdirs = load_jetdirs()
    @p let
        innerjoin(
            (_=VGmatches, J=jetdirs),
            by_key(x -> x.R.names.J2000, :J2000),
            cardinality=(0:1, 0:1,)
        )
        mapinsert(VG_J_pa=Circular.center_angle(_.VG.pa - _.J.pa))
    end
end

# ── load_and_fit: allsrc + lightcurves + MC flare ───────────────────

const NPARTICLES = 500

function load_and_fit(root)
    allsrc = load_allsrc(root)
    glightcurves = load_glightcurves(root)

    merged = @p let
        innerjoin((__=allsrc, PHOTg=glightcurves), by_key(x->x.G.source_id, x->x.source_id))
        @set __.PHOTg = __.PHOTg.measurements
    end

    seed!(1)

    @p let
        merged

        mapinsert(pm_J_pa=function(r)
            pa = atan(U.value(r.G.pm)...) - r.J.pa
            mod(pa, -π/2..(2π-π/2))
        end)
        map(@insert _.J.nvec = SVector(sincos(_.J.pa)))
        filter(U.nσ(_.VG.vec) > 3 && U.nσ(_.G.pm) > 0)

        # MC flare estimation: PM-only in 2D
        mapinsert(mc=function(r)
            nvec = r.J.nvec

            # MC particles for Gaia astrometric parameters
            VG_2d = U.uconvert(MCM.Particles{Any,NPARTICLES}, ustrip(u"mas", r.VG.vec))
            μ = U.uconvert(MCM.Particles{Any,NPARTICLES}, ustrip(u"mas/yr", r.G.pm))

            # MC particles for reciprocal flux q = 1/S
            q = [1 / U.uconvert(MCM.Particles{Any,NPARTICLES}, ustrip(p.flux)) for p in r.PHOTg]

            # Deterministic basis components
            t = [yeardecimal(p.time) for p in r.PHOTg]
            Δt = t .- GAIA_DR3_EPOCH
            Nt = length(t)

            # q(t) linear trend, referenced to GAIA_DR3_EPOCH
            avgt = mean(t)
            Δt_c = t .- avgt
            mean_q = mean(q)
            slope_q = dot(Δt_c, q .- mean_q) / dot(Δt_c, Δt_c)
            q_ref = mean_q + slope_q * (GAIA_DR3_EPOCH - avgt)

            # Covariance matrices
            Σ_VG = ustrip.(u"mas^2", U.uncertainty(r.VG.vec).cov)
            Σ_μ  = ustrip.(u"mas^2/yr^2", U.uncertainty(r.G.pm).cov)

            # PM-only
            y_pm  = [VG_2d[1], VG_2d[2], μ[1], μ[2]]
            Z_pm  = build_Z_pm(Nt, Δt)
            Σi_pm = build_Sigma_inv_pm(Σ_VG, Σ_μ)
            res_pm = estimate_flare_2d(y_pm, q, Z_pm, Σi_pm)

            # 1D projection along jet (for summary plot)
            x_flare_pm_1d = dot(res_pm.x_flare, nvec)

            return (;
                pm = res_pm,
                x_flare_pm_1d,
                q_slope = U.Value(slope_q),
                q_ref = U.Value(q_ref),
            )
        end)
    end
end
