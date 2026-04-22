# Packages for the plot itself
using CairoMakie
using MakieExtra
using IntervalSets: ±
using Statistics: median
using Unitful: NoUnits

const ROOT = joinpath(@__DIR__, "..")

include(joinpath(@__DIR__, "prepare_data.jl"))
include(joinpath(@__DIR__, "theme.jl"))

# --- Load data ---

enriched = load_and_fit(ROOT)

# --- Filter sources ---

sources = @p let
    enriched
    filter(U.nσ(_.VG.vec) > 3)
    filter(Circular.distance(U.value(_.VG_J_pa), 0u"°") < 45u"°")
    filter(!isnan(_.pm_J_pa))
    filter(U.nσ(_.G.pm) > 3)
    filter(HIGHLIGHTED_SOURCES_PRED(_.R.names.J2000))
    filter(!isempty(_.PHOTg))
    sort(by=__ -> __.PHOTg |> extrema(_.flux) |> __[2] / __[1], rev=true)
    collect
end

println("Generating lightcurve+position plots for $(length(sources)) sources...")

# --- Generate per-source figures ---

outdir = joinpath(ROOT, "figs/gaia_lc")
mkpath(outdir)

for r in sources
    j2000 = r.R.names.J2000
    fig = Figure()

    lc = @p r.PHOTg sort(by=_.time)
    nflux = @p lc map(_.flux) median

    # --- Row 1: Gaia lightcurve ---
    ax1 = Axis(fig[1, 1], width=300, height=200,
        xticks=2010:2030, xminorticksvisible=true, xminorticks=2010:(1 // 12):2030,
        yscale=SymLog(1), yticks=BaseMulTicks([1,2,5]), yminorticksvisible=true, yminorticks=BaseMulTicks(1:9),
        xlabel="Date", ylabel=rich("Gaia", font=:italic) * " G-band flux (cnt/s)")
    fplt = FPlot(lc, @o(yeardecimal(_.time)), @o(ustrip(u"s^-1", _.flux));
        color=:black, markersize=8)
    multiplot!(ax1, (scatterlines, rangebars), fplt, label="Measurements")
    lines!(ax1, range(2010..2020, 50),
        t -> ustrip(u"s^-1", nflux) * r.mc.q_ref / (r.mc.q_ref + r.mc.q_slope * (t - GAIA_DR3_EPOCH));
        to_xy_attrs(autolimits=false)...,
        label="Linear trend in 1/S(t)")
    axislegend(ax1, j2000, merge=true, position=(0.5, 1), backgroundcolor=(:white, 0.5), framecolor=:transparent)

    # --- Row 2: predicted Gaia position along the jet ---
    ax2 = Axis(fig[2, 1], height=200,
        xticks=2010:2030, xminorticksvisible=true, xminorticks=2010:(1 // 12):2030,
        xlabel="Date", ylabel=rich("Gaia", font=:italic) * " offset along the jet (mas)")
    let
        v = dot(r.G.pm, r.J.nvec)  # uncertain
        xavg = dot(r.VG.vec, r.J.nvec)  # uncertain

        # Flare model: position along jet from MC
        x_flare_1d = dot(MCM.pmean.(r.mc.pm.x_flare), r.J.nvec)
        A_1d = dot(MCM.pmean.(r.mc.pm.A), r.J.nvec)

        # Proper motion line with 1σ error band
        x_of_t(t) = xavg + v * (t - GAIA_DR3_EPOCH) * u"yr"
        x_mas(t) = let x = x_of_t(t); ustrip(u"mas", U.value(x)) ± ustrip(u"mas", U.uncertainty(x)) end
        band!(ax2, x_mas; alpha=0.2, label="Proper motion (measured)")
        lines!(ax2, t -> ustrip(u"mas", U.value(x_of_t(t))), label="Proper motion (measured)")

        fplt = FPlot(lc, @o(yeardecimal(_.time)), @o(x_flare_1d + A_1d / ustrip(_.flux));
            color=:black)
        multiplot!(ax2, (scatter, lines), fplt, label="Per-epoch (prediction)")
    end
    axislegend(ax2, merge=true, position=(r.mc.q_slope > 0 ? 1 : 0, 0), backgroundcolor=(:white, 0.5), framecolor=:transparent)

    resize_to_layout!()
    outpath = joinpath(outdir, "$j2000.pdf")
    save(outpath, fig)
    println("  $j2000 -> $outpath")
end

println("Done.")
