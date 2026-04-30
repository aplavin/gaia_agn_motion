using CairoMakie
using MakieExtra
using MakieExtra: @lift
using MakieBake

const ROOT = joinpath(@__DIR__, "..")

include(joinpath(@__DIR__, "cartoon_common.jl"))

# --- Gaussian flare light curve ---

const F0 = 0.2   # baseline flux (constant component)
const A = 1.0     # flare amplitude
const T_PEAK = 2.5
const SIGMA_T = 1.2
const X_CONST = 1.8  # position of constant jet emission centroid

function flare_flux(t; t_peak=T_PEAK, σ=SIGMA_T, A=A, F0=F0)
    return F0 + A * exp(-0.5 * ((t - t_peak) / σ)^2)
end

# --- Reactive drawing: flare glow ---

function draw_flare_glow_reactive!(ax, x_obs, n_rings_obs, max_radius_obs, max_alpha_obs)
    MAX_RINGS = 15
    gc = Makie.Colors.red(GLOW_COLOR), Makie.Colors.green(GLOW_COLOR), Makie.Colors.blue(GLOW_COLOR)
    for i in MAX_RINGS:-1:1
        ring_color = @lift begin
            nr = $n_rings_obs
            ma = $max_alpha_obs
            if i > nr
                RGBAf(0, 0, 0, 0)
            else
                alpha = ma * (1 - (i - 1) / nr)^2
                RGBAf(gc[1], gc[2], gc[3], alpha)
            end
        end
        ring_circle = @lift Circle(Point2f($x_obs, 0), i > $n_rings_obs ? 0f0 : Float32($max_radius_obs * i / $n_rings_obs))
        poly!(ax, ring_circle; color=ring_color, strokewidth=0)
    end
    core_circle = @lift Circle(Point2f($x_obs, 0), Float32($max_radius_obs * 0.25))
    poly!(ax, core_circle; color=@lift(RGBAf(1.0, 1.0, 0.85, $max_alpha_obs)), strokewidth=0)
end

# --- Reactive drawing: centroid marker ---

function draw_centroid_marker_reactive!(ax, x_obs, color_obs; size=0.12)
    inner_r = size * 0.65
    fill_color = @lift Makie.Colors.weighted_color_mean(0.5, $color_obs, Makie.Colors.colorant"white")
    poly!(ax, @lift(Circle(Point2f($x_obs, 0), size)); color=fill_color, strokecolor=color_obs, strokewidth=2.5)
    poly!(ax, @lift(Circle(Point2f($x_obs, 0), inner_r)); color=:white, strokewidth=0)
end

# --- Gaia panel background (reactive PSF + static jet outline) ---

function draw_gaia_background!(ax, psf_center_x, flux_obs)
    nx, ny = 401, 401
    σ = 0.9
    xs_base = range(-4, 4, nx)
    ys = range(-2, 2, ny)
    z_base = [exp(-0.5 * ((x / σ)^2 + (y / σ)^2)) for x in xs_base, y in ys]

    psf_xs = @lift range(-4 + $psf_center_x, 4 + $psf_center_x, nx)
    z = @lift ($flux_obs / S_TOTAL) .* z_base
    cmap = RGBAf.(Makie.to_colormap(:afmhot), 0.7)
    hm = contourf!(ax, psf_xs, ys, z, colormap=cmap, levels=range(0, 1.2, 15))
    translate!(hm, 0, 0, -100)

    draw_jet_outline!(ax; origin=(-0.5, 0), scale=0.33)

    textglow!(ax, (0, 0.5), text=rich(rich("Gaia", font=:italic), " view:\nconvolved with PSF"),
          align=(:left, :center), color=:black, offset=(2, 0), space=:relative,
          glowcolor=(:white, 0.5), glowwidth=2)
end

# =====================================================================
# Main: build figure with Observables, then bake
# =====================================================================

const X_FLARE = 0.5  # fixed flare position

params = Observable((t=T_PEAK, flare_frac=0.5))

# --- Derived observables ---
const S_TOTAL = F0 + A
F0_eff = @lift (1 - $params.flare_frac) * S_TOTAL
A_eff = @lift $params.flare_frac * S_TOTAL
flux = @lift $F0_eff + $A_eff * exp(-0.5 * (($params.t - T_PEAK) / SIGMA_T)^2)
flux_frac = @lift clamp(($flux - $F0_eff) / max($A_eff, 1e-10), 0, 1)

centroid_x = @lift (($flux - $F0_eff) * X_FLARE + $F0_eff * X_CONST) / $flux

n_rings = @lift round(Int, 6 + 9 * $flux_frac)
max_radius = @lift 0.2 + 0.3 * $flux_frac
max_alpha = @lift 0.5 + 0.4 * $flux_frac

state_color = @lift Makie.Colors.weighted_color_mean(Float64($flux_frac), COLOR_BRIGHT, COLOR_FAINT)

gaia_centroid_x = @lift -0.5 + $centroid_x * 0.33

# --- Figure ---
fig = Figure(fontsize=15, size=(500, 550))

# --- Panel 1: Light curve ---
ts = range(T_PEAK - 4, T_PEAK + 4, 500)
ax_lc = Axis(fig[1, 1],
    xautolimitmargin=(0, 0),
    limits=(extrema(ts), (0, S_TOTAL * 1.1)),
    xlabel="Time", ylabel="Flux",
    width=400, height=120,
    topspinevisible=false, rightspinevisible=false,
    leftspinecolor=:gray80, bottomspinecolor=:gray80,
    xgridvisible=false, ygridvisible=false,
    xticklabelsvisible=false, yticklabelsvisible=false,
    xticksvisible=false, yticksvisible=false)

hlines!(ax_lc, [S_TOTAL]; color=:gray70, linewidth=1, linestyle=:dash)
text!(ax_lc, last(ts), S_TOTAL; text="S" * subscript("const") * " + S" * subscript("flare"), align=(:right, :top), offset=(-4, -2), color=:gray50)

hlines!(ax_lc, @lift([$F0_eff]); color=:gray70, linewidth=1, linestyle=:dash)
text!(ax_lc, @lift(Point2f(last(ts), $F0_eff));
    text="S" * subscript("const"), align=(:right, :bottom), offset=(-4, 2), color=:gray50)

lc_flux = @lift [$F0_eff + $A_eff * exp(-0.5 * ((t - T_PEAK) / SIGMA_T)^2) for t in ts]
lines!(ax_lc, ts, lc_flux, color=:gray30, linewidth=2)
text!(ax_lc, 0.06, 0.95, text="Optical lightcurve",
      align=(:left, :top), offset=(2, -2), color=:black, space=:relative)

scatter!(ax_lc, @lift([$params.t]), @lift([$flux]);
    color=state_color, markersize=14, strokewidth=1.5, strokecolor=:white)

# --- Panel 2: Physical jet ---
ax_phys = Axis(fig[2, 1], aspect=DataAspect(),
    width=400, height=300,
    limits=(-0.25, 2.8, -0.9, 0.65))
hidedecorations!(ax_phys)
hidespines!(ax_phys)

draw_jet!(ax_phys)

draw_flare_glow_reactive!(ax_phys, Observable(X_FLARE), n_rings, max_radius, max_alpha)

draw_centroid_marker_reactive!(ax_phys, centroid_x, state_color; size=0.14)

lines!(ax_phys, [X_FLARE, X_FLARE], [-0.55, -0.2];
    color=:gray40, linewidth=1, linestyle=:dot)
text!(ax_phys, X_FLARE, -0.6;
    text=rich("x", subscript("flare")), align=(:center, :top), color=:black)

lines!(ax_phys, [X_CONST, X_CONST], [-0.55, 0]; color=:gray40, linewidth=1, linestyle=:dot)
text!(ax_phys, X_CONST, -0.6; text=rich("x", subscript("const")),
    align=(:center, :top), color=:black)

lines!(ax_phys, @lift([$centroid_x, $centroid_x]), [-0.55, -0.2];
    color=state_color, linewidth=1, linestyle=:dot)
text!(ax_phys, @lift(Point2f($centroid_x, -0.7));
    text=rich("Gaia", font=:italic) * " centroid",
    align=(:center, :top), color=state_color)

# --- Panel 3: Gaia PSF view ---
ax_gaia = Axis(fig[3, 1], aspect=DataAspect(),
    width=400, height=150,
    limits=(-1.5, 0.8, -0.3, 0.3))
hidedecorations!(ax_gaia)
hidespines!(ax_gaia)

draw_gaia_background!(ax_gaia, gaia_centroid_x, flux)

draw_centroid_marker_reactive!(ax_gaia, gaia_centroid_x, state_color; size=0.055)

# --- Layout ---
rowgap!(fig.layout, 5)
resize_to_layout!(fig)

# --- Bake ---
outdir = joinpath(ROOT, "figs/interactive/corejet_cartoon")
bake_html(
    params => (
        (@o _.flare_frac) => 0.1:0.1:0.9,
        (@o _.t) => range(T_PEAK - 4, T_PEAK + 4, 20),
    );
    blocks=[ax_lc, ax_phys, ax_gaia],
    outdir,
)
write(joinpath(outdir, "layout.js"), """
const HEADER = 'Schematic of how <span style="color: orange;">optical flares</span> drive systematic <i>Gaia</i> position shifts';
const TITLE = 'Gaia centroid motion';
const LAYOUT = ["A S", "B S", "C S"];
const MAXWIDTH = '900px';
const AUTOLOOP = {
  control: 't',
  seconds: 4,
  range: [-100, 100]
};
""")

@info "Saved to $outdir"
