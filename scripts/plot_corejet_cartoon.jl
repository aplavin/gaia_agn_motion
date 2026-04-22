using CairoMakie
using MakieExtra
using IntervalSets

const ROOT = joinpath(@__DIR__, "..")

include(joinpath(@__DIR__, "cartoon_common.jl"))

COLOR_MID    = Makie.Colors.weighted_color_mean(0.5, COLOR_FAINT, COLOR_BRIGHT)
STATE_COLORS = [COLOR_BRIGHT, COLOR_MID, COLOR_FAINT]

# --- Gaussian flare light curve ---

function flare_flux(t; t_peak=2.5, sigma=1.2, A=1.0, F0=0.2)
    return F0 + A * exp(-0.5 * ((t - t_peak) / sigma)^2)
end

FLARE_PEAK = 2.5
# Three states on the decay side: bright (peak) -> intermediate -> faint (baseline)
T_BRIGHT  = FLARE_PEAK       # at peak
T_MID     = 4.0              # midway down
T_FAINT   = 5.5              # near baseline

# --- Drawing functions ---

function draw_physical_panel!(ax, state_color, n_rings, max_radius, max_alpha,
                              centroid_x, flare_x; show_labels=false, show_centroid_label=false)
    draw_jet!(ax)
    draw_flare_glow!(ax, flare_x, 0.0, GLOW_COLOR, n_rings, max_radius, max_alpha)
    draw_centroid_marker!(ax, centroid_x, 0.0, state_color; size=0.14)

    if show_centroid_label
        lines!(ax, [centroid_x, centroid_x], [-0.55, -0.2], color=state_color, linewidth=1, linestyle=:dot)
        text!(ax, centroid_x, -0.6, text=rich("Gaia", font=:italic) * " centroid",
              align=(:center, :top), color=state_color)
    end

    if show_labels
        lines!(ax, [flare_x, flare_x], [-0.55, -0.2], color=:gray40, linewidth=1, linestyle=:dot)
        text!(ax, flare_x, -0.6, text=rich("x", subscript("var")),
              align=(:center, :top), color=:black)

        r_const_x = 1.8
        lines!(ax, [r_const_x, r_const_x], [-0.55, 0], color=:gray40, linewidth=1, linestyle=:dot)
        text!(ax, r_const_x, -0.6, text=rich("x", subscript("const")),
              align=(:center, :top), color=:black)
    end
end

function draw_lightcurve_panel!(ax)
    ts = range(0, 7, 500)
    flux = flare_flux.(ts)

    lines!(ax, ts, flux, color=:gray30, linewidth=2)

    for (t, c) in zip([T_BRIGHT, T_MID, T_FAINT], [COLOR_BRIGHT, COLOR_MID, COLOR_FAINT])
        scatter!(ax, [t], [flare_flux(t)], color=c, markersize=14, strokewidth=1.5, strokecolor=:white)
    end

    text!(ax, 0, 1, text="Optical lightcurve",
          align=(:left, :top), offset=(2, -2), color=:black, space=:relative)
end

function draw_gaia_view_panel!(ax, centroid_xs, centroid_colors=STATE_COLORS)
    # 2D Gaussian PSF as heatmap
    nx, ny = 201, 201
    sigma = 0.9
    xs = range(-1.2, 1.2, nx)
    ys = range(-0.9, 0.9, ny)
    z = [exp(-0.5 * ((x / sigma)^2 + (y / sigma)^2)) for x in xs, y in ys]

    cmap = RGBAf.(Makie.to_colormap(:afmhot), 0.7)
    hm = contourf!(ax, xs, ys, z, colormap=cmap, levels=0:0.1:1.5)
    translate!(hm, 0, 0, -100)

    draw_jet_outline!(ax; origin=(-0.5, 0), scale=0.33)

    # Centroid markers
    for (x, c) in zip(centroid_xs, centroid_colors)
        draw_centroid_marker!(ax, x, 0.0, c; size=0.055)
    end

    # Connect centroids with a thin line
    lines!(ax, [centroid_xs[1], centroid_xs[end]], [0.0, 0.0], color=:gray50, linewidth=1, linestyle=:dash)

    # Double-headed arrows with labels
    arrow_y = -0.24
    mid_x = (centroid_xs[1] + centroid_xs[end]) / 2

    arrowlines!([Point2f(mid_x, arrow_y), Point2f(centroid_xs[1] - 0.02, arrow_y)];
        color=[:gray, COLOR_BRIGHT], linewidth=1.5, markersize=10, arrowstyle="->")
    arrowlines!([Point2f(mid_x, arrow_y), Point2f(centroid_xs[end] + 0.02, arrow_y)];
        color=[:gray, COLOR_FAINT], linewidth=1.5, markersize=10, arrowstyle="->")

    text!(ax, (centroid_xs[1] - 0.05, arrow_y), text=rich("Centroid motion:", color=:black) * " brightening",
          align=(:right, :center), color=COLOR_BRIGHT)
    text!(ax, (centroid_xs[end] + 0.05, arrow_y), text="fading",
          align=(:left, :center), color=COLOR_FAINT)

    lims = ax.limits[]
    text!(ax, (lims[1], lims[4]), text=rich(rich("Gaia", font=:italic), " view:\nconvolved with PSF"),
          align=(:left, :top), color=:black, offset=(2, -2))
end

function draw_connectors!(fig, ax_lc, ax_states, ax_gaia, times, colors, gaia_centroid_xs)
    colorbuffer(fig)

    for (i, (t, color)) in enumerate(zip(times, colors))
        p_top = MakieExtra.fullproject(ax_lc, (t, flare_flux(t)))
        p_bot = MakieExtra.fullproject(ax_states[i], (1.4, 0.75))
        p = lines!(fig.scene, [p_top, p_bot], color=(color, 0.5), linestyle=:dash, linewidth=1.5)
        translate!(p, 0, 0, 1000)

        p_phys_bot = MakieExtra.fullproject(ax_states[i], (1.4, -0.6))
        p_gaia_top = MakieExtra.fullproject(ax_gaia, (gaia_centroid_xs[i], 0.18))
        p2 = lines!(fig.scene, [p_phys_bot, p_gaia_top], color=(color, 0.5), linestyle=:dash, linewidth=1.5)
        translate!(p2, 0, 0, 1000)
    end
end

# --- Main ---

fig = Figure(fontsize=15)

# --- Row 1: Light curve ---
ax_lc = Axis(fig[1, 1:3],
    xlabel="Time", ylabel="Flux",
    height=120,
    topspinevisible=false, rightspinevisible=false,
    leftspinecolor=:gray80, bottomspinecolor=:gray80,
    xgridvisible=false, ygridvisible=false)
ax_lc.xticklabelsvisible = false
ax_lc.yticklabelsvisible = false
ax_lc.xticksvisible = false
ax_lc.yticksvisible = false
draw_lightcurve_panel!(ax_lc)

# --- Row 2: Three physical state panels ---
flare_x = 0.5
states = [
    (color=COLOR_BRIGHT, n_rings=15, max_r=0.50, alpha=0.9,  centroid_x=0.8),
    (color=COLOR_MID,    n_rings=10, max_r=0.35, alpha=0.7,  centroid_x=1.3),
    (color=COLOR_FAINT,  n_rings=6,  max_r=0.20, alpha=0.5,  centroid_x=1.8),
]

ax_states = Axis[]
for (i, s) in enumerate(states)
    ax = Axis(fig[2, i], aspect=DataAspect(),
              height=130,
              limits=(-0.25, 2.8, -1, 0.85))
    hidedecorations!(ax)
    hidespines!(ax)
    draw_physical_panel!(ax, s.color, s.n_rings, s.max_r, s.alpha,
                         s.centroid_x, flare_x;
                         show_labels=(i == 2), show_centroid_label=(i == 1))
    push!(ax_states, ax)
end
colgap!(fig.layout, 8)

# --- Row 3: Gaia's view ---
ax_gaia = Axis(fig[3, 1:3], aspect=DataAspect(),
               height=120,
               limits=(-1.2, 1.2, -0.3, 0.2))
hidedecorations!(ax_gaia)
hidespines!(ax_gaia)

gaia_centroid_xs = [-0.25, 0.0, 0.25]
draw_gaia_view_panel!(ax_gaia, gaia_centroid_xs)

# --- Layout tuning ---
rowgap!(fig.layout, 0)
resize_to_layout!(fig)

# --- Cross-panel connectors ---
draw_connectors!(fig, ax_lc, ax_states, ax_gaia, [T_BRIGHT, T_MID, T_FAINT], [COLOR_BRIGHT, COLOR_MID, COLOR_FAINT], gaia_centroid_xs)

outpath = joinpath(ROOT, "figs/corejet_cartoon.pdf")
mkpath(dirname(outpath))
save(outpath, fig)
println("Saved to $outpath")
