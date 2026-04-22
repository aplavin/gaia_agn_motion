# Shared helpers for cartoon/sketch scripts.
# Included by: plot_corejet_cartoon.jl, plot_corejet_cartoon_interactive.jl, plot_vlbi_gaia_sketch.jl

# --- Colors ---
const GLOW_COLOR   = Makie.Colors.colorant"#FFD700"
const JET_STROKE   = Makie.Colors.colorant"#666666"
const _DIVERGING_CMAP = Makie.ColorSchemes.diverging_gwr_55_95_c38_n256
const COLOR_FAINT  = get(_DIVERGING_CMAP, 1)   # red  — faint / fading
const COLOR_BRIGHT = get(_DIVERGING_CMAP, 0)    # green — bright / rising

# --- Quadratic bezier sampling ---
qbezier(p0, p1, p2; n=100) =
    [Point2f((1-t)^2 .* p0 .+ 2(1-t)*t .* p1 .+ t^2 .* p2) for t in range(0, 1, n)]

# --- Jet geometry ---
const JET_LEN = 2.8
const JET_CP  = JET_LEN * 0.375
jet_upper(x) = 0.3 + 0.15 * sqrt(x)
jet_lower(x) = -0.3 - 0.15 * sqrt(x)

# --- Draw the full jet shape (5 gradient fill layers + boundary curves + core) ---
function draw_jet!(ax)
    for j in 0:4
        alpha = (40 - j * 8) / 255
        shrink = j * 0.04

        bottom = qbezier(
            (0.0, -0.1 - shrink),
            (JET_CP, jet_lower(JET_CP) + shrink),
            (JET_LEN, jet_lower(JET_LEN) + shrink),
        )
        top_curve = qbezier(
            (JET_LEN, jet_upper(JET_LEN) - shrink),
            (JET_CP, jet_upper(JET_CP) - shrink),
            (0.0, 0.1 + shrink),
        )

        pts = vcat(
            [Point2f(0, 0.1 + shrink), Point2f(0, -0.1 - shrink)],
            bottom,
            top_curve,
        )
        poly!(ax, pts, color=RGBAf(200/255, 150/255, 100/255, alpha), strokewidth=0)
    end

    lines!(ax, qbezier((0, 0.3), (JET_CP, jet_upper(JET_CP)), (JET_LEN, jet_upper(JET_LEN))),
        color=JET_STROKE, linewidth=1)
    lines!(ax, qbezier((0, -0.3), (JET_CP, jet_lower(JET_CP)), (JET_LEN, jet_lower(JET_LEN))),
        color=JET_STROKE, linewidth=1)

    poly!(ax, Circle(Point2f(0, 0), 0.25f0), color=:black, strokewidth=0)
end

# --- Draw jet outline only (boundary curves + core), at given origin and scale ---
function draw_jet_outline!(ax; origin=Point2f(0, 0), scale=1.0)
    ox, oy = origin
    lines!(ax, qbezier(
            (ox, oy + jet_upper(0) * scale),
            (ox + JET_CP * scale, oy + jet_upper(JET_CP) * scale),
            (ox + JET_LEN * scale, oy + jet_upper(JET_LEN) * scale)),
        color=JET_STROKE, linewidth=1)
    lines!(ax, qbezier(
            (ox, oy + jet_lower(0) * scale),
            (ox + JET_CP * scale, oy + jet_lower(JET_CP) * scale),
            (ox + JET_LEN * scale, oy + jet_lower(JET_LEN) * scale)),
        color=JET_STROKE, linewidth=1)
    poly!(ax, Circle(Point2f(ox, oy), 0.25f0 * scale), color=:black, strokewidth=0)
end

# --- Flare glow (concentric semi-transparent rings) ---
function draw_flare_glow!(ax, x, y, color, n_rings, max_radius, max_alpha)
    gc = Makie.Colors.red(color), Makie.Colors.green(color), Makie.Colors.blue(color)
    for i in n_rings:-1:1
        r = max_radius * i / n_rings
        alpha = max_alpha * (1 - (i - 1) / n_rings)^2
        poly!(ax, Circle(Point2f(x, y), r),
              color=RGBAf(gc..., alpha), strokewidth=0)
    end
    poly!(ax, Circle(Point2f(x, y), max_radius * 0.25),
          color=RGBAf(1.0, 1.0, 0.85, max_alpha), strokewidth=0)
end

# --- Centroid marker (nested circles, Gaia style) ---
function draw_centroid_marker!(ax, x, y, color; size=0.12)
    fill_color = Makie.Colors.weighted_color_mean(0.5, color, Makie.Colors.colorant"white")
    poly!(ax, Circle(Point2f(x, y), size), color=fill_color, strokecolor=color, strokewidth=2.5)
    poly!(ax, Circle(Point2f(x, y), size * 0.65), color=:white, strokewidth=0)
end
