using CairoMakie

include(joinpath(@__DIR__, "cartoon_common.jl"))

function jet_base!(ax, vlbi_x, gaia_x)
    vlbi_color = Makie.Colors.colorant"#0066cc"
    gaia_color = Makie.Colors.colorant"#ff9900"

    draw_jet!(ax)

    # VLBI marker
    vlbi_fill = Makie.Colors.weighted_color_mean(0.4, vlbi_color, colorant"white")
    poly!(ax, Circle(Point2f(vlbi_x, 0), 0.14f0), color=vlbi_fill, strokecolor=vlbi_color, strokewidth=2)
    lines!(ax, [vlbi_x, vlbi_x], [-0.6, 0.0], color=vlbi_color, linewidth=1.5, linestyle=:dot)
    text!(ax, vlbi_x, -0.75, text="VLBI", align=(:center, :center), color=vlbi_color)

    # Gaia marker
    poly!(ax, Circle(Point2f(gaia_x, 0), 0.18f0), color=colorant"#ffcc00", strokecolor=gaia_color, strokewidth=2.5)
    poly!(ax, Circle(Point2f(gaia_x, 0), 0.12f0), color=:white, strokewidth=0)
    lines!(ax, [gaia_x, gaia_x], [-0.6, 0.0], color=gaia_color, linewidth=1.5, linestyle=:dot)
    text!(ax, gaia_x, -0.75, text="Gaia", align=(:center, :center), color=gaia_color, font=:italic)
end

function jet_panel_vlbi_gaia!(ax, vlbi_x, gaia_x, arrow_color)
    jet_base!(ax, vlbi_x, gaia_x)

    # Arrow from VLBI to Gaia
    dx = gaia_x - vlbi_x
    arrows2d!(ax, [vlbi_x], [0.0], [dx], [0.0],
        color=arrow_color, shaftwidth=3, tipwidth=15)
end

function scale_panel_vlbi_gaia!(ax, vlbi_x, gaia_x)
    jet_base!(ax, vlbi_x, gaia_x)

    # Black dotted line from black hole down
    lines!(ax, [0, 0], [-0.6, 0.0], color=:black, linewidth=1.5, linestyle=:dot)

    # Bracket between black hole and VLBI dotted lines
    bracket_y = -0.55
    bar_color = :gray30
    lines!(ax, [0, vlbi_x], [bracket_y, bracket_y], color=bar_color, linewidth=1.5)
    text!(ax, vlbi_x / 2, bracket_y, text="< 1 mas", align=(:center, :bottom), color=bar_color)
end
