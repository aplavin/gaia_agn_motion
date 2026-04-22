# Packages for the plot itself
using CairoMakie
using MakieExtra
using IntervalSets: ±
using Statistics

const ROOT = joinpath(@__DIR__, "..")

include(joinpath(@__DIR__, "prepare_data.jl"))
include(joinpath(@__DIR__, "theme.jl"))
include(joinpath(@__DIR__, "plot_vlbi_gaia_sketch.jl"))

# -- Load data --

enriched = load_and_fit(ROOT)

# -- Compute derived fields --

data = @p let
    enriched
    filter(U.nσ(_.VG.vec) > 3 && Circular.distance(U.value(_.VG_J_pa), 0u"°") < 45u"°")
    mapinsert(
        x1_pm = _.mc.x_flare_pm_1d,
        xvg = dot(_.VG.vec, _.J.nvec) |> u"mas" |> ustrip,
        σ2d_pm = flare_major_std(_.mc.pm.x_flare),
    )
end

# -- Filter: 2D major-axis σ < 0.5 mas --

plotdata = @p let
    data
    filter(_.σ2d_pm < 0.5)
    sort(by=MCM.pmean(_.x1_pm))
    @insert __.ix = 1:length(__)
end

println("Flare-location horizontal plot (PM-only, σ<0.5): $(length(plotdata)) sources")

# -- Plot --

const PX_PER_POINT = 8

# -- Layout parameters --
VLBI_X = 1.1          # sketch x-coordinate of VLBI position (= data x=0)
GAIA_X = 2.4          # sketch x-coordinate of Gaia marker
DATA_XLIMS = (-1.1, 3.1)  # x-limits on the bottom axis (mas)

flare_color = Makie.ColorSchemes.seaborn_bright[4]
color = flare_color
xs = [0.7, 1.4]

mkpath(joinpath(ROOT, "figs"))

# -- Top panel: jet schematic --

fig_top = Figure(fontsize=15, size=(480, 200))
ax_top = Axis(fig_top[1, 1], aspect=DataAspect())
hidedecorations!(ax_top)
hidespines!(ax_top)
jet_base!(ax_top, VLBI_X, GAIA_X)
let fx = first(xs) * 1.1
    lines!(ax_top, [fx, fx], [0.0, 0.6], color=flare_color, linewidth=1.5, linestyle=:dot)
    scatter!(ax_top, [fx], [0.0]; color=flare_color, marker=:star5, markersize=19)
    text!(ax_top, fx, 0.65; text="Optical flare", align=(:center, :bottom), color=flare_color)
end
# Red bracket from B.H. (x=0) to 1.5*VLBI distance
let bh_x = 0.0, end_x = 1.5 * VLBI_X, bracket_y = 0.55
    lines!(ax_top, [bh_x, bh_x], [0.0, bracket_y], color=flare_color, linewidth=1.5, linestyle=:dot)
    lines!(ax_top, [end_x, end_x], [0.0, bracket_y], color=flare_color, linewidth=1.5, linestyle=:dot)
    lines!(ax_top, [bh_x, end_x], [bracket_y, bracket_y], color=flare_color, linewidth=1.5)
end
scatter!(ax_top, [0], [1.1]; color=:transparent)
resize_to_layout!(fig_top)
outpath_top = joinpath(ROOT, "figs/flare_location_h_top.pdf")
save(outpath_top, fig_top)
println("Saved figure to $outpath_top")

fig_bot = Figure(fontsize=15, size=(480, 400))
fplt2 = FPlot(plotdata,
    AxFunc(label="Offset along the jet, from the VLBI position (mas)", ticks=-10:1:10, @o _.x1_pm),
    AxFunc(@o _.ix))
(ax2, _), = multiplot(fig_bot[1, 1], (axplot(scatter) => (;marker=:star5, markersize=11), rangebars => (;direction=:x, linewidth=1)), fplt2; color, label="Optical flare")
fplt_vg2 = FPlot(plotdata,
    AxFunc(@o _.xvg),
    AxFunc(@o _.ix + 0.25))
multiplot!(ax2, (scatter, rangebars => (;direction=:x, linewidth=1)), fplt_vg2, color=Cycled(2), label=rich("Gaia", font=:italic) * " centroid")
ax2.height = length(plotdata) * PX_PER_POINT
ax2.ylabel = "Blazars, sorted by flare location"
ax2.yticklabelsvisible = false
ax2.yticksvisible = false
axislegend(ax2, position=:rb, merge=true)
vlines!(ax2, 0, color=colorant"#0066cc", linestyle=:dot)
xlims!(ax2, DATA_XLIMS...)
resize_to_layout!(fig_bot)
outpath_bot = joinpath(ROOT, "figs/flare_location_h_bot.pdf")
save(outpath_bot, fig_bot)
println("Saved figure to $outpath_bot")
