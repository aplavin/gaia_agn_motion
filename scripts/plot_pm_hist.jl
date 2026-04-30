# Packages for the plot itself
using CairoMakie
using MakieExtra

const ROOT = joinpath(@__DIR__, "..")

include(joinpath(@__DIR__, "prepare_data.jl"))
include(joinpath(@__DIR__, "theme.jl"))

# Load data
enriched = load_and_fit(ROOT)

# Bins for angular histogram
bins = range(-π/2..(2π-π/2), 14+1)
colormap = Makie.ColorSchemes.diverging_gwr_55_95_c38_n256

pm_angle_from_jet(r) = r.pm_J_pa

function filter_and_bin(enriched, bins; nσ_VG, nσ_pm, angle_filt=Returns(true), pm_angle=pm_angle_from_jet)
    data = @p enriched filter(U.nσ(_.VG.vec) > nσ_VG)
    data = @p data filter(angle_filt(Circular.distance(U.value(_.VG_J_pa), 0u"°")))
    @p data filter(!isnan(_.pm_J_pa)) filter(U.nσ(_.G.pm) > nσ_pm) mapinsert(x=pm_angle(_)) mapinsert(bin=something(findlast(≤(_.x), bins))) mapinsert(xbin=bins[_.bin] + step(bins)/2) sort(by=-U.nσ(_.mc.q_slope))
end

function plot_pm_hist!(pos, wbins, bins, colormap; title=nothing, xlabel=rich("Gaia", font=:italic) * " proper motion angle from the jet")
    isempty(wbins) && return
    ax = Axis(pos;
        xticks=Makie.AngularTicks(rad2deg(1), "°", [(0,4)]),
        xlabel,
        ylabel="Count", xautolimitmargin=(0,0), yautolimitmargin=(0,0.05),
        title=isnothing(title) ? "" : title)
    plt = barplot!(ax, wbins.xbin, fill(1, 1:length(wbins)); stack=wbins.bin,
        gap=0.03,
        color=(@p wbins map(-U.nσ(_.mc.q_slope) * sign(_.mc.q_slope))),
        colorrange=(-1, 1) .* 200, colorscale=SymLog(5), colormap=Reverse(colormap))
    translate!(plt, 0,0,-100)
    stephist!(ax, wbins.x, bins=bins, color=:black, linewidth=1)
    hlines!(ax, 0, color=:black, linewidth=2)
    hidespines!(ax)
    return plt
end

# --- Main figure (nσ > 3, angle < 45°) ---

wbins = filter_and_bin(enriched, bins; nσ_VG=3, nσ_pm=3, angle_filt=<(45u"°"))
fig = Figure()
plt = plot_pm_hist!(fig[1,1], wbins, bins, colormap)
Colorbar(fig[1,2][1,1], plt, tickformat=xs -> string.(xs) .* "σ", height=200, tellheight=false)
Label(fig[1,2][1,2][1,1], rich("↑ Rising"; color=get(colormap, 0)); halign=:left, tellheight=false)
Label(fig[1,2][1,2][3,1], rich("↓ Fading"; color=get(colormap, 1)); halign=:left, tellheight=false)
Label(fig[1,2][1,2][2,1], rich("Gaia", font=:italic) * " flux trend", rotation=π/2)

outpath = joinpath(ROOT, "figs/pm_hist_color.pdf")
mkpath(dirname(outpath))
save(outpath, fig)
println("Saved figure to $outpath")
