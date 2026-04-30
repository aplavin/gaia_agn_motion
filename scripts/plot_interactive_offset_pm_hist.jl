# Packages for the plot itself
import CairoMakie
using MakieExtra
using MakieBake
using IntervalSets: ±
using Statistics

const ROOT = joinpath(@__DIR__, "..")

include(joinpath(@__DIR__, "prepare_data.jl"))
include(joinpath(@__DIR__, "theme.jl"))

# Load data: offset histogram needs all sources, enriched for PM histogram
allsrc = load_allsrc(ROOT)
enriched = load_and_fit(ROOT)

# --- Constants ---

# Offset histogram bins & colors (from plot_offset_hist.jl)
offset_bins = deg2rad(1) * (-90:15:270)
offset_centers = [(offset_bins[i] + offset_bins[i+1])/2 for i in 1:length(offset_bins)-1]
seaborn_bright = Makie.ColorSchemes.colorschemes[:seaborn_bright]
color_blue = seaborn_bright[1]
color_blue_alpha = Makie.coloralpha(color_blue, 0.75)
color_gray_alpha = Makie.coloralpha(Makie.Colors.Gray(0), 0.15)

# PM histogram bins & colormap (from plot_pm_hist.jl)
pm_bins = range(-π/2..(2π-π/2), 14+1)
colormap = Makie.ColorSchemes.diverging_gwr_55_95_c38_n256

# --- Reused functions from plot_pm_hist.jl ---

pm_angle_from_jet(r) = r.pm_J_pa

function filter_and_bin(enriched, bins; nσ_VG, nσ_pm, angle_filt=Returns(true), pm_angle=pm_angle_from_jet)
    data = @p enriched filter(U.nσ(_.VG.vec) > nσ_VG)
    data = @p data filter(angle_filt(Circular.distance(U.value(_.VG_J_pa), 0u"°")))
    @p data filter(!isnan(_.pm_J_pa)) filter(U.nσ(_.G.pm) > nσ_pm) mapinsert(x=pm_angle(_)) mapinsert(bin=something(findlast(≤(_.x), bins))) mapinsert(xbin=bins[_.bin] + step(bins)/2) sort(by=-U.nσ(_.mc.q_slope))
end

# --- Dynamic bar colors for offset histogram ---

function offset_bar_colors(bins, max_VG_angle)
    s = step(bins)
    map(1:length(bins)-1) do i
        center = bins[i] + s/2
        Circular.distance(center, 0) <= max_VG_angle ? color_blue_alpha : color_gray_alpha
    end
end

# --- Figure setup ---

params = Observable((nσ_VG=3, nσ_pm=3, max_VG_angle=45u"°"))

fig = Figure(fontsize=20, figure_padding=50)
ax_left = Axis(fig[1,1];
    width=350, height=300,
    yautolimitmargin=(0, 0.1),
    xticks=Makie.AngularTicks(rad2deg(1), "°", [(0,4)]),
    xlabel="VLBI ➡️ Gaia offset direction\nrelative to the jet",
    ylabel="Number of AGNs",
    limits=(-π/2..3π/2, (0, nothing)))
ax_right = Axis(fig[1,2][1,1];
    width=350, height=300,
    xticks=Makie.AngularTicks(rad2deg(1), "°", [(0,4)]),
    xlabel=rich("Gaia", font=:italic) * " proper motion direction\nrelative to the jet",
    ylabel="Count",
    xautolimitmargin=(0,0), yautolimitmargin=(0, 0.05))
Colorbar(fig[1,2][2,1]; colorrange=(-1, 1) .* 200, scale=SymLog(5),
    colormap=Reverse(colormap), tickformat=xs -> string.(Int.(xs)) .* "σ",
    ticks=Int[-100, -10, 0, 10, 100],
    label="Optical flux trend:",
    vertical=false,
    width=300, tellwidth=false)
rowgap!(fig[1,2], 40)
hidespines!(ax_right)

# --- Left panel: offset histogram (reactive) ---

offset_result = @lift let
    data = @p allsrc filter(U.nσ(_.VG.vec) > $params.nσ_VG)
    angles = @p data map(mod(U.value(_.VG_J_pa), -π/2..3π/2))
    max_rad = ustrip(u"rad", $params.max_VG_angle)
    counts = [count(a -> offset_bins[i] <= a < offset_bins[i+1], angles) for i in 1:length(offset_bins)-1]
    colors = offset_bar_colors(offset_bins, max_rad)
    n_in = count(x -> x ∈ (0 ± max_rad), angles)
    (; counts, colors, n_in, angles)
end

barplot!(ax_left, offset_centers, @lift($offset_result.counts);
    color=@lift($offset_result.colors), gap=0)
stephist!(ax_left, @lift($offset_result.angles); bins=offset_bins, color=:black)
text!(ax_left, (0, 0); align=(:left, :center), offset=(0, 25), rotation=90u"°",
    text=@lift("$($offset_result.n_in) AGNs"), color=:white)

# --- Right panel: PM histogram (reactive) ---

pm_result = @lift let
    wbins = filter_and_bin(enriched, pm_bins; nσ_VG=$params.nσ_VG, nσ_pm=$params.nσ_pm, angle_filt=<($params.max_VG_angle))
    if isempty(wbins)
        (; xbin=Float64[], heights=Int[], bin=Int[], color=Float64[], x=Float64[])
    else
        (; xbin=collect(Float64, wbins.xbin),
           heights=fill(1, length(wbins)),
           bin=collect(Int, wbins.bin),
           color=Float64[(-U.nσ(r.mc.q_slope) * sign(r.mc.q_slope)) for r in wbins],
           x=collect(Float64, wbins.x))
    end |> StructArray
end

barplot!(ax_right, (@lift FPlot($pm_result, (@o _.xbin), (@o _.heights), stack=(@o _.bin), color=(@o _.color))),
    gap=0.03,
    colorrange=(-1, 1) .* 200, colorscale=SymLog(5), colormap=Reverse(colormap))
stephist!(ax_right, @lift($pm_result.x); bins=pm_bins, color=:black, linewidth=1)
hlines!(ax_right, 0; color=:black, linewidth=2)
on(params) do _; reset_limits!(ax_right); end

# --- Label block: filtering funnel counts ---

label_text = @lift let
    n_offset = $offset_result.n_in
    data_e = @p enriched filter(U.nσ(_.VG.vec) > $params.nσ_VG)
    data_e = @p data_e filter(Circular.distance(U.value(_.VG_J_pa), 0u"°") < $params.max_VG_angle)
    n_phot = length(data_e)
    data_pm = @p data_e filter(!isnan(_.pm_J_pa)) filter(U.nσ(_.G.pm) > $params.nσ_pm)
    n_pm = length(data_pm)
    n_phot_str = rich("$n_phot", color=:gray60)
    n_phot_str * " / " * rich("$n_offset", color=color_blue) * " with epoch photometry\n$n_pm / " * n_phot_str * " with proper motion > $($params.nσ_pm)σ"
end
Label(fig[2, 1], label_text; halign=:left)

colgap!(fig.layout, 50)
resize_to_layout!()

# --- Bake to HTML ---

outdir = joinpath(ROOT, "figs/interactive/offset_pm_hist")
rm(outdir; force=true, recursive=true)
bake_html(
    params => (
        (@o _.nσ_VG) => [3, 4, 5, 7, 10],
        (@o _.nσ_pm) => [0, 1, 2, 3, 4, 5, 7, 10],
        (@o ustrip(u"°", _.max_VG_angle)) => [15, 30, 45, 60, 90, 180],
    );
    blocks=[ax_left, fig[1,2], fig[2, 1]],
    outdir
)
write(joinpath(outdir, "layout.js"), """
const HEADER = 'Population of <i>Gaia</i> astrometric offsets and apparent proper motions';
const TITLE = 'Gaia offsets and motions';
const LAYOUT = ["A B S", ". C ."];
const MAXWIDTH = '1400px';
""")
println("Saved interactive plot to $outdir/index.html")
