# Packages for the plot itself
using CairoMakie
using MakieExtra
using IntervalSets: ±
using PyFormattedStrings
using Statistics

const ROOT = joinpath(@__DIR__, "..")

include(joinpath(@__DIR__, "prepare_data.jl"))
include(joinpath(@__DIR__, "theme.jl"))

# Load data from scratch --- all sources with jet directions, not just those with light curves
allsrc = load_allsrc(ROOT)

# --- Figure: V-G offset angle from the jet histogram ---

nσ = 3
data = @p allsrc filter(U.nσ(_.VG.vec) > nσ)
angles = @p data map(mod(rad2deg(U.value(_.VG_J_pa)), -90..270))

seaborn_bright = Makie.ColorSchemes.colorschemes[:seaborn_bright]
color1 = seaborn_bright[1]  # blue
color2 = seaborn_bright[2]  # orange
color1_dark = @set convert(Makie.Colors.LCHuv, color1).l *= 0.6
color2_dark = @set convert(Makie.Colors.LCHuv, color2).l *= 0.7

bins = -90:15:270
# 24 bins: 3 gray + 6 blue + 6 gray + 6 orange + 3 gray
barcolors = [
    fill(Makie.coloralpha(Makie.Colors.Gray(0), 0.15), 3);
    fill(Makie.coloralpha(color1, 0.75), 6);
    fill(Makie.coloralpha(Makie.Colors.Gray(0), 0.15), 3+3);
    fill(Makie.coloralpha(Makie.Colors.Gray(0), 0.15), 6);
    fill(Makie.coloralpha(Makie.Colors.Gray(0), 0.15), 3);
]

fp = FPlot(data, (@o mod(rad2deg(U.value(_.VG_J_pa)), -90..270)),
    bins=Ref(bins),
    axis=(yautolimitmargin=(0, 0.1), xticks=-360:45:360,
          xlabel="V-G offset angle from the jet, Ψ (°)",
          ylabel="Number of AGNs",
          limits=(-90..270, (0, nothing))))

include(joinpath(@__DIR__, "plot_vlbi_gaia_sketch.jl"))

update_theme!(fontsize=17)
let
    fig = Figure(size=(450, 500))

    axh,_ = axplot(hist)(fig[1,1], fp, color=barcolors)
    stephist!(fp, color=:black)

    text!((0, 30), align=(:center, :bottom), text=f"{count(x -> x ∈ (0 ± 45), angles):d} AGNs", color=:white)

    # Inset jet sketch drawn natively via jet_panel_vlbi_gaia!
    inset_rect = Rect((200, 385), 1.8 .* (120, 50))
    ax_inset = Axis(fig, bbox=inset_rect, aspect=DataAspect(), limits=(nothing, (-1, nothing)))
    hidedecorations!(ax_inset)
    hidespines!(ax_inset)
    jet_panel_vlbi_gaia!(ax_inset, 0.6, 2.2, Makie.Colors.colorant"#023EFF")

    colorbuffer(fig)
    r1 = MakieExtra.fullproject(axh, (0, 100))
    r2 = MakieExtra.fullproject(ax_inset, (1.6, -0.15))
    p = lines!(fig.scene, [r1, r2], color=color1, linestyle=:dash); translate!(p, 0,0,1000)

    outpath = joinpath(ROOT, "figs/offset_hist.pdf")
    mkpath(dirname(outpath))
    save(outpath, fig)
    println("Saved figure to $outpath")
end
