# Packages for the plot itself
using CairoMakie
using MakieExtra
using IntervalSets: ±
using VLBIFiles
using AxisKeysExtra
using PyFormattedStrings
using Statistics: median

const ROOT = joinpath(@__DIR__, "..")

include(joinpath(@__DIR__, "prepare_data.jl"))
include(joinpath(@__DIR__, "theme.jl"))

Base.get(c::Makie.Reverse, x) = get(c.data, 1 - x)

# --- Load data ---

enriched = load_and_fit(ROOT)
imgdir = joinpath(ROOT, "data/vlbi_images")

# --- Filter sources ---

sources = @p enriched filter(
    U.nσ(_.VG.vec) > 3 &&
    Circular.distance(U.value(_.VG_J_pa), 0u"°") < 45u"°"
) filter(!isnan(_.pm_J_pa)) filter(U.nσ(_.G.pm) > 3) filter(HIGHLIGHTED_SOURCES_PRED(_.R.names.J2000))

println("Generating images for $(length(sources)) sources...")

# --- Generate per-source figures ---

seaborn_bright = Makie.ColorSchemes.colorschemes[:seaborn_bright]
color1 = seaborn_bright[1]  # blue

outdir = joinpath(ROOT, "figs/image_motion")
mkpath(outdir)

colors = Makie.wong_colors()
flux_colormap = Makie.ColorSchemes.diverging_gwr_55_95_c38_n256

for r in sources
    j2000 = r.R.names.J2000

    # Load VLBI FITS image as KeyedArray (axes in mas)
    fits_file = only(filter(f -> startswith(f, j2000), readdir(imgdir)))
    vlbi_img = VLBIFiles.load(KeyedArray, joinpath(imgdir, fits_file))

    # Compute axis limits (all in mas)
    vg_mas = ustrip.(u"mas", U.value(r.VG.vec))
    pm_mas = ustrip.(u"mas/yr", U.value(r.G.pm))
    maxVG = norm(vg_mas) + 2 * norm(pm_mas)
    lim = max(1.4maxVG, 3.0)

    colormap = Makie.ColorSchemes.afmhot

    fig = Figure()

    # --- Row 1: Gaia lightcurve ---
    symlog = SymLog(5)
    val = -U.nσ(r.mc.q_slope) * sign(r.mc.q_slope)
    cval = (symlog(U.value(val)) - symlog(-200)) / (symlog(200) - symlog(-200))
    src_color = get(Reverse(flux_colormap), cval)

    ax_lc = Axis(fig[1, 1]; limits=(nothing, nothing), xticks=2010:2030,
        width=300, height=300, xminorticksvisible=true, xminorticks=2010:(1 // 12):2030,
        yscale=SymLog(1), yticks=BaseMulTicks([1,2,5]), yminorticksvisible=true, yminorticks=BaseMulTicks(1:9),
        xlabel="Date", ylabel=rich("Gaia", font=:italic) * " G-band flux (cnt/s)")
    text!((0,1), space=:relative, align=(:left, :top), offset=(2, -2), color=:black, text="$j2000, " * rich("Gaia", font=:italic) * " lightcurve")
    nflux = @p r.PHOTg map(_.flux) median
    fplt = FPlot(r.PHOTg, @o(yeardecimal(_.time)), @o(ustrip(u"s^-1", _.flux));
        color=:black, markersize=5, linewidth=1)
    multiplot!(ax_lc, (scatterlines, rangebars), fplt, label="Measurements")
    lines!(ax_lc, range(2010..2020, 50),
        t -> ustrip(u"s^-1", nflux) * r.mc.q_ref / (r.mc.q_ref + r.mc.q_slope * (t - GAIA_DR3_EPOCH));
        to_xy_attrs(autolimits=false)..., linewidth=3, linestyle=:dash, color=src_color,
        label="Trend: " * rich(r.mc.q_slope > 0 ? "fading" : "rising", color=src_color))
    Legend(fig[2, 1], ax_lc, merge=true)

    # --- VLBI image panel ---
    sz = 300
    ax_img = Axis(fig[1, 2]; limits=(0 ± lim, 0 ± lim), width=sz, height=sz,
        backgroundcolor=get(colormap, 0))
    plt = image!(ax_img, vlbi_img; colorscale=SymLog(1e-3; vmin=0), colormap)
    translate!(plt, 0, 0, -100)
    hidedecorations!(ax_img)
    scalebar!(1u"mas"; color=:white, position=(0.15, 0.9))
    text!(ax_img, (0, 1), space=:relative, align=(:left, :top),
        offset=(2, -2), color=:white, text=f"{j2000}, VLBI image")

    scatter!(ax_img, (0, 0); color=color1, markersize=20,
        marker=marker_lw(:cross, 0.2), label="VLBI position\n(shown at peak)")

    arrowlines!(ax_img, Point2f[(0, 0), 1.3maxVG .* sincos(r.J.pa)]; color=color1,
        to_xy_attrs(autolimits=false)..., linestyle=:dash, label="VLBI jet direction")
    poly!(ax_img, U.boundary(ustrip(u"mas", r.VG.vec); mul=2);
        color=(:black, 0.07), strokewidth=1, strokecolor=:black, label=rich("Gaia", font=:italic) * " position\n(2σ uncertainty)")
    poly!(ax_img, U.boundary(ustrip(u"mas", r.VG.vec); mul=2.1);
        color=(:black, 0.0), strokewidth=1, strokecolor=:white)

    nyrs = 1.5
    arrowlines!(ax_img, [vg_mas - nyrs * pm_mas, vg_mas + nyrs * pm_mas]; color=colors[3],
        arrowstyle="-|>", linewidth=2, markersize=15,
        label=rich("Gaia", font=:italic) * " proper motion\n(over 3 years)")

    Legend(fig[2, end], content(fig[1, end]); merge=true, tellheight=true, orientation=:horizontal, nbanks=2)

    rowgap!(fig.layout, 0)
    resize_to_layout!()
    save(joinpath(outdir, "$j2000.pdf"), fig)
    println("  $j2000")
end

println("Done. Saved figures to $outdir")
