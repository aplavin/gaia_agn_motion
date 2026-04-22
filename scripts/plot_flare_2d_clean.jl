# Clean per-source flare localization plots (PM-only)

# Packages for the plot itself
using CairoMakie
using MakieExtra
using VLBIFiles
using AxisKeysExtra
using PyFormattedStrings
using Statistics: mean, cov

const ROOT = joinpath(@__DIR__, "..")

include(joinpath(@__DIR__, "prepare_data.jl"))
include(joinpath(@__DIR__, "theme.jl"))

# -- Load data --

enriched = load_and_fit(ROOT)
imgdir = joinpath(ROOT, "data/vlbi_images")

# -- Filter sources --

sources = @p enriched filter(
    U.nσ(_.VG.vec) > 3 &&
    Circular.distance(U.value(_.VG_J_pa), 0u"°") < 45u"°" &&
    flare_major_std(_.mc.pm.x_flare) < 0.5
)

println("Found $(length(sources)) sources for clean 2D flare plots")

# -- Generate per-source figures --

seaborn_bright = Makie.ColorSchemes.colorschemes[:seaborn_bright]
vlbi_color = seaborn_bright[1]

outdir = joinpath(ROOT, "figs/flare_2d_clean")
mkpath(outdir)

colors = Makie.wong_colors()

for r in sources
    j2000 = r.R.names.J2000

    # Load VLBI FITS image as KeyedArray (axes in mas)
    fits_files = filter(f -> startswith(f, j2000), readdir(imgdir))
    if isempty(fits_files)
        println("  $j2000: no FITS file, skipping")
        continue
    end
    vlbi_img = VLBIFiles.load(KeyedArray, joinpath(imgdir, only(fits_files)))

    # -- Extract astrometric parameters --
    vg_mas = ustrip.(u"mas", U.value(r.VG.vec))
    pm_mas = ustrip.(u"mas/yr", U.value(r.G.pm))

    # Time range from lightcurve
    lc_times = @p r.PHOTg map(yeardecimal(_.time))
    t_start = minimum(lc_times)
    t_end = maximum(lc_times)

    # Flare extent for axis limits
    xf = r.mc.pm.x_flare
    flare_extent = max(abs(MCM.pmean(xf[1])), abs(MCM.pmean(xf[2]))) +
        2 * max(MCM.pstd(xf[1]), MCM.pstd(xf[2]))

    maxVG = max(norm(vg_mas) + 2 * norm(pm_mas), flare_extent + 0.5)
    lim = max(1.4maxVG, 3.0)

    # -- Figure --
    colormap = Makie.ColorSchemes.afmhot
    fig = Figure()

    # -- Col 1: Gaia lightcurve --
    ax_lc = Axis(fig[1, 1]; limits=(nothing, nothing), xticks=2010:2030,
        width=300, height=300, xminorticksvisible=true, xminorticks=2010:(1 // 12):2030,
        yscale=SymLog(1), yticks=BaseMulTicks([1,2,5]), yminorticksvisible=true, yminorticks=BaseMulTicks(1:9),
        xlabel="Date", ylabel="G-band flux (cnt/s)")
    text!((0,1), space=:relative, align=(:left, :top), offset=(2, -2), color=:black,
        text=j2000 * ", " * rich("Gaia", font=:italic) * " lightcurve")

    fplt = FPlot(r.PHOTg, @o(yeardecimal(_.time)), @o(ustrip(u"s^-1", _.flux));
        color=:black, markersize=5, linewidth=1)
    multiplot!(ax_lc, (scatterlines, rangebars), fplt, label="Measurements")

    # Flux trend line from MC q(t) slope
    t_dense = collect(range(t_start, t_end, length=500))
    S_pm_dense = map(t_dense) do ti
        q = r.mc.q_ref + r.mc.q_slope * (ti - GAIA_DR3_EPOCH)
        1 / q
    end

    lines!(ax_lc, t_dense, S_pm_dense;
        to_xy_attrs(autolimits=false)..., linewidth=2.5, linestyle=:dash, color=colors[3],
        label="Flux trend")

    # -- Col 2: VLBI image + PM arrow + flare --
    sz = 300
    ax_img = Axis(fig[1, 2]; limits=(0 ± lim, 0 ± lim), width=sz, height=sz,
        backgroundcolor=get(colormap, 0))
    plt = image!(ax_img, vlbi_img; colorscale=SymLog(1e-3; vmin=0), colormap)
    translate!(plt, 0, 0, -100)
    hidedecorations!(ax_img)
    scalebar!(1u"mas"; color=:white, position=(0.15, 0.9))
    text!(ax_img, (0, 1), space=:relative, align=(:left, :top),
        offset=(2, -2), color=:white, text=j2000 * ", VLBI image")

    # VLBI position
    scatter!(ax_img, (0, 0); color=vlbi_color, markersize=20,
        marker=marker_lw(:cross, 0.2), label="VLBI position\n(shown at peak)")

    # Jet direction
    arrowlines!(ax_img, [(0, 0), 1.3maxVG .* sincos(r.J.pa)]; color=vlbi_color,
        to_xy_attrs(autolimits=false)..., linestyle=:dash, label="VLBI jet direction")

    # Gaia centroid position
    poly!(ax_img, U.boundary(ustrip(u"mas", r.VG.vec); mul=2);
        color=(:black, 0.07), strokewidth=1, strokecolor=:black,
        label=rich("Gaia", font=:italic) * " position\n(2σ uncertainty)")
    poly!(ax_img, U.boundary(ustrip(u"mas", r.VG.vec); mul=2.1);
        color=(:black, 0.0), strokewidth=1, strokecolor=:white)

    # Apparent proper motion arrow
    nyrs = 1.5
    arrowlines!(ax_img, [vg_mas - nyrs * pm_mas, vg_mas + nyrs * pm_mas]; color=colors[3],
        arrowstyle="->", linewidth=2, markersize=15,
        label=rich("Gaia", font=:italic) * " proper motion\n(over 3 years)")

    # Flare position with uncertainty ellipse
    flare_color = Makie.ColorSchemes.seaborn_bright[4]
    flare_val = U.Value(r.mc.pm.x_flare)
    poly!(ax_img, U.boundary(flare_val; mul=2);
        color=(flare_color, 0.15), strokewidth=2, strokecolor=flare_color,
        to_xy_attrs(autolimits=false)...,
        label="Optical flare location\n(2σ uncertainty)")
    scatter!(ax_img, U.value(flare_val); color=flare_color, markersize=12,
        marker=:star5, to_xy_attrs(autolimits=false)..., label="Optical flare location\n(2σ uncertainty)")

    # -- Legends --
    Legend(fig[2, 1], ax_lc; merge=true, tellheight=true)
    Legend(fig[2, 2], content(fig[1, 2]); merge=true, tellheight=true, orientation=:vertical, nbanks=2)

    rowgap!(fig.layout, 0)
    resize_to_layout!()
    save(joinpath(outdir, "$j2000.pdf"), fig)
    println("  $j2000")
end

println("Done. Saved figures to $outdir")
