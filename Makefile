JULIA ?= julia --project=scripts

SCRIPTS = scripts/prepare_data.jl
VLBI_IMAGES = data/vlbi_images

.PHONY: all clean

all: scripts/.instantiated \
     figs/offset_hist.pdf figs/corejet_cartoon.pdf \
     figs/image_motion/J1058+0133.pdf figs/image_motion/J1824+5651.pdf \
     figs/pm_hist_color.pdf \
     figs/flare_location_h_top.pdf figs/flare_location_h_bot.pdf \
     figs/flare_2d_clean \
     figs/gaia_lc/J1058+0133.pdf figs/gaia_lc/J1824+5651.pdf \
     figs/interactive/corejet_cartoon figs/interactive/offset_pm_hist

scripts/.instantiated: scripts/Project.toml scripts/Manifest.toml
	$(JULIA) -e 'using Pkg; Pkg.instantiate()'
	touch $@

$(VLBI_IMAGES): scripts/download_vlbi_images.jl
	$(JULIA) $<

# -- Figures --

figs/offset_hist.pdf: scripts/plot_offset_hist.jl $(SCRIPTS) scripts/theme.jl scripts/plot_vlbi_gaia_sketch.jl scripts/cartoon_common.jl
	$(JULIA) $<

figs/corejet_cartoon.pdf: scripts/plot_corejet_cartoon.jl scripts/cartoon_common.jl
	$(JULIA) $<

figs/image_motion/J1058+0133.pdf figs/image_motion/J1824+5651.pdf &: scripts/plot_image_motion.jl $(SCRIPTS) scripts/theme.jl $(VLBI_IMAGES)
	$(JULIA) $<

figs/pm_hist_color.pdf: scripts/plot_pm_hist.jl $(SCRIPTS) scripts/theme.jl
	$(JULIA) $<

figs/flare_location_h_top.pdf figs/flare_location_h_bot.pdf &: scripts/plot_flare_location_h.jl $(SCRIPTS) scripts/theme.jl scripts/plot_vlbi_gaia_sketch.jl scripts/cartoon_common.jl
	$(JULIA) $<

figs/flare_2d_clean: scripts/plot_flare_2d_clean.jl $(SCRIPTS) scripts/theme.jl $(VLBI_IMAGES)
	$(JULIA) $<

figs/gaia_lc/J1058+0133.pdf figs/gaia_lc/J1824+5651.pdf &: scripts/plot_gaia_lc.jl $(SCRIPTS) scripts/theme.jl
	$(JULIA) $<

# -- Interactive figures --

figs/interactive/corejet_cartoon: scripts/plot_corejet_cartoon_interactive.jl scripts/cartoon_common.jl
	$(JULIA) $<

figs/interactive/offset_pm_hist: scripts/plot_interactive_offset_pm_hist.jl $(SCRIPTS) scripts/theme.jl
	$(JULIA) $<

clean:
	rm -rf data/vlbi_images figs/

crop:
	find figs -name '*.pdf' -exec sh -c 'uvx pdfcropmargins -p 0 -o "$$1.tmp" "$$1" && mv "$$1.tmp" "$$1"' _ {} \;
