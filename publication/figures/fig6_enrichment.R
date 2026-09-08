#!/usr/bin/env Rscript
# Figure 6: the enrichment map -- colour is magnitude, outline is significance.
#
#   Rscript figures/fig6_enrichment.R [GENESET_NAME]
#
# Needs enrichment_reference.rds (build_stats_cache.R) in the repository root.
#
# Fill is the within-window deviation, not the raw score, because that is the
# quantity the test is actually about. Painting the raw score instead produces a
# figure that argues with itself: the cerebellum carries its highest absolute
# ssGSEA in the two prenatal windows, but it is only *distinctive* among
# structures in infancy and adulthood, so a raw-score fill shows a blue parcel
# ringed as significant. Colour = how far above this window's other structures,
# outline = whether that survives both tests.
#
# Figures 3-5 already show absolute magnitude; this one is not a second copy.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggseg); library(cowplot); library(GSVA)
})

find_root <- function() {
  for (p in c(".", "..", "../..", "../../..")) {
    if (file.exists(file.path(p, "enrichment_reference.rds"))) return(normalizePath(p))
  }
  stop("enrichment_reference.rds not found; run build_stats_cache.R first")
}
ROOT <- find_root()
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(HERE, "enrichment_stats.R"))

AGES <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)",
          "Infant (n = 14)", "Adult (n = 8)")
# Diverging, centred on zero: the fill is a signed deviation, so a sequential
# ramp would imply an origin the quantity does not have.
PAL  <- c("#1A318B", "#4F71BE", "#F2F2F2", "#D1498C", "#7A0845")
SIG  <- "#00E676"   # outline colour for an enriched structure

# Same four views and facetting as plumber.R:build_brain_grid, with significance
# carried by the outline of the very same layer that carries magnitude.
#
# A separate outline layer was the obvious approach and does not work: ggseg
# appends its own scale_fill_manual whenever `fill` is absent from a layer's
# mapping, so a fill = NA overlay drags a discrete fill scale into a plot that
# already has a continuous one. Mapping colour and linewidth inside the existing
# geom_brain sidesteps that, and has the side benefit of needing no change to
# build_brain_grid's signature in plumber.R -- only its aes().
build_brain_grid <- function(data, fill_var, fill_scale) {
  panel <- function(atlas, pos, strip, legend) {
    ggplot(data) +
      geom_brain(atlas = atlas, position = position_brain(pos),
                 mapping = aes(fill = .data[[fill_var]],
                               colour = enriched, linewidth = enriched)) +
      fill_scale +
      scale_colour_manual(values = c(`TRUE` = SIG, `FALSE` = "grey25"),
                          na.value = "grey25", guide = "none") +
      scale_linewidth_manual(values = c(`TRUE` = 1.2, `FALSE` = 0.3),
                             na.value = 0.3, guide = "none") +
      facet_wrap(~broad_age, ncol = 5) +
      theme_void() +
      theme(strip.text        = if (strip) element_text(face = "bold") else element_blank(),
            strip.background  = element_blank(),
            legend.position   = legend,
            plot.margin       = margin(t = 20, r = 5, b = 0, l = 5),
            legend.text       = element_text(size = 8),
            legend.title      = element_text(size = 8),
            plot.background   = element_rect(fill = "white", colour = NA))
  }
  plot_grid(panel(ggseg::dk(),   "right lateral", TRUE,  "none"),
            panel(ggseg::dk(),   "right medial",  FALSE, "none"),
            panel(ggseg::aseg(), "sagittal",      FALSE, "none"),
            panel(ggseg::aseg(), "coronal_1",     FALSE, "bottom"),
            nrow = 4, align = "v", rel_heights = c(1, 0.89, 1.10, 1.4))
}

geneset <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(geneset)) geneset <- "GOBP_CEREBELLAR_CORTEX_DEVELOPMENT"

message("loading reference from ", ROOT)
ref <- readRDS(file.path(ROOT, "enrichment_reference.rds"))
map <- read.csv(file.path(ROOT, "mapeamento_regioes.csv"), stringsAsFactors = FALSE)

message("testing ", geneset)
self <- match(geneset, rownames(ref$fit$R))
st <- enrich(ref, scores_ontology(geneset, ref, ROOT), ref$eff[[geneset]],
             B = 10000, self = self)

n_sig <- sum(st$enriched)
message("  ", n_sig, " of ", sum(st$tested), " tested cells enriched at q < 0.05")
if (n_sig) print(st %>% filter(enriched) %>%
                 transmute(broad_age, structure_original, n_donors,
                           q_competitive = signif(q_competitive, 3),
                           q_spatial = signif(q_spatial, 3)) %>% as.data.frame())

# Painted outward for display only -- a parcel inherits its structure's verdict.
df <- paint_outward(st, map) %>%
  mutate(broad_age = factor(broad_age, levels = AGES),
         # scores live on the raw scale; put the deviation back on the one the
         # rest of the paper's legends use
         delta_within_window = delta_within_window / ref$K_PRE)

lim <- max(abs(df$delta_within_window), na.rm = TRUE)
sc <- scale_fill_gradientn(colors = PAL, limits = c(-lim, lim), n.breaks = 5,
                           name = "ssGSEA score,\ndeviation from\nwindow mean",
                           na.value = "grey85",
                           guide = guide_colourbar(barwidth = 12, barheight = 0.8,
                                                   title.position = "left",
                                                   title.vjust = 0.9))

out <- file.path(HERE, "fig6_enrichment.pdf")
ggsave(out, build_brain_grid(df, "delta_within_window", sc),
       width = 11, height = 9, device = cairo_pdf)
message("wrote ", out)
