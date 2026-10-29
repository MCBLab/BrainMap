#!/usr/bin/env Rscript
# Figure 6: the enrichment map -- colour is magnitude, an asterisk is significance.
#
#   Rscript figures/fig6_enrichment.R [GENESET_NAME]
#
# Needs dados_otimizados.rds and genesets_list.rds (preparaDados.R) in the
# repository root. An asterisk marks a structure whose within-window deviation
# beats 10,000 expression-matched random gene sets at q < 0.05 and survives
# leaving any single donor out (enrichment_stats.R).
#
# Fill is the within-window deviation, not the raw score, because that is the
# quantity the test is actually about. Painting the raw score instead produces a
# figure that argues with itself: the cerebellum carries its highest absolute
# ssGSEA in the two prenatal windows, but it is only *distinctive* among
# structures in infancy and adulthood, so a raw-score fill shows a blue parcel
# starred as significant. Colour = how far above this window's other structures,
# asterisk = whether that survives both tests.
#
# Figures 3-5 already show absolute magnitude; this one is not a second copy.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggseg); library(cowplot)
})

find_root <- function() {
  for (p in c(".", "..", "../..", "../../..")) {
    if (file.exists(file.path(p, "dados_otimizados.rds"))) return(normalizePath(p))
  }
  stop("dados_otimizados.rds not found; run preparaDados.R first")
}
ROOT <- find_root()
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(HERE, "enrichment_stats.R"))

AGES <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)",
          "Infant (n = 14)", "Adult (n = 8)")
# Diverging, centred on zero: the fill is a signed deviation, so a sequential
# ramp would imply an origin the quantity does not have.
PAL  <- c("#1A318B", "#4F71BE", "#F2F2F2", "#D1498C", "#7A0845")
STAR <- list(label = "*", size = 9)   # marker for an enriched parcel

# Where to put the asterisk: a point guaranteed to lie INSIDE the parcel's
# largest polygon (sf::st_point_on_surface; the area-weighted centroid falls
# outside thin, curved parcels such as superior parietal in the medial view),
# computed from the geometry ggseg actually draws for that view -- one branch per
# ggseg version. Hemispheres drawn on top of each other (aseg's sagittal
# cerebellum) collapse to one point.
parcel_centroids <- function(atlas, pos) {
  medir <- function(pg, region, hemi, view) {
    bb <- sf::st_bbox(pg)
    pt <- tryCatch(sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(pg))),
                   error = function(e) sf::st_coordinates(sf::st_centroid(pg)))
    tibble(region = region, hemi = hemi, view = view, area = as.numeric(sf::st_area(pg)),
           x = pt[1, 1], y = pt[1, 2],
           xmin = bb[["xmin"]], xmax = bb[["xmax"]], ymin = bb[["ymin"]], ymax = bb[["ymax"]])
  }
  if (exists("prepare_polygon_atlas", envir = asNamespace("ggseg"), inherits = FALSE)) {
    # ggseg >= 2.2 draws flat polygons (x, y)
    plano <- ggseg:::prepare_polygon_atlas(atlas, position = position_brain(pos)) %>%
      filter(!is.na(region)) %>%
      mutate(chave = paste(region, hemi, .feature_id, subgroup, sep = "|"))
    aneis <- bind_rows(lapply(split(plano, plano$chave), function(d) {
      xy <- cbind(d$x, d$y)
      if (any(xy[1, ] != xy[nrow(xy), ])) xy <- rbind(xy, xy[1, ])
      medir(sf::st_sfc(sf::st_polygon(list(xy))), d$region[1], d$hemi[1], d$view[1])
    }))
  } else {
    # ggseg 2.1.x draws repositioned sf geometry
    rep <- ggseg::reposition_brain(as.data.frame(atlas), position = pos)
    rep <- rep[!is.na(rep$region), ]
    geo <- sf::st_geometry(rep)
    aneis <- bind_rows(lapply(seq_len(nrow(rep)), function(i) {
      polys <- suppressWarnings(sf::st_cast(geo[i], "POLYGON"))
      bind_rows(lapply(seq_along(polys), function(k) medir(polys[k], rep$region[i], rep$hemi[i], rep$view[i])))
    }))
  }
  # a parcel that is only a sliver in this view (precentral seen from the medial
  # side: 0.05% of the panel) is starred where it is actually drawn, not here
  view_area <- (max(aneis$xmax) - min(aneis$xmin)) * (max(aneis$ymax) - min(aneis$ymin))
  aneis %>%
    group_by(region, hemi, view) %>%
    slice_max(area, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    filter(area >= 0.0025 * view_area) %>%
    group_by(region, view) %>%
    # keep both hemispheres only when their bounding boxes are disjoint
    filter(n() == 1 | max(xmin) > min(xmax) | max(ymin) > min(ymax) | area == max(area)) %>%
    ungroup() %>%
    select(region, hemi, view, x, y)
}


# Same four views and facetting as plumber.R:build_brain_grid. Significance is a
# separate text layer -- an asterisk on every enriched parcel -- so the brain
# layer itself is unchanged from the rest of the paper's maps.
#
# The asterisk layer carries its own data (one row per enriched parcel per
# window), so facet_wrap places it in the right window, and it sets no fill, so
# ggseg's discrete fill scale is never triggered.
build_brain_grid <- function(data, fill_var, fill_scale) {
  panel <- function(atlas, pos, strip, legend, star_size = STAR$size) {
    stars <- data %>%
      filter(enriched) %>%
      distinct(broad_age, region) %>%
      inner_join(parcel_centroids(atlas, pos), by = "region",
                 relationship = "many-to-many")
    # ggseg >= 2.2 no longer reads the plot's data: it must be handed to
    # geom_brain, grouped by the facet variable so every window gets a full atlas.
    ggplot() +
      geom_brain(data = group_by(data, broad_age),
                 atlas = atlas, position = position_brain(pos),
                 mapping = aes(fill = .data[[fill_var]]),
                 colour = "grey25", linewidth = 0.3) +
      # thin dark halo under a white glyph, so the mark reads on any fill
      geom_text(data = stars, mapping = aes(x = x, y = y), inherit.aes = FALSE,
                label = STAR$label, size = star_size * 1.2, colour = "grey20",
                fontface = "bold", vjust = 0.3) +
      geom_text(data = stars, mapping = aes(x = x, y = y), inherit.aes = FALSE,
                label = STAR$label, size = star_size, colour = "white",
                fontface = "bold", vjust = 0.3) +
      fill_scale +
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
            # small, adjacent parcels: smaller glyph so asterisks do not overlap
            panel(ggseg::aseg(), "coronal_1",     FALSE, "bottom", star_size = 6),
            nrow = 4, align = "v", rel_heights = c(1, 0.89, 1.10, 1.4))
}

geneset <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(geneset)) geneset <- "GOBP_CEREBELLAR_CORTEX_DEVELOPMENT"

message("loading data from ", ROOT)
dados_app <- readRDS(file.path(ROOT, "dados_otimizados.rds"))
genesets  <- readRDS(file.path(ROOT, "genesets_list.rds"))
map <- read.csv(file.path(ROOT, "mapeamento_regioes.csv"), stringsAsFactors = FALSE)
ref <- cells_from_meta(dados_app$col_meta)
# null genes come from the MSigDB universe, the kind of gene a curated list has
universe <- intersect(unique(unlist(genesets, use.names = FALSE)), rownames(dados_app$expression_matrix))
rt  <- rank_tables(dados_app$expression_matrix, universe = universe)
if (is.null(genesets[[geneset]])) stop("gene set not in genesets_list.rds: ", geneset)

message("testing ", geneset)
# k = 10,000 for the figure: a finer p floor than the API's 4,000, at ~2 s.
st <- enrich(ref, rt, genesets[[geneset]], k = 10000)

n_sig <- sum(st$enriched)
message("  ", n_sig, " of ", sum(st$tested), " tested cells enriched at q < 0.05 ",
        "(", sum(st$pass_competitive & !st$enriched), " passed BH but not the donor guard)")
if (n_sig) print(st %>% filter(enriched) %>%
                 transmute(broad_age, structure_original, n_donors,
                           q_competitive = signif(q_competitive, 3),
                           p_loo = signif(p_loo, 3)) %>% as.data.frame())

# Structures are mapped onto atlas parcels for display; a parcel is starred when
# any structure feeding it is enriched (see enrichment_stats.R:paint_outward).
df <- paint_outward(st, map) %>%
  mutate(broad_age = factor(broad_age, levels = AGES),
         # scores live on the raw scale; put the deviation back on the one the
         # rest of the paper's legends use (K_PRE is display-only)
         delta_within_window = delta_within_window / K_PRE)

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
