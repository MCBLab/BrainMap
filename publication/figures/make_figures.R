#!/usr/bin/env Rscript
# Publication figures for the VEGABrain manuscript.
#
# Renders locally with the same ggseg code path the API serves
# (see ../../plumber.R:build_brain_grid), so what the paper shows is what
# the tool produces. Writes PDF directly -- no SVG conversion step, and no
# dependency on the deployed API being warm.
#
#   Rscript figures/make_figures.R            # all figures
#   Rscript figures/make_figures.R fig4        # one figure
#
# Needs the prepared data objects in the repository root:
#   dados_otimizados.rds, ontologia_micro.rds, ontologia_macro.rds
# built by ../../preparaDados.R.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggseg); library(cowplot); library(GSVA)
})

# Repository root: walk up from the working directory until the prepared data
# objects appear, so the script runs from figures/, from publication/ or from
# the repository root.
find_root <- function() {
  for (p in c(".", "..", "../..", "../../..")) {
    if (file.exists(file.path(p, "dados_otimizados.rds"))) return(normalizePath(p))
  }
  stop("dados_otimizados.rds not found; run preparaDados.R first (see README.md)")
}
ROOT <- find_root()
OUT  <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                                 value = TRUE)[1])), "."),
                      mustWork = FALSE)
if (is.na(OUT) || !dir.exists(OUT)) OUT <- "."

AGES <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)",
          "Infant (n = 8)", "Adult (n = 14)")
PAL  <- c("#1A318B", "#4F71BE", "#C2B4D6", "#D1498C", "#7A0845")

message("loading data from ", ROOT)
dados_app       <- readRDS(file.path(ROOT, "dados_otimizados.rds"))
ontologia_micro <- readRDS(file.path(ROOT, "ontologia_micro.rds"))
ontologia_macro <- readRDS(file.path(ROOT, "ontologia_macro.rds"))
matrix_dados    <- as.matrix(dados_app$expression_matrix)

meta_limpo <- dados_app$col_meta %>%
  rename(region = structure_mapped) %>%
  filter(!is.na(region), !is.na(broad_age)) %>%
  select(column_num, region, broad_age, macro_region)

regioes_por_macro <- meta_limpo %>%
  filter(!is.na(macro_region)) %>%
  distinct(broad_age, region, macro_region)

# Mirrors build_brain_grid() in plumber.R: four stacked views, faceted by age.
build_brain_grid <- function(data, fill_var, fill_scale, caption_text = NULL) {
  panel <- function(atlas, pos, strip, legend, caption) {
    ggplot(data) +
      geom_brain(atlas = atlas, position = position_brain(pos),
                 mapping = aes(fill = .data[[fill_var]]), color = "black", size = 0.5) +
      fill_scale +
      facet_wrap(~broad_age, ncol = 5) +
      labs(caption = caption) +
      theme_void() +
      theme(strip.text      = if (strip) element_text(face = "bold") else element_blank(),
            strip.background = element_blank(),
            legend.position = legend,
            plot.margin     = margin(t = 20, r = 5, b = 0, l = 5),
            legend.key.height = unit(2, "cm"),
            plot.background = element_rect(fill = "white", colour = NA),
            plot.caption    = element_text(color = "firebrick", face = "bold",
                                           size = 11, hjust = 0.5, margin = margin(t = 15)))
  }
  plot_grid(panel(ggseg::dk(),   "right lateral", TRUE,  "none",   NULL),
            panel(ggseg::dk(),   "right medial",  FALSE, "none",   NULL),
            panel(ggseg::aseg(), "sagittal",      FALSE, "none",   NULL),
            panel(ggseg::aseg(), "coronal_1",     FALSE, "bottom", caption_text),
            nrow = 4, align = "v", rel_heights = c(1, 0.89, 1.10, 1.4))
}

to_scale <- function(df, value_col, escala) {
  if (escala == "macro") {
    df %>% filter(!is.na(macro_region)) %>%
      group_by(broad_age, macro_region) %>%
      mutate(v = mean(.data[[value_col]], na.rm = TRUE)) %>% ungroup() %>%
      distinct(broad_age, region, v)
  } else {
    df %>% group_by(broad_age, region) %>%
      summarise(v = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")
  }
}

# ---- Figure 3: single gene ------------------------------------------------
fig_gene <- function(gene = "SOX10", escala = "micro", out = "fig3_gene_sox10.pdf") {
  df <- data.frame(column_num = as.numeric(colnames(matrix_dados)),
                   Expressao  = as.numeric(matrix_dados[gene, ])) %>%
    inner_join(meta_limpo, by = "column_num") %>%
    mutate(Expressao = pmin(Expressao, 6),
           broad_age = factor(broad_age, levels = AGES)) %>%
    to_scale("Expressao", escala)
  sc <- scale_fill_gradientn(colors = PAL, values = scales::rescale(c(0, 1.5, 3, 4.5, 6)),
                             limits = c(0, 6), name = "Log2 Expr", na.value = "darkgray")
  ggsave(file.path(OUT, out), build_brain_grid(df, "v", sc), width = 11, height = 9, device = cairo_pdf)
  message("wrote ", file.path(OUT, out))
}

# ---- Figure 4: precomputed ontology score ---------------------------------
fig_ontology <- function(geneset = "GOBP_FOREBRAIN_GENERATION_OF_NEURONS",
                         escala = "micro", out = "fig4_ontology.pdf") {
  df <- if (escala == "macro") {
    ontologia_macro %>% filter(GeneSet == geneset) %>%
      inner_join(regioes_por_macro, by = c("broad_age", "macro_region"),
                 relationship = "many-to-many") %>%
      select(broad_age, region, v = Score_Medio_ssGSEA)
  } else {
    ontologia_micro %>% filter(GeneSet == geneset) %>%
      select(broad_age, region, v = Score_Medio_ssGSEA)
  }
  stopifnot(nrow(df) > 0)
  df <- df %>% mutate(broad_age = factor(broad_age, levels = AGES))
  sc <- scale_fill_gradientn(colors = PAL, name = "ssGSEA score", na.value = "darkgray")
  ggsave(file.path(OUT, out), build_brain_grid(df, "v", sc), width = 11, height = 9, device = cairo_pdf)
  message("wrote ", file.path(OUT, out))
}

# ---- Figure 5: on-demand ssGSEA for a user gene list ----------------------
fig_genelist <- function(genes = c("GAD1", "GAD2", "SLC32A1", "GABRB3"),
                         escala = "micro", out = "fig5_genelist.pdf") {
  valid <- intersect(genes, rownames(matrix_dados))
  stopifnot(length(valid) > 0)
  res <- GSVA::gsva(GSVA::ssgseaParam(exprData = matrix_dados,
                                      geneSets = list(UserSet = valid)), verbose = FALSE)
  df <- data.frame(column_num   = as.numeric(colnames(res)),
                   Score_ssGSEA = as.numeric(res["UserSet", ])) %>%
    inner_join(meta_limpo, by = "column_num") %>%
    mutate(broad_age = factor(broad_age, levels = AGES)) %>%
    to_scale("Score_ssGSEA", escala)
  sc <- scale_fill_gradientn(colors = PAL, name = "Custom ssGSEA", na.value = "darkgray")
  ggsave(file.path(OUT, out), build_brain_grid(df, "v", sc), width = 11, height = 9, device = cairo_pdf)
  message("wrote ", file.path(OUT, out))
}

args <- commandArgs(trailingOnly = TRUE)
want <- function(x) length(args) == 0 || x %in% args
if (want("fig3")) fig_gene()
if (want("fig4")) fig_ontology()
if (want("fig5")) fig_genelist()
message("done. Figure 2 (interface) is a screen capture, not generated here.")
