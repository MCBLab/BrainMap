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
#   dados_otimizados.rds, ontologyssGSEA.csv, mapeamento_regioes.csv
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
dados_app    <- readRDS(file.path(ROOT, "dados_otimizados.rds"))
matrix_dados <- as.matrix(dados_app$expression_matrix)

# Sample -> atlas region is rebuilt from mapeamento_regioes.csv the way
# preparaDados.R builds it, instead of read off col_meta$structure_mapped.
# The .rds in the repository root predates the last mapping fix on two
# parcels: it still names one "bankssts", which the DK atlas does not carry
# and which therefore renders grey, and it still paints the three cerebellum
# samples onto "corpus callosum". Going back to the CSV keeps the figures on
# the same mapping the deployed API serves.
mapeamento <- read.csv(file.path(ROOT, "mapeamento_regioes.csv"),
                       stringsAsFactors = FALSE)

meta_limpo <- dados_app$col_meta %>%
  distinct(column_num, broad_age, structure_original) %>%
  inner_join(mapeamento %>% select(structure_name, region, macro_region = macro_regions),
             by = c("structure_original" = "structure_name"),
             relationship = "many-to-many") %>%
  filter(!is.na(region), !is.na(broad_age)) %>%
  distinct(column_num, region, broad_age, macro_region)

regioes_por_macro <- meta_limpo %>%
  filter(!is.na(macro_region)) %>%
  distinct(broad_age, region, macro_region)

# One gene set's 524 sample-level scores, pulled out of the 9.7M-row export.
# Figure 4 aggregates these here rather than reading ontologia_micro.rds for
# the reason above: the precomputed tables carry the stale region names, and
# aggregating from samples puts all three figures on one code path.
scores_por_amostra <- function(geneset) {
  csv <- file.path(ROOT, "ontologyssGSEA.csv")
  if (!file.exists(csv))
    stop("ontologyssGSEA.csv not found; run preparaDados.R first (see README.md)")
  linhas <- system2("awk", c("-F,", shQuote(sprintf('$1 == "%s" { print $2 "," $3 }', geneset)),
                             shQuote(csv)), stdout = TRUE)
  if (length(linhas) == 0) stop("gene set not in ontologyssGSEA.csv: ", geneset)
  read.csv(text = paste(c("column_num,Score_ssGSEA", linhas), collapse = "\n"))
}

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

# ---- Figure 4: ontology score ---------------------------------------------
fig_ontology <- function(geneset = "GOBP_FOREBRAIN_GENERATION_OF_NEURONS",
                         escala = "micro", out = "fig4_ontology.pdf") {
  df <- scores_por_amostra(geneset) %>%
    inner_join(meta_limpo, by = "column_num") %>%
    mutate(broad_age = factor(broad_age, levels = AGES)) %>%
    to_scale("Score_ssGSEA", escala)
  stopifnot(nrow(df) > 0)
  sc <- scale_fill_gradientn(colors = PAL, name = "ssGSEA score", na.value = "darkgray")
  ggsave(file.path(OUT, out), build_brain_grid(df, "v", sc), width = 11, height = 9, device = cairo_pdf)
  message("wrote ", file.path(OUT, out))
}

# ---- Figure 5: on-demand ssGSEA for a user gene list ----------------------
# The demonstration list is the dorsal/glutamatergic specification programme
# that Fischer et al. (2021, Alcohol Clin Exp Res 45:979-995) report as
# upregulated by chronic intermittent ethanol in human pluripotent stem cells
# differentiated to first-trimester-equivalent cortical neurons. Only genes
# whose direction that paper states in the text are included. It stands in for
# a real differential-expression result rather than a marker panel, and its
# developmental window is experimentally explicit, which is what makes the map
# checkable: these genes should peak in the first trimester.
FAS_DORSAL <- c("EMX2", "LHX1", "LHX2", "LHX5", "LHX9", "OTX2", "FEZF2",
                "NEUROD1", "NEUROD2", "NEUROD6", "NEUROG2", "TBR1", "EOMES")

fig_genelist <- function(genes = FAS_DORSAL,
                         escala = "micro", out = "fig5_genelist.pdf") {
  valid <- intersect(genes, rownames(matrix_dados))
  stopifnot(length(valid) > 0)
  if (length(valid) < length(genes))
    message("not in the matrix: ", paste(setdiff(genes, valid), collapse = ", "))
  res <- GSVA::gsva(GSVA::ssgseaParam(exprData = matrix_dados,
                                      geneSets = list(UserSet = valid)), verbose = FALSE)
  df <- data.frame(column_num   = as.numeric(colnames(res)),
                   Score_ssGSEA = as.numeric(res["UserSet", ])) %>%
    inner_join(meta_limpo, by = "column_num") %>%
    mutate(broad_age = factor(broad_age, levels = AGES)) %>%
    to_scale("Score_ssGSEA", escala)
  # Per-window means back the face-validity paragraph in Section 4.2.
  print(df %>% group_by(broad_age) %>%
          summarise(mean_score = mean(v), .groups = "drop") %>% as.data.frame())
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
