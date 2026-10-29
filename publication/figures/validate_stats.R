#!/usr/bin/env Rscript
# Calibration, power and controls for the enrichment test.
#
#   Rscript figures/validate_stats.R [n_null] [k]
#
# A test that reports significance is only worth having if its p-values mean what
# they say, so this is not optional tooling -- the numbers it prints are what
# Section 4 of the manuscript and IMPLEMENTATION.md quote. Needs
# dados_otimizados.rds and genesets_list.rds in the repository root; GSVA is
# loaded only for check 0. ~15 min at the defaults.

suppressPackageStartupMessages({ library(dplyr) })

find_root <- function() {
  for (p in c(".", "..", "../..", "../../..")) {
    if (file.exists(file.path(p, "dados_otimizados.rds"))) return(normalizePath(p))
  }
  stop("dados_otimizados.rds not found; run preparaDados.R first")
}
ROOT <- find_root()
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(HERE, "enrichment_stats.R"))

args   <- commandArgs(trailingOnly = TRUE)
N_NULL <- if (length(args) >= 1) as.integer(args[1]) else 300
K      <- if (length(args) >= 2) as.integer(args[2]) else 4000
set.seed(1)

dados_app <- readRDS(file.path(ROOT, "dados_otimizados.rds"))
m   <- as.matrix(dados_app$expression_matrix)
gs  <- readRDS(file.path(ROOT, "genesets_list.rds"))
ref <- cells_from_meta(dados_app$col_meta)
message("cells: ", nrow(ref$cells), ", tested (>= ", MIN_DONORS, " donors): ",
        sum(ref$cells$n_donors >= MIN_DONORS))

# Random lists, and the null sets, are drawn from the union of MSigDB symbols, as
# a user's list would be, not from all 47,808 matrix rows: the other ~30k rows
# are regionally noisier at the same mean expression and, as null, made the test
# conservative (3.2% at nominal 5%).
universe <- intersect(unique(unlist(gs, use.names = FALSE)), rownames(m))
message("MSigDB universe: ", length(universe), " symbols")
t0 <- Sys.time(); rt <- rank_tables(m, universe = universe)
message("rank tables: ", round(as.numeric(Sys.time() - t0, units = "secs")), " s")

# cortical = drawn by the cortical atlas; everything else is subcortical
map <- read.csv(file.path(ROOT, "mapeamento_regioes.csv"), stringsAsFactors = FALSE)
cortical_structures <- unique(map$structure_name[map$region %in% ggseg::dk()$core$region])
classify <- function(s) ifelse(s %in% cortical_structures, "cortical", "subcortical")

rates <- function(p) c(`p<=0.01` = mean(p <= 0.01, na.rm = TRUE),
                       `p<=0.05` = mean(p <= 0.05, na.rm = TRUE),
                       `p<=0.10` = mean(p <= 0.10, na.rm = TRUE))
run <- function(genes, ...) enrich(ref, rt, genes, k = K, ...)

# ------------------------------------------------------------ 0. score identity
# The whole design rests on ssgsea_from_ranks() being GSVA's ssGSEA, not an
# approximation of it. A GSVA upgrade that changes ranking or filtering shows up
# here first.
message("\n==> 0. ssgsea_from_ranks() vs GSVA::gsva(ssgseaParam(normalize = FALSE))")
suppressPackageStartupMessages(library(GSVA))
chk <- list(GOBP_CEREBELLAR_CORTEX_DEVELOPMENT = gs[["GOBP_CEREBELLAR_CORTEX_DEVELOPMENT"]],
            HALLMARK_OXIDATIVE_PHOSPHORYLATION = gs[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]],
            KEGG_DOPAMINERGIC_SYNAPSE          = gs[["KEGG_DOPAMINERGIC_SYNAPSE"]],
            RANDOM_15  = sample(universe, 15),
            RANDOM_50  = sample(universe, 50),
            RANDOM_200 = sample(universe, 200))
chk <- Filter(Negate(is.null), chk)                       # a name absent from this library is skipped
chk <- lapply(chk, function(g) intersect(unique(g), rownames(m)))
gv  <- GSVA::gsva(GSVA::ssgseaParam(m, chk, minSize = 1, normalize = FALSE), verbose = FALSE)
for (nm in names(chk)) {
  lin <- ssgsea_from_ranks(rt, membership(match(chk[[nm]], rt$genes), rt$n))[1, ]
  cat(sprintf("   %-38s n=%4d  max rel. diff %.1e\n", nm, length(chk[[nm]]),
              max(abs(lin - gv[nm, ]) / abs(gv[nm, ]))))
}

# ----------------------------------------------------- 1-2. calibration on null
# The null is "no more regionally structured than random genes of the same size
# and expression level", so a random list IS the null and its p-values must be
# uniform: ~5% of tested cells at p <= 0.05, in every size class. Then the
# question a reader cares about: how often does a null list still produce a
# figure with at least one asterisk, before and after the donor guard.
message("\n==> 1-2. calibration on ", N_NULL, " random lists per size")
for (n in c(15, 50, 200)) {
  cellp <- c(); fam_bh <- c(); fam_guard <- c()
  for (b in seq_len(N_NULL)) {
    st <- run(sample(universe, n), seed = b) %>% filter(tested)
    cellp <- c(cellp, st$p_competitive)
    fam_bh <- c(fam_bh, any(st$pass_competitive)); fam_guard <- c(fam_guard, any(st$enriched))
  }
  cat(sprintf("   size %3d  cell-level: p<=0.01 %.4f  p<=0.05 %.4f  p<=0.10 %.4f   |  family (>=1 hit): BH %.3f  BH+guard %.3f\n",
              n, rates(cellp)[1], rates(cellp)[2], rates(cellp)[3], mean(fam_bh), mean(fam_guard)))
}

# ----------------------------------------------------------------- 3. power
# What fraction of real, curated sets is called somewhere, how concentrated the
# hits are, what the guard removes, and where the hits land. Hippocampus and
# cerebellar cortex are the most distinct tissues from the cortical mean and
# collect more hits than any cortical area -- a property of the tissue, printed
# here so it is never a surprise.
message("\n==> 3. power on 200 real MSigDB sets (effective size 15-300)")
eff  <- vapply(gs, function(g) length(intersect(unique(g), rownames(m))), integer(1))
pick <- sample(names(gs)[eff >= 15 & eff <= 300], 200)
res  <- bind_rows(lapply(pick, function(nm) run(gs[[nm]]) %>% filter(tested) %>% mutate(set = nm)))
hits <- res %>% filter(enriched)
cat(sprintf("   sets with >= 1 enriched cell: %d / 200 (%.0f%%)\n", n_distinct(hits$set), 100 * n_distinct(hits$set) / 200))
cat(sprintf("   cells passing BH: %d; after donor guard: %d (guard removed %.0f%%)\n",
            sum(res$pass_competitive), nrow(hits), 100 * (1 - nrow(hits) / sum(res$pass_competitive))))
cat("   enriched cells per set (among sets with hits):\n")
print(table(cut(table(hits$set), c(0, 1, 2, 4, 8, 16, 70))))
cls <- table(classify(hits$structure_original))
cat(sprintf("   enriched cells by class: cortical %d, subcortical %d  (tested cells: cortical %d, subcortical %d)\n",
            cls["cortical"], cls["subcortical"],
            sum(classify(ref$cells$structure_original[ref$cells$n_donors >= MIN_DONORS]) == "cortical"),
            sum(classify(ref$cells$structure_original[ref$cells$n_donors >= MIN_DONORS]) == "subcortical")))
cat("   hits per structure:\n")
print(hits %>% count(structure_original, sort = TRUE) %>% mutate(class = classify(structure_original)) %>% as.data.frame())

# -------------------------------------------------------------- 4. the guard
message("\n==> 4. donor guard on HALLMARK_OXIDATIVE_PHOSPHORYLATION")
st <- run(gs[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]])
cat(sprintf("   %d cells pass BH, %d survive leaving any single donor out\n", sum(st$pass_competitive), sum(st$enriched)))
print(st %>% filter(pass_competitive) %>%
        transmute(broad_age, structure = substr(structure_original, 1, 32), n_donors,
                  q = signif(q_competitive, 2), p_loo = signif(p_loo, 2), enriched) %>% as.data.frame())

# ---------------------------------------------------------------- 5. controls
message("\n==> 5. controls (anatomically named sets)")
for (nm in c("GOBP_CEREBELLAR_CORTEX_DEVELOPMENT", "GOBP_STRIATUM_DEVELOPMENT",
             "GOBP_THALAMUS_DEVELOPMENT", "GOBP_HIPPOCAMPUS_DEVELOPMENT",
             "GOBP_CEREBRAL_CORTEX_NEURON_DIFFERENTIATION", "GOBP_FOREBRAIN_GENERATION_OF_NEURONS",
             "GOBP_MYELIN_ASSEMBLY")) {
  if (is.null(gs[[nm]])) next
  st <- run(gs[[nm]])
  hit <- st %>% filter(enriched) %>% count(structure_original, sort = TRUE)
  cat(sprintf("   %-44s %2d cells  %s\n", substr(nm, 1, 44), sum(st$enriched),
              if (nrow(hit)) paste0(sub(" \\(.*", "", hit$structure_original), " x", hit$n, collapse = ", ") else "-"))
}

# ---------------------------------------------------------- 6. weighted path
# Same null family in the weighted score's own currency. Random genes with
# random log2FC are the null, so the same uniformity is expected.
message("\n==> 6. weighted signatures: 100 random lists of 50 genes, weights ~ N(0, 1)")
wref <- weighted_reference(m, universe = universe)
cellp <- c(); fam <- c()
for (b in 1:100) {
  g <- sample(universe, 50)
  st <- enrich_weighted(ref, wref, g, rnorm(50), k = K, seed = b) %>% filter(tested)
  cellp <- c(cellp, st$p_competitive); fam <- c(fam, any(st$enriched))
}
cat(sprintf("   cell-level: p<=0.01 %.4f  p<=0.05 %.4f  p<=0.10 %.4f   |  family (>=1 hit) %.3f\n",
            rates(cellp)[1], rates(cellp)[2], rates(cellp)[3], mean(fam)))
message("\ndone.")
