#!/usr/bin/env Rscript
# Calibration and controls for the enrichment test.
#
#   Rscript figures/validate_stats.R [n_null] [B]
#
# A test that reports significance is only worth having if its p-values mean what
# they say, so this is not optional tooling -- the numbers it prints are what
# Section 4 of the manuscript quotes.

suppressPackageStartupMessages({ library(dplyr); library(GSVA) })

find_root <- function() {
  for (p in c(".", "..", "../..", "../../..")) {
    if (file.exists(file.path(p, "enrichment_reference.rds"))) return(normalizePath(p))
  }
  stop("enrichment_reference.rds not found; run build_stats_cache.R first")
}
ROOT <- find_root()
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(HERE, "enrichment_stats.R"))

args   <- commandArgs(trailingOnly = TRUE)
N_NULL <- if (length(args) >= 1) as.integer(args[1]) else 300
B      <- if (length(args) >= 2) as.integer(args[2]) else 2000
set.seed(1)

ref <- readRDS(file.path(ROOT, "enrichment_reference.rds"))
dados_app <- readRDS(file.path(ROOT, "dados_otimizados.rds"))
m  <- as.matrix(dados_app$expression_matrix)
gs <- readRDS(file.path(ROOT, "genesets_list.rds"))

# Random sets are drawn from the union of MSigDB symbols, not from all 47,808
# matrix rows. A set of arbitrary rows scores far lower than any curated set
# (mean raw ES 1,635 vs 11,981 at size 13) because most rows are barely
# expressed, and calibrating against those would flatter the test.
universe <- intersect(unique(unlist(gs, use.names = FALSE)), rownames(m))
message("MSigDB universe: ", length(universe), " symbols")

rates <- function(p) c(`p<=0.01` = mean(p <= 0.01, na.rm = TRUE),
                       `p<=0.05` = mean(p <= 0.05, na.rm = TRUE),
                       `p<=0.10` = mean(p <= 0.10, na.rm = TRUE))

# ------------------------------------------------------- 1a. competitive, own null
# The competitive p is a rank within the reference population, so its null is
# "this set is an ordinary member of that population". A real set, compared to
# its size-matched peers with itself excluded, is exactly that null -- and the
# rate should be nominal and, crucially, FLAT across size bands. A trend here
# means the size adjustment is not doing its job.
#
# Do not be tempted to build this null by scrambling structure labels instead.
# Scrambling removes a set's regional signal but also strips the systematic
# offset every real set carries in a given cell, so a scrambled query is no
# longer exchangeable with the reference and the test reads ~1.5x
# anti-conservative -- an artifact of the null, not of the statistic.
message("\n==> 1a. competitive calibration (real sets vs their own population)")
idx <- sample(nrow(ref$fit$R), min(500, nrow(ref$fit$R)))
tested0 <- ref$cells$n_donors >= MIN_DONORS
pc0 <- t(vapply(idx, function(i)
  p_competitive(ref$fit$R, ref$log_eff, ref$fit$R[i, ], ref$eff[i], self = i),
  numeric(ncol(ref$fit$R))))
cat("   observed vs nominal:\n"); print(round(rates(pc0[, tested0]), 4))
cat("   by size band (want flat ~0.05):\n")
print(round(tapply(seq_along(idx), cut(ref$eff[idx], c(0, 10, 25, 50, 100, 200, Inf)),
                   function(i) mean(pc0[i, tested0] <= 0.05, na.rm = TRUE)), 4))

# ------------------------------------------------- 1b. competitive, user gene list
# A user's list is not drawn from MSigDB, so its exchangeability with the
# reference is an assumption rather than a fact. Random lists probe the far end
# of that: they score well below curated sets and come out conservative, which
# is the safe direction but is worth quoting as a limitation.
message("\n==> 1b. competitive calibration (", N_NULL, " random gene lists)")
sizes <- sample(ref$eff, N_NULL, replace = TRUE)
rand  <- lapply(sizes, function(n) sample(universe, n))
names(rand) <- sprintf("RANDOM_%04d", seq_len(N_NULL))

res <- GSVA::gsva(GSVA::ssgseaParam(exprData = m, geneSets = rand,
                                    minSize = 1, normalize = FALSE), verbose = FALSE)
Aq <- cell_means(res, ref$meta)
Aq <- Aq[, match(ref$cells$cell, attr(Aq, "cells")$cell), drop = FALSE]
Dq <- within_window_delta(Aq, ref$cells)

pc <- t(vapply(seq_len(nrow(Dq)), function(i) {
  rq <- size_adjust_query(Dq[i, ], log(sizes[i]), ref$fit)
  p_competitive(ref$fit$R, ref$log_eff, rq, sizes[i])
}, numeric(ncol(Dq))))

tested <- ref$cells$n_donors >= MIN_DONORS
cat("   observed vs nominal (expect conservative):\n")
print(round(rates(pc[, tested]), 4))
cat("   by query size band:\n")
print(round(tapply(seq_along(sizes), cut(sizes, c(0, 10, 25, 50, 100, 200, Inf)),
                   function(i) mean(pc[i, tested] <= 0.05, na.rm = TRUE)), 4))

# ---------------------------------------------------------------- 2. spatial
# Null here is built by scrambling which structure each sample came from, within
# donor, and then running *real* gene sets. That destroys spatial signal while
# keeping every other property of real data, which a random gene set would not.
message("\n==> 2. spatial calibration (", N_NULL, " real gene sets, labels scrambled)")
scrambled <- ref$meta %>% group_by(broad_age, donor_id) %>%
  mutate(structure_original = structure_original[sample(n())]) %>% ungroup()

pick <- sample(rownames(ref$fit$R), N_NULL)
ps <- lapply(pick, function(nm) {
  x <- scores_ontology(nm, ref, ROOT)
  p_spatial(x, scrambled, B = B)$p_spatial
})
cat("   observed vs nominal:\n")
print(round(rates(unlist(ps)), 4))

# ---------------------------------------------------------------- 3. family
# What a reader actually cares about: how often does a figure with no real signal
# still show at least one significant structure.
message("\n==> 3. family-level false positive rate (>=1 enriched cell)")
fam <- vapply(seq_len(min(N_NULL, 150)), function(i) {
  x <- stats::setNames(as.numeric(res[i, ]), colnames(res))
  st <- enrich(ref, x, sizes[i], B = B)
  any(st$enriched)
}, logical(1))
cat(sprintf("   %.3f of null gene sets produce at least one enriched cell (target <= 0.05)\n",
            mean(fam)))

# ---------------------------------------------------------------- 4. controls
message("\n==> 4. controls")
for (nm in c("GOBP_CEREBELLAR_CORTEX_DEVELOPMENT", "GOBP_FOREBRAIN_GENERATION_OF_NEURONS",
             "GOBP_MYELIN_ASSEMBLY")) {
  st <- enrich(ref, scores_ontology(nm, ref, ROOT), ref$eff[[nm]], B = 10000,
               self = match(nm, rownames(ref$fit$R)))
  hit <- st %>% filter(enriched) %>% pull(structure_original) %>% unique()
  cat(sprintf("   %-40s %d cells  %s\n", substr(nm, 1, 40), sum(st$enriched),
              if (length(hit)) paste(hit, collapse = ", ") else "-"))
}
message("\ndone.")
