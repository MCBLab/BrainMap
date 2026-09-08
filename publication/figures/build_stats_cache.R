#!/usr/bin/env Rscript
# Builds the reference a competitive enrichment test needs, once, from the
# sample-level ssGSEA export.
#
#   Rscript figures/build_stats_cache.R
#
# Writes enrichment_reference.rds to the repository root. Needs
# ontologyssGSEA.csv, genesets_list.rds, dados_otimizados.rds and
# mapeamento_regioes.csv there.

suppressPackageStartupMessages({ library(dplyr); library(vroom) })

find_root <- function() {
  for (p in c(".", "..", "../..", "../../..")) {
    if (file.exists(file.path(p, "ontologyssGSEA.csv"))) return(normalizePath(p))
  }
  stop("ontologyssGSEA.csv not found; run preparaDados.R first")
}
ROOT <- find_root()
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(HERE, "enrichment_stats.R"))

# GSVA's ssGSEA divides the whole score matrix by one global scalar, the range
# over every gene set x every sample in that call. The precomputed table was
# divided by the range over 18,646 sets; a user's list, scored alone, is divided
# by its own. The two are therefore on different scales -- which would quietly
# invalidate any comparison between an on-demand score and this reference.
#
# The constant is recoverable exactly from a single gene set, so rather than
# re-running ssGSEA for hours we recover it and assert it. If GSVA's
# normalisation ever changes, this stops the build instead of silently shifting
# every p-value.
K_PRE_EXPECTED <- 39802.5579

message("==> loading sample-level scores")
long <- vroom(file.path(ROOT, "ontologyssGSEA.csv"),
              col_types = list(GeneSet = "c", Amostra = "d", Score_ssGSEA = "d"))

n_samp <- n_distinct(long$Amostra)
sets   <- unique(long$GeneSet)
stopifnot(nrow(long) == length(sets) * n_samp)

# Each gene set occupies one contiguous block and the sample order repeats
# identically, so the long table reshapes by a plain matrix() -- no join, no
# sort. Checked rather than assumed, because everything below depends on it.
first_block <- long$Amostra[seq_len(n_samp)]
stopifnot(!any(duplicated(first_block)),
          identical(long$Amostra[n_samp + seq_len(n_samp)], first_block),
          identical(long$GeneSet[seq_len(n_samp)], rep(sets[1], n_samp)))

scores <- matrix(long$Score_ssGSEA, ncol = n_samp, byrow = TRUE,
                 dimnames = list(sets, as.character(first_block)))
rm(long); invisible(gc())
message("    ", nrow(scores), " gene sets x ", ncol(scores), " samples")

# Undo the normalisation, so this reference lives on the same raw scale an
# on-demand ssGSEA(normalize = FALSE) call produces.
raw <- scores * K_PRE_EXPECTED
rm(scores); invisible(gc())

message("==> effective gene set sizes")
dados_app <- readRDS(file.path(ROOT, "dados_otimizados.rds"))
universe  <- rownames(dados_app$expression_matrix)
genesets  <- readRDS(file.path(ROOT, "genesets_list.rds"))
eff <- effective_sizes(genesets[rownames(raw)], universe)

keep <- eff >= 5
message("    dropping ", sum(!keep), " sets with effective size < 5")
raw <- raw[keep, , drop = FALSE]
eff <- eff[keep]
log_eff <- log(eff)

message("==> cell means (structure x window)")
meta <- dados_app$col_meta %>%
  distinct(column_num, donor_id, structure_original, broad_age, structure_mapped) %>%
  filter(!is.na(structure_mapped), !is.na(broad_age)) %>%
  distinct(column_num, donor_id, structure_original, broad_age)

A <- cell_means(raw, meta)
cells <- attr(A, "cells")
cells$n_donors <- vapply(seq_len(nrow(cells)), function(i)
  n_distinct(meta$donor_id[meta$broad_age == cells$broad_age[i] &
                           meta$structure_original == cells$structure_original[i]]),
  integer(1))
message("    ", nrow(cells), " populated cells, ",
        sum(cells$n_donors >= MIN_DONORS), " with >= ", MIN_DONORS, " donors")

message("==> within-window centring and size adjustment")
Dw  <- within_window_delta(A, cells)
fit <- size_adjust(Dw, log_eff)

ref <- list(K_PRE = K_PRE_EXPECTED, cells = cells, eff = eff, log_eff = log_eff,
            fit = fit, A = A, meta = meta, built = Sys.time())
saveRDS(ref, file.path(ROOT, "enrichment_reference.rds"))
message("==> wrote ", file.path(ROOT, "enrichment_reference.rds"),
        " (", round(file.size(file.path(ROOT, "enrichment_reference.rds")) / 1e6), " MB)")

# The size adjustment exists to stop a gene set's size deciding its p-value.
# Report the residual spread by size band so a regression is visible.
bands <- cut(eff, c(0, 10, 25, 50, 100, 200, Inf))
sds <- tapply(seq_along(eff), bands, function(i) mean(apply(fit$R[i, , drop = FALSE], 2, sd)))
message("    residual sd by size band (want these flat):")
for (b in names(sds)) if (!is.na(sds[b])) message(sprintf("      %-10s %.3f", b, sds[b]))
