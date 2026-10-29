#!/usr/bin/env Rscript
# Per-structure enrichment testing for VEGABrain.
#
# The maps paint a mean score per region and window. These functions add the
# missing question: is that score *high enough to mean anything* in this
# structure, in this window. One competitive test, plus a donor guard:
#
#   competitive  in this cell, is the set's within-window deviation higher than
#                that of random gene sets of the same size and expression level?
#   guard        does that survive dropping any single donor?
#
# The unit of inference is the BrainSpan structure, never the ggseg parcel. The
# atlas draws 29 parcels but BrainSpan only sampled 22 structures, so several
# parcels carry byte-identical values (accumbens/caudate/putamen are all
# "striatum"). Testing parcels would count one measurement as up to four
# independent findings. Results are painted outward onto parcels only at the very
# end, for display.
#
# History, so the reasoning is not lost. An earlier design ranked the query
# against the 4,000 nearest-size MSigDB sets and gated it with a within-donor
# label-permutation test. Both halves were wrong for the question asked:
#   - Ranks within a population of curated sets are uniform by construction, so
#     at most ~1% of sets can ever be called per cell, and a set had to be more
#     regionally specific than 99.9% of MSigDB -- itself brain-structured. Only
#     2% of real MSigDB sets produced any enriched cell.
#   - The permutation test compared a structure with its donor's average over the
#     other structures, which is dominated by cortex (11 of 16). Subcortical
#     tissue differs from that average for *every* gene set: on arbitrary sets it
#     passed 32% of subcortical cells and 4% of cortical ones. It filtered cortex
#     out instead of checking donor consistency.
# The random-gene null below is the standard competitive null of gene set
# analysis, is calibrated (see validate_stats.R), calls about a third of real
# MSigDB sets somewhere, and its cortical controls find cortex.
#
# Statistics are pure functions. The only inputs are the expression matrix and
# its sample metadata (dados_otimizados.rds); there is no precomputed reference.

suppressPackageStartupMessages({ library(dplyr) })

# Cells with fewer donors than this are not tested at all. n donors == n samples
# in every cell (BrainSpan has one sample per donor per structure), so this is
# literally "fewer than three people".
MIN_DONORS <- 3

# The 3rd trimester's donors contribute 1, 1, 3, 15 and 15 samples, so its median
# structure rests on two people. It is flagged, not suppressed.
UNDERPOWERED_WINDOWS <- "3rd trimester (n = 5)"

# GSVA's ssGSEA divides the whole score matrix by one global scalar, the range
# over every gene set x every sample in that call. The precomputed ontology
# tables were divided by the range over 18,646 sets, and this is that constant.
# It is display-only: it puts a raw score on the legend scale the rest of the
# paper's maps use. No p-value depends on it.
K_PRE <- 39802.5579


# ---- cells ----------------------------------------------------------------
# One row per sample with the four columns the statistics need, and the list of
# populated (window, structure) cells with donor counts. The structure_mapped
# filter is what restricts the design to structures the atlas can draw: 87 cells,
# 70 with >= MIN_DONORS donors.
cells_from_meta <- function(col_meta) {
  meta <- col_meta %>%
    filter(!is.na(structure_mapped), !is.na(broad_age)) %>%
    distinct(column_num, donor_id, structure_original, broad_age) %>%
    mutate(cell = paste(broad_age, structure_original, sep = " || "))
  cells <- meta %>%
    group_by(cell, broad_age, structure_original) %>%
    summarise(n_donors = n_distinct(donor_id), .groups = "drop") %>%
    as.data.frame()
  list(meta = meta, cells = cells)
}


# ---- cell means -----------------------------------------------------------
# scores: matrix [gene sets x samples], raw (unnormalised) scores.
# meta:   one row per sample, columns column_num, structure_original, broad_age.
# Returns [gene sets x cells] with a `cells` attribute describing the columns.
cell_means <- function(scores, meta) {
  meta <- meta %>%
    filter(!is.na(structure_original), !is.na(broad_age)) %>%
    mutate(cell = paste(broad_age, structure_original, sep = " || "))
  keep <- meta$column_num[meta$column_num %in% as.numeric(colnames(scores))]
  meta <- meta[match(keep, meta$column_num), ]
  idx  <- match(as.character(meta$column_num), colnames(scores))

  # One matrix multiply rather than a split-apply.
  ind <- model.matrix(~ 0 + factor(meta$cell))
  colnames(ind) <- levels(factor(meta$cell))
  A <- scores[, idx, drop = FALSE] %*% sweep(ind, 2, colSums(ind), "/")

  cells <- meta %>%
    distinct(cell, broad_age, structure_original) %>%
    arrange(match(cell, colnames(A)))
  structure(A, cells = as.data.frame(cells))
}


# ---- within-window centring ----------------------------------------------
# Subtracting each window's own mean over structures is what makes the test
# spatial. Centre globally instead and any gene set with a developmental trend
# lights up every structure of its peak window -- a true statement about time
# and a false one about place.
within_window_delta <- function(A, cells) {
  Dw <- A
  for (w in unique(cells$broad_age)) {
    j <- which(cells$broad_age == w)
    Dw[, j] <- A[, j, drop = FALSE] - rowMeans(A[, j, drop = FALSE])
  }
  Dw
}


# ---- ssGSEA from rank tables ----------------------------------------------
# GSVA's ssGSEA statistic (GSVA:::.fastRndWalk) is linear in per-gene terms:
# with per-sample ranks R (ties averaged), weights Ra = R^alpha and C[g] = n -
# position of g in decreasing order + 1,
#   ES = sum_{g in S} Ra*C / sum_{g in S} Ra - (n(n+1)/2 - sum_{g in S} C) / (n - |S|).
# So thousands of gene sets are three sparse-matrix products over tables that
# are computed once. This reproduces GSVA::gsva(ssgseaParam(normalize = FALSE))
# to relative error 1e-6 (validate_stats.R check 0), on the unfiltered matrix.
#
# bins: mean-expression bin of every gene, for expression-matched null sets.
# universe: the genes null sets may be drawn from (default: all rows). Pass the
# union of MSigDB symbols so the null is made of the same kind of gene a curated
# or user list is made of; the matrix's other ~30k rows (non-coding, poorly
# annotated) are regionally noisier at the same mean expression and, used as
# null, make the test conservative (3.2% vs 5% at size 15). Bin edges are
# quantiles of the universe, so a query gene outside it still gets a bin.
rank_tables <- function(matrix_dados, alpha = 0.25, n_bins = 50, universe = NULL) {
  m  <- as.matrix(matrix_dados)
  n  <- nrow(m)
  Rk <- apply(m, 2, rank, ties.method = "average")
  Cm <- apply(Rk, 2, function(r) {
    o <- order(r, decreasing = TRUE)
    cc <- integer(n); cc[o] <- n - seq_len(n) + 1L
    cc
  })
  Ra <- Rk^alpha
  dimnames(Ra) <- dimnames(Cm) <- dimnames(m)
  mu <- rowMeans(m)
  in_universe <- if (is.null(universe)) rep(TRUE, n) else rownames(m) %in% universe
  edges <- quantile(mu[in_universe], probs = seq(0, 1, length.out = n_bins + 1))
  edges[1] <- -Inf; edges[n_bins + 1] <- Inf
  bins <- as.integer(cut(mu, edges, labels = FALSE, include.lowest = TRUE))
  list(Ra = Ra, Cm = Cm, n = n, genes = rownames(m), bins = bins,
       pool = which(in_universe))
}

# W: sparse 0/1 membership matrix [sets x genes] (Matrix::sparseMatrix), columns
# in rt$genes order. Returns raw ES [sets x samples].
ssgsea_from_ranks <- function(rt, W) {
  W@x[] <- 1                                   # a gene drawn twice counts once
  sz <- Matrix::rowSums(W)
  n  <- rt$n
  ES <- as.matrix(W %*% (rt$Ra * rt$Cm)) / as.matrix(W %*% rt$Ra) -
        (n * (n + 1) / 2 - as.matrix(W %*% rt$Cm)) / (n - sz)
  dimnames(ES) <- list(rownames(W), colnames(rt$Ra))
  ES
}

# idx: gene indices, either one set (a vector) or k sets (a [k x |S|] matrix).
membership <- function(idx, n) {
  k <- if (is.matrix(idx)) nrow(idx) else 1L
  Matrix::sparseMatrix(i = rep(seq_len(k), times = length(idx) / k),
                       j = as.integer(idx), x = 1, dims = c(k, n))
}


# ---- random gene sets, expression-matched ---------------------------------
# k null sets: each query gene replaced by a gene drawn from the same
# mean-expression bin. The matching is load-bearing. Roughly half the matrix is
# barely detected; a null drawn uniformly from it would be noise, narrower in
# every cell than any real signature, and would call almost anything enriched.
# pool: indices the null may use (default all). Returns an integer matrix [k x |S|].
draw_matched <- function(bins, idx_query, k = 4000, seed = 1, pool = seq_along(bins)) {
  if (!is.null(seed)) set.seed(seed)
  pools <- split(pool, bins[pool])
  draws <- vapply(idx_query, function(i) {
    pool <- pools[[as.character(bins[i])]]
    pool[sample.int(length(pool), k, replace = TRUE)]
  }, integer(k))
  matrix(draws, nrow = k)
}

null_sets <- function(rt, idx_query, k = 4000, seed = 1) {
  membership(draw_matched(rt$bins, idx_query, k, seed, pool = rt$pool), rt$n)
}


# ---- weighted signatures --------------------------------------------------
# A log2FC-weighted signature, score_j = sum_i w_i z_ij / sum_i |w_i|, is not an
# ssGSEA score, so it gets the same null in its own currency: random
# expression-matched gene sets carrying the query's weights.

# Z-scores of every gene with variance, plus its mean-expression bin. ~200 MB
# for the 47,808 x 524 matrix; build once and keep.
weighted_reference <- function(matrix_dados, n_bins = 50, universe = NULL) {
  m  <- as.matrix(matrix_dados)
  mu <- rowMeans(m)
  sd <- sqrt(rowSums((m - mu)^2) / (ncol(m) - 1))
  ok <- is.finite(sd) & sd > 0
  mu <- mu[ok]
  in_universe <- if (is.null(universe)) rep(TRUE, sum(ok)) else rownames(m)[ok] %in% universe
  edges <- quantile(mu[in_universe], probs = seq(0, 1, length.out = n_bins + 1))
  edges[1] <- -Inf; edges[n_bins + 1] <- Inf
  list(Z    = (m[ok, , drop = FALSE] - mu) / sd[ok],
       bins = as.integer(cut(mu, edges, labels = FALSE, include.lowest = TRUE)),
       pool = which(in_universe))
}

# k null sets, one row each, scored on every sample. Weights keep their position
# in the signature; only the gene identities are drawn.
weighted_null_scores <- function(wref, idx_query, w, k = 4000, seed = 1) {
  draws <- draw_matched(wref$bins, idx_query, k, seed, pool = wref$pool)
  W <- Matrix::sparseMatrix(i = rep(seq_len(k), times = length(idx_query)),
                            j = as.integer(draws), x = rep(w, each = k),
                            dims = c(k, nrow(wref$Z)))
  S <- as.matrix(W %*% wref$Z) / sum(abs(w))
  dimnames(S) <- list(NULL, colnames(wref$Z))
  S
}


# ---- competitive p and BH -------------------------------------------------
# BH over the whole family of tested cells for this gene set (~70), not per
# window. Five per-window families would be more powerful but would run five
# error budgets across one figure.
combine_and_adjust <- function(df, alpha = 0.05) {
  df %>%
    mutate(
      tested = .data$n_donors >= MIN_DONORS,
      q_competitive = ifelse(tested, p.adjust(ifelse(tested, p_competitive, NA), "BH"), NA),
      pass_competitive = !is.na(q_competitive) & q_competitive < alpha,
      underpowered_window = .data$broad_age %in% UNDERPOWERED_WINDOWS)
}


# ---- donor guard: leave one donor out -------------------------------------
# "Is this one donor's tissue?" is answered directly: drop each donor of the cell
# in turn, recompute the cell mean and the window mean for the query and for
# every null set, and re-rank. A cell is robust if it stays below alpha every
# time. This is a guard on cells that already passed BH, so it uses the nominal
# alpha; it is not part of the multiplicity family. Unlike a within-donor label
# permutation, it is tissue-neutral: the null sets share whatever baseline the
# structure has.
#
# S: [1 + k, samples], row 1 the query. Returns one row per requested cell.
loo_guard <- function(S, meta, which_cells, alpha = 0.05) {
  A  <- cell_means(S, meta)
  cl <- attr(A, "cells")
  k  <- nrow(S) - 1
  meta <- meta %>% mutate(cell = paste(broad_age, structure_original, sep = " || "))
  out <- lapply(which_cells, function(cell) {
    jc <- match(cell, cl$cell)
    jw <- which(cl$broad_age == cl$broad_age[jc])
    wm <- rowMeans(A[, jw, drop = FALSE])
    mc <- meta[meta$cell == cell, ]
    ns <- nrow(mc)
    p_d <- vapply(split(as.character(mc$column_num), mc$donor_id), function(cols) {
      Aj <- (A[, jc] * ns - rowSums(S[, cols, drop = FALSE])) / (ns - length(cols))
      d  <- Aj - (wm + (Aj - A[, jc]) / length(jw))
      (1 + sum(d[-1] >= d[1])) / (k + 1)
    }, numeric(1))
    data.frame(cell = cell, p_loo = max(p_d), stringsAsFactors = FALSE)
  })
  bind_rows(out) %>% mutate(robust = p_loo < alpha)
}


# ---- shared tail ----------------------------------------------------------
# S: [1 + k, samples], row 1 the query. One row per cell of ref$cells.
test_cells <- function(S, ref, alpha = 0.05) {
  k  <- nrow(S) - 1
  A  <- cell_means(S, ref$meta)
  cl <- attr(A, "cells")
  Dw <- within_window_delta(A, cl)
  j  <- match(ref$cells$cell, cl$cell)
  Dq <- Dw[1, j]
  Dn <- Dw[-1, j, drop = FALSE]

  st <- ref$cells %>%
    mutate(mean_score          = A[1, j],
           delta_within_window = Dq,
           p_competitive       = (1 + colSums(Dn >= rep(Dq, each = k))) / (k + 1)) %>%
    combine_and_adjust(alpha = alpha) %>%
    mutate(p_loo = NA_real_, robust = NA)

  hit <- which(st$pass_competitive)
  if (length(hit)) {
    g <- loo_guard(S, ref$meta, st$cell[hit], alpha = alpha)
    st$p_loo[hit]  <- g$p_loo[match(st$cell[hit], g$cell)]
    st$robust[hit] <- g$robust[match(st$cell[hit], g$cell)]
  }
  st %>% mutate(enriched = pass_competitive & !is.na(robust) & robust)
}


# ---- top level ------------------------------------------------------------
# Gene list (ssGSEA). genes: symbols; those absent from the matrix are dropped.
# The query's raw per-sample scores are attached as attr "query" so a caller
# that also paints the map does not score the list twice.
enrich <- function(ref, rt, genes, k = 4000, alpha = 0.05, seed = 1) {
  idx <- match(unique(genes), rt$genes)
  idx <- idx[!is.na(idx)]
  if (length(idx) < 5) stop("need at least 5 genes present in the matrix; got ", length(idx))
  W <- rbind(membership(idx, rt$n), null_sets(rt, idx, k, seed))
  S <- ssgsea_from_ranks(rt, W)
  st <- test_cells(S, ref, alpha = alpha)
  attr(st, "query") <- S[1, ]
  st
}

# Weighted signature. Genes without variance are dropped, as score_ponderado()
# in plumber.R drops them.
enrich_weighted <- function(ref, wref, genes, w, k = 4000, alpha = 0.05, seed = 1) {
  idx  <- match(genes, rownames(wref$Z))
  keep <- !is.na(idx)
  idx <- idx[keep]; w <- w[keep]
  if (length(idx) == 0) stop("none of the genes varies across samples")
  if (sum(abs(w)) == 0) stop("all weights are zero")

  x <- as.numeric(crossprod(wref$Z[idx, , drop = FALSE], w)) / sum(abs(w))
  names(x) <- colnames(wref$Z)
  S <- rbind(query = x, weighted_null_scores(wref, idx, w, k = k, seed = seed))
  st <- test_cells(S, ref, alpha = alpha)
  attr(st, "query") <- x
  st
}


# ---- painting outward -----------------------------------------------------
# A parcel is enriched if ANY structure feeding it is enriched. Scores can be
# averaged across contributing structures; p-values cannot, so this is a rule and
# not an aggregation. n_structures records where the ambiguity is: under the
# current mapping Thalamus and Cerebellum each receive two structures.
paint_outward <- function(stats_df, mapeamento) {
  stats_df %>%
    inner_join(mapeamento %>% select(structure_name, region),
               by = c("structure_original" = "structure_name"),
               relationship = "many-to-many") %>%
    group_by(broad_age, region) %>%
    summarise(across(any_of(c("mean_score", "delta_within_window")), mean),
              enriched      = any(enriched),
              # untested cells are all-NA: min() would warn, Inf becomes NA below
              q_competitive = suppressWarnings(min(q_competitive, na.rm = TRUE)),
              p_loo         = suppressWarnings(min(p_loo, na.rm = TRUE)),
              robust        = any(robust, na.rm = TRUE),
              n_structures  = n(),
              tested        = any(tested),
              .groups = "drop") %>%
    mutate(across(c(q_competitive, p_loo), ~ ifelse(is.finite(.x), .x, NA_real_)))
}
