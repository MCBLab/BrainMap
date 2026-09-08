#!/usr/bin/env Rscript
# Per-structure enrichment testing for VEGABrain.
#
# The maps paint a mean ssGSEA score per region and window. These functions add
# the missing question: is that score *high enough to mean anything* in this
# structure, in this window. Two tests, and a cell has to pass both.
#
#   competitive  is this gene set unusual among gene sets of the same size, here?
#   spatial      is that consistent across donors, or is it one donor's tissue?
#
# The unit of inference is the BrainSpan structure, never the ggseg parcel. The
# atlas draws 29 parcels but BrainSpan only sampled 22 structures, so several
# parcels carry byte-identical values (accumbens/caudate/putamen are all
# "striatum"). Testing parcels would count one measurement as up to four
# independent findings. Results are painted outward onto parcels only at the very
# end, for display.
#
# Statistics are pure functions; the two score fetchers at the bottom are the
# only I/O. build_stats_cache.R builds the reference these read.

suppressPackageStartupMessages({ library(dplyr) })

# Cells with fewer donors than this are not tested at all. n donors == n samples
# in every cell (BrainSpan has one sample per donor per structure), so this is
# literally "fewer than three people".
MIN_DONORS <- 3

# The 3rd trimester's donors contribute 1, 1, 3, 15 and 15 samples, so its median
# structure rests on two people and its smallest attainable p is 1.5e-3. It is
# not null -- it rejects above the null rate -- but it must not be presented as
# though it carried the same weight as the windows either side of it.
UNDERPOWERED_WINDOWS <- "3rd trimester (n = 5)"


# ---- gene set size --------------------------------------------------------
# Effective size, not length(): 2,894 of the 18,650 sets repeat symbols
# (BIOCARTA_41BB_PATHWAY opens with ATF2 five times), and a symbol the matrix
# does not carry contributes nothing to a score. Both inflate length() and would
# push a set into the wrong size band.
effective_sizes <- function(genesets, universe) {
  vapply(genesets, function(g) length(intersect(unique(g), universe)), integer(1))
}


# ---- cell means -----------------------------------------------------------
# scores: matrix [gene sets x samples], raw (unnormalised) ssGSEA.
# meta:   one row per sample, columns column_num, structure_original, broad_age.
# Returns [gene sets x cells] with a `cells` attribute describing the columns.
cell_means <- function(scores, meta) {
  meta <- meta %>%
    filter(!is.na(structure_original), !is.na(broad_age)) %>%
    mutate(cell = paste(broad_age, structure_original, sep = " || "))
  keep <- meta$column_num[meta$column_num %in% as.numeric(colnames(scores))]
  meta <- meta[match(keep, meta$column_num), ]
  idx  <- match(as.character(meta$column_num), colnames(scores))

  # One matrix multiply rather than a split-apply: 18,646 x 524 by 524 x n_cells.
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
# lights up every structure of its peak window -- GOBP_MYELIN_ASSEMBLY fires on
# all eight adult structures at p = 0.002, which is a true statement about time
# and a false one about place.
within_window_delta <- function(A, cells) {
  Dw <- A
  for (w in unique(cells$broad_age)) {
    j <- which(cells$broad_age == w)
    Dw[, j] <- A[, j, drop = FALSE] - rowMeans(A[, j, drop = FALSE])
  }
  Dw
}


# ---- size adjustment ------------------------------------------------------
# ssGSEA spread falls ~2.2x from the smallest gene sets to the largest, so an
# unadjusted score is mostly a statement about set size. A mean+sd model does not
# fix it -- the tail is far thinner than Gaussian and unevenly so across sizes
# (P(z > 1.645) ranges 0.003 to 0.034 where it should be 0.05), which is why the
# adjustment is only a pre-whitening step and the p-value below is an empirical
# rank, never a normal quantile.
size_adjust <- function(Dw, log_eff, f = 0.3) {
  # Fits are kept on the grid of distinct sizes as well as applied, so a gene
  # list that was never in the reference can be placed on the same scale by
  # interpolation instead of refitting 87 lowess curves per query.
  grid <- sort(unique(log_eff))
  R   <- Dw
  loc_grid <- matrix(NA_real_, length(grid), ncol(Dw))
  scl_grid <- matrix(NA_real_, length(grid), ncol(Dw))

  for (j in seq_len(ncol(Dw))) {
    lo  <- stats::lowess(log_eff, Dw[, j], f = f)
    loc <- suppressWarnings(stats::approx(lo$x, lo$y, xout = log_eff, rule = 2)$y)
    ls  <- stats::lowess(log_eff, abs(Dw[, j] - loc), f = f)
    scl <- suppressWarnings(stats::approx(ls$x, ls$y, xout = log_eff, rule = 2)$y)
    floor_scl <- stats::median(scl[scl > 0], na.rm = TRUE)
    scl[!is.finite(scl) | scl <= 0] <- floor_scl
    R[, j] <- (Dw[, j] - loc) / scl
    loc_grid[, j] <- suppressWarnings(stats::approx(lo$x, lo$y, xout = grid, rule = 2)$y)
    sg <- suppressWarnings(stats::approx(ls$x, ls$y, xout = grid, rule = 2)$y)
    sg[!is.finite(sg) | sg <= 0] <- floor_scl
    scl_grid[, j] <- sg
  }
  list(R = R, grid = grid, loc = loc_grid, scl = scl_grid)
}


# Place a query's within-window deltas on the reference's adjusted scale.
size_adjust_query <- function(Dw_query, log_eff_query, fit) {
  vapply(seq_along(Dw_query), function(j) {
    loc <- suppressWarnings(stats::approx(fit$grid, fit$loc[, j], xout = log_eff_query, rule = 2)$y)
    scl <- suppressWarnings(stats::approx(fit$grid, fit$scl[, j], xout = log_eff_query, rule = 2)$y)
    (Dw_query[j] - loc) / scl
  }, numeric(1))
}


# ---- competitive p --------------------------------------------------------
# Rank the query against the k reference sets nearest it in log effective size.
#
# k is a power/resolution trade-off with a hard floor: the smallest p this can
# return is 1/(k+1), and after BH over ~70 cells the smallest reachable q is
# 70/(k+1). At k = 1000 that is 0.070 -- above any sane threshold, so nothing
# could ever be called significant. k = 4000 puts the floor at 0.0175.
# `self` is the query's own row in the reference, when it has one. An ontology
# term is a member of the library it is being compared against, and leaving it in
# means it counts as a set that beat itself -- inflating every ontology p-value
# by exactly one count while a user gene list, absent from the library, pays no
# such penalty. Dropping it makes the two modes comparable and the +1 below the
# ordinary permutation correction rather than a double one.
p_competitive <- function(R_ref, log_eff_ref, r_query, eff_query, k = 4000, self = NULL) {
  ord <- order(abs(log_eff_ref - log(eff_query)))
  if (!is.null(self)) ord <- ord[ord != self]
  band <- ord[seq_len(min(k, length(ord)))]
  vapply(seq_along(r_query), function(j) {
    (1 + sum(R_ref[band, j] >= r_query[j])) / (length(band) + 1)
  }, numeric(1))
}


# ---- spatial p ------------------------------------------------------------
# Permute structure labels WITHIN each donor. Donor identity explains up to 88%
# of within-window variance, and the design is unbalanced -- dorsal thalamus is
# 1st-trimester-only, cerebellum drops out after the 2nd -- so shuffling labels
# freely across a window would let "this donor ran hot" masquerade as "this
# structure is enriched". Donor-centring first, then permuting inside the donor,
# tests the only thing the data can speak to: whether a structure sits above its
# own donor's average, repeatedly, across different donors.
#
# x: named numeric, one raw score per sample (names are column_num).
p_spatial <- function(x, meta, B = 10000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  meta <- meta %>%
    filter(!is.na(structure_original), !is.na(broad_age),
           as.character(column_num) %in% names(x))
  out <- list()

  for (w in unique(meta$broad_age)) {
    mw <- meta[meta$broad_age == w, ]
    # A donor sampled in only one structure has zero residual after centring and
    # nothing left to permute; keeping them only dilutes.
    multi <- names(which(tapply(mw$structure_original, mw$donor_id,
                                function(s) length(unique(s))) > 1))
    mw <- mw[mw$donor_id %in% multi, ]
    if (nrow(mw) == 0) next

    v <- x[as.character(mw$column_num)]
    r <- v - ave(v, mw$donor_id, FUN = mean)
    s <- factor(mw$structure_original)
    d <- factor(mw$donor_id)

    M <- model.matrix(~ 0 + s)
    M <- sweep(M, 2, colSums(M), "/")          # column means
    obs <- as.numeric(crossprod(M, r))

    # B permutations at once: shuffle row order within donor, reuse M.
    perm <- replicate(B, ave(seq_along(d), d, FUN = sample))
    null <- crossprod(M, matrix(r[perm], nrow = length(r)))

    out[[w]] <- data.frame(
      broad_age          = w,
      structure_original = levels(s),
      T_obs              = obs,
      p_spatial          = (1 + rowSums(null >= obs)) / (B + 1),
      stringsAsFactors   = FALSE)
  }
  bind_rows(out)
}


# ---- combine --------------------------------------------------------------
# BH over the whole family of tested cells for this gene set (~70), not per
# window. Five per-window families would be more powerful but would run five
# error budgets across one figure, taking the chance of at least one false cell
# from 3% to 11%.
#
# The two tests are combined by intersection, not by pooling p-values. They are
# dependent through the same cell means, and Fisher/Stouffer would let one very
# small p carry a cell that fails the other test, which is exactly what the
# second test exists to prevent.
combine_and_adjust <- function(df, alpha = 0.05) {
  df %>%
    mutate(
      tested = .data$n_donors >= MIN_DONORS,
      q_competitive = ifelse(tested, p.adjust(ifelse(tested, p_competitive, NA), "BH"), NA),
      q_spatial     = ifelse(tested, p.adjust(ifelse(tested, p_spatial,     NA), "BH"), NA),
      enriched = !is.na(q_competitive) & !is.na(q_spatial) &
                 q_competitive < alpha & q_spatial < alpha,
      underpowered_window = .data$broad_age %in% UNDERPOWERED_WINDOWS)
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
    summarise(mean_ssGSEA = mean(mean_ssGSEA),
              delta_within_window = mean(delta_within_window),
              enriched     = any(enriched),
              q_competitive = min(q_competitive, na.rm = TRUE),
              q_spatial     = min(q_spatial,     na.rm = TRUE),
              n_structures  = n(),
              tested        = any(tested),
              .groups = "drop") %>%
    mutate(across(c(q_competitive, q_spatial), ~ ifelse(is.finite(.x), .x, NA_real_)))
}


# ---- top level ------------------------------------------------------------
# One gene set in, one row per (structure, window) out.
#   ref      enrichment_reference.rds
#   x        named numeric, raw per-sample ssGSEA (names are column_num)
#   eff_size the query's effective size
enrich <- function(ref, x, eff_size, B = 10000, alpha = 0.05, k = 4000, seed = 1,
                   self = NULL) {
  qm <- matrix(x, nrow = 1, dimnames = list("query", names(x)))
  Aq <- cell_means(qm, ref$meta)
  Aq <- as.numeric(Aq)[match(ref$cells$cell, attr(Aq, "cells")$cell)]

  Dq <- Aq
  for (w in unique(ref$cells$broad_age)) {
    j <- which(ref$cells$broad_age == w)
    Dq[j] <- Aq[j] - mean(Aq[j])
  }
  rq <- size_adjust_query(Dq, log(eff_size), ref$fit)

  sp <- p_spatial(x, ref$meta, B = B, seed = seed)

  ref$cells %>%
    mutate(mean_ssGSEA_raw      = Aq,
           mean_ssGSEA          = Aq / ref$K_PRE,
           delta_within_window  = Dq,
           size_adjusted_effect = rq,
           p_competitive        = p_competitive(ref$fit$R, ref$log_eff, rq, eff_size,
                                                k = k, self = self)) %>%
    left_join(sp[, c("broad_age", "structure_original", "p_spatial")],
              by = c("broad_age", "structure_original")) %>%
    # A window where every donor sat in one structure yields no permutable
    # residual; such cells cannot be spatially supported, so they fail rather
    # than pass by default.
    mutate(p_spatial = ifelse(is.na(p_spatial), 1, p_spatial)) %>%
    combine_and_adjust(alpha = alpha)
}


# ---- raw per-sample scores for a query ------------------------------------
# Precomputed set: pull its 524-row block and undo the normalisation, so it
# lands on the same raw scale as the reference.
scores_ontology <- function(geneset, ref, root = ".") {
  txt <- system2("grep", c(shQuote(paste0("^", geneset, ",")),
                           shQuote(file.path(root, "ontologyssGSEA.csv"))), stdout = TRUE)
  if (length(txt) == 0) stop("gene set not in ontologyssGSEA.csv: ", geneset)
  d <- utils::read.csv(text = paste(c("GeneSet,Amostra,Score_ssGSEA", txt), collapse = "\n"))
  stats::setNames(d$Score_ssGSEA * ref$K_PRE, as.character(d$Amostra))
}

# User gene list: normalize = FALSE is what puts an on-demand score on the
# reference's scale. With normalize = TRUE, GSVA would divide by the range of
# this call alone -- one gene set -- and the result would not be comparable to
# anything precomputed.
scores_genelist <- function(genes, matrix_dados) {
  valid <- intersect(unique(genes), rownames(matrix_dados))
  if (length(valid) < 5) stop("need at least 5 genes present in the matrix; got ", length(valid))
  res <- GSVA::gsva(GSVA::ssgseaParam(exprData = matrix_dados,
                                      geneSets = list(UserSet = valid),
                                      minSize = 1, normalize = FALSE), verbose = FALSE)
  list(x = stats::setNames(as.numeric(res["UserSet", ]), colnames(res)),
       eff = length(valid), missing = setdiff(unique(genes), valid))
}
