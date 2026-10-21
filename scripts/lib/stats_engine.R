# ============================================================================
# stats_engine.R
#
# Vectorized Spearman correlation, p-value and OLS residualization routines
# used throughout the pipeline, together with shared statistical helpers
# for the immune-adjustment, clustering, meta-analysis and matched-null
# analyses. validate_stats_engine() checks the vectorized routines against
# cor.test(method = "spearman", exact = FALSE) and per-row lm()
# residualization on real data before a stage proceeds, so that matrix
# computation over millions of miRNA-gene-cancer combinations reproduces
# the scalar textbook definitions within floating-point tolerance.
# ============================================================================

# rank_rows() and safe_rank() are the ranking entry points for Spearman
# computation on fully observed data (one_vs_many_spearman() assumes its
# inputs were ranked by one of them). Both stop on any non-finite value,
# because base R's rank() with the default na.last = TRUE assigns NA the
# largest rank, as if a missing value were the most extreme observation.
# Inputs that may contain missing values are handled by
# one_vs_many_spearman_pairwise() below.
rank_rows <- function(m) {
  n_bad <- sum(!is.finite(m))
  if (n_bad > 0) {
    stop(sprintf(
      "rank_rows(): input matrix contains %d non-finite (NA/NaN/Inf) value(s) out of %d. rank()'s default na.last=TRUE would silently assign each missing value the LARGEST rank rather than excluding it, corrupting the Spearman correlation for any row containing one. Refusing to rank; the caller must resolve or exclude missing values first.",
      n_bad, length(m)), call. = FALSE)
  }
  t(apply(m, 1, rank))
}

# Single-vector counterpart of rank_rows() (e.g. a CYT vector, a
# residualized CYT vector, or one miRNA's expression row).
safe_rank <- function(x, label = "vector") {
  n_bad <- sum(!is.finite(x))
  if (n_bad > 0) {
    stop(sprintf(
      "safe_rank(): %s contains %d non-finite (NA/NaN/Inf) value(s) out of %d. rank()'s default na.last=TRUE would silently assign each missing value the LARGEST rank rather than excluding it, corrupting the resulting Spearman correlation. Refusing to rank; the caller must resolve or exclude missing values first.",
      label, n_bad, length(x)), call. = FALSE)
  }
  rank(x)
}

# Implements the asymptotic t-approximation used by
# cor.test(method = "spearman", exact = FALSE): t = r * sqrt((n - 2) / (1 - r^2)),
# p = 2 * pt(-|t|, df = n - 2). Agreement with cor.test() to < 1e-9 across
# sample sizes, effect sizes, signs and tied ranks is tested in
# tests/test_spearman_pvalue_equivalence.R.
#
# As |r| approaches 1 the denominator 1 - r^2 approaches zero and the
# statistic becomes numerically unstable: machine-epsilon differences in r
# produce order-of-magnitude differences in p-values that are already far
# below any significance threshold (e.g. 1.06e-62 vs 1.70e-61 for n = 10,
# |r| = 1). This boundary is handled explicitly below.
spearman_pvalue <- function(r, n) {
  one_minus_r2 <- 1 - r^2
  # At or near |r| = 1 the denominator is numerically zero; a perfect
  # monotonic relationship is the strongest possible evidence against the
  # null, so the smallest representable non-negative p-value is returned.
  ifelse(one_minus_r2 <= .Machine$double.eps^0.5, 0,
         2 * pt(-abs(r * sqrt((n - 2) / one_minus_r2)), df = n - 2))
}

one_vs_many_spearman <- function(x_ranked, Y_ranked) {
  n <- length(x_ranked)
  r <- as.numeric(cor(x_ranked, t(Y_ranked)))
  list(cor = r, p = spearman_pvalue(r, n), n = n)
}

residualize_matrix <- function(Y, covariate) {
  if (!all(is.finite(covariate))) {
    stop(sprintf("residualize_matrix(): covariate contains %d non-finite value(s) out of %d. The caller must restrict to samples with a finite covariate value before calling this function -- this has always been this function's contract (every existing caller pre-filters its index set on the covariate); it is Y, not covariate, that may legitimately contain missing values per-row (see the pairwise-complete fallback below).",
                 sum(!is.finite(covariate)), length(covariate)), call. = FALSE)
  }
  xc <- covariate - mean(covariate)
  ssq <- sum(xc^2)
  if (!is.finite(ssq) || ssq < .Machine$double.eps^0.5) {
    stop(sprintf("residualize_matrix(): covariate has (near-)zero variance (sum of centered squares = %.3e over %d values); OLS residualization is undefined. This must fail loudly rather than silently divide by zero and propagate NaN/Inf into downstream ranks and correlations (see reports/AUDIT_CYCLE_2_STATISTICAL_VALIDITY.md, finding C2-003).",
                 ssq, length(covariate)), call. = FALSE)
  }

  out <- matrix(NA_real_, nrow = nrow(Y), ncol = ncol(Y), dimnames = dimnames(Y))
  row_bad <- rowSums(!is.finite(Y)) > 0
  clean_idx <- which(!row_bad)
  affected_idx <- which(row_bad)

  # Fully finite rows: closed-form vectorized computation; one matrix
  # product gives every row's OLS slope and intercept against the covariate.
  if (length(clean_idx) > 0) {
    Yc <- Y[clean_idx, , drop = FALSE]
    b <- as.numeric(Yc %*% xc) / ssq
    a <- rowMeans(Yc) - b * mean(covariate)
    out[clean_idx, ] <- Yc - outer(b, covariate) - a
  }

  # Rows with at least one non-finite value are residualized individually
  # on that row's finite samples (an OLS fit on the row's complete-case
  # subset, as in reference_residuals() below). The closed-form path above
  # computes one slope per row from all columns, so a single missing value
  # would otherwise make every residual in that row missing.
  for (i in affected_idx) {
    y <- Y[i, ]
    ok <- is.finite(y)
    # sum(ok) < 3 is the minimum for a one-covariate OLS fit to be defined
    # (intercept, slope and at least one residual degree of freedom). The
    # prespecified scientific minimum sample size
    # (analysis_config$discovery$min_samples_per_cancer) is applied once,
    # at the correlation step (one_vs_many_spearman_pairwise()).
    if (sum(ok) < 3 || length(unique(covariate[ok])) < 2) next  # left NA (non-estimable): too few finite samples, or a constant covariate within the row's complete-case subset
    fit <- tryCatch(lm(y[ok] ~ covariate[ok]), error = function(e) NULL)
    if (!is.null(fit)) out[i, ok] <- residuals(fit)
  }
  out
}

spearman_pair <- function(x, y, min_n = 3) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < min_n || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(list(rho = NA_real_, p = NA_real_, n = length(x)))
  }
  r <- suppressWarnings(cor(rank(x), rank(y)))
  list(rho = r, p = spearman_pvalue(r, length(x)), n = length(x))
}

# Pairwise-complete Spearman correlation of one reference vector `x` with
# every row of a matrix `Y`. For each row the analysis set is the samples
# finite in both x and that row, and both are re-ranked within that subset;
# ranking once on the full sample set and deleting missing positions
# afterwards would not reproduce pairwise-complete ranks. This function is
# used wherever `x` or `Y` may contain non-finite values; the strict
# rank_rows() + safe_rank() + one_vs_many_spearman() path refuses such
# input.
#
# Implementation: samples with a non-finite value in `x` are excluded from
# every row's subset, so they are dropped once up front (an exact
# reduction). Rows of Y that are then fully finite use the vectorized path
# (rank_rows() + one_vs_many_spearman()); the remaining rows use a per-row
# call to spearman_pair(), which re-ranks on that row's complete-case
# subset (validated against stats::cor.test() in
# tests/test_pairwise_complete_spearman.R).
#
# `min_n` is the prespecified minimum paired sample size; production callers
# pass analysis_config$discovery$min_samples_per_cancer explicitly. A row
# whose paired n is below `min_n`, or whose paired x or y has zero variance
# within the subset, is returned as non-estimable (cor = NA, p = NA,
# non_estimable = TRUE). Every row of Y produces exactly one result, in the
# order of Y's rows.
#
# Returns a list of four vectors of length nrow(Y): $cor, $p, $n (the paired
# sample size used for that row) and $non_estimable (logical).
one_vs_many_spearman_pairwise <- function(x, Y, min_n = 3) {
  if (is.null(dim(Y))) Y <- matrix(Y, nrow = 1, dimnames = list(NULL, names(Y)))
  n_rows <- nrow(Y)
  out_cor <- rep(NA_real_, n_rows)
  out_p <- rep(NA_real_, n_rows)
  out_n <- rep(0L, n_rows)
  out_ne <- rep(TRUE, n_rows)

  x_ok <- is.finite(x)
  if (sum(x_ok) == 0) {
    return(list(cor = out_cor, p = out_p, n = out_n, non_estimable = out_ne))
  }
  x2 <- x[x_ok]
  Y2 <- Y[, x_ok, drop = FALSE]
  n_avail <- length(x2)

  row_bad <- rowSums(!is.finite(Y2)) > 0
  clean_idx <- which(!row_bad)
  affected_idx <- which(row_bad)

  if (length(clean_idx) > 0) {
    if (n_avail >= min_n && length(unique(x2)) >= 2) {
      ranked_x <- rank(x2)
      ranked_Yc <- rank_rows(Y2[clean_idx, , drop = FALSE])
      res_c <- suppressWarnings(one_vs_many_spearman(ranked_x, ranked_Yc))
      # A row that is constant within the subset yields NaN from cor();
      # it is recorded explicitly as non-estimable.
      bad_c <- !is.finite(res_c$cor)
      out_cor[clean_idx[!bad_c]] <- res_c$cor[!bad_c]
      out_p[clean_idx[!bad_c]] <- res_c$p[!bad_c]
      out_n[clean_idx] <- n_avail
      out_ne[clean_idx] <- bad_c
    } else {
      out_n[clean_idx] <- n_avail
      out_ne[clean_idx] <- TRUE
    }
  }

  for (ri in affected_idx) {
    sp <- spearman_pair(x2, Y2[ri, ], min_n = min_n)
    out_cor[ri] <- sp$rho
    out_p[ri] <- sp$p
    out_n[ri] <- sp$n
    out_ne[ri] <- is.na(sp$rho)
  }

  list(cor = out_cor, p = out_p, n = out_n, non_estimable = out_ne)
}

# Residualizes each column of `vars_df` on all columns of `covariate_df`
# jointly, one lm() fit per column, on the shared complete-case subset. This
# is the multi-covariate counterpart of residualize_matrix(); at the scale
# it is used (per-triplet, per-cancer covariate sets) lm() is itself the
# reference calculation.
residualize_columns_multi <- function(vars_df, covariate_df) {
  ok <- complete.cases(vars_df) & complete.cases(covariate_df)
  resid <- matrix(NA_real_, nrow = nrow(vars_df), ncol = ncol(vars_df), dimnames = list(NULL, colnames(vars_df)))
  if (sum(ok) < 10) return(resid)
  for (j in seq_len(ncol(vars_df))) {
    fit <- tryCatch(lm(vars_df[ok, j] ~ as.matrix(covariate_df[ok, , drop = FALSE])), error = function(e) NULL)
    if (!is.null(fit)) resid[ok, j] <- residuals(fit)
  }
  resid
}

# --- Stage 12 immune-adjustment helpers (shared by
# scripts/12_immune_composition_analysis.R and the synthetic tests in
# tests/test_immune_adjustment_scenarios.R) --------------------------------

# Validates the population names observed in estimator output (e.g.
# rownames(MCPcounter.estimate(...))) against the names expected from an
# independently validated source (the signature table's "Cell population"
# grouping), never against a hard-coded count or list. Both arguments are
# treated as sets of names, so order is irrelevant and downstream code can
# rely on name-keyed matching. Discrepancies are reported by name.
#
# Returns a list: $non_empty (every observed name has at least one
# non-whitespace character), $unique (no duplicated observed name),
# $missing (expected names absent from observed), $unexpected (observed
# names absent from expected), $ok (all of the above hold), $n_observed,
# $n_expected.
validate_population_schema <- function(observed, expected) {
  non_empty <- length(observed) > 0 && all(nzchar(trimws(observed)))
  is_unique <- !any(duplicated(observed))
  missing <- sort(setdiff(expected, observed))
  unexpected <- sort(setdiff(observed, expected))
  list(ok = non_empty && is_unique && length(missing) == 0 && length(unexpected) == 0,
       non_empty = non_empty, unique = is_unique,
       missing = missing, unexpected = unexpected,
       n_observed = length(observed), n_expected = length(expected))
}

# Zero-variance and collinearity diagnostics for a two-covariate adjustment
# (CPE plus one immune-population score). For two covariates the VIF equals
# 1 / (1 - r^2) between them. Returns a list with $collinearity_r,
# $collinearity_vif, $zero_variance_population and $aliased_covariate
# (TRUE when |r| > 0.999, the aliasing threshold used by Stage 12).
immune_covariate_diagnostics <- function(pop_v, cpe_v) {
  pop_v_finite <- pop_v[is.finite(pop_v) & is.finite(cpe_v)]
  cpe_v_finite <- cpe_v[is.finite(pop_v) & is.finite(cpe_v)]
  zero_variance_population <- length(unique(pop_v_finite)) < 2
  collinearity_r <- if (!zero_variance_population && length(unique(cpe_v_finite)) >= 2 && length(pop_v_finite) >= 3)
    suppressWarnings(cor(cpe_v_finite, pop_v_finite)) else NA_real_
  collinearity_vif <- if (!is.na(collinearity_r)) 1 / (1 - collinearity_r^2) else NA_real_
  aliased_covariate <- !is.na(collinearity_r) && abs(collinearity_r) > 0.999
  list(collinearity_r = collinearity_r, collinearity_vif = collinearity_vif,
       zero_variance_population = zero_variance_population, aliased_covariate = aliased_covariate,
       n_finite = length(pop_v_finite))
}

# Builds the gene -> overlapping-populations exclusion map. split() is used
# because setNames() followed by [[ ]] returns only the first match for a
# duplicated name, whereas a gene may overlap several populations.
build_gene_exclusion_map <- function(gene_overlap) {
  split(gene_overlap$overlapping_population, gene_overlap$primary_gene)
}

# Computes BH-FDR once per unique biological hypothesis and broadcasts the
# result back to repeated triplet rows. `assoc_subset` contains one arm and
# one immune population; `identity_col` is "mirna_id" for miRNA-CYT and
# "gene" for gene-CYT. Repeated rows for one biological hypothesis must
# carry the same p-value.
dedup_shared_arm_fdr <- function(assoc_subset, identity_col) {
  required <- c(identity_col, "cancer", "p")
  missing_cols <- setdiff(required, names(assoc_subset))
  if (length(missing_cols) > 0) {
    stop(sprintf("dedup_shared_arm_fdr(): missing required column(s): %s",
                 paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  keys <- data.table::data.table(
    identity = assoc_subset[[identity_col]],
    cancer = assoc_subset$cancer,
    p = assoc_subset$p
  )
  inconsistent <- keys[, .(n_distinct_p = data.table::uniqueN(p, na.rm = FALSE)),
                       by = .(identity, cancer)][n_distinct_p > 1]
  if (nrow(inconsistent) > 0) {
    stop(sprintf(
      "dedup_shared_arm_fdr(): %d repeated (%s, cancer) hypothesis key(s) carry inconsistent p-values",
      nrow(inconsistent), identity_col
    ), call. = FALSE)
  }

  key_dt <- unique(keys)
  key_dt[, fdr_key := p.adjust(p, "BH"), by = cancer]
  lookup_key <- paste(key_dt$identity, key_dt$cancer, sep = "\r")
  row_key <- paste(keys$identity, keys$cancer, sep = "\r")
  key_dt$fdr_key[match(row_key, lookup_key)]
}

dedup_mirna_cyt_fdr <- function(assoc_subset) {
  dedup_shared_arm_fdr(assoc_subset, "mirna_id")
}

dedup_gene_cyt_fdr <- function(assoc_subset) {
  dedup_shared_arm_fdr(assoc_subset, "gene")
}

reference_cor_test <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(list(estimate = NA_real_, p.value = NA_real_, n = length(x)))
  }
  out <- tryCatch(suppressWarnings(cor.test(x, y, method = method, exact = FALSE)), error = function(e) NULL)
  if (is.null(out)) return(list(estimate = NA_real_, p.value = NA_real_, n = length(x)))
  list(estimate = unname(out$estimate)[1], p.value = out$p.value, n = length(x))
}

reference_residuals <- function(y, covariate) {
  ok <- is.finite(y) & is.finite(covariate)
  out <- rep(NA_real_, length(y))
  if (sum(ok) < 3 || length(unique(covariate[ok])) < 2) return(out)
  fit <- tryCatch(lm(y[ok] ~ covariate[ok]), error = function(e) NULL)
  if (is.null(fit)) return(out)
  out[ok] <- residuals(fit)
  out
}

validate_stats_engine <- function(x_matrix, y_reference, covariate = NULL, n_checks = 15, seed = 1) {
  set.seed(seed)
  rows_to_check <- sample(seq_len(nrow(x_matrix)), min(n_checks, nrow(x_matrix)))
  max_cor_diff <- 0
  max_p_diff <- 0
  for (i in rows_to_check) {
    ref <- reference_cor_test(x_matrix[i, ], y_reference, "spearman")
    # Exercise the production ranking and correlation path
    # (rank_rows()/one_vs_many_spearman(), as used by Stages 06, 08 and 10)
    # rather than an inline reimplementation. Restricted to the finite
    # (x, y) pairs used by reference_cor_test(), because rank_rows() stops
    # on non-finite input.
    ok <- is.finite(x_matrix[i, ]) & is.finite(y_reference)
    if (sum(ok) >= 3) {
      x_ranked <- rank_rows(x_matrix[i, ok, drop = FALSE])
      y_ranked <- safe_rank(y_reference[ok])
      prod <- one_vs_many_spearman(y_ranked, x_ranked)
      r <- prod$cor[1]
      p <- prod$p[1]
      if (is.finite(ref$estimate)) {
        max_cor_diff <- max(max_cor_diff, abs(ref$estimate - r))
        max_p_diff <- max(max_p_diff, abs(ref$p.value - p))
      }
    }
  }
  resid_diff <- NA_real_
  if (!is.null(covariate)) {
    one_row <- x_matrix[rows_to_check[1], , drop = FALSE]
    ref_resid <- reference_residuals(one_row[1, ], covariate)
    vec_resid <- residualize_matrix(one_row, covariate)[1, ]
    resid_diff <- max(abs(ref_resid - vec_resid))
  }
  list(max_cor_diff = max_cor_diff, max_p_diff = max_p_diff, max_resid_diff = resid_diff)
}

# --- Stage 18 clustering null-calibration helpers -----------------------

# Adjusted Rand Index between two clusterings. Used by the stability,
# seed-sensitivity and scaling-sensitivity checks in
# scripts/17_clustering_and_state_analysis.R and by
# run_null_cluster_calibration() below, so observed and null stability
# statistics use the same formula.
adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  sum_comb <- function(x) sum(choose(x, 2))
  a_sums <- rowSums(tab); b_sums <- colSums(tab)
  index <- sum_comb(tab)
  expected <- sum_comb(a_sums) * sum_comb(b_sums) / choose(n, 2)
  max_index <- 0.5 * (sum_comb(a_sums) + sum_comb(b_sums))
  if (max_index == expected) return(1)
  (index - expected) / (max_index - expected)
}

# Independently permutes each row (feature) of `mat` across its columns
# (samples). Each feature's marginal distribution and the matrix dimensions
# are preserved exactly, while cross-feature, sample-level covariance (the
# structure multivariate clustering depends on) is destroyed. This is the
# standard permutation null for cluster-significance testing.
permute_features_independently <- function(mat) {
  out <- t(apply(mat, 1, sample))
  dimnames(out) <- dimnames(mat)
  out
}

# Finite-sample empirical p-value with the Davison & Hinkley (1997)
# correction, for a one-sided test of an unusually large observed
# statistic (silhouette and stability ARI are both larger-is-more-clustered
# statistics). The result is never exactly zero.
empirical_p_upper <- function(observed, null_dist) {
  nd <- null_dist[is.finite(null_dist)]
  if (!is.finite(observed) || length(nd) == 0) return(NA_real_)
  (1 + sum(nd >= observed)) / (1 + length(nd))
}

# Standard error of the Fisher z-transform of a Spearman correlation, using
# the Fieller, Hartley & Pearson (1957) variance correction:
# Var(z_s) ~= correction / (n - 3), so SE = sqrt(correction / (n - 3)).
# Used by scripts/13_cross_cancer_meta_analysis.R for every REML
# meta-analysis model. `correction` is normally
# analysis_config$meta_analysis$fisher_z_spearman_correction (1.06) and is
# an explicit argument so the function has no hidden configuration
# dependency. Returns NA for n <= 3, where the variance is undefined.
fisher_z_se_spearman <- function(n, correction) {
  ifelse(is.na(n) | n <= 3, NA_real_, sqrt(correction / (n - 3)))
}

# --- Stage 19 matched-null target-specificity helpers --------------------
# Shared by scripts/18_matched_null_target_specificity_test.R and the
# known-answer tests in tests/test_matched_null_helpers.R.

# Raw-nominal-p three-arm pass count, applied identically to a real
# (miRNA, gene) pair and, called once per null gene by the caller, to each
# matched-null gene: counts how many of `cancers` satisfy |rho| >= rho_thr
# and raw p < p_thr in all three arms (miRNA-CYT, miRNA-gene, gene-CYT),
# using the CPE-adjusted evidence tables. `mirgene`, `genecyt` and `mircyt`
# must be keyed as in scripts/18_matched_null_target_specificity_test.R
# (setkey(mirgene, mirna_id, gene, cancer), etc.); this function only reads.
three_arm_raw_p_pass_count <- function(mir, gene_x, cancers, mirgene, genecyt, mircyt, rho_thr, p_thr) {
  mg <- mirgene[.(mir, gene_x, cancers), on = c("mirna_id", "gene", "cancer")]
  gc <- genecyt[.(gene_x, cancers), on = c("gene", "cancer")]
  mc <- mircyt[.(mir, cancers), on = c("mirna_id", "cancer")]
  # Cardinality guard: a duplicated key for one of `cancers` would make the
  # X[Y, on =] join return more than one row for that cancer and misalign
  # the pass_mc/pass_mg/pass_gc vectors through recycling. Upstream
  # cardinality assertions (Stages 10 and 14) already guarantee unique keys,
  # so this check is expected never to trigger.
  if (nrow(mg) != length(cancers) || nrow(gc) != length(cancers) || nrow(mc) != length(cancers)) {
    stop(sprintf(
      "three_arm_raw_p_pass_count(): evidence-table join returned an unexpected row count for mirna_id='%s', gene='%s' (expected %d, one per cancer; got mirgene=%d, genecyt=%d, mircyt=%d) -- indicates a duplicated key in one of the evidence tables, which would otherwise silently misalign the pass/fail vectors.",
      mir, gene_x, length(cancers), nrow(mg), nrow(gc), nrow(mc)
    ), call. = FALSE)
  }
  pass_mc <- !is.na(mc$p_mirna_cyt_adj) & abs(mc$cor_mirna_cyt_adj) >= rho_thr & mc$p_mirna_cyt_adj < p_thr
  pass_mg <- !is.na(mg$p_mir_gene_adj) & abs(mg$cor_mir_gene_adj) >= rho_thr & mg$p_mir_gene_adj < p_thr
  pass_gc <- !is.na(gc$p_gene_cyt_adj) & abs(gc$cor_gene_cyt_adj) >= rho_thr & gc$p_gene_cyt_adj < p_thr
  sum(pass_mc & pass_mg & pass_gc, na.rm = TRUE)
}

# Builds the matched-null caliper for one real gene: the standardized
# Euclidean distance from `real_gene` to every gene in `candidate_pool`,
# ordered nearest first and truncated to `caliper_size` (or fewer if the
# pool is smaller). `gene_stats_z` must be keyed by `gene` and z-scored on
# every column in `match_vars`. The random null draw itself is made by the
# caller under its own set.seed(); this function only determines the
# eligible caliper pool.
build_caliper <- function(real_gene, candidate_pool, gene_stats_z, match_vars, caliper_size) {
  real_z <- gene_stats_z[.(real_gene), on = "gene"]
  cand_z <- gene_stats_z[.(candidate_pool), on = "gene"]
  dist <- sqrt(rowSums(sweep(as.matrix(cand_z[, ..match_vars]), 2, as.numeric(real_z[, ..match_vars]), "-")^2))
  ord <- order(dist)
  candidate_pool[ord][seq_len(min(caliper_size, length(candidate_pool)))]
}

# Runs the Stage 18 k-selection procedure (maximum mean silhouette over
# k_range, via cluster::silhouette()) and the stability-ARI computation on
# one matrix, returning the chosen k, its silhouette and its stability ARI.
# Shared by the observed computation and every null iteration.
select_k_and_stability <- function(mat_scaled, k_range, n_init, n_boot_stability, boot_nstart = 10, min_boot_n = 10) {
  X <- t(mat_scaled)
  dist_X <- dist(X)
  best_sil <- -Inf; best_km <- NULL; best_k <- NA_integer_
  for (k in k_range) {
    km <- kmeans(X, centers = k, nstart = n_init, iter.max = 50)
    sil <- mean(cluster::silhouette(km$cluster, dist_X)[, "sil_width"])
    if (sil > best_sil) { best_sil <- sil; best_km <- km; best_k <- k }
  }
  n_s <- nrow(X)
  ari_vals <- rep(NA_real_, n_boot_stability)
  if (n_s >= min_boot_n) {
    for (b in seq_len(n_boot_stability)) {
      idx <- sample(seq_len(n_s), n_s, replace = TRUE)
      km_b <- tryCatch(kmeans(X[idx, , drop = FALSE], centers = best_k, nstart = boot_nstart, iter.max = 50), error = function(e) NULL)
      if (!is.null(km_b)) ari_vals[b] <- adjusted_rand_index(best_km$cluster[idx], km_b$cluster)
    }
  }
  list(k = best_k, max_silhouette = best_sil, stability_ari = mean(ari_vals, na.rm = TRUE), km = best_km)
}

# Null-cluster-structure calibration: B independent permutation-null
# matrices (permute_features_independently(), one draw per iteration), each
# processed with the same within-cancer scaling (within_cancer_scale_fn,
# e.g. within_cancer_scale() from scripts/17_clustering_and_state_analysis.R)
# and the same k-selection and stability procedure
# (select_k_and_stability()) as the observed data, so that observed and
# null results differ only in the presence of multivariate structure.
# Returns the null distributions of maximum silhouette and stability ARI
# (length B each); the caller computes empirical_p_upper() against the
# observed values. Tested in tests/test_cluster_null_calibration.R.
run_null_cluster_calibration <- function(raw_mat, cancer_of, within_cancer_scale_fn, k_range, n_init, B, seed,
                                          n_boot_stability = 20, boot_nstart = 5) {
  set.seed(seed)
  null_max_silhouette <- numeric(B)
  null_stability <- numeric(B)
  for (b in seq_len(B)) {
    null_raw <- permute_features_independently(raw_mat)
    null_scaled <- within_cancer_scale_fn(null_raw, cancer_of)
    res_b <- select_k_and_stability(null_scaled, k_range, n_init, n_boot_stability, boot_nstart)
    null_max_silhouette[b] <- res_b$max_silhouette
    null_stability[b] <- res_b$stability_ari
  }
  list(null_max_silhouette = null_max_silhouette, null_stability = null_stability, B = B)
}

# Checks whether a kmeans() solution is invariant to a substantially larger
# iter.max, i.e. whether "Quick-TRANSfer stage steps exceeded maximum"
# warnings affected the returned best-of-nstart partition or only
# restarts that were not selected. The warning reflects exhaustion of the
# Hartigan-Wong quick-transfer step budget (iter.max * nrow(X)) on a single
# random start; it occurs predominantly on the permutation-null matrices of
# the Stage 18 calibration, which lack cross-feature structure.
#
# If raising iter.max leaves the result unchanged (identical assignment,
# ARI == 1 and tot.withinss equal within a tight relative tolerance), the
# original fit is accepted. Otherwise the function reports
# invariant = FALSE, so an unconverged fit is never accepted silently.
kmeans_convergence_invariant <- function(X, centers, nstart, seed, iter.max = 50,
                                          algorithm = "Hartigan-Wong",
                                          iter.max_check_mult = 10) {
  n_warn_base <- 0
  set.seed(seed)
  km_base <- withCallingHandlers(
    kmeans(X, centers = centers, nstart = nstart, iter.max = iter.max, algorithm = algorithm),
    warning = function(w) {
      if (grepl("Quick-TRANSfer", conditionMessage(w))) n_warn_base <<- n_warn_base + 1
      invokeRestart("muffleWarning")
    }
  )
  n_warn_check <- 0
  set.seed(seed)
  km_check <- withCallingHandlers(
    kmeans(X, centers = centers, nstart = nstart, iter.max = iter.max * iter.max_check_mult, algorithm = algorithm),
    warning = function(w) {
      if (grepl("Quick-TRANSfer", conditionMessage(w))) n_warn_check <<- n_warn_check + 1
      invokeRestart("muffleWarning")
    }
  )
  ari <- adjusted_rand_index(km_base$cluster, km_check$cluster)
  withinss_rel_diff <- abs(km_base$tot.withinss - km_check$tot.withinss) / max(km_check$tot.withinss, .Machine$double.eps)
  list(
    n_warnings_base = n_warn_base, n_warnings_check = n_warn_check,
    tot_withinss_base = km_base$tot.withinss, tot_withinss_check = km_check$tot.withinss,
    withinss_rel_diff = withinss_rel_diff, ari_base_vs_check = ari,
    invariant = isTRUE(ari > 0.999) && isTRUE(withinss_rel_diff < 1e-6)
  )
}
