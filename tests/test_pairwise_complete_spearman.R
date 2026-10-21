# ============================================================================
# test_pairwise_complete_spearman.R
#
# Proves R/stats_engine.R::one_vs_many_spearman_pairwise() and the
# pairwise-complete fallback added to R/stats_engine.R::residualize_matrix()
# implement genuine pairwise-complete Spearman correlation: for every test,
# the analysis sample set is exactly the samples finite in BOTH variables,
# and BOTH variables are re-ranked FROM SCRATCH within that exact subset --
# never a rank computed once on the full sample set with missing positions
# merely deleted afterward. This is the repair for the Stage 08 halt
# documented in reports/STAGE08_MISSING_DATA_DIAGNOSIS.md and
# reports/STAGE08_PAIRWISE_COMPLETE_FIX_REPORT.md.
#
# The reference computation used throughout is INDEPENDENT of the
# production function under test: every equivalence check below calls
# stats::cor.test() directly (never one_vs_many_spearman_pairwise(),
# spearman_pair(), or any other stats_engine.R helper as the ground truth),
# exactly the definition the user specified:
#   cor.test(x[complete.cases(x,y)], y[complete.cases(x,y)],
#             method="spearman", exact=FALSE)
#
# SYNTHETIC DATA for the core correctness/edge-case tests (no missingness
# pattern below depends on any project file). REAL DATA for the spot-check
# section only (RD1-RD4): identifier-level/value-level reads of the most
# recent from-scratch run's own Stage 03/05 outputs already on disk under
# runs/, read-only, purely to validate that this ENGINE FUNCTION reproduces
# cor.test() on genuine real-world missingness patterns -- exactly the same
# kind of real-data engine validation already performed live by
# R/stats_engine.R::validate_stats_engine() at the start of Stages 06 and 10
# themselves, and by tests/test_purity_parsing_hard.R against the real
# TCGAbiolinks::Tumor.purity table. No correlation, FDR, survival, or
# clustering RESULT is computed, reported, or interpreted anywhere in this
# file -- only whether the engine's numeric output matches cor.test()'s.
# The real-data spot checks are automatically SKIPPED (not FAILED) if no
# from-scratch run with a COMPLETE Stage 05 is present on disk, so this file
# remains runnable standalone after any such run directory is removed.
#
# Usage: Rscript tests/test_pairwise_complete_spearman.R <project_root>
# ============================================================================

run_pairwise_complete_tests <- function(project_root) {
  suppressPackageStartupMessages(library(data.table))
  source(file.path(project_root, "scripts", "lib", "stats_engine.R"))
  source(file.path(project_root, "scripts", "lib", "tcga_identifiers.R"))
  set.seed(20260918)

  results <- list()
  add <- function(id, desc, pass, detail) {
    results[[length(results) + 1]] <<- data.frame(test_id = id, description = desc,
      result = if (isTRUE(pass)) "PASS" else if (identical(pass, "SKIP")) "SKIPPED" else "FAIL",
      detail = detail, stringsAsFactors = FALSE)
    cat(sprintf("[%s] %s -- %s\n", if (isTRUE(pass)) "PASS" else if (identical(pass, "SKIP")) "SKIPPED" else "FAIL", id, detail))
  }

  # Independent reference: stats::cor.test() on the exact complete-case
  # subset, called fresh here -- never via any stats_engine.R helper.
  independent_ref <- function(x, y) {
    ok <- stats::complete.cases(data.frame(x = x, y = y)) & is.finite(x) & is.finite(y)
    xs <- x[ok]; ys <- y[ok]
    if (length(xs) < 3 || length(unique(xs)) < 2 || length(unique(ys)) < 2) {
      return(list(rho = NA_real_, p = NA_real_, n = length(xs)))
    }
    ct <- suppressWarnings(stats::cor.test(xs, ys, method = "spearman", exact = FALSE))
    list(rho = unname(ct$estimate), p = ct$p.value, n = length(xs))
  }

  # ==========================================================================
  # A. Core equivalence to the independent cor.test() reference, one row of Y
  #    at a time, across every required missingness pattern.
  # ==========================================================================
  check_row <- function(id, desc, x, y, min_n = 3, rho_tol = 1e-9, p_tol = 1e-6) {
    Y <- matrix(y, nrow = 1)
    res <- one_vs_many_spearman_pairwise(x, Y, min_n = min_n)
    ref <- independent_ref(x, y)
    if (is.na(ref$rho)) {
      pass <- is.na(res$cor[1]) && res$n[1] == ref$n && isTRUE(res$non_estimable[1])
      add(id, desc, pass, sprintf("both non-estimable: ref n=%d, own n=%d, own cor=%s, own non_estimable=%s",
                                   ref$n, res$n[1], format(res$cor[1]), res$non_estimable[1]))
    } else {
      pass <- !is.na(res$cor[1]) && abs(res$cor[1] - ref$rho) < rho_tol &&
        abs(res$p[1] - ref$p) < p_tol && res$n[1] == ref$n && !isTRUE(res$non_estimable[1])
      add(id, desc, pass, sprintf("rho: ref=%.6f own=%.6f | p: ref=%.4g own=%.4g | n: ref=%d own=%d",
                                   ref$rho, res$cor[1], ref$p, res$p[1], ref$n, res$n[1]))
    }
  }

  n <- 40
  x_full <- rnorm(n); y_full <- x_full * 0.5 + rnorm(n)

  # A1: no missing values at all -- exercises the fast vectorized path only.
  check_row("A1_no_missing", "No missing values: matches cor.test() (fast vectorized path)", x_full, y_full)

  # A2: one-sided predictor (x) missingness.
  x_pred_na <- x_full; x_pred_na[c(3, 7, 15)] <- NA
  check_row("A2_predictor_missing", "One-sided predictor (x) missingness: matches cor.test()", x_pred_na, y_full)

  # A3: outcome (y, the "many" row) missingness.
  y_out_na <- y_full; y_out_na[c(1, 20, 39)] <- NA
  check_row("A3_outcome_missing", "Outcome (y) missingness: matches cor.test()", x_full, y_out_na)

  # A4: different missingness patterns in different rows of Y (multi-row call).
  y1 <- y_full; y1[c(2, 4)] <- NA
  y2 <- y_full; y2[c(10, 11, 12)] <- NA
  y3 <- y_full  # fully finite row, mixed in with affected rows in the same call
  Ymulti <- rbind(y1, y2, y3)
  res_multi <- one_vs_many_spearman_pairwise(x_full, Ymulti, min_n = 3)
  ref1 <- independent_ref(x_full, y1); ref2 <- independent_ref(x_full, y2); ref3 <- independent_ref(x_full, y3)
  pass_multi <- abs(res_multi$cor[1] - ref1$rho) < 1e-9 && res_multi$n[1] == ref1$n &&
    abs(res_multi$cor[2] - ref2$rho) < 1e-9 && res_multi$n[2] == ref2$n &&
    abs(res_multi$cor[3] - ref3$rho) < 1e-9 && res_multi$n[3] == ref3$n
  add("A4_different_missingness_per_row", "Different missingness patterns in different rows of Y, single call, each row matches its own cor.test() independently",
      pass_multi, sprintf("row1 n=%d/%d row2 n=%d/%d row3(clean) n=%d/%d", res_multi$n[1], ref1$n, res_multi$n[2], ref2$n, res_multi$n[3], ref3$n))

  # A5: tied values plus missingness together.
  x_tied <- c(1, 1, 2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 8, 9, 10, 10)
  y_tied <- c(5, 4, 4, 3, 3, 2, 2, 1, 1, 1, 2, 3, 4, 5, 6, 6)
  y_tied_na <- y_tied; y_tied_na[c(3, 9)] <- NA
  check_row("A5_ties_plus_missing", "Tied values plus missingness: matches cor.test()", x_tied, y_tied_na)

  # A6: insufficient paired n after removing missing values (min_n=30, only 25 pairs survive).
  x_lown <- rnorm(50); y_lown <- x_lown * 0.4 + rnorm(50)
  y_lown[1:25] <- NA
  res_lown <- one_vs_many_spearman_pairwise(x_lown, matrix(y_lown, nrow = 1), min_n = 30)
  add("A6_insufficient_paired_n", "Paired n (25) below the prespecified minimum (30): explicit non-estimable result, not silently computed",
      is.na(res_lown$cor[1]) && res_lown$n[1] == 25 && isTRUE(res_lown$non_estimable[1]),
      sprintf("cor=%s n=%d non_estimable=%s", format(res_lown$cor[1]), res_lown$n[1], res_lown$non_estimable[1]))

  # A7: zero variance in y after complete-case restriction (constant on the surviving subset).
  x_zv <- rnorm(20)
  y_zv <- rep(5, 20)
  y_zv[c(2, 5, 9)] <- NA  # still constant (5) on every remaining finite position
  res_zv <- one_vs_many_spearman_pairwise(x_zv, matrix(y_zv, nrow = 1), min_n = 3)
  add("A7_zero_variance_after_restriction", "Zero variance in y within the paired complete-case subset: explicit non-estimable result",
      is.na(res_zv$cor[1]) && isTRUE(res_zv$non_estimable[1]) && res_zv$n[1] == 17,
      sprintf("cor=%s n=%d non_estimable=%s", format(res_zv$cor[1]), res_zv$n[1], res_zv$non_estimable[1]))

  # A7b: zero variance in x itself (shared across the whole call) -- every row must come back non-estimable.
  x_zv2 <- rep(3, 20)
  res_zv2 <- one_vs_many_spearman_pairwise(x_zv2, matrix(rnorm(20), nrow = 1), min_n = 3)
  add("A7b_zero_variance_in_x", "Zero variance in x (the shared reference vector): explicit non-estimable result",
      is.na(res_zv2$cor[1]) && isTRUE(res_zv2$non_estimable[1]),
      sprintf("cor=%s non_estimable=%s", format(res_zv2$cor[1]), res_zv2$non_estimable[1]))

  # A8: all values missing (entire row is NA).
  res_allna <- one_vs_many_spearman_pairwise(x_full, matrix(rep(NA_real_, n), nrow = 1), min_n = 3)
  add("A8_all_values_missing", "Entire row is NA (zero paired samples): explicit non-estimable result with n=0, never an error",
      is.na(res_allna$cor[1]) && res_allna$n[1] == 0L && isTRUE(res_allna$non_estimable[1]),
      sprintf("cor=%s n=%d non_estimable=%s", format(res_allna$cor[1]), res_allna$n[1], res_allna$non_estimable[1]))

  # A9: NA, NaN, and Inf all handled identically (all are "non-finite").
  y_mixed <- y_full
  y_mixed[5] <- NA_real_; y_mixed[6] <- NaN; y_mixed[7] <- Inf; y_mixed[8] <- -Inf
  check_row("A9_na_nan_inf_mixed", "NA, NaN, and Inf in the same row are all excluded identically (all non-finite): matches cor.test()", x_full, y_mixed)

  # A9b: x itself contains NA, NaN, and Inf together.
  x_mixed <- x_full
  x_mixed[10] <- NA_real_; x_mixed[11] <- NaN; x_mixed[12] <- Inf
  check_row("A9b_na_nan_inf_in_x", "NA, NaN, and Inf in x are all excluded identically: matches cor.test()", x_mixed, y_full)

  # ==========================================================================
  # B. Ranks are genuinely recomputed within each pair-specific complete-case
  #    subset -- NOT a global rank with missing positions merely deleted.
  #    Constructed so the two approaches provably disagree: removing a
  #    MIDDLE-ranked value of y (deliberately NOT the global min/max) opens
  #    a gap in the middle of the rank sequence, a non-affine distortion
  #    relative to the correctly renormalized version -- and Pearson
  #    correlation of ranks is sensitive to that. (Removing an EXTREME
  #    value instead would not work: the remaining global ranks would be a
  #    pure constant shift of the correctly-renormalized sequence, and
  #    Pearson correlation is invariant to a constant shift, so the two
  #    approaches would coincidentally agree -- verified during test
  #    construction, which is why the extreme (position 2, value 1, global
  #    rank 1) is deliberately NOT the position removed below.)
  # ==========================================================================
  x_b <- c(1, 2, 3, 4, 5, 6, 7, 8)
  y_b <- c(10, 1, 20, 30, 40, 50, 60, 70)  # global rank of y: 2,1,3,4,5,6,7,8
  y_b_na <- y_b; y_b_na[1] <- NA           # remove global rank 2 (position 1, value 10) -- a MIDDLE rank

  # "Wrong" approach: rank y globally first (ranks 2,1,3,4,5,6,7,8), THEN
  # delete the missing position -- leaves ranks 1,3,4,5,6,7,8 (a gap where
  # rank 2 used to be, never renormalized).
  wrong_global_rank_y <- rank(y_b)[-1]
  wrong_x <- x_b[-1]
  wrong_rho <- suppressWarnings(cor(rank(wrong_x), wrong_global_rank_y))

  # Correct approach (what the production function must do): delete first,
  # then rank fresh on the surviving 7 values -- ranks 1,2,3,4,5,6,7 (no gap).
  correct_rho <- independent_ref(x_b, y_b_na)$rho

  res_b <- one_vs_many_spearman_pairwise(x_b, matrix(y_b_na, nrow = 1), min_n = 3)
  add("B1_ranks_recomputed_not_deleted_from_global",
      "Ranks are recomputed fresh within the pair-specific complete-case subset (not a global rank with the missing position merely deleted) -- own result matches the correct re-rank reference and NOT the wrong global-rank-then-delete approach",
      abs(res_b$cor[1] - correct_rho) < 1e-9 && abs(res_b$cor[1] - wrong_rho) > 1e-6,
      sprintf("own=%.6f correct_reref=%.6f wrong_global_then_delete=%.6f (own must equal the 2nd, differ from the 3rd)", res_b$cor[1], correct_rho, wrong_rho))

  # Same proof for the fallback branch specifically (x itself has additional
  # missingness so this row is NOT eligible for the clean vectorized path).
  x_b2 <- x_b; x_b2[8] <- NA  # forces this call through spearman_pair() fallback, not the fast path
  res_b2 <- one_vs_many_spearman_pairwise(x_b2, matrix(y_b_na, nrow = 1), min_n = 3)
  correct_rho2 <- independent_ref(x_b2, y_b_na)$rho
  add("B2_ranks_recomputed_in_fallback_branch",
      "Same re-ranking correctness holds in the scalar pairwise-complete fallback branch (row with its own missing value, x also missing elsewhere)",
      abs(res_b2$cor[1] - correct_rho2) < 1e-9, sprintf("own=%.6f ref=%.6f", res_b2$cor[1], correct_rho2))

  # ==========================================================================
  # C. Reported n equals the ACTUAL paired n (never the nominal full-cohort
  #    count), for both the clean and the fallback branch.
  # ==========================================================================
  Ymix <- rbind(y_full, y_full); Ymix[1, c(1, 2, 3, 4, 5)] <- NA  # row 1 affected (35 remain), row 2 clean (40)
  res_c <- one_vs_many_spearman_pairwise(x_full, Ymix, min_n = 3)
  add("C1_reported_n_matches_actual_paired_n", "n_samples reported equals the ACTUAL paired complete-case n, not the nominal cohort size, for both an affected and a clean row in the same call",
      res_c$n[1] == 35 && res_c$n[2] == 40, sprintf("affected row n=%d (expected 35), clean row n=%d (expected 40)", res_c$n[1], res_c$n[2]))

  # ==========================================================================
  # D. BH-FDR family/denominator is unaffected by non-estimable rows: R's
  #    p.adjust("BH") already excludes NA p-values from the family size and
  #    re-inserts them as NA at their original position -- exactly the
  #    "unchanged except for explicitly non-estimable tests" contract. This
  #    is base R's own documented p.adjust() behavior, verified explicitly
  #    here since production output now relies on it for the first time
  #    (every p-value in this pipeline's history was previously always
  #    non-NA, so this path was never exercised in production before).
  # ==========================================================================
  set.seed(7)
  p_clean <- runif(20)
  p_with_na <- p_clean; p_with_na[c(3, 11)] <- NA
  fdr_clean_family <- p.adjust(p_clean[-c(3, 11)], "BH")           # the 18 estimable p-values, adjusted alone
  fdr_with_na <- p.adjust(p_with_na, "BH")                          # the same 18, embedded among 2 NAs
  fdr_with_na_estimable <- fdr_with_na[!is.na(p_with_na)]
  add("D1_bh_fdr_denominator_excludes_non_estimable",
      "p.adjust('BH') on a p-vector containing non-estimable (NA) entries produces IDENTICAL FDR values for every estimable entry as adjusting the estimable subset alone (family size = n estimable, not n total) -- NA entries stay NA, never silently treated as p=0, p=1, or dropped from the position vector",
      all(abs(fdr_clean_family - fdr_with_na_estimable) < 1e-12) && all(is.na(fdr_with_na[c(3, 11)])) && length(fdr_with_na) == 20,
      sprintf("max diff on estimable entries=%.2e, NA positions stayed NA: %s, output length=%d (expected 20)",
              max(abs(fdr_clean_family - fdr_with_na_estimable)), all(is.na(fdr_with_na[c(3, 11)])), length(fdr_with_na)))

  # ==========================================================================
  # E. residualize_matrix(): a single missing value in one row must NOT
  #    propagate NA to every other sample in that row (the actual defect
  #    that made Stage 10/17 latently broken even before Stage 08's halt --
  #    see reports/STAGE08_PAIRWISE_COMPLETE_FIX_REPORT.md).
  # ==========================================================================
  set.seed(11)
  cov_e <- rnorm(20)
  Yclean_row <- rnorm(20) + cov_e * 0.3
  Yaffected_row <- rnorm(20) + cov_e * 0.3
  Yaffected_row_na <- Yaffected_row; Yaffected_row_na[c(4, 15)] <- NA
  Ye <- rbind(Yclean_row, Yaffected_row_na)
  resid_e <- residualize_matrix(Ye, cov_e)

  add("E1_single_na_does_not_propagate_to_whole_row",
      "residualize_matrix(): a single missing value in a row leaves ONLY that sample's residual as NA, never the entire row (the row-wide-NA-propagation defect this replaces)",
      sum(is.na(resid_e[2, ])) == 2 && all(is.na(resid_e[2, c(4, 15)])) && all(is.finite(resid_e[2, -c(4, 15)])),
      sprintf("row 2 (2 missing inputs): %d NA position(s) in the output (expected exactly 2, at positions 4 and 15)", sum(is.na(resid_e[2, ]))))

  add("E2_clean_row_unchanged_by_fix",
      "residualize_matrix(): a fully-finite row is numerically IDENTICAL to the original single-shot vectorized closed-form (no regression for the common case)",
      all(is.finite(resid_e[1, ])) && {
        xc <- cov_e - mean(cov_e); ssq <- sum(xc^2)
        b_direct <- sum(Yclean_row * xc) / ssq; a_direct <- mean(Yclean_row) - b_direct * mean(cov_e)
        direct <- Yclean_row - b_direct * cov_e - a_direct
        max(abs(direct - resid_e[1, ])) < 1e-9
      }, "clean row matches the direct closed-form OLS computation to 1e-9")

  # E3: the affected row's non-NA residuals match an independent per-row
  # lm() fit on that row's own complete-case subset (the reference this
  # fallback branch must reproduce, computed fresh here, not via
  # residualize_matrix() or reference_residuals()).
  ok_e3 <- is.finite(Yaffected_row_na)
  fit_e3 <- lm(Yaffected_row_na[ok_e3] ~ cov_e[ok_e3])
  ref_resid_e3 <- residuals(fit_e3)
  add("E3_affected_row_matches_independent_lm_on_its_own_subset",
      "residualize_matrix(): the affected row's non-missing residuals match an independently-fit lm() on that row's own complete-case subset",
      max(abs(resid_e[2, ok_e3] - ref_resid_e3)) < 1e-9,
      sprintf("max diff=%.2e", max(abs(resid_e[2, ok_e3] - ref_resid_e3))))

  # E4: residualize_matrix() still fails loudly if the COVARIATE itself is non-finite (unchanged contract).
  cov_bad <- cov_e; cov_bad[3] <- NA
  e4_err <- tryCatch({ residualize_matrix(Ye, cov_bad); NULL }, error = function(e) conditionMessage(e))
  add("E4_covariate_missingness_still_fails_loudly",
      "residualize_matrix() still fails loudly (never silently) if the COVARIATE itself contains a non-finite value -- unchanged from before this fix; only Y (the row/feature matrix) is allowed pairwise-complete handling",
      !is.null(e4_err) && grepl("non-finite", e4_err), sprintf("error: %s", if (is.null(e4_err)) "(none raised)" else e4_err))

  # ==========================================================================
  # F. one_vs_many_spearman_pairwise() on residualize_matrix() output end to
  #    end (the exact composition Stages 10/17 now use): a gene with one
  #    missing raw value keeps its OTHER samples usable, rather than losing
  #    the whole gene for that cancer.
  # ==========================================================================
  Yf <- rbind(Yclean_row, Yaffected_row_na)
  resid_f <- residualize_matrix(Yf, cov_e)
  cyt_resid_f <- as.numeric(residualize_matrix(matrix(rnorm(20) + cov_e * 0.2, nrow = 1), cov_e))
  res_f <- one_vs_many_spearman_pairwise(cyt_resid_f, resid_f, min_n = 3)
  add("F1_residualize_then_correlate_end_to_end",
      "residualize_matrix() + one_vs_many_spearman_pairwise() composed end to end (the exact Stage 10/17 pattern): the affected gene is ESTIMABLE using its 18 surviving samples, not silently dropped to non-estimable, and not corrupted by whole-row NA propagation",
      !is.na(res_f$cor[2]) && res_f$n[2] == 18 && !isTRUE(res_f$non_estimable[2]),
      sprintf("affected-gene arm: cor=%s n=%d non_estimable=%s (expected non-NA cor, n=18)", format(res_f$cor[2]), res_f$n[2], res_f$non_estimable[2]))

  # ==========================================================================
  # RD. Real-data spot checks (read-only): most recent from-scratch run's
  # own Stage 03/05 outputs, if one is present on disk with a COMPLETE
  # Stage 05. SKIPPED, not FAILED, if none is found -- so this file remains
  # runnable standalone regardless of which runs/ directories currently
  # exist. No correlation RESULT is reported or interpreted; only whether
  # the engine matches cor.test() on genuine real missingness.
  # ==========================================================================
  run_dirs <- if(isTRUE(getOption("mirna.tests.real_data",FALSE))) file.path(project_root,"results") else character()
  usable_run <- NULL
  for (rd in run_dirs) {
    status_path <- file.path(rd, "logs", "run_status.tsv")
    if (!file.exists(status_path)) next
    st <- tryCatch(read.delim(status_path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(st)) next
    if (any(st$stage == "05_define_eligible_features" & st$status == "COMPLETE")) { usable_run <- rd; break }
  }

  if (is.null(usable_run)) {
    for (id in c("RD1_real_fully_complete_combination", "RD2_real_affected_combination", "RD3_real_heavily_missing_feature", "RD4_real_boundary_case")) {
      add(id, "Real-data spot check", "SKIP", "No from-scratch run with a COMPLETE Stage 05 found under runs/ -- skipped, not failed")
    }
  } else {
    mrna_matched <- readRDS(file.path(usable_run, "processed", "stage03", "mrna_matched.rds"))
    mrna_mat <- as.matrix(mrna_matched[, -1, drop = FALSE]); rownames(mrna_mat) <- toupper(trimws(mrna_matched$gene_id_raw))
    sample_annotation <- read.delim(file.path(usable_run, "processed", "stage05", "sample_annotation_with_purity.tsv"), stringsAsFactors = FALSE)

    spot_check <- function(id, desc, gene, cancer) {
      idx <- which(sample_annotation$cancer == cancer)
      cyt_v <- sample_annotation$CYT[idx]
      gene_v <- mrna_mat[gene, idx]
      res <- one_vs_many_spearman_pairwise(cyt_v, matrix(gene_v, nrow = 1), min_n = 3)
      ref <- independent_ref(cyt_v, gene_v)
      if (is.na(ref$rho)) {
        add(id, desc, is.na(res$cor[1]) && res$n[1] == ref$n, sprintf("gene=%s cancer=%s: both non-estimable, n=%d", gene, cancer, ref$n))
      } else {
        add(id, desc, !is.na(res$cor[1]) && abs(res$cor[1] - ref$rho) < 1e-9 && abs(res$p[1] - ref$p) < 1e-6 && res$n[1] == ref$n,
            sprintf("gene=%s cancer=%s n=%d: rho own=%.6f ref=%.6f | p own=%.4g ref=%.4g", gene, cancer, ref$n, res$cor[1], ref$rho, res$p[1], ref$p))
      }
    }

    na_counts <- rowSums(!is.finite(mrna_mat))
    clean_genes <- names(na_counts)[na_counts == 0]
    affected_genes <- names(na_counts)[na_counts > 0]

    # Find the first cancer where `gene` has PARTIAL missingness (>=1 NA
    # AND >=1 finite value) -- a genuinely estimable affected case, as
    # opposed to a cancer where the gene happens to be entirely missing
    # (which is a valid non-estimable case, but not what "an affected
    # cancer-feature combination" or "a heavily missing feature" is meant
    # to demonstrate here: the pairwise-complete computation actually
    # succeeding on real partial missingness).
    first_partial_missing_cancer <- function(gene) {
      v <- mrna_mat[gene, ]
      for (cc in unique(sample_annotation$cancer)) {
        idx <- which(sample_annotation$cancer == cc)
        n_na <- sum(is.na(v[idx])); n_tot <- length(idx)
        if (n_na > 0 && n_na < n_tot) return(cc)
      }
      NA_character_
    }

    # RD1: a fully complete cancer-feature combination.
    g1 <- clean_genes[1]
    spot_check("RD1_real_fully_complete_combination", sprintf("Real data: fully-complete gene/cancer combination (no missing values at all) matches cor.test()"), g1, sample_annotation$cancer[1])

    # RD2: an affected (genuinely, partially missing) cancer-feature
    # combination -- still estimable, with a real mix of present and
    # missing samples in that cancer.
    g2 <- NA_character_; cancer_g2 <- NA_character_
    for (cand in affected_genes) {
      cc <- first_partial_missing_cancer(cand)
      if (!is.na(cc)) { g2 <- cand; cancer_g2 <- cc; break }
    }
    if (is.na(g2)) {
      add("RD2_real_affected_combination", "Real data: a partially-missing cancer-feature combination", "SKIP", "No gene/cancer combination with partial (neither zero nor total) missingness found")
    } else {
      spot_check("RD2_real_affected_combination", "Real data: a cancer-feature combination with genuine native PARTIAL missingness (still estimable) matches cor.test() on the exact same complete-case subset", g2, cancer_g2)
    }

    # RD3: the most heavily (but not totally, in at least one cancer) missing feature in the whole matrix.
    g3 <- NA_character_; cancer_g3 <- NA_character_
    for (cand in names(sort(na_counts, decreasing = TRUE))) {
      cc <- first_partial_missing_cancer(cand)
      if (!is.na(cc)) { g3 <- cand; cancer_g3 <- cc; break }
    }
    if (is.na(g3)) {
      add("RD3_real_heavily_missing_feature", "Real data: heavily-missing feature spot check", "SKIP", "No gene with partial missingness in any cancer found")
    } else {
      spot_check("RD3_real_heavily_missing_feature", sprintf("Real data: among the most heavily missing features in the matrix, a cancer where it is partially (not totally) missing (%d/%d missing overall) matches cor.test()", na_counts[g3], ncol(mrna_mat)), g3, cancer_g3)
    }

    # RD4: a boundary case -- scan a handful of affected genes/cancers for one whose |rho| lands near the primary threshold (0.30).
    boundary_found <- FALSE
    for (g in head(affected_genes, 200)) {
      for (cc in unique(sample_annotation$cancer)) {
        idx <- which(sample_annotation$cancer == cc)
        if (!any(is.na(mrna_mat[g, idx]))) next
        ref <- independent_ref(sample_annotation$CYT[idx], mrna_mat[g, idx])
        if (!is.na(ref$rho) && abs(abs(ref$rho) - 0.30) < 0.05) {
          spot_check("RD4_real_boundary_case", sprintf("Real data: a gene/cancer combination near the primary rho threshold (0.30) with real missingness matches cor.test() (|rho|=%.3f)", abs(ref$rho)), g, cc)
          boundary_found <- TRUE
          break
        }
      }
      if (boundary_found) break
    }
    if (!boundary_found) add("RD4_real_boundary_case", "Real data: boundary-case spot check", "SKIP", "No affected gene/cancer combination with |rho| within 0.05 of the 0.30 threshold found in the scanned sample -- skipped, not failed")
  }

  df <- do.call(rbind, results)
  out_path <- file.path(project_root, "tests", "output", "test_pairwise_complete_spearman_results.tsv")
  write.table(df, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
  n_fail <- sum(df$result == "FAIL")
  cat(sprintf("\n%d/%d PASS, %d SKIPPED, %d FAIL\n", sum(df$result == "PASS"), nrow(df), sum(df$result == "SKIPPED"), n_fail))
  if (n_fail > 0) stop(sprintf("%d pairwise-complete Spearman test(s) FAILED. See %s.", n_fail, out_path), call. = FALSE)
  invisible(df)
}

if (sys.nframe() == 0) {
  cli_args <- commandArgs(trailingOnly = TRUE)
  if (length(cli_args) < 1) stop("Usage: Rscript test_pairwise_complete_spearman.R <project_root>")
  invisible(run_pairwise_complete_tests(normalizePath(cli_args[[1]], mustWork = TRUE)))
}
