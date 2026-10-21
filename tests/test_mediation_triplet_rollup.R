# ============================================================================
# test_mediation_triplet_rollup.R
#
# Known-answer, boundary, missing-value, duplicated-key, and direction-
# consistency tests for
# scripts/audit/mediation_triplet_rollup.R::compute_mediation_triplet_rollup(),
# the deterministic aggregation written to resolve Blocker 2 of
# reports/MANUSCRIPT_AND_SUPPLEMENTARY_UPDATE_ACTION_PLAN.md. Calls the
# REAL production function against constructed synthetic data with a KNOWN
# hand-computed answer -- never a grep-based source-text check.
#
# SYNTHETIC DATA ONLY.
#
# Usage: Rscript tests/test_mediation_triplet_rollup.R <project_root>
# ============================================================================

run_mediation_triplet_rollup_tests <- function(project_root) {
  suppressPackageStartupMessages(library(data.table))
  source(file.path(project_root, "scripts", "lib", "mediation_triplet_rollup.R"), local = TRUE)

  results <- list()
  add <- function(id, desc, pass, detail) {
    results[[length(results) + 1]] <<- data.frame(test_id = id, description = desc,
      result = if (identical(pass,"SKIP")) "SKIP" else if (isTRUE(pass)) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE)
    cat(sprintf("[%s] %s -- %s\n", if (identical(pass,"SKIP")) "SKIP" else if (isTRUE(pass)) "PASS" else "FAIL", id, detail))
  }

  mk_row <- function(mirna_id, gene, cancer, raw_bootstrap_p, bh_fdr, prop_mediated,
                      prop_mediated_meaningful, point_indirect) {
    data.table(mirna_id = mirna_id, gene = gene, cancer = cancer,
               raw_bootstrap_p = raw_bootstrap_p, bh_fdr = bh_fdr, prop_mediated = prop_mediated,
               prop_mediated_meaningful = prop_mediated_meaningful, point_indirect = point_indirect)
  }

  # ==========================================================================
  # K1: hand-computed known answer, one triplet, 4 contexts, mixed
  # significance -- every field checked against a value worked out by hand.
  #   raw_bootstrap_p: 0.001, 0.02, 0.04, 0.20  -> nominal (<0.05): 3/4 = 75%
  #   bh_fdr:          0.005, 0.03, 0.06, 0.40  -> <0.05: 2/4=50% (0.005,0.03 qualify; 0.06 does NOT); <0.01: 1/4=25%
  #   prop_mediated:   0.10, 0.30, 0.50, 0.90   -> min=0.10, median=0.40, max=0.90, range=0.80
  #   meaningful:      T, T, F, F               -> 2/4 = 50%
  #   point_indirect:  all positive             -> direction_consistent = TRUE
  # ==========================================================================
  dt1 <- rbind(
    mk_row("hsa-mir-X", "GENEA", "AAA", 0.001, 0.005, 0.10, TRUE,  0.5),
    mk_row("hsa-mir-X", "GENEA", "BBB", 0.02,  0.03,  0.30, TRUE,  0.2),
    mk_row("hsa-mir-X", "GENEA", "CCC", 0.04,  0.06,  0.50, FALSE, 0.1),
    mk_row("hsa-mir-X", "GENEA", "DDD", 0.20,  0.40,  0.90, FALSE, 0.05)
  )
  r1 <- compute_mediation_triplet_rollup(dt1)
  add("K1_n_contexts", "Known-answer: 4 contexts counted for the single triplet", r1$n_contexts == 4, sprintf("got %d, expected 4", r1$n_contexts))
  add("K1_nominal_sig", "Known-answer: 3/4=75% nominal significance (raw_bootstrap_p<0.05)", r1$n_nominal_sig == 3 && abs(r1$pct_nominal_sig - 75) < 1e-9, sprintf("n=%d pct=%.4f", r1$n_nominal_sig, r1$pct_nominal_sig))
  add("K1_fdr05", "Known-answer: 2/4=50% at BH-FDR<0.05 (0.06 does not qualify)", r1$n_fdr05 == 2 && abs(r1$pct_fdr05 - 50) < 1e-9, sprintf("n=%d pct=%.4f", r1$n_fdr05, r1$pct_fdr05))
  add("K1_fdr01", "Known-answer: 1/4=25% at BH-FDR<0.01", r1$n_fdr01 == 1 && abs(r1$pct_fdr01 - 25) < 1e-9, sprintf("n=%d pct=%.4f", r1$n_fdr01, r1$pct_fdr01))
  add("K1_min_median_fdr", "Known-answer: min_bh_fdr=0.005, median_bh_fdr=(0.03+0.06)/2=0.045", abs(r1$min_bh_fdr - 0.005) < 1e-9 && abs(r1$median_bh_fdr - 0.045) < 1e-9, sprintf("min=%.6f median=%.6f", r1$min_bh_fdr, r1$median_bh_fdr))
  add("K1_prop_mediated_stats", "Known-answer: min=0.10, median=0.40, max=0.90, range=0.80", abs(r1$min_prop_mediated-0.10)<1e-9 && abs(r1$median_prop_mediated-0.40)<1e-9 && abs(r1$max_prop_mediated-0.90)<1e-9 && abs(r1$range_prop_mediated-0.80)<1e-9,
      sprintf("min=%.4f median=%.4f max=%.4f range=%.4f", r1$min_prop_mediated, r1$median_prop_mediated, r1$max_prop_mediated, r1$range_prop_mediated))
  add("K1_meaningful", "Known-answer: 2/4=50% prop_mediated_meaningful", r1$n_prop_mediated_meaningful == 2 && abs(r1$pct_prop_mediated_meaningful - 50) < 1e-9, sprintf("n=%d pct=%.4f", r1$n_prop_mediated_meaningful, r1$pct_prop_mediated_meaningful))
  add("K1_direction_consistent_true", "Known-answer: all 4 point_indirect positive -> direction_consistent_across_contexts=TRUE", isTRUE(r1$direction_consistent_across_contexts) && r1$n_distinct_indirect_directions == 1, sprintf("consistent=%s n_distinct=%d", r1$direction_consistent_across_contexts, r1$n_distinct_indirect_directions))

  # ==========================================================================
  # K2: direction INCONSISTENCY -- known answer, one triplet, mixed-sign
  # point_indirect must be flagged direction_consistent_across_contexts=FALSE.
  # ==========================================================================
  dt2 <- rbind(
    mk_row("hsa-mir-Y", "GENEB", "AAA", 0.01, 0.02, 0.20, TRUE, 0.3),
    mk_row("hsa-mir-Y", "GENEB", "BBB", 0.01, 0.02, 0.20, TRUE, -0.3)
  )
  r2 <- compute_mediation_triplet_rollup(dt2)
  add("K2_direction_inconsistent", "Known-answer: one positive + one negative point_indirect -> direction_consistent_across_contexts=FALSE, n_distinct=2", !isTRUE(r2$direction_consistent_across_contexts) && r2$n_distinct_indirect_directions == 2, sprintf("consistent=%s n_distinct=%d", r2$direction_consistent_across_contexts, r2$n_distinct_indirect_directions))

  # ==========================================================================
  # K3: multiple triplets in one input -- confirms grouping is per-(mirna,gene),
  # not pooled, and confirms the `supporting_cancers` concatenation is exact.
  # ==========================================================================
  dt3 <- rbind(dt1, dt2,
    mk_row("hsa-mir-Z", "GENEC", "ZZZ", 0.5, 0.6, 1.0, FALSE, 0.01))
  r3 <- compute_mediation_triplet_rollup(dt3)
  add("K3_three_distinct_triplets", "Three distinct (mirna,gene) groups in the output for a 7-row, 3-triplet input", nrow(r3) == 3, sprintf("got %d rows (expected 3)", nrow(r3)))
  add("K3_group_row_counts_isolated", "Each group's n_contexts reflects ONLY its own rows, not pooled across triplets (hsa-mir-X=4, hsa-mir-Y=2, hsa-mir-Z=1)",
      r3[mirna_id == "hsa-mir-X"]$n_contexts == 4 && r3[mirna_id == "hsa-mir-Y"]$n_contexts == 2 && r3[mirna_id == "hsa-mir-Z"]$n_contexts == 1,
      paste(sprintf("%s=%d", r3$mirna_id, r3$n_contexts), collapse = ", "))
  add("K3_supporting_cancers_string", "supporting_cancers is the sorted, semicolon-joined list of that triplet's own cancer contexts",
      r3[mirna_id == "hsa-mir-X"]$supporting_cancers == "AAA;BBB;CCC;DDD",
      r3[mirna_id == "hsa-mir-X"]$supporting_cancers)

  # ==========================================================================
  # B1: single-context triplet (n=1) -- min==median==max, range=0, boundary case.
  # ==========================================================================
  dt_single <- mk_row("hsa-mir-S", "GENES", "ONLY", 0.03, 0.04, 0.5, TRUE, 0.2)
  r_single <- compute_mediation_triplet_rollup(dt_single)
  add("B1_single_context_boundary", "A triplet with exactly 1 supporting context: min=median=max=prop_mediated, range=0, pct fields are 0 or 100 (never NaN/NA)",
      r_single$n_contexts == 1 && r_single$min_prop_mediated == 0.5 && r_single$median_prop_mediated == 0.5 && r_single$max_prop_mediated == 0.5 && r_single$range_prop_mediated == 0 && r_single$pct_nominal_sig == 100,
      sprintf("n=%d min=%.2f median=%.2f max=%.2f range=%.2f pct_nominal=%.1f", r_single$n_contexts, r_single$min_prop_mediated, r_single$median_prop_mediated, r_single$max_prop_mediated, r_single$range_prop_mediated, r_single$pct_nominal_sig))

  # ==========================================================================
  # M1: missing value (NA) in a required column must FAIL LOUDLY, never be
  # silently dropped or silently coerced.
  # ==========================================================================
  dt_na <- rbind(dt1, mk_row("hsa-mir-X", "GENEA", "EEE", NA_real_, 0.02, 0.3, TRUE, 0.1))
  r_na <- tryCatch(compute_mediation_triplet_rollup(dt_na), error = function(e) e)
  add("M1_na_raises_loud_error", "A row with NA in a required column (raw_bootstrap_p) causes compute_mediation_triplet_rollup() to stop() with an informative message, not a silent drop",
      inherits(r_na, "error") && grepl("NA", conditionMessage(r_na)),
      if (inherits(r_na, "error")) substr(conditionMessage(r_na), 1, 150) else "NOT REJECTED (defect)")

  # ==========================================================================
  # D1: duplicated (mirna_id, gene, cancer) key must FAIL LOUDLY, never
  # silently double-count or arbitrarily pick one row.
  # ==========================================================================
  dt_dup <- rbind(dt1, dt1[1])
  r_dup <- tryCatch(compute_mediation_triplet_rollup(dt_dup), error = function(e) e)
  add("D1_duplicate_key_raises_loud_error", "A duplicated (mirna_id, gene, cancer) key causes an informative stop(), not silent double-counting",
      inherits(r_dup, "error") && grepl("duplicat", conditionMessage(r_dup), ignore.case = TRUE),
      if (inherits(r_dup, "error")) substr(conditionMessage(r_dup), 1, 150) else "NOT REJECTED (defect)")

  # ==========================================================================
  # E1: empty input must FAIL LOUDLY, never return a silent zero-row result
  # that could be mistaken for "zero triplets found."
  # ==========================================================================
  r_empty <- tryCatch(compute_mediation_triplet_rollup(dt1[0]), error = function(e) e)
  add("E1_empty_input_raises_loud_error", "A zero-row input causes an informative stop(), not a silently-empty result", inherits(r_empty, "error"),
      if (inherits(r_empty, "error")) substr(conditionMessage(r_empty), 1, 150) else "NOT REJECTED (defect)")

  # ==========================================================================
  # C1: missing required column must FAIL LOUDLY.
  # ==========================================================================
  dt_missing_col <- dt1[, !"bh_fdr", with = FALSE]
  r_missing_col <- tryCatch(compute_mediation_triplet_rollup(dt_missing_col), error = function(e) e)
  add("C1_missing_column_raises_loud_error", "An input missing a required column (bh_fdr) causes an informative stop() naming the missing column",
      inherits(r_missing_col, "error") && grepl("bh_fdr", conditionMessage(r_missing_col)),
      if (inherits(r_missing_col, "error")) substr(conditionMessage(r_missing_col), 1, 150) else "NOT REJECTED (defect)")

  # ==========================================================================
  # R1: real-data spot check (read-only): if a completed run with Stage 15's
  # own output exists, confirm the production aggregation's row count and
  # totals match an INDEPENDENTLY, freshly re-derived aggregation computed
  # here via a completely separate code path (base R aggregate(), not
  # data.table grouping) -- SKIPPED, not FAILED, if none is found.
  # ==========================================================================
  run_dirs <- if(isTRUE(getOption("mirna.tests.real_data",FALSE))) file.path(project_root,"results") else character()
  run_dirs <- run_dirs[order(run_dirs, decreasing = TRUE)]
  usable_run <- NULL
  for (rd in run_dirs) {
    p <- file.path(rd, "tables", "downstream", "mediation_primary_support_subset.tsv")
    if (file.exists(p)) { usable_run <- rd; break }
  }
  if (is.null(usable_run)) {
    add("R1_real_data_matches_independent_reaggregation", "Real-run spot check", "SKIP", "Optional real-data spot check not requested or completed result absent; synthetic checks are reported separately")
  } else {
    real_dt <- fread(file.path(usable_run, "tables", "downstream", "mediation_primary_support_subset.tsv"))
    prod <- compute_mediation_triplet_rollup(real_dt)

    # Independent re-derivation via base R aggregate(), deliberately NOT
    # reusing any of compute_mediation_triplet_rollup()'s own code.
    key <- paste(real_dt$mirna_id, real_dt$gene, sep = "\r")
    n_contexts_indep <- table(key)
    n_fdr01_indep <- tapply(real_dt$bh_fdr < 0.01, key, sum)
    total_triplets_indep <- length(unique(key))

    ok <- nrow(prod) == total_triplets_indep &&
      all(sort(prod$n_contexts) == sort(as.integer(n_contexts_indep))) &&
      sum(prod$n_fdr01) == sum(n_fdr01_indep)
    add("R1_real_data_matches_independent_reaggregation", "Production aggregation on the real run's own mediation_primary_support_subset.tsv matches an independently-coded re-aggregation (base R table()/tapply(), not data.table) for triplet count, per-triplet n_contexts, and total n_fdr01",
        ok, sprintf("production: %d triplets, sum(n_fdr01)=%d; independent: %d triplets, sum(n_fdr01)=%d",
                     nrow(prod), sum(prod$n_fdr01), total_triplets_indep, sum(n_fdr01_indep)))
  }

  df <- do.call(rbind, results)
  out_path <- file.path(project_root, "tests", "output", "test_mediation_triplet_rollup_results.tsv")
  write.table(df, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
  n_fail <- sum(df$result == "FAIL")
  n_skip <- sum(df$result == "SKIP")
  cat(sprintf("\n%d/%d PASS, %d FAIL, %d SKIP\n", sum(df$result == "PASS"), nrow(df), n_fail, n_skip))
  if (n_fail > 0) stop(sprintf("%d test(s) FAILED. See %s.", n_fail, out_path), call. = FALSE)
  invisible(df)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop("Usage: Rscript test_mediation_triplet_rollup.R <project_root>")
  invisible(run_mediation_triplet_rollup_tests(normalizePath(args[[1]], mustWork = TRUE)))
}
