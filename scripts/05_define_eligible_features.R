# ============================================================================
# 05_define_eligible_features.R
#
# Purpose: determine which cancer types are eligible for (a) unadjusted
#          per-cancer discovery, (b) CPE-adjusted discovery, and (c)
#          LUMP-adjusted discovery, and attach tumor-purity (CPE, LUMP)
#          values to every matched sample.
# Inputs:  <run_root>/processed/stage03/sample_annotation.tsv (same-run);
#          <run_root>/processed/stage04/cyt_scores.tsv (same-run);
#          TCGAbiolinks::Tumor.purity (immutable external reference,
#          bundled inside the installed TCGAbiolinks package)
# Outputs: <run_root>/processed/stage05/sample_annotation_with_purity.tsv,
#          eligible_cancers_unadjusted.tsv, eligible_cancers_cpe.tsv,
#          eligible_cancers_lump.tsv
# Method:  a cancer type is eligible for unadjusted analysis if it has at
#          least analysis_config$discovery$min_samples_per_cancer matched
#          samples, and additionally eligible for CPE-adjusted (or
#          LUMP-adjusted) analysis only if that many samples also have
#          non-missing CPE (or LUMP) and CYT.
#          Purity values are matched to samples by the 15-character TCGA
#          barcode (Project-TSS-Participant-SampleType, e.g.
#          "TCGA-BL-A0C8-01"; see R/tcga_identifiers.R::
#          tcga_sample_barcode_15()), the granularity of
#          TCGAbiolinks::Tumor.purity. Because this granularity excludes the
#          vial letter, the purity table contains more than one row for some
#          identifiers (different vials of one participant and sample type);
#          these are detected, reported and resolved by the deterministic
#          rule described under "Duplicate-barcode handling" below.
#          TCGAbiolinks::Tumor.purity$CPE and $LUMP are factors with a comma
#          decimal separator and are parsed through their character labels
#          (R/tcga_identifiers.R::parse_purity_value()).
# Depends on: 03_harmonize_identifiers.R, 04_compute_validate_cyt.R
# ============================================================================

run_stage05 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "05_define_eligible_features")
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "tcga_identifiers.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "sample_exclusions.R"), local = TRUE)
  suppressPackageStartupMessages(library(data.table))
  script_path <- file.path(ctx$data_root, "scripts", "05_define_eligible_features.R")
  exclusions_path <- file.path(ctx$data_root, "config", "sample_exclusions.tsv")
  exclusions <- load_sample_exclusions(exclusions_path)

  sample_annotation <- read_run_artifact(ctx, file.path("processed", "stage03", "sample_annotation.tsv"), "03_harmonize_identifiers",
                                          function(p) read.delim(p, stringsAsFactors = FALSE))
  cyt_scores <- read_run_artifact(ctx, file.path("processed", "stage04", "cyt_scores.tsv"), "04_compute_validate_cyt",
                                   function(p) read.delim(p, stringsAsFactors = FALSE))
  input_paths <- file.path(ctx$run_root, c(file.path("processed", "stage03", "sample_annotation.tsv"),
                                            file.path("processed", "stage04", "cyt_scores.tsv")))

  run_stage_assert(ctx, log_fn, identical(sample_annotation$sample_id, cyt_scores$sample_id) ||
                      setequal(sample_annotation$sample_id, cyt_scores$sample_id),
                    "Sample sets in stage03 and stage04 outputs match")

  sample_annotation <- merge(sample_annotation, cyt_scores[, c("sample_id", "CYT")], by = "sample_id")
  sample_annotation$sample_id_15 <- tcga_sample_barcode_15(sample_annotation$sample_id)

  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    run_stage_fail(ctx, log_fn, "TCGAbiolinks is required to obtain CPE/LUMP tumor-purity estimates but is not installed.")
  }
  data("Tumor.purity", package = "TCGAbiolinks", envir = environment())
  purity <- as.data.frame(Tumor.purity, stringsAsFactors = FALSE)
  purity$Sample.ID <- standardize_tcga_id(purity$Sample.ID)
  # sample_id_15 is not unique in TCGAbiolinks::Tumor.purity: the same
  # participant and sample type can have more than one vial (distinguished
  # by the letter at position 16). In the installed table, 67 sample_id_15
  # groups among 9,364 rows are duplicated, and about half of them carry
  # differing CPE/LUMP values between vials. Deduplication to the
  # 15-character identifier is therefore performed by the explicit rule
  # below.
  purity$sample_id_15 <- tcga_sample_barcode_15(purity$Sample.ID)
  # TCGAbiolinks::Tumor.purity$CPE and $LUMP are factors whose labels use a
  # comma decimal separator (e.g. "0,9246"). as.numeric() applied directly to
  # a factor returns the internal level code, so values are parsed with
  # parse_purity_value() (R/tcga_identifiers.R), which converts through the
  # character labels and normalizes the separator (tests in
  # tests/test_purity_value_parsing.R).
  raw_factors <- list(CPE = purity[["CPE"]], LUMP = purity[["LUMP"]])
  for (nm in c("CPE", "LUMP")) purity[[nm]] <- parse_purity_value(purity[[nm]])

  # Integrity check against level-code leakage, out-of-range values and
  # value equality with the raw factor's level codes. It must run before
  # any deduplication or exclusion, while `purity` and `raw_factors` are
  # still full-length and positionally aligned.
  for (nm in c("CPE", "LUMP")) {
    corruption <- check_purity_value_corruption(purity[[nm]], raw_factor = raw_factors[[nm]])
    run_stage_assert(ctx, log_fn, nrow(corruption) == 0,
                      "%s purity values show no sign of factor-level-code corruption, out-of-range values, or level-order rank equivalence%s",
                      nm, if (nrow(corruption) > 0) paste0(" -- ", paste(sprintf("[%s] %s", corruption$check, corruption$detail), collapse = "; ")) else "")
  }

  # --- Prespecified sample exclusions, applied after the integrity check
  #     above (which requires full-length, positionally aligned vectors).
  #     Excluded samples are already absent from `sample_annotation`
  #     because Stage 03 removed them before matching; filtering `purity`
  #     here and asserting the final result below makes that guarantee
  #     explicit. -----------------------------------------------------------
  excl_purity <- drop_excluded_sample_rows(purity, "Sample.ID", exclusions, "purity",
                                            applicable_scopes = c("all_layers", "purity_only"))
  n_purity_before_excl <- nrow(purity)
  purity <- excl_purity$df
  if (excl_purity$n_removed > 0) {
    log_fn("Prespecified exclusion: removed %d purity-table row(s) (%s); %d -> %d rows.",
           excl_purity$n_removed, paste(excl_purity$removed_ids, collapse = ", "), n_purity_before_excl, nrow(purity))
  }

  # ==========================================================================
  # Duplicate-barcode handling: report the count, check agreement between
  # vials, and apply a deterministic, documented resolution rule.
  # ==========================================================================
  dup_ids <- unique(purity$sample_id_15[duplicated(purity$sample_id_15)])
  n_dup_groups <- length(dup_ids)
  conflict_rows <- list()
  n_conflict_cpe <- 0L; n_conflict_lump <- 0L
  for (id in dup_ids) {
    sub <- purity[purity$sample_id_15 == id, , drop = FALSE]
    cpe_vals <- unique(sub$CPE[!is.na(sub$CPE)])
    lump_vals <- unique(sub$LUMP[!is.na(sub$LUMP)])
    cpe_conflict <- length(cpe_vals) > 1
    lump_conflict <- length(lump_vals) > 1
    if (cpe_conflict) n_conflict_cpe <- n_conflict_cpe + 1L
    if (lump_conflict) n_conflict_lump <- n_conflict_lump + 1L
    if (cpe_conflict || lump_conflict) {
      conflict_rows[[length(conflict_rows) + 1]] <- data.frame(
        sample_id_15 = id, n_rows = nrow(sub), sample_ids = paste(sub$Sample.ID, collapse = ";"),
        cpe_values = paste(sub$CPE, collapse = ";"), lump_values = paste(sub$LUMP, collapse = ";"),
        cpe_conflict = cpe_conflict, lump_conflict = lump_conflict, stringsAsFactors = FALSE
      )
    }
  }
  conflict_dt <- if (length(conflict_rows) > 0) do.call(rbind, conflict_rows) else
    data.frame(sample_id_15 = character(0), n_rows = integer(0), sample_ids = character(0),
               cpe_values = character(0), lump_values = character(0), cpe_conflict = logical(0), lump_conflict = logical(0))
  write_run_artifact(ctx, conflict_dt, file.path("tables", "audit", "purity_duplicate_barcode_conflicts.tsv"),
                      "05_define_eligible_features", script_path, input_paths,
                      function(x, p) write.table(x, p, sep = "\t", row.names = FALSE, quote = FALSE))
  log_fn("Duplicate sample_id_15 groups: %d (of %d total rows, %d unique). CPE conflicts: %d. LUMP conflicts: %d. Full detail: tables/audit/purity_duplicate_barcode_conflicts.tsv.",
         n_dup_groups, nrow(purity), length(unique(purity$sample_id_15)), n_conflict_cpe, n_conflict_lump)

  # Per-metric, per-group resolution without vial preference.
  #
  # Rule (prespecified primary analysis):
  #   - For a sample_id_15 whose CPE values conflict across vials (more
  #     than one distinct non-missing value), CPE is set to missing (NA);
  #     values are neither averaged nor selected by vial. LUMP is resolved
  #     independently by the same rule, so a sample can have a resolvable
  #     CPE and an unresolvable LUMP, or vice versa.
  #   - For a sample_id_15 without a conflict for a given metric (values
  #     agree, or only one vial has a non-missing value), the single
  #     available value is kept.
  #   - A sample whose CPE (or LUMP) is set to missing is excluded only
  #     from analyses that require that purity metric as a covariate
  #     (Stage 10 filters on `!is.na(sample_annotation[[metric_name]])`);
  #     its expression, clinical, survival and CYT data remain available to
  #     every unadjusted analysis.
  #   - Stage 11 computes a sensitivity analysis in which these samples are
  #     excluded from every layer and compares it with this primary
  #     analysis; it does not replace the primary analysis.
  # The rule is implemented once in
  # R/tcga_identifiers.R::resolve_purity_duplicate_value() and unit-tested
  # in tests/test_winsorization_purity_conflicts_and_path_containment.R.
  purity_dt <- as.data.table(purity)[, c("sample_id_15", "CPE", "LUMP")]
  conflict_cpe_ids <- if (nrow(conflict_dt) > 0) conflict_dt$sample_id_15[conflict_dt$cpe_conflict] else character(0)
  conflict_lump_ids <- if (nrow(conflict_dt) > 0) conflict_dt$sample_id_15[conflict_dt$lump_conflict] else character(0)
  purity <- as.data.frame(purity_dt[, .(
    CPE = resolve_purity_duplicate_value(CPE, any(sample_id_15 %in% conflict_cpe_ids)),
    LUMP = resolve_purity_duplicate_value(LUMP, any(sample_id_15 %in% conflict_lump_ids))
  ), by = sample_id_15])
  n_cpe_forced_missing <- length(conflict_cpe_ids)
  n_lump_forced_missing <- length(conflict_lump_ids)
  log_fn("Purity conflict resolution (Option 1, no arbitrary vial tie-break): %d sample(s) had CPE forced to explicitly missing (unresolvable vial conflict); %d sample(s) had LUMP forced to explicitly missing. Affected identifiers: tables/audit/purity_duplicate_barcode_conflicts.tsv. These samples remain fully available to every analysis that does not require the affected purity metric as a covariate.",
         n_cpe_forced_missing, n_lump_forced_missing)

  sample_annotation <- merge(sample_annotation, purity, by = "sample_id_15", all.x = TRUE)
  log_fn("Samples with non-missing CPE: %d / %d. Samples with non-missing LUMP: %d / %d.",
         sum(!is.na(sample_annotation$CPE)), nrow(sample_annotation),
         sum(!is.na(sample_annotation$LUMP)), nrow(sample_annotation))

  # Flag, carried in sample_annotation_with_purity.tsv, marking samples in
  # the matched cohort whose CPE and/or LUMP was set to missing by the
  # conflict rule, so that downstream stages (including the Stage 11
  # full-exclusion sensitivity analysis) can identify them directly.
  sample_annotation$purity_cpe_forced_missing_conflict <- sample_annotation$sample_id_15 %in% conflict_cpe_ids
  sample_annotation$purity_lump_forced_missing_conflict <- sample_annotation$sample_id_15 %in% conflict_lump_ids

  # Per-sample detail table (identifier, cancer, conflicting values,
  # resolution): joins the duplicate-barcode conflict census with the
  # matched-cohort cancer assignment, restricted to samples in the matched
  # cohort (conflict_dt is computed from the raw purity table and has no
  # cancer column).
  cohort_conflict_ids <- unique(sample_annotation$sample_id_15[
    sample_annotation$purity_cpe_forced_missing_conflict | sample_annotation$purity_lump_forced_missing_conflict])
  conflict_detail <- merge(
    conflict_dt[conflict_dt$sample_id_15 %in% cohort_conflict_ids, , drop = FALSE],
    unique(sample_annotation[sample_annotation$sample_id_15 %in% cohort_conflict_ids,
                              c("sample_id_15", "sample_id", "cancer")]),
    by = "sample_id_15", all.x = TRUE)
  write_run_artifact(ctx, conflict_detail, file.path("tables", "audit", "purity_conflict_resolution_matched_cohort.tsv"),
                      "05_define_eligible_features", script_path, input_paths,
                      function(x, p) write.table(x, p, sep = "\t", row.names = FALSE, quote = FALSE))
  log_fn("Purity conflict resolution, matched-cohort detail (%d sample(s), spanning %d cancer type(s)): tables/audit/purity_conflict_resolution_matched_cohort.tsv.",
         length(cohort_conflict_ids), length(unique(conflict_detail$cancer)))

  # Verifies that no prespecified-excluded sample_id_15 is present in the
  # final merged sample_annotation, whether it was removed by Stage 03 or
  # by the purity-table filter above.
  excl_scope_all_or_purity <- exclusions$sample_id_15[exclusions$scope %in% c("all_layers", "purity_only")]
  run_stage_assert(ctx, log_fn, !any(sample_annotation$sample_id_15 %in% excl_scope_all_or_purity),
                    "No prespecified-excluded sample is present in the final merged sample_annotation_with_purity table (%d exclusion(s) checked)",
                    length(excl_scope_all_or_purity))
  # A purity fraction is always in [0, 1]; a factor level code would not be.
  run_stage_assert(ctx, log_fn,
    all(sample_annotation$CPE[!is.na(sample_annotation$CPE)] >= 0 & sample_annotation$CPE[!is.na(sample_annotation$CPE)] <= 1) &&
    all(sample_annotation$LUMP[!is.na(sample_annotation$LUMP)] >= 0 & sample_annotation$LUMP[!is.na(sample_annotation$LUMP)] <= 1),
    "Every non-missing CPE and LUMP value lies in the valid [0,1] purity-fraction range")

  min_n <- analysis_config$discovery$min_samples_per_cancer
  cancer_counts <- as.data.table(sample_annotation)[, .N, by = cancer]
  eligible_unadjusted <- sort(cancer_counts[N >= min_n, cancer])

  cpe_counts <- as.data.table(sample_annotation)[!is.na(CPE) & !is.na(CYT), .N, by = cancer]
  eligible_cpe <- sort(cpe_counts[N >= min_n, cancer])
  lump_counts <- as.data.table(sample_annotation)[!is.na(LUMP) & !is.na(CYT), .N, by = cancer]
  eligible_lump <- sort(lump_counts[N >= min_n, cancer])

  run_stage_assert(ctx, log_fn, length(eligible_unadjusted) >= 1, "At least one cancer type is eligible for unadjusted analysis")
  run_stage_assert(ctx, log_fn, all(eligible_cpe %in% eligible_unadjusted), "Every CPE-eligible cancer is also unadjusted-eligible")
  log_fn("Eligible cancer types: unadjusted=%d, CPE-adjusted=%d, LUMP-adjusted=%d",
         length(eligible_unadjusted), length(eligible_cpe), length(eligible_lump))

  write_fn <- function(x, p) write.table(x, p, sep = "\t", row.names = FALSE, quote = FALSE)
  write_run_artifact(ctx, sample_annotation, file.path("processed", "stage05", "sample_annotation_with_purity.tsv"),
                      "05_define_eligible_features", script_path, input_paths, write_fn)
  write_run_artifact(ctx, data.frame(cancer = eligible_unadjusted), file.path("processed", "stage05", "eligible_cancers_unadjusted.tsv"),
                      "05_define_eligible_features", script_path, input_paths, write_fn)
  write_run_artifact(ctx, data.frame(cancer = eligible_cpe), file.path("processed", "stage05", "eligible_cancers_cpe.tsv"),
                      "05_define_eligible_features", script_path, input_paths, write_fn)
  write_run_artifact(ctx, data.frame(cancer = eligible_lump), file.path("processed", "stage05", "eligible_cancers_lump.tsv"),
                      "05_define_eligible_features", script_path, input_paths, write_fn)

  # --- Complete the sample-exclusion application summary (Stage 03 wrote
  #     the miRNA/mRNA/clinical rows; this appends the purity row and
  #     re-writes the same path under this stage's attribution) -----------
  prior_summary <- read_run_artifact(ctx, file.path("tables", "audit", "sample_exclusion_expression_summary.tsv"),
                                      "03_harmonize_identifiers", function(p) read.delim(p, stringsAsFactors = FALSE))
  purity_summary_row <- if (excl_purity$n_removed == 0) {
    data.frame(layer = "purity", n_matched_records = 0L, removed_identifiers = "",
               n_remaining = n_purity_before_excl, reason = "(no exclusion matched this layer)", stringsAsFactors = FALSE)
  } else {
    data.frame(layer = "purity", n_matched_records = excl_purity$n_removed,
               removed_identifiers = paste(excl_purity$removed_ids, collapse = "; "),
               n_remaining = n_purity_before_excl - excl_purity$n_removed,
               reason = paste(unique(excl_purity$audit$reason), collapse = " | "), stringsAsFactors = FALSE)
  }
  full_summary <- rbind(prior_summary, purity_summary_row)
  write_run_artifact(ctx, full_summary, file.path("tables", "audit", "sample_exclusion_application_summary.tsv"),
                      "05_define_eligible_features", script_path, c(exclusions_path, input_paths), write_fn)
  log_fn("Sample-exclusion application summary completed with the purity layer (%d total layer row(s)).", nrow(full_summary))

  record_stage_status(ctx, "05_define_eligible_features", "COMPLETE")
  log_fn("Stage complete.")
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=5 or source scripts/run.R in RStudio.")
}
