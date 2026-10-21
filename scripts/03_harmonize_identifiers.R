# ============================================================================
# 03_harmonize_identifiers.R
#
# Purpose: apply prespecified sample-level exclusions
#          (config/sample_exclusions.tsv); standardize TCGA barcodes,
#          restrict to primary tumor samples, apply expression preprocessing
#          (declared-scale validation, duplicate and aliquot barcode
#          classification and resolution, zero-variance removal), and match
#          samples across the miRNA, mRNA and clinical layers by exact
#          barcode.
# Inputs:  <run_root>/processed/stage02/raw_mirna.rds, raw_mrna.rds,
#          raw_survival.rds (same-run, provenance-verified);
#          config/sample_exclusions.tsv (version-controlled configuration;
#          see R/sample_exclusions.R)
# Outputs: <run_root>/processed/stage03/mirna_matched.rds, mrna_matched.rds,
#          survival_matched.rds, sample_annotation.tsv,
#          normalization_qc.tsv; tables/audit/expression_duplicate_column_
#          resolutions.tsv, tables/audit/expression_aliquot_resolutions.tsv,
#          tables/audit/clinical_duplicate_row_resolutions.tsv,
#          tables/audit/prespecified_sample_exclusions.tsv,
#          tables/audit/sample_exclusion_application_summary.tsv (miRNA,
#          mRNA and clinical rows; Stage 05 appends the purity row)
# Method:  prespecified exclusions (loaded and validated by
#          R/sample_exclusions.R::load_sample_exclusions(); no sample
#          identity is hard-coded) are applied to the raw miRNA/mRNA
#          matrices and the raw clinical table before duplicate-barcode
#          classification, so an excluded sample never enters that logic.
#          Each source's declared scale (config/analysis_config.R,
#          $expression$mirna_scale and $expression$mrna_scale) is validated
#          against the observed range. Duplicate barcodes are classified as
#          technical_duplicate (resolved automatically) or
#          conflicting_record (the stage stops and requires a documented
#          decision), for expression columns and clinical rows alike.
#          Aliquot duplicates (one 15-character sample with more than one
#          physical vial or portion) are resolved by a deterministic rule.
#          Gene symbols appearing on more than one row are deduplicated
#          (scripts/lib/expression_preprocessing.R::deduplicate_matrix_rows(), highest
#          variance retained), so every row-name lookup by gene symbol is
#          unambiguous. Winsorization is not applied here; it is applied
#          downstream by the stages whose methods are value-sensitive
#          (Stages 10, 12, 15, 16, 17, 18).
# Depends on: 02_ingest_and_provenance.R
# ============================================================================

run_stage03 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "03_harmonize_identifiers")
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "tcga_identifiers.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "expression_preprocessing.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "sample_exclusions.R"), local = TRUE)
  suppressPackageStartupMessages(library(data.table))

  script_path <- file.path(ctx$data_root, "scripts", "03_harmonize_identifiers.R")
  in02 <- file.path("processed", "stage02")
  raw_mirna <- read_run_artifact(ctx, file.path(in02, "raw_mirna.rds"), "02_ingest_and_provenance", readRDS)
  raw_mrna <- read_run_artifact(ctx, file.path(in02, "raw_mrna.rds"), "02_ingest_and_provenance", readRDS)
  raw_survival <- read_run_artifact(ctx, file.path(in02, "raw_survival.rds"), "02_ingest_and_provenance", readRDS)
  input_paths_02 <- file.path(ctx$run_root, in02, c("raw_mirna.rds", "raw_mrna.rds", "raw_survival.rds"))
  fwrite_fn <- function(x, p) write.table(x, p, sep = "\t", row.names = FALSE, quote = FALSE)

  # --- Prespecified sample exclusions (config/sample_exclusions.tsv),
  #     applied before any duplicate-barcode classification ----------------
  exclusions_path <- file.path(ctx$data_root, "config", "sample_exclusions.tsv")
  exclusions <- load_sample_exclusions(exclusions_path)
  log_fn("Loaded %d prespecified sample exclusion(s) from config/sample_exclusions.tsv.", nrow(exclusions))
  if (nrow(exclusions) > 0) {
    write_run_artifact(ctx, exclusions, file.path("tables", "audit", "prespecified_sample_exclusions.tsv"),
                        "03_harmonize_identifiers", script_path, exclusions_path, fwrite_fn)
  }

  # --- Clinical/survival table: standardize, restrict to primary tumor,
  #     classify and resolve duplicate rows ---------------------------------
  colnames(raw_survival)[colnames(raw_survival) == "sample"] <- "sample_id"
  colnames(raw_survival)[colnames(raw_survival) == "cancer type abbreviation"] <- "cancer_abbrev"
  raw_survival$sample_id <- standardize_tcga_id(raw_survival$sample_id)
  raw_survival <- raw_survival[is_primary_tumor_barcode(raw_survival$sample_id), , drop = FALSE]

  n_clinical_before_excl <- nrow(raw_survival)
  excl_clinical <- drop_excluded_sample_rows(raw_survival, "sample_id", exclusions, "clinical",
                                              applicable_scopes = c("all_layers", "clinical_only"))
  raw_survival <- excl_clinical$df
  if (excl_clinical$n_removed > 0) {
    log_fn("Prespecified exclusion: removed %d clinical/survival row(s) (%s); %d -> %d rows.",
           excl_clinical$n_removed, paste(excl_clinical$removed_ids, collapse = ", "), n_clinical_before_excl, nrow(raw_survival))
  }

  clinical_dup_report <- classify_duplicate_rows(raw_survival, "sample_id")
  clinical_conflicting <- clinical_dup_report[clinical_dup_report$classification == "conflicting_record", , drop = FALSE]
  run_stage_assert(ctx, log_fn, nrow(clinical_conflicting) == 0,
    "No conflicting duplicate clinical/survival rows (0 unresolved out of %d duplicate identifier group(s) found; see tables/audit/clinical_duplicate_row_resolutions.tsv). If this fails, a human must inspect the conflicting rows and decide which record is authoritative before this stage can proceed -- it is never auto-resolved.",
    nrow(clinical_dup_report))
  if (nrow(clinical_dup_report) > 0) {
    raw_survival <- resolve_technical_duplicate_rows(raw_survival, "sample_id", clinical_dup_report)
    write_run_artifact(ctx, clinical_dup_report, file.path("tables", "audit", "clinical_duplicate_row_resolutions.tsv"),
                        "03_harmonize_identifiers", script_path, input_paths_02[3], fwrite_fn)
    log_fn("Clinical/survival table: %d duplicate sample_id group(s) found, all technical_duplicate (identical across all fields), resolved to one row each; see tables/audit/clinical_duplicate_row_resolutions.tsv.", nrow(clinical_dup_report))
  }
  log_fn("Survival table after primary-tumor restriction and duplicate-row resolution: %d samples", nrow(raw_survival))

  # --- Expression matrices: prespecified sample exclusions, then
  #     declared-scale validation, duplicate/aliquot barcode classification
  #     and resolution, zero-variance removal -------------------------------
  n_mirna_cols_before_excl <- ncol(raw_mirna) - 1
  excl_mirna <- drop_excluded_sample_columns(raw_mirna, exclusions, "miRNA",
                                              applicable_scopes = c("all_layers", "expression_only"))
  raw_mirna <- excl_mirna$df
  if (excl_mirna$n_removed > 0) {
    log_fn("Prespecified exclusion: removed %d miRNA-matrix column(s) (%s); %d -> %d sample columns.",
           excl_mirna$n_removed, paste(excl_mirna$removed_names, collapse = ", "), n_mirna_cols_before_excl, ncol(raw_mirna) - 1)
  }

  n_mrna_cols_before_excl <- ncol(raw_mrna) - 1
  excl_mrna <- drop_excluded_sample_columns(raw_mrna, exclusions, "mRNA",
                                             applicable_scopes = c("all_layers", "expression_only"))
  raw_mrna <- excl_mrna$df
  if (excl_mrna$n_removed > 0) {
    log_fn("Prespecified exclusion: removed %d mRNA-matrix column(s) (%s); %d -> %d sample columns.",
           excl_mrna$n_removed, paste(excl_mrna$removed_names, collapse = ", "), n_mrna_cols_before_excl, ncol(raw_mrna) - 1)
  }

  cfg_expr <- analysis_config$expression
  log_fn("Preprocessing miRNA matrix (declared scale: %s) ...", cfg_expr$mirna_scale)
  mirna_res <- preprocess_expression_matrix(raw_mirna, "mirna_id",
                                             declared_scale = cfg_expr$mirna_scale, pseudo_count = cfg_expr$pseudo_count,
                                             scale_sanity_max = cfg_expr$scale_sanity_max, remove_zero_variance = TRUE)
  run_stage_assert(ctx, log_fn, isTRUE(mirna_res$ok),
    "miRNA matrix preprocessing completed without an unresolved conflict (%s)",
    if (!isTRUE(mirna_res$ok)) mirna_res$message else "ok")
  run_stage_assert(ctx, log_fn, isTRUE(mirna_res$scale_check$consistent),
    "miRNA declared scale is consistent with the observed empirical range: %s", mirna_res$scale_check$detail)
  run_stage_assert(ctx, log_fn, mirna_res$n_infinite == 0 && mirna_res$n_nan == 0,
    "No infinite or NaN values in the miRNA matrix after declared-scale transformation (n_infinite=%d, n_nan=%d)", mirna_res$n_infinite, mirna_res$n_nan)
  mirna <- mirna_res$data
  log_fn("miRNA: transformed_by_pipeline=%s, %d features x %d samples; %d value(s) newly NA from non-numeric coercion; %d zero-variance feature(s) removed.",
         mirna_res$transformed_by_pipeline, nrow(mirna), ncol(mirna) - 1, mirna_res$n_newly_na_from_coercion, mirna_res$n_zero_variance_removed)

  log_fn("Preprocessing mRNA matrix (declared scale: %s) ...", cfg_expr$mrna_scale)
  mrna_res <- preprocess_expression_matrix(raw_mrna, "gene_id_raw",
                                            declared_scale = cfg_expr$mrna_scale, pseudo_count = cfg_expr$pseudo_count,
                                            scale_sanity_max = cfg_expr$scale_sanity_max, remove_zero_variance = TRUE)
  run_stage_assert(ctx, log_fn, isTRUE(mrna_res$ok),
    "mRNA matrix preprocessing completed without an unresolved conflict (%s)",
    if (!isTRUE(mrna_res$ok)) mrna_res$message else "ok")
  run_stage_assert(ctx, log_fn, isTRUE(mrna_res$scale_check$consistent),
    "mRNA declared scale is consistent with the observed empirical range: %s", mrna_res$scale_check$detail)
  run_stage_assert(ctx, log_fn, mrna_res$n_infinite == 0 && mrna_res$n_nan == 0,
    "No infinite or NaN values in the mRNA matrix after declared-scale transformation (n_infinite=%d, n_nan=%d)", mrna_res$n_infinite, mrna_res$n_nan)
  mrna <- mrna_res$data
  log_fn("mRNA: transformed_by_pipeline=%s, %d features x %d samples (before gene-symbol dedup); %d value(s) newly NA from non-numeric coercion; %d zero-variance feature(s) removed.",
         mrna_res$transformed_by_pipeline, nrow(mrna), ncol(mrna) - 1, mrna_res$n_newly_na_from_coercion, mrna_res$n_zero_variance_removed)

  dup_col_report <- rbind(
    if (nrow(mirna_res$duplicate_column_report)) cbind(layer = "miRNA", mirna_res$duplicate_column_report) else NULL,
    if (nrow(mrna_res$duplicate_column_report)) cbind(layer = "mRNA", mrna_res$duplicate_column_report) else NULL
  )
  if (!is.null(dup_col_report) && nrow(dup_col_report) > 0) {
    write_run_artifact(ctx, dup_col_report, file.path("tables", "audit", "expression_duplicate_column_resolutions.tsv"),
                        "03_harmonize_identifiers", script_path, input_paths_02[1:2], fwrite_fn)
    log_fn("%d exact-duplicate standardized-barcode column group(s) found across both matrices (all technical_duplicate; any conflicting_record would have already stopped this stage above); see tables/audit/expression_duplicate_column_resolutions.tsv.", nrow(dup_col_report))
  }
  aliquot_report <- rbind(
    if (nrow(mirna_res$aliquot_report)) cbind(layer = "miRNA", mirna_res$aliquot_report) else NULL,
    if (nrow(mrna_res$aliquot_report)) cbind(layer = "mRNA", mrna_res$aliquot_report) else NULL
  )
  if (!is.null(aliquot_report) && nrow(aliquot_report) > 0) {
    write_run_artifact(ctx, aliquot_report, file.path("tables", "audit", "expression_aliquot_resolutions.tsv"),
                        "03_harmonize_identifiers", script_path, input_paths_02[1:2], fwrite_fn)
    log_fn("%d aliquot-duplicate 15-char-sample group(s) found (>1 physical vial/portion for the same participant+sample-type); resolved to the most-complete aliquot per group (fewest missing values, alphabetical tie-break); see tables/audit/expression_aliquot_resolutions.tsv.", nrow(aliquot_report))
  }

  n_before_dedup <- nrow(mrna)
  mrna_symbols <- toupper(trimws(mrna$gene_id_raw))
  dup_symbols <- unique(mrna_symbols[duplicated(mrna_symbols)])
  if (length(dup_symbols) > 0) {
    mrna_mat_tmp <- as.matrix(mrna[, -1, drop = FALSE])
    keep <- !logical(nrow(mrna))
    for (s in dup_symbols) {
      idx <- which(mrna_symbols == s)
      keep[setdiff(idx, highest_variance_row(mrna_mat_tmp, idx))] <- FALSE
    }
    mrna <- mrna[keep, , drop = FALSE]
    log_fn("mRNA: %d duplicate gene symbol(s) resolved by highest-variance row (%s); %d -> %d rows",
           length(dup_symbols), paste(dup_symbols, collapse = ","), n_before_dedup, nrow(mrna))
  }
  run_stage_assert(ctx, log_fn, !any(duplicated(toupper(trimws(mrna$gene_id_raw)))),
                    "No duplicate gene symbols remain in the mRNA matrix after deduplication")

  common_ids <- Reduce(intersect, list(
    raw_survival$sample_id,
    setdiff(colnames(mirna), "mirna_id"),
    setdiff(colnames(mrna), "gene_id_raw")
  ))
  common_ids <- sort(unique(common_ids))
  run_stage_assert(ctx, log_fn, length(common_ids) > 0, "At least one sample is common to all three layers")
  log_fn("Samples common to miRNA, mRNA, and survival layers: %d", length(common_ids))

  survival_matched <- raw_survival[match(common_ids, raw_survival$sample_id), , drop = FALSE]
  mirna_matched <- mirna[, c("mirna_id", common_ids), drop = FALSE]
  mrna_matched <- mrna[, c("gene_id_raw", common_ids), drop = FALSE]

  run_stage_assert(ctx, log_fn, identical(colnames(mirna_matched)[-1], common_ids), "miRNA matrix column order matches common_ids")
  run_stage_assert(ctx, log_fn, identical(colnames(mrna_matched)[-1], common_ids), "mRNA matrix column order matches common_ids")
  run_stage_assert(ctx, log_fn, identical(survival_matched$sample_id, common_ids), "Survival table row order matches common_ids")
  run_stage_assert(ctx, log_fn, !any(duplicated(common_ids)), "No duplicate sample identifiers in the matched cohort")

  write_run_artifact(ctx, mirna_matched, file.path("processed", "stage03", "mirna_matched.rds"),
                      "03_harmonize_identifiers", script_path, input_paths_02, saveRDS)
  write_run_artifact(ctx, mrna_matched, file.path("processed", "stage03", "mrna_matched.rds"),
                      "03_harmonize_identifiers", script_path, input_paths_02, saveRDS)
  write_run_artifact(ctx, survival_matched, file.path("processed", "stage03", "survival_matched.rds"),
                      "03_harmonize_identifiers", script_path, input_paths_02, saveRDS)

  cancer_final <- ifelse(!is.na(survival_matched$cancer_abbrev) & survival_matched$cancer_abbrev != "",
                          survival_matched$cancer_abbrev, "UNKNOWN")
  sample_annotation <- data.frame(sample_id = common_ids, cancer = cancer_final, stringsAsFactors = FALSE)
  write_run_artifact(ctx, sample_annotation, file.path("processed", "stage03", "sample_annotation.tsv"),
                      "03_harmonize_identifiers", script_path, input_paths_02, fwrite_fn)

  qc <- data.frame(
    layer = c("miRNA", "mRNA"),
    declared_scale = c(cfg_expr$mirna_scale, cfg_expr$mrna_scale),
    scale_check_consistent = c(mirna_res$scale_check$consistent, mrna_res$scale_check$consistent),
    transformed_by_pipeline = c(mirna_res$transformed_by_pipeline, mrna_res$transformed_by_pipeline),
    n_features = c(nrow(mirna_matched), nrow(mrna_matched)),
    n_samples = length(common_ids),
    n_newly_na_from_coercion = c(mirna_res$n_newly_na_from_coercion, mrna_res$n_newly_na_from_coercion),
    n_zero_variance_removed = c(mirna_res$n_zero_variance_removed, mrna_res$n_zero_variance_removed),
    winsorized = FALSE,
    stringsAsFactors = FALSE
  )
  write_run_artifact(ctx, qc, file.path("processed", "stage03", "normalization_qc.tsv"),
                      "03_harmonize_identifiers", script_path, input_paths_02, fwrite_fn)

  # --- Sample-exclusion application summary (miRNA/mRNA/clinical rows;
  #     Stage 05 reads it back, appends the purity-layer row and re-writes
  #     the same path under its own stage attribution) ---------------------
  exclusion_summary_rows <- function(layer, n_before, excl_res) {
    if (excl_res$n_removed == 0) {
      return(data.frame(layer = layer, n_matched_records = 0L, removed_identifiers = "",
                         n_remaining = n_before, reason = "(no exclusion matched this layer)", stringsAsFactors = FALSE))
    }
    ids_removed <- if (!is.null(excl_res$removed_names)) excl_res$removed_names else excl_res$removed_ids
    reasons <- unique(excl_res$audit$reason)
    data.frame(layer = layer, n_matched_records = excl_res$n_removed,
               removed_identifiers = paste(ids_removed, collapse = "; "),
               n_remaining = n_before - excl_res$n_removed,
               reason = paste(reasons, collapse = " | "), stringsAsFactors = FALSE)
  }
  exclusion_summary <- rbind(
    exclusion_summary_rows("miRNA", n_mirna_cols_before_excl, excl_mirna),
    exclusion_summary_rows("mRNA", n_mrna_cols_before_excl, excl_mrna),
    exclusion_summary_rows("clinical", n_clinical_before_excl, excl_clinical)
  )
  write_run_artifact(ctx, exclusion_summary, file.path("tables", "audit", "sample_exclusion_expression_summary.tsv"),
                      "03_harmonize_identifiers", script_path,
                      c(exclusions_path, input_paths_02), fwrite_fn)
  log_fn("Sample-exclusion application summary written (miRNA/mRNA/clinical); Stage 05 will append the purity layer.")

  record_stage_status(ctx, "03_harmonize_identifiers", "COMPLETE")
  log_fn("Stage complete. Matched cohort: %d samples, %d miRNA features, %d mRNA features.",
         length(common_ids), nrow(mirna_matched), nrow(mrna_matched))
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=3 or source scripts/run.R in RStudio.")
}
