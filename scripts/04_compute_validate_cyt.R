# ============================================================================
# 04_compute_validate_cyt.R
#
# Purpose: compute the cytolytic activity (CYT) score for every matched
#          sample as the arithmetic mean of log2-scale GZMA and PRF1
#          expression, and validate its distribution and completeness.
# Inputs:  <run_root>/processed/stage03/mrna_matched.rds (same-run)
# Outputs: <run_root>/processed/stage04/cyt_scores.tsv
# Method:  CYT_i = (GZMA_i + PRF1_i) / 2 on the log2-scale matched mRNA
#          matrix from Stage 03, i.e. the log2 of the geometric mean of the
#          two linear-scale values. Absolute CYT values are therefore not
#          directly comparable with studies that use the linear-scale
#          geometric mean (Rooney et al. 2015), but the sample ranking used
#          by the rank-based analyses of this pipeline is preserved.
# Depends on: 03_harmonize_identifiers.R
# ============================================================================

run_stage04 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "04_compute_validate_cyt")
  source(file.path(ctx$data_root, "scripts", "lib", "expression_preprocessing.R"), local = TRUE)
  script_path <- file.path(ctx$data_root, "scripts", "04_compute_validate_cyt.R")

  mrna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mrna_matched.rds"), "03_harmonize_identifiers", readRDS)
  input_path <- file.path(ctx$run_root, "processed", "stage03", "mrna_matched.rds")

  mrna_mat <- as.matrix(mrna_matched[, -1, drop = FALSE])
  rownames(mrna_mat) <- mrna_matched$gene_id_raw

  gzma <- extract_gene_row(mrna_mat, "GZMA")
  prf1 <- extract_gene_row(mrna_mat, "PRF1")
  cyt <- (gzma + prf1) / 2

  cyt_scores <- data.frame(sample_id = colnames(mrna_mat), CYT = cyt, GZMA = gzma, PRF1 = prf1, stringsAsFactors = FALSE)

  run_stage_assert(ctx, log_fn, all(is.finite(cyt_scores$CYT)), "CYT is finite for every matched sample")
  run_stage_assert(ctx, log_fn, nrow(cyt_scores) == ncol(mrna_mat), "CYT computed for every matched sample (%d)", ncol(mrna_mat))
  log_fn("CYT range: [%.3f, %.3f], mean=%.3f, n=%d", min(cyt), max(cyt), mean(cyt), length(cyt))

  write_run_artifact(ctx, cyt_scores, file.path("processed", "stage04", "cyt_scores.tsv"),
                      "04_compute_validate_cyt", script_path, input_path,
                      function(x, p) write.table(x, p, sep = "\t", row.names = FALSE, quote = FALSE))

  record_stage_status(ctx, "04_compute_validate_cyt", "COMPLETE")
  log_fn("Stage complete.")
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=4 or source scripts/run.R in RStudio.")
}
