# ============================================================================
# 01_configure_and_validate_inputs.R
#
# Purpose: load the analysis configuration and verify, before any
#          computation, that every required raw and reference input is
#          present, non-empty and listed with a matching SHA-256 in
#          config/FROM_SCRATCH_ALLOWED_INPUTS.tsv. Later stages read raw
#          inputs only after re-verifying them against the same allowlist,
#          and read upstream results only as same-run, provenance-verified
#          artifacts through scripts/lib/runtime.R.
# Inputs:  config/analysis_config.R; the required paths in
#          analysis_config$paths; config/FROM_SCRATCH_ALLOWED_INPUTS.tsv.
# Outputs: <run_root>/logs/01_configure_and_validate_inputs.log;
#          <run_root>/processed/stage01/run_manifest.tsv (recorded with
#          provenance; it has no upstream scientific input, so input_paths
#          is empty).
# Method:  deterministic file-existence, non-zero-size and allowlist-hash
#          checks; no statistical computation.
# Depends on: none (first stage).
# ============================================================================

run_stage01 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "01_configure_and_validate_inputs")
  log_fn("Data root: %s", ctx$data_root)
  log_fn("Run root (isolated): %s", ctx$run_root)
  log_fn("Run ID: %s", ctx$run_id)

  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  validate_analysis_config(analysis_config)
  log_fn("Configuration loaded and validated: primary threshold=%.2f, sensitivity=%s, purity primary=%s",
         analysis_config$thresholds$primary,
         paste(analysis_config$thresholds$sensitivity, collapse = ","),
         analysis_config$purity$primary_metric)

  # Required inputs are those consumed by the downstream stages, verified
  # here so that a missing or corrupted file is detected before any
  # long-running stage starts. mirtarbase_mti, mirtarbase_strong,
  # drug_gene_interactions and raw_gene_id are declared in
  # analysis_config$paths but are not read by any stage. mir_family_info is
  # verified against the allowlist by Stages 07 and 20 immediately before
  # they read it.
  required <- analysis_config$paths[c("raw_mrna", "raw_mirna", "raw_survival", "targetscan",
                                       "mcpcounter_genes", "mcpcounter_signature_manifest")]
  manifest_rows <- list()
  missing <- character(0)
  for (nm in names(required)) {
    full_path <- file.path(ctx$data_root, required[[nm]])
    exists_ok <- file.exists(full_path)
    size_ok <- exists_ok && file.info(full_path)$size > 0
    allowlist_ok <- FALSE
    allowlist_err <- NA_character_
    if (exists_ok && size_ok) {
      allowlist_result <- tryCatch({
        verify_allowed_raw_input(ctx$data_root, full_path)
        TRUE
      }, error = function(e) { allowlist_err <<- conditionMessage(e); FALSE })
      allowlist_ok <- isTRUE(allowlist_result)
    }
    if (!exists_ok || !size_ok) missing <- c(missing, sprintf("%s (%s): file missing or empty", nm, full_path))
    if (exists_ok && size_ok && !allowlist_ok) missing <- c(missing, sprintf("%s (%s): %s", nm, full_path, allowlist_err))
    manifest_rows[[nm]] <- data.frame(
      input_name = nm, path = full_path, exists = exists_ok,
      size_bytes = if (exists_ok) file.info(full_path)$size else NA_integer_,
      allowlist_verified = allowlist_ok,
      stringsAsFactors = FALSE
    )
  }
  manifest <- do.call(rbind, manifest_rows)

  if (length(missing) > 0) {
    run_stage_fail(ctx, log_fn, "Required input(s) missing, empty, or not allowlist-verified: %s", paste(missing, collapse = "; "))
  }
  for (i in seq_len(nrow(manifest))) {
    log_fn("Input OK (allowlist-verified): %s (%d bytes)", manifest$path[i], manifest$size_bytes[i])
  }

  manifest$run_id <- ctx$run_id
  manifest$config_version <- analysis_config$config_version
  write_run_artifact(ctx, manifest, file.path("processed", "stage01", "run_manifest.tsv"),
                      stage = "01_configure_and_validate_inputs",
                      script_path = file.path(ctx$data_root, "scripts", "01_configure_and_validate_inputs.R"),
                      input_paths = character(0),
                      write_fn = function(x, p) write.table(x, p, sep = "\t", row.names = FALSE, quote = FALSE))

  record_stage_status(ctx, "01_configure_and_validate_inputs", "COMPLETE")
  log_fn("Stage complete. %d required inputs verified against the from-scratch allowlist.", nrow(manifest))
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=1 or source scripts/run.R in RStudio.")
}
