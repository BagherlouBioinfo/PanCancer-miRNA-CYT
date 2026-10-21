# ============================================================================
# 02_ingest_and_provenance.R
#
# Purpose: read the raw miRNA expression, mRNA expression and
#          clinical/survival tables and record a SHA-256 hash and dimension
#          summary for each before any transformation. The three files are
#          re-verified against the allowlist immediately before parsing.
# Inputs:  data/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena;
#          data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena;
#          data/Survival_SupplementalTable_S1_20171025_xena_sp
# Outputs: <run_root>/processed/stage02/input_hashes.tsv, raw_mirna.rds,
#          raw_mrna.rds, raw_survival.rds
# Method:  file hashing (digest::digest, SHA-256); tab-delimited ingestion
#          (data.table::fread).
# Depends on: 01_configure_and_validate_inputs.R
# ============================================================================

run_stage02 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "02_ingest_and_provenance")
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  suppressPackageStartupMessages({ library(data.table); library(digest) })

  paths <- list(
    raw_mirna = file.path(ctx$data_root, analysis_config$paths$raw_mirna),
    raw_mrna = file.path(ctx$data_root, analysis_config$paths$raw_mrna),
    raw_survival = file.path(ctx$data_root, analysis_config$paths$raw_survival)
  )
  for (nm in names(paths)) verify_allowed_raw_input(ctx$data_root, paths[[nm]])
  log_fn("All 3 raw inputs re-verified against the from-scratch allowlist immediately before reading.")

  hash_rows <- list()
  for (nm in names(paths)) {
    log_fn("Hashing %s ...", paths[[nm]])
    h <- digest(paths[[nm]], algo = "sha256", file = TRUE)
    hash_rows[[nm]] <- data.frame(input = nm, path = paths[[nm]], sha256 = h,
                                   size_bytes = file.info(paths[[nm]])$size, stringsAsFactors = FALSE)
    log_fn("SHA-256(%s) = %s", nm, h)
  }
  write_run_artifact(ctx, do.call(rbind, hash_rows), file.path("processed", "stage02", "input_hashes.tsv"),
                      stage = "02_ingest_and_provenance",
                      script_path = file.path(ctx$data_root, "scripts", "02_ingest_and_provenance.R"),
                      input_paths = unlist(paths),
                      write_fn = function(x, p) write.table(x, p, sep = "\t", row.names = FALSE, quote = FALSE))

  log_fn("Reading raw miRNA expression matrix ...")
  raw_mirna <- fread(paths$raw_mirna, data.table = FALSE)
  log_fn("Raw miRNA matrix: %d rows x %d columns", nrow(raw_mirna), ncol(raw_mirna))

  log_fn("Reading raw mRNA expression matrix (this can take a moment for a >1GB file) ...")
  raw_mrna <- fread(paths$raw_mrna, data.table = FALSE)
  log_fn("Raw mRNA matrix: %d rows x %d columns", nrow(raw_mrna), ncol(raw_mrna))

  log_fn("Reading raw survival/clinical table ...")
  raw_survival <- fread(paths$raw_survival, data.table = FALSE)
  log_fn("Raw survival table: %d rows x %d columns", nrow(raw_survival), ncol(raw_survival))

  run_stage_assert(ctx, log_fn, nrow(raw_mirna) > 0 && nrow(raw_mrna) > 0 && nrow(raw_survival) > 0,
                    "All three raw inputs contain at least one row")

  write_run_artifact(ctx, raw_mirna, file.path("processed", "stage02", "raw_mirna.rds"),
                      stage = "02_ingest_and_provenance", script_path = file.path(ctx$data_root, "scripts", "02_ingest_and_provenance.R"),
                      input_paths = paths$raw_mirna, write_fn = saveRDS)
  write_run_artifact(ctx, raw_mrna, file.path("processed", "stage02", "raw_mrna.rds"),
                      stage = "02_ingest_and_provenance", script_path = file.path(ctx$data_root, "scripts", "02_ingest_and_provenance.R"),
                      input_paths = paths$raw_mrna, write_fn = saveRDS)
  write_run_artifact(ctx, raw_survival, file.path("processed", "stage02", "raw_survival.rds"),
                      stage = "02_ingest_and_provenance", script_path = file.path(ctx$data_root, "scripts", "02_ingest_and_provenance.R"),
                      input_paths = paths$raw_survival, write_fn = saveRDS)

  record_stage_status(ctx, "02_ingest_and_provenance", "COMPLETE")
  log_fn("Stage complete.")
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=2 or source scripts/run.R in RStudio.")
}
