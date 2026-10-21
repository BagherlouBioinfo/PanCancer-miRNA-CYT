# ============================================================================
# 09_recurrence_analysis.R
#
# Purpose: aggregate the per-cancer same-cancer three-arm pass indicator
#          into a triplet-level recurrence count, and apply the
#          prespecified qualification rule (recurrence in at least
#          analysis_config$discovery$min_recurrence cancer types).
# Inputs:  <run_root>/processed/stage08/three_arm_threshold_*.tsv (same-run)
# Outputs: <run_root>/tables/primary/triplets_rho030_unadjusted.tsv;
#          tables/sensitivity/triplets_rho025_unadjusted.tsv,
#          triplets_rho020_unadjusted.tsv;
#          tables/audit/cancer_level_evidence_unadjusted_rho*.tsv
# Method:  strict_recurrence(m, g) = number of unique cancer types in which
#          three_arm_pass(m, g, c) is TRUE; a triplet qualifies if
#          strict_recurrence >= min_recurrence. Recurrence is a direct row
#          count on the jointly evaluated long table from Stage 08, so it is
#          never assembled from arms observed in different cancers. The
#          number of qualifying triplets is an output of this computation and
#          is not compared with any expected value.
# Depends on: 08_same_cancer_three_arm_testing.R
# ============================================================================

run_stage09 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "09_recurrence_analysis")
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  suppressPackageStartupMessages(library(data.table))
  script_path <- file.path(ctx$data_root, "scripts", "09_recurrence_analysis.R")

  summarize_triplets <- function(long_dt) {
    agg <- long_dt[, .(
      n_cancers_tested = uniqueN(cancer),
      strict_recurrence = sum(three_arm_pass, na.rm = TRUE),
      supporting_cancers = paste(sort(unique(cancer[three_arm_pass %in% TRUE])), collapse = "; ")
    ), by = .(mirna_id, gene)]
    agg[, qualifies := strict_recurrence >= analysis_config$discovery$min_recurrence]
    agg
  }

  fwrite_fn <- function(x, p) fwrite(x, p, sep = "\t")
  for (T in analysis_config$thresholds$all) {
    fname <- sprintf("three_arm_threshold_%.2f.tsv", T)
    rel_path <- file.path("processed", "stage08", fname)
    long_dt <- read_run_artifact(ctx, rel_path, "08_same_cancer_three_arm_testing", fread)
    input_path <- file.path(ctx$run_root, rel_path)
    long_dt[, run_id := ctx$run_id]; long_dt[, threshold := T]; long_dt[, purity_setting := "unadjusted"]
    long_dt[, generating_script := "08_same_cancer_three_arm_testing.R;09_recurrence_analysis.R"]

    summary_dt <- summarize_triplets(long_dt)

    n_listed <- ifelse(summary_dt$supporting_cancers == "", 0L, lengths(strsplit(summary_dt$supporting_cancers, "; ")))
    run_stage_assert(ctx, log_fn, all(summary_dt$strict_recurrence == n_listed),
                      "[threshold %.2f] strict_recurrence equals the count of unique listed supporting cancers", T)
    qualifying <- summary_dt[qualifies == TRUE]
    if (nrow(qualifying) > 0) {
      check_rows <- long_dt[paste(mirna_id, gene) %in% paste(qualifying$mirna_id, qualifying$gene) & three_arm_pass %in% TRUE]
      run_stage_assert(ctx, log_fn, all(check_rows$pass_mirna_cyt & check_rows$pass_mirna_gene & check_rows$pass_gene_cyt),
                        "[threshold %.2f] Every counted supporting-cancer row independently satisfies all three arm conditions", T)
    }

    summary_dt$run_id <- ctx$run_id
    summary_dt$threshold <- T
    summary_dt$purity_setting <- "unadjusted"
    summary_dt$generating_script <- "09_recurrence_analysis.R"

    thr_tag <- sprintf("rho%03d", round(T * 100))
    out_table <- sprintf("triplets_%s_unadjusted.tsv", thr_tag)
    out_subdir <- if (T == analysis_config$thresholds$primary) "primary" else "sensitivity"
    write_run_artifact(ctx, summary_dt, file.path("tables", out_subdir, out_table),
                        "09_recurrence_analysis", script_path, input_path, fwrite_fn)
    write_run_artifact(ctx, long_dt, file.path("tables", "audit", sprintf("cancer_level_evidence_unadjusted_%s.tsv", thr_tag)),
                        "09_recurrence_analysis", script_path, input_path, fwrite_fn)

    log_fn("Threshold %.2f: %d qualifying triplets (%d unique miRNAs, %d unique genes) written to tables/%s/%s",
           T, sum(summary_dt$qualifies), uniqueN(qualifying$mirna_id), uniqueN(qualifying$gene),
           out_subdir, out_table)
  }

  record_stage_status(ctx, "09_recurrence_analysis", "COMPLETE")
  log_fn("Stage complete.")
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=9 or source scripts/run.R in RStudio.")
}
