# ============================================================================
# 06_mirna_cyt_discovery.R
#
# Purpose: compute the Spearman correlation and Benjamini-Hochberg-adjusted
#          significance of every matched miRNA against CYT, independently
#          within every eligible cancer type, and derive a threshold-
#          specific recurrent-candidate miRNA list for each of the three
#          prespecified effect-size thresholds.
# Inputs:  <run_root>/processed/stage03/mirna_matched.rds;
#          <run_root>/processed/stage05/sample_annotation_with_purity.tsv,
#          eligible_cancers_unadjusted.tsv (all same-run)
# Outputs: <run_root>/processed/stage06/master_mirna_cyt.tsv (threshold-invariant raw
#          correlations/p-values/FDR, one row per miRNA x cancer);
#          candidates_threshold_0.30.tsv, 0.25.tsv, 0.20.tsv
# Method:  pairwise-complete Spearman correlation via
#          scripts/lib/stats_engine.R::one_vs_many_spearman_pairwise() (rank transform
#          and Pearson correlation of ranks within each miRNA's complete-case
#          subset; validated against cor.test(method = "spearman",
#          exact = FALSE) in tests/test_pairwise_complete_spearman.R).
#          Benjamini-Hochberg FDR is computed within each cancer type across
#          all eligible miRNAs. This is the complete hypothesis family for
#          the miRNA-CYT arm: it is fixed by data eligibility alone, before
#          any candidate selection, and is independent of later TargetScan
#          or threshold restrictions. Stages 08 and 10 reuse this FDR value
#          for the miRNA-CYT arm and never recompute it on a
#          threshold-restricted subset. Each threshold's candidate list
#          (miRNAs with |rho| >= T and FDR < 0.05 in at least
#          analysis_config$discovery$min_recurrence cancers) is derived
#          independently from this shared table; no list is derived by
#          filtering another threshold's list.
# Depends on: 03_harmonize_identifiers.R, 05_define_eligible_features.R
# ============================================================================

run_stage06 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "06_mirna_cyt_discovery")
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "stats_engine.R"), local = TRUE)
  suppressPackageStartupMessages(library(data.table))
  script_path <- file.path(ctx$data_root, "scripts", "06_mirna_cyt_discovery.R")

  mirna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mirna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mirna_mat <- as.matrix(mirna_matched[, -1, drop = FALSE])
  rownames(mirna_mat) <- mirna_matched$mirna_id

  sample_annotation <- read_run_artifact(ctx, file.path("processed", "stage05", "sample_annotation_with_purity.tsv"), "05_define_eligible_features",
                                          function(p) read.delim(p, stringsAsFactors = FALSE))
  eligible_cancers <- read_run_artifact(ctx, file.path("processed", "stage05", "eligible_cancers_unadjusted.tsv"), "05_define_eligible_features",
                                         function(p) read.delim(p, stringsAsFactors = FALSE))$cancer
  input_paths <- file.path(ctx$run_root, c(file.path("processed", "stage03", "mirna_matched.rds"),
                                            file.path("processed", "stage05", "sample_annotation_with_purity.tsv"),
                                            file.path("processed", "stage05", "eligible_cancers_unadjusted.tsv")))

  run_stage_assert(ctx, log_fn, identical(sample_annotation$sample_id, colnames(mirna_mat)),
                    "Sample order in purity-annotated table matches the miRNA matrix columns")

  log_fn("Validating vectorized statistics engine on real data before large-scale computation ...")
  first_cancer_idx <- which(sample_annotation$cancer == eligible_cancers[1])
  validation <- validate_stats_engine(mirna_mat[, first_cancer_idx, drop = FALSE],
                                       sample_annotation$CYT[first_cancer_idx], seed = analysis_config$seeds$master_seed)
  run_stage_assert(ctx, log_fn, validation$max_cor_diff < 1e-9 && validation$max_p_diff < 1e-6,
                    "Vectorized Spearman correlation/p-value matches the reference cor.test implementation (max cor diff=%.2e, max p diff=%.2e)",
                    validation$max_cor_diff, validation$max_p_diff)

  log_fn("Computing master miRNA-CYT correlation table across %d eligible cancers x %d miRNAs ...",
         length(eligible_cancers), nrow(mirna_mat))
  rows <- vector("list", length(eligible_cancers))
  for (i in seq_along(eligible_cancers)) {
    cc <- eligible_cancers[i]
    idx <- which(sample_annotation$cancer == cc)
    # Pairwise-complete correlation: the matched miRNA matrix is not
    # guaranteed to be finite for every cancer and miRNA, and
    # rank_rows()/safe_rank() refuse non-finite input.
    # one_vs_many_spearman_pairwise() uses the vectorized path for fully
    # finite rows and a per-row pairwise-complete computation otherwise.
    res <- one_vs_many_spearman_pairwise(sample_annotation$CYT[idx], mirna_mat[, idx, drop = FALSE],
                                          min_n = analysis_config$discovery$min_samples_per_cancer)
    rows[[i]] <- data.table(cancer = cc, mirna_id = rownames(mirna_mat), n_samples = res$n,
                             cor_mirna_cyt = res$cor, p_mirna_cyt = res$p)
  }
  master <- rbindlist(rows)
  master[, fdr_mirna_cyt := p.adjust(p_mirna_cyt, "BH"), by = cancer]

  run_stage_assert(ctx, log_fn, nrow(master) == nrow(mirna_mat) * length(eligible_cancers),
                    "Master table has %d rows (%d miRNAs x %d cancers)", nrow(mirna_mat) * length(eligible_cancers),
                    nrow(mirna_mat), length(eligible_cancers))
  run_stage_assert(ctx, log_fn, all(master$fdr_mirna_cyt[!is.na(master$fdr_mirna_cyt)] >= 0 &
                                       master$fdr_mirna_cyt[!is.na(master$fdr_mirna_cyt)] <= 1), "All FDR values lie in [0,1]")
  run_stage_assert(ctx, log_fn, all(abs(master$cor_mirna_cyt[!is.na(master$cor_mirna_cyt)]) <= 1), "All correlations lie in [-1,1]")

  write_fn <- function(x, p) write.table(x, p, sep = "\t", row.names = FALSE, quote = FALSE)
  write_run_artifact(ctx, master, file.path("processed", "stage06", "master_mirna_cyt.tsv"),
                      "06_mirna_cyt_discovery", script_path, input_paths, write_fn)

  candidate_lists <- list()
  for (T in analysis_config$thresholds$all) {
    m <- copy(master)
    m[, pass := abs(cor_mirna_cyt) >= T & fdr_mirna_cyt < analysis_config$discovery$fdr_cutoff]
    agg <- m[, .(n_pass = sum(pass, na.rm = TRUE),
                 cancers_passing = paste(sort(cancer[pass %in% TRUE]), collapse = "; ")),
             by = mirna_id]
    cand <- agg[n_pass >= analysis_config$discovery$min_recurrence]
    setorder(cand, -n_pass, mirna_id)
    fname <- sprintf("candidates_threshold_%.2f.tsv", T)
    write_run_artifact(ctx, cand, file.path("processed", "stage06", fname),
                        "06_mirna_cyt_discovery", script_path, input_paths, write_fn)
    candidate_lists[[as.character(T)]] <- cand$mirna_id
    log_fn("Threshold %.2f: %d candidate miRNAs (independently derived from the shared master table).", T, nrow(cand))
  }

  sorted_thresholds <- sort(analysis_config$thresholds$all, decreasing = TRUE)
  for (i in seq_len(length(sorted_thresholds) - 1)) {
    tighter <- as.character(sorted_thresholds[i]); looser <- as.character(sorted_thresholds[i + 1])
    run_stage_assert(ctx, log_fn, all(candidate_lists[[tighter]] %in% candidate_lists[[looser]]),
                      "Candidate nesting holds: threshold %s candidates are a subset of threshold %s candidates",
                      tighter, looser)
  }

  record_stage_status(ctx, "06_mirna_cyt_discovery", "COMPLETE")
  log_fn("Stage complete.")
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=6 or source scripts/run.R in RStudio.")
}
