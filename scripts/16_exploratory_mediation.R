# ============================================================================
# 17_exploratory_mediation.R
#
# Purpose: exploratory, cancer-specific bootstrap mediation analysis across
#          every evaluable primary-triplet-cancer context of the current run
#          ("Option B": every context with adequate data and valid CPE, not
#          only the contexts that pass the three-arm discovery rule). For
#          each context, tests whether the miRNA's association with CYT is
#          statistically consistent with mediation through its target gene,
#          using CPE-residualized, within-cancer models. This is a
#          post-discovery, exploratory statistical characterization, not a
#          causal mediation analysis: the triplets were selected in the
#          discovery stages using the same TCGA data, temporal ordering of
#          miRNA, gene and CYT expression is unavailable in cross-sectional
#          RNA-seq data, and unmeasured confounding remains possible.
#          Contexts are not selected on any mediation result.
# Inputs:  the current run's primary triplet set (scripts/lib/run_discovery_source.R);
#          <run_root>/processed/stage03/mirna_matched.rds, mrna_matched.rds
#          (same-run); processed/stage04/cyt_scores.tsv (same-run);
#          processed/stage05/sample_annotation_with_purity.tsv (same-run);
#          tables/downstream/meta_analysis_input_all_cancers.tsv (Stage 14,
#          same-run; its CPE eligibility defines the Option B universe).
# Outputs: <run_root>/tables/downstream/mediation_all_contexts.tsv (one row
#          per context, including the bootstrap diagnostics n_successful,
#          n_failed and runtime_sec);
#          mediation_primary_support_subset.tsv (the primary-supporting
#          contexts, flagged rather than re-estimated);
#          processed/stage15_mediation/iteration_level_indirect_effects.tsv
#          (every bootstrap draw of the indirect effect for every context);
#          reports/31_exploratory_mediation_report.md.
# Method:  for each context, miRNA, gene and CYT are each residualized on
#          CPE within the cancer (as in the discovery stages); a-path:
#          gene_resid ~ mirna_resid; b-path: CYT_resid ~ mirna_resid +
#          gene_resid; indirect = a * b; direct (c') = the mirna_resid
#          coefficient of the b-path model; total = indirect + direct.
#          analysis_config$resampling$mediation_n_boot within-cancer
#          bootstrap resamples (with replacement) are drawn per context and
#          both models are refitted in every resample. Percentile confidence
#          intervals use analysis_config$resampling$mediation_ci_level. The
#          bootstrap p-value is twice the proportion of resamples on the
#          opposite side of zero from the point estimate, with the finite-
#          resample correction described below, capped at 1. BH-FDR is
#          applied once across all evaluable contexts (one prespecified
#          family).
#
#          Every context is processed in full on every execution; results
#          are accumulated in memory and written once at the end through
#          write_run_artifact(). There is no checkpoint or resume mechanism,
#          consistent with the isolated-run design.
# Depends on: 03_harmonize_identifiers.R, 04_compute_validate_cyt.R,
#             05_define_eligible_features.R, 09_recurrence_analysis.R,
#             10_purity_adjusted_analysis.R, 14_cross_cancer_meta_analysis.R
#             (all same-run).
# ============================================================================

run_stage17 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "17_exploratory_mediation")
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "run_discovery_source.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "expression_preprocessing.R"), local = TRUE)
  suppressPackageStartupMessages(library(data.table))
  script_path <- file.path(ctx$data_root, "scripts", "16_exploratory_mediation.R")
  fwrite_fn <- function(x, p) fwrite(x, p, sep = "\t")

  triplets <- load_run_primary_triplets(ctx)
  meta_input <- read_run_artifact(ctx, file.path("tables", "downstream", "meta_analysis_input_all_cancers.tsv"),
                                   "14_cross_cancer_meta_analysis", fread)

  # Option B universe: every (mirna, gene, cancer) context in which all three arms have a valid
  # CPE-adjusted estimate (the eligibility definition of Stage 14). primary_support flags the
  # current run's three-arm-passing subset.
  cpe_eligible <- meta_input[adjustment == "CPE", .(mirna_id, gene, cancer, arm, eligible, primary_support)]
  wide_elig <- dcast(cpe_eligible, mirna_id + gene + cancer + primary_support ~ arm, value.var = "eligible")
  universe <- wide_elig[`miRNA-CYT` == TRUE & `miRNA-gene` == TRUE & `gene-CYT` == TRUE,
                         .(mirna_id, gene, cancer, primary_support)]
  # The universe size follows from the run's discovery output and the CPE-eligibility rule; it is
  # logged and not compared with any fixed value.
  run_stage_assert(ctx, log_fn, nrow(universe) > 0, "Option B universe is non-empty")

  n_boot <- analysis_config$resampling$mediation_n_boot
  ci_level <- analysis_config$resampling$mediation_ci_level
  seed <- analysis_config$seeds$bootstrap
  log_fn("Option B universe: %d contexts (%d this-run primary-support). n_boot=%d, ci_level=%.2f, seed=%d.",
         nrow(universe), sum(universe$primary_support), n_boot, ci_level, seed)

  mirna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mirna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mirna_mat <- as.matrix(mirna_matched[, -1, drop = FALSE]); rownames(mirna_mat) <- mirna_matched$mirna_id
  mrna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mrna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mrna_mat <- as.matrix(mrna_matched[, -1, drop = FALSE]); rownames(mrna_mat) <- toupper(trimws(mrna_matched$gene_id_raw))
  # Per-feature winsorization before the mediation path regressions, which
  # are OLS fits on expression values and sensitive to extreme samples.
  mirna_mat <- winsorize_matrix_rows(mirna_mat)
  mrna_mat <- winsorize_matrix_rows(mrna_mat)
  cyt_scores <- read_run_artifact(ctx, file.path("processed", "stage04", "cyt_scores.tsv"), "04_compute_validate_cyt", fread)
  sample_annot <- read_run_artifact(ctx, file.path("processed", "stage05", "sample_annotation_with_purity.tsv"), "05_define_eligible_features", fread)

  build_cancer_df <- function(mir, gene, cancer_code) {
    idx <- which(sample_annot$cancer == cancer_code & !is.na(sample_annot$CPE) & !is.na(sample_annot$CYT))
    d <- data.frame(
      sample_id = sample_annot$sample_id[idx], CPE = sample_annot$CPE[idx],
      mirna_expr = as.numeric(mirna_mat[mir, idx]), gene_expr = as.numeric(mrna_mat[gene, idx]),
      CYT = cyt_scores$CYT[match(sample_annot$sample_id[idx], cyt_scores$sample_id)]
    )
    d[complete.cases(d), ]
  }

  cpe_residualized_paths <- function(d) {
    m_resid <- residuals(lm(mirna_expr ~ CPE, data = d))
    g_resid <- residuals(lm(gene_expr ~ CPE, data = d))
    y_resid <- residuals(lm(CYT ~ CPE, data = d))
    a_fit <- tryCatch(lm(g_resid ~ m_resid), error = function(e) NULL)
    b_fit <- tryCatch(lm(y_resid ~ m_resid + g_resid), error = function(e) NULL)
    if (is.null(a_fit) || is.null(b_fit)) return(list(converged = FALSE, a = NA, b = NA, c_prime = NA, indirect = NA, total = NA))
    a <- unname(coef(a_fit)["m_resid"]); b <- unname(coef(b_fit)["g_resid"]); c_prime <- unname(coef(b_fit)["m_resid"])
    list(converged = is.finite(a) && is.finite(b) && is.finite(c_prime), a = a, b = b, c_prime = c_prime,
         indirect = a * b, total = a * b + c_prime)
  }

  universe[, context_key := paste(mirna_id, gene, cancer, sep = "||")]
  log_fn("%d contexts to process.", nrow(universe))

  iter_row_list <- vector("list", nrow(universe))
  summary_row_list <- vector("list", nrow(universe))

  for (i in seq_len(nrow(universe))) {
    mir <- universe$mirna_id[i]; tgt_gene <- universe$gene[i]; cc <- universe$cancer[i]; psupport <- universe$primary_support[i]
    t0 <- Sys.time()
    d <- build_cancer_df(mir, tgt_gene, cc)
    n_total <- nrow(d)
    point <- cpe_residualized_paths(d)

    set.seed(seed + i)
    iter_rows <- vector("list", n_boot)
    n_success <- 0L
    for (b in seq_len(n_boot)) {
      boot_idx <- sample(seq_len(n_total), n_total, replace = TRUE)
      est <- cpe_residualized_paths(d[boot_idx, ])
      if (isTRUE(est$converged)) n_success <- n_success + 1L
      iter_rows[[b]] <- data.table(mirna_id = mir, gene = tgt_gene, cancer = cc, iteration = b,
                                    indirect = est$indirect, a = est$a, b_coef = est$b, c_prime = est$c_prime,
                                    total = est$total, converged = est$converged)
    }
    iter_dt <- rbindlist(iter_rows)
    iter_row_list[[i]] <- iter_dt

    conv_dt <- iter_dt[converged == TRUE]
    n_failed <- n_boot - n_success
    if (nrow(conv_dt) > 1) {
      ci <- quantile(conv_dt$indirect, c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2), na.rm = TRUE)
      point_sign <- sign(point$indirect)
      # Finite-bootstrap correction (Davison & Hinkley 1997): (count + 1) /
      # (B + 1) rather than count / B, so that the bootstrap p-value is never
      # exactly zero when no resample falls on the opposite side of zero.
      opposite_count <- sum(sign(conv_dt$indirect) != point_sign & conv_dt$indirect != 0)
      raw_p <- min(1, 2 * (opposite_count + 1) / (nrow(conv_dt) + 1))
      n_same_sign <- sum(sign(conv_dt$indirect) == point_sign)
      pct_dir_consistent <- round(100 * n_same_sign / nrow(conv_dt), 1)
      prop_mediated <- if (!is.na(point$total) && abs(point$total) > 1e-6 && sign(point$indirect) == sign(point$total)) point$indirect / point$total else NA_real_
      prop_meaningful <- !is.na(prop_mediated)
    } else {
      ci <- c(NA_real_, NA_real_); raw_p <- NA_real_; n_same_sign <- NA_integer_; pct_dir_consistent <- NA_real_
      prop_mediated <- NA_real_; prop_meaningful <- FALSE
    }
    runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    summary_row_list[[i]] <- data.table(
      mirna_id = mir, gene = tgt_gene, cancer = cc, primary_support = psupport, n = n_total,
      n_iterations = n_boot, n_successful = n_success, n_failed = n_failed,
      point_a = point$a, point_b = point$b, point_c_prime = point$c_prime, point_indirect = point$indirect, point_total = point$total,
      ci_lower = unname(ci[1]), ci_upper = unname(ci[2]), prop_mediated = prop_mediated, prop_mediated_meaningful = prop_meaningful,
      raw_bootstrap_p = raw_p, n_same_sign_as_point = n_same_sign, pct_direction_consistent = pct_dir_consistent,
      runtime_sec = round(runtime, 3), seed = seed + i
    )

    if (i %% 20 == 0 || i == nrow(universe)) log_fn("  Completed %d/%d contexts (%s|%s|%s, %.1fs).", i, nrow(universe), mir, tgt_gene, cc, runtime)
  }

  base_inputs <- file.path(ctx$run_root, c(
    file.path("tables", "primary", run_primary_triplet_filename(ctx)),
    file.path("tables", "downstream", "meta_analysis_input_all_cancers.tsv"),
    file.path("processed", "stage03", "mirna_matched.rds"), file.path("processed", "stage03", "mrna_matched.rds"),
    file.path("processed", "stage04", "cyt_scores.tsv"), file.path("processed", "stage05", "sample_annotation_with_purity.tsv")
  ))

  iter_dt_all <- rbindlist(iter_row_list, use.names = TRUE)
  write_run_artifact(ctx, iter_dt_all, file.path("processed", "stage15_mediation", "iteration_level_indirect_effects.tsv"),
                      "17_exploratory_mediation", script_path, base_inputs, fwrite_fn)

  full_summary <- rbindlist(summary_row_list, use.names = TRUE)
  full_summary[, bh_fdr := p.adjust(raw_bootstrap_p, method = "BH")]
  write_run_artifact(ctx, full_summary, file.path("tables", "downstream", "mediation_all_contexts.tsv"),
                      "17_exploratory_mediation", script_path, base_inputs, fwrite_fn)

  primary_subset <- full_summary[primary_support == TRUE]
  write_run_artifact(ctx, primary_subset, file.path("tables", "downstream", "mediation_primary_support_subset.tsv"),
                      "17_exploratory_mediation", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "mediation_all_contexts.tsv"), fwrite_fn)

  n_unreliable <- sum(full_summary$n_successful / full_summary$n_iterations < 0.95)
  log_fn("mediation_all_contexts.tsv written (%d rows, %d flagged unreliable due to <95%% successful iterations).",
         nrow(full_summary), n_unreliable)
  log_fn("mediation_primary_support_subset.tsv written (%d rows).", nrow(primary_subset))

  # ==========================================================================
  # REPORT 31
  # ==========================================================================
  min_reliable_prop <- 0.95
  full_summary[, reliable := (n_successful / n_iterations) >= min_reliable_prop]
  n_sig_fdr <- sum(full_summary$bh_fdr < 0.05 & full_summary$reliable, na.rm = TRUE)
  # Membership is matched on the exact (mirna_id, gene, cancer) key rather
  # than on independent mirna_id and gene tests. primary_subset is a row
  # subset of full_summary, so its bh_fdr, n_successful and n_iterations
  # columns already refer to the same context.
  primary_subset[, reliable := (n_successful / n_iterations) >= min_reliable_prop]
  n_sig_fdr_primary <- sum(primary_subset$bh_fdr < 0.05 & primary_subset$reliable, na.rm = TRUE)

  report31 <- c(
    "# Stage 17: Exploratory Cancer-Specific Mediation (Option B -- All Evaluable Contexts)",
    "",
    sprintf("Generated by `scripts/16_exploratory_mediation.R`, run `%s`.", ctx$run_id),
    "",
    "## Design",
    "",
    sprintf("Option B universe: all %d evaluable (miRNA, gene, cancer) contexts with adequate matched data and valid CPE (n>=%d, same eligibility already verified in Stage 14), not only the %d contexts that already pass this run's own three-arm discovery threshold. The %d primary-support contexts are flagged (`primary_support` column) as a descriptive subset, not used as the inferential universe.",
            nrow(universe), analysis_config$discovery$min_samples_per_cancer, sum(universe$primary_support), sum(universe$primary_support)),
    sprintf("Each context: CPE-residualized within-cancer mediator/outcome models (a-path: gene-residual ~ miRNA-residual; b-path: CYT-residual ~ miRNA-residual + gene-residual), %d bootstrap iterations resampling subjects with replacement WITHIN that cancer only, refitting both models every iteration. Percentile %.0f%% CI. Raw bootstrap inference measure: twice the proportion of iterations on the opposite side of zero from the point estimate, capped at 1. BH-FDR: single flat family across all %d contexts, decided prospectively, before any context's result was inspected.",
            n_boot, ci_level * 100, nrow(universe)),
    "No cancer was pooled with any other at any point.",
    "",
    "## Completion and reliability",
    "",
    sprintf("- Contexts processed: %d of %d.", nrow(full_summary), nrow(universe)),
    sprintf("- Contexts with >=%.0f%% successful (converged) iterations (reliable): %d.", min_reliable_prop * 100, sum(full_summary$reliable)),
    sprintf("- Contexts below the %.0f%% reliability threshold (labeled unreliable, not significant/nonsignificant): %d.", min_reliable_prop * 100, sum(!full_summary$reliable)),
    sprintf("- BH-FDR < 0.05 among reliable contexts: %d of %d evaluated (%d of these are among the %d primary-support contexts).", n_sig_fdr, sum(full_summary$reliable), n_sig_fdr_primary, sum(universe$primary_support)),
    "",
    "## Selection-bias and interpretive limitations",
    "",
    sprintf("- The %d this-run primary triplets underlying every context here were selected during discovery (same-cancer three-arm recurrence >= %d); this residual selection is not removed by broadening the mediation universe to all %d contexts.",
            uniqueN(triplets[, .(mirna_id, gene)]), analysis_config$discovery$min_recurrence, nrow(universe)),
    "- Cancer contexts within each triplet were **not** selected based on any mediation result -- eligibility depends only on sample-size and CPE-availability, decided before this script was written.",
    "- This analysis is post-discovery and exploratory. It characterizes statistical consistency with a mediation pattern; it does not establish causal mediation.",
    "- The same TCGA cross-sectional cohort was used for discovery, meta-analysis, survival, and this mediation analysis -- there is no independent replication cohort anywhere in this project.",
    "- Temporal ordering between miRNA expression, target-gene expression, and CYT is unavailable in cross-sectional RNA-seq data; the a-path/b-path directionality is a modeling assumption, not an observed sequence.",
    "- Unmeasured confounding (beyond tumor purity, which is explicitly adjusted for) remains possible throughout.",
    "- No result in this report or its output tables should be read as a confirmed or causal mediation claim.",
    "",
    "## Outputs",
    "",
    "- `tables/downstream/mediation_all_contexts.tsv` (includes bootstrap diagnostic columns: n_successful, n_failed, runtime_sec), `mediation_primary_support_subset.tsv`",
    "- `processed/stage15_mediation/iteration_level_indirect_effects.tsv` (all iteration-level draws)",
    "- `reports/31_exploratory_mediation_report.md` (this file)"
  )
  write_run_artifact(ctx, report31, file.path("reports", "31_exploratory_mediation_report.md"),
                      "17_exploratory_mediation", script_path, character(0), function(x, p) writeLines(x, p))
  log_fn("reports/31_exploratory_mediation_report.md written.")

  record_stage_status(ctx, "17_exploratory_mediation", "COMPLETE")
  log_fn("Stage 17 complete.")
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=16 or source scripts/run.R in RStudio.")
}
