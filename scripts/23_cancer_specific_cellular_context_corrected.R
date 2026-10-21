# ============================================================================
# 24_cancer_specific_cellular_context_corrected.R
#
# Cancer-specific cellular-composition context analysis with conjunction
# inference. This stage refines the Stage 23 analysis. It reads Stage 23's
# summary table (same-run) only for the historical Stage 23 reconciliation
# and recomputes all other quantities from the same-run matched-cohort
# artifacts; Stage 23 outputs are not modified. The Stage 24 result is the
# primary evidence funnel, which starts from every conjunction-supported
# combination in the full corrected universe (R/cell_context_funnel.R).
#
# Methods, relative to Stage 23:
#   1. Purity (CPE) is an explicit covariate inside the interaction model,
#      rank(Y) ~ rank(X) + Z_pop + rank(X):Z_pop + purity, rather than a
#      residualization applied to X and Y beforehand, so the MCPcounter
#      population score is also purity-adjusted.
#   2. Inference for the interaction coefficient uses HC3
#      heteroscedasticity-consistent standard errors
#      (R/cell_context_helpers.R::fit_cell_context_interaction_hc3()).
#      Classical OLS standard errors and p-values are retained only as a
#      diagnostic comparison. Every fit records the maximum leverage (hat
#      value) and any warning raised during HC3 estimation; fits with a
#      warning or maximum leverage >= 0.99 are flagged leverage_unstable and
#      excluded from the conjunction family.
#   3. The "all three arms" claim is tested as a conjunction
#      (intersection-union) test: p_conjunction = max(p1, p2, p3) of the
#      three arms' raw HC3 p-values, BH-corrected once across the family of
#      evaluable conjunction tests (rather than requiring three separately
#      adjusted arm-level FDRs). Benjamini-Yekutieli adjustment is reported
#      as a dependence-robust sensitivity analysis.
#   4. Every (gene, population) pair is classified as direct_tautology,
#      CYT_biological_overlap or nonoverlapping_population; results in
#      CYT-overlapping populations are reported as expected immune-context
#      sensitivity, not as independent mechanistic evidence.
#   5. Population correlations and a collinearity sensitivity refit (VIF,
#      condition number, shared lymphocyte-infiltration PC1 axis) for every
#      conjunction-supported finding. A finding without a correlated
#      population is recorded as not applicable, not as demonstrated
#      robustness.
#   6. Fixed-seed permutation negative controls on a documented subset of
#      combinations (see config/analysis_config.R); synthetic known-answer
#      scenarios are tested in tests/test_cell_context_conjunction_helpers.R.
#   7. Per-cancer power and sample-size reporting, including whether the
#      cancer-type distribution of findings is explained by cohort size or
#      population-score variance.
#   8. A prespecified nonlinear (spline) sensitivity check for
#      conjunction-supported candidates only (descriptive).
#   9. LUMP-purity sensitivity refit of the complete three-arm conjunction
#      for cancers with valid LUMP values, and a completeness-gated AJCC
#      stage sensitivity refit restricted to cancers with adequate stage
#      completeness (cancers without stage data are excluded from that
#      check, not imputed).
#
# Inputs (all same-run, provenance-verified): Stage 23's
#   tables/downstream/triplet_cell_context_summary.tsv (historical Stage 23
#   reconciliation only); tables/primary/<primary triplets>;
#   processed/stage03/{mirna,mrna}_matched.rds;
#   processed/stage05/sample_annotation_with_purity.tsv,
#   eligible_cancers_cpe.tsv; processed/stage11_immune/
#   sample_level_immune_estimates.tsv; processed/stage03/survival_matched.rds
#   (AJCC stage sensitivity only); reference/mcpcounter_signatures/genes.txt.
# Outputs: tables/downstream/stage24_corrected_arm_level.tsv,
#   stage24_corrected_conjunction.tsv,
#   stage24_nonoverlapping_population_results.tsv,
#   stage24_cyt_overlap_results.tsv, stage24_population_correlations.tsv,
#   stage24_collinearity_sensitivity.tsv, stage24_negative_controls.tsv,
#   stage24_power_by_cancer.tsv, stage24_nonlinearity_sensitivity.tsv,
#   stage24_lump_and_stage_sensitivity.tsv,
#   stage24_primary_evidence_funnel.tsv (primary result funnel, from every
#   conjunction-supported combination), stage24_historical_reconciliation.tsv
#   (Stage 23 all-three-arms set followed through the same checks; comparison
#   only), stage24_final_robust_nonoverlapping.tsv (row-level final result),
#   stage24_final_triplet_summary.tsv; figures/cell_context_corrected/*.pdf,
#   *.png; reports/STAGE24_CELLULAR_CONTEXT_METHOD_AND_RESULTS.md.
#   Stage 23 outputs keep their own names (triplet_cell_context_*.tsv).
# Depends on (same-run, enforced by read_run_artifact()): 03, 05, 09, 10,
#   12, 23.
# ============================================================================

run_stage24 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "24_cancer_specific_cellular_context_corrected")
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "run_discovery_source.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "expression_preprocessing.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "stats_engine.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "cell_context_helpers.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "cell_context_funnel.R"), local = environment(run_stage24))
  suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(sandwich)
  })
  script_path <- file.path(ctx$data_root, "scripts", "23_cancer_specific_cellular_context_corrected.R")
  fwrite_fn <- function(x, p) fwrite(x, p, sep = "\t")
  cfg <- analysis_config$cell_context_corrected
  min_n <- cfg$min_n
  conj_alpha <- cfg$conjunction_fdr_cutoff
  cyt_overlap_pops <- cfg$cyt_overlap_populations
  corr_threshold <- cfg$population_correlation_threshold

  fit_hc3_checked <- function(predictor, outcome, population, purity, min_n) {
    warning_messages <- character(0)
    fit <- withCallingHandlers(
      fit_cell_context_interaction_hc3(
        predictor, outcome, population, purity, min_n = min_n
      ),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )

    complete <- is.finite(predictor) & is.finite(outcome) &
      is.finite(population) & is.finite(purity)
    max_hat <- NA_real_
    if (sum(complete) >= min_n && stats::sd(population[complete]) > 0) {
      diagnostic_data <- data.frame(
        r_out = rank(outcome[complete], ties.method = "average"),
        r_pred = rank(predictor[complete], ties.method = "average"),
        pop_z = as.numeric(scale(population[complete])),
        purity = purity[complete]
      )
      diagnostic_fit <- tryCatch(
        lm(r_out ~ r_pred * pop_z + purity, data = diagnostic_data),
        error = function(e) NULL
      )
      if (!is.null(diagnostic_fit)) {
        max_hat <- tryCatch(max(hatvalues(diagnostic_fit), na.rm = TRUE),
                            error = function(e) NA_real_)
      }
    }

    hc3_warning <- length(warning_messages) > 0
    leverage_unstable <- isTRUE(hc3_warning) || isTRUE(is.finite(max_hat) && max_hat >= 0.99)
    fit$max_hat <- max_hat
    fit$hc3_warning <- hc3_warning
    fit$hc3_warning_message <- if (length(warning_messages)) {
      paste(unique(warning_messages), collapse = " | ")
    } else {
      ""
    }
    fit$leverage_unstable <- leverage_unstable
    fit
  }

  # --------------------------------------------------------------------------
  # Same-run inputs (sourced as in Stage 23).
  # --------------------------------------------------------------------------
  triplets <- load_run_primary_triplets(ctx)
  run_stage_assert(ctx, log_fn, all(c("mirna_id", "gene") %in% names(triplets)) && nrow(triplets) > 0,
                    "This run's primary triplet table has mirna_id/gene columns and at least 1 row (%d triplets loaded)", nrow(triplets))

  sig_path <- file.path(ctx$data_root, "reference", "mcpcounter_signatures", "genes.txt")
  verify_allowed_raw_input(ctx$data_root, sig_path)
  sig_raw <- read.table(sig_path, sep = "\t", stringsAsFactors = FALSE, header = TRUE, colClasses = "character", check.names = FALSE)
  signature_map <- data.frame(gene = toupper(trimws(sig_raw[["HUGO symbols"]])), population = sig_raw[["Cell population"]], stringsAsFactors = FALSE)
  populations <- sort(unique(signature_map$population))

  mirna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mirna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mirna_mat <- as.matrix(mirna_matched[, -1, drop = FALSE]); rownames(mirna_mat) <- mirna_matched$mirna_id
  mrna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mrna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mrna_mat <- as.matrix(mrna_matched[, -1, drop = FALSE]); rownames(mrna_mat) <- toupper(trimws(mrna_matched$gene_id_raw))
  mrna_mat <- deduplicate_matrix_rows(mrna_mat)

  sample_annotation <- read_run_artifact(ctx, file.path("processed", "stage05", "sample_annotation_with_purity.tsv"), "05_define_eligible_features",
                                          function(p) read.delim(p, stringsAsFactors = FALSE))
  eligible_cancers <- sort(read_run_artifact(ctx, file.path("processed", "stage05", "eligible_cancers_cpe.tsv"), "05_define_eligible_features",
                                              function(p) read.delim(p, stringsAsFactors = FALSE))$cancer)

  cancer_group_mirna <- sample_annotation$cancer[match(colnames(mirna_mat), sample_annotation$sample_id)]
  cancer_group_mrna <- sample_annotation$cancer[match(colnames(mrna_mat), sample_annotation$sample_id)]
  mirna_mat <- winsorize_matrix_rows_by_group(mirna_mat, cancer_group_mirna)
  mrna_mat <- winsorize_matrix_rows_by_group(mrna_mat, cancer_group_mrna)

  immune_est <- read_run_artifact(ctx, file.path("processed", "stage11_immune", "sample_level_immune_estimates.tsv"), "12_immune_composition_analysis", fread)
  run_stage_assert(ctx, log_fn, all(populations %in% names(immune_est)),
                    "Every MCPcounter population in this run's own signature file has a matching column in Stage 12's sample-level immune estimates")

  # Stage 23 output, read for the historical Stage 23 reconciliation only;
  # it does not define the Stage 24 result.
  stage23_summary <- read_run_artifact(ctx, file.path("tables", "downstream", "triplet_cell_context_summary.tsv"), "23_cancer_specific_cellular_context", fread)
  historical_stage23_supported <- stage23_summary[all_three_arms_fdr_significant == TRUE]
  log_fn("Historical Stage 23 all-three-arms-significant set (unmodified, reconciliation only): %d combinations.", nrow(historical_stage23_supported))

  n_triplets <- nrow(triplets); n_cancers <- length(eligible_cancers); n_populations <- length(populations)
  log_fn("Corrected cancer-specific cellular-context analysis: %d primary triplets x up to %d eligible cancers x %d MCPcounter populations x 3 arms, HC3 inference.",
         n_triplets, n_cancers, n_populations)

  # ==========================================================================
  # PART 1: arm-level inference for every combination.
  # ==========================================================================
  arm_rows <- vector("list", n_triplets * n_cancers * n_populations * 3)
  aidx <- 0L
  for (ti in seq_len(n_triplets)) {
    mir <- triplets$mirna_id[ti]; gene_x <- triplets$gene[ti]
    if (!(mir %in% rownames(mirna_mat)) || !(gene_x %in% rownames(mrna_mat))) next
    for (cc in eligible_cancers) {
      idx <- which(sample_annotation$cancer == cc & !is.na(sample_annotation$CPE) & !is.na(sample_annotation$CYT))
      idx_samples <- sample_annotation$sample_id[idx]
      if (length(idx) < min_n) {
        for (pop in populations) for (arm in c("mirna_cyt", "gene_cyt", "mirna_gene")) {
          aidx <- aidx + 1L
          overlap_cat <- classify_overlap_category(gene_x, pop, signature_map, cyt_overlap_pops)
          arm_rows[[aidx]] <- data.table(mirna_id = mir, gene = gene_x, cancer = cc, population = pop, arm = arm,
                                          overlap_category = overlap_cat, n = length(idx), estimable = FALSE, non_evaluable_reason = "insufficient_n",
                                          beta_main = NA_real_, beta_interaction = NA_real_, se_interaction_hc3 = NA_real_,
                                          p_interaction_hc3 = NA_real_, ci_lower_hc3 = NA_real_, ci_upper_hc3 = NA_real_,
                                          beta_interaction_ols = NA_real_, se_interaction_ols = NA_real_, p_interaction_ols = NA_real_,
                                          slope_low = NA_real_, slope_high = NA_real_, df_residual = NA_integer_,
                                          max_hat = NA_real_, hc3_warning = FALSE, hc3_warning_message = "", leverage_unstable = FALSE)
        }
        next
      }
      purity_vec <- sample_annotation$CPE[idx]
      mir_v <- as.numeric(mirna_mat[mir, idx]); gene_v <- as.numeric(mrna_mat[gene_x, idx]); cyt_v <- sample_annotation$CYT[idx]
      immune_match <- match(idx_samples, immune_est$sample_id)

      for (pop in populations) {
        pop_vec <- immune_est[[pop]][immune_match]
        overlap_cat <- classify_overlap_category(gene_x, pop, signature_map, cyt_overlap_pops)
        arm_pairs <- list(mirna_cyt = list(predictor = mir_v, outcome = cyt_v),
                           gene_cyt = list(predictor = gene_v, outcome = cyt_v),
                           mirna_gene = list(predictor = mir_v, outcome = gene_v))
        for (arm in names(arm_pairs)) {
          aidx <- aidx + 1L
          if (identical(overlap_cat, "direct_tautology") && arm != "mirna_cyt") {
            arm_rows[[aidx]] <- data.table(mirna_id = mir, gene = gene_x, cancer = cc, population = pop, arm = arm,
                                            overlap_category = overlap_cat, n = length(idx), estimable = FALSE, non_evaluable_reason = "direct_tautology",
                                            beta_main = NA_real_, beta_interaction = NA_real_, se_interaction_hc3 = NA_real_,
                                            p_interaction_hc3 = NA_real_, ci_lower_hc3 = NA_real_, ci_upper_hc3 = NA_real_,
                                            beta_interaction_ols = NA_real_, se_interaction_ols = NA_real_, p_interaction_ols = NA_real_,
                                            slope_low = NA_real_, slope_high = NA_real_, df_residual = NA_integer_,
                                            max_hat = NA_real_, hc3_warning = FALSE, hc3_warning_message = "", leverage_unstable = FALSE)
            next
          }
          pred <- arm_pairs[[arm]]$predictor; out <- arm_pairs[[arm]]$outcome
          fit <- fit_hc3_checked(pred, out, pop_vec, purity_vec, min_n = min_n)
          arm_rows[[aidx]] <- data.table(mirna_id = mir, gene = gene_x, cancer = cc, population = pop, arm = arm,
                                          overlap_category = overlap_cat, n = fit$n, estimable = fit$estimable, non_evaluable_reason = fit$non_evaluable_reason,
                                          beta_main = fit$beta_main, beta_interaction = fit$beta_interaction, se_interaction_hc3 = fit$se_interaction_hc3,
                                          p_interaction_hc3 = fit$p_interaction_hc3, ci_lower_hc3 = fit$ci_lower_hc3, ci_upper_hc3 = fit$ci_upper_hc3,
                                          beta_interaction_ols = fit$beta_interaction_ols, se_interaction_ols = fit$se_interaction_ols, p_interaction_ols = fit$p_interaction_ols,
                                          slope_low = fit$slope_low, slope_high = fit$slope_high, df_residual = fit$df_residual,
                                          max_hat = fit$max_hat, hc3_warning = fit$hc3_warning,
                                          hc3_warning_message = fit$hc3_warning_message,
                                          leverage_unstable = fit$leverage_unstable)
        }
      }
    }
    if (ti %% 5 == 0) log_fn("  Corrected arm-level fits: processed %d/%d triplets.", ti, n_triplets)
  }
  arm_dt <- rbindlist(arm_rows[seq_len(aidx)])
  run_stage_assert(ctx, log_fn, nrow(arm_dt) > 0, "At least one corrected arm-level row was produced")

  # Arm-specific BH: one global family per arm across all eligible
  # (triplet, cancer, population) tests for that arm. direct_tautology arms
  # are not evaluable and are excluded; CYT-overlap arms are included,
  # because the overlap label is an interpretive distinction rather than an
  # eligibility criterion.
  arm_dt[, fdr_hc3 := NA_real_]
  for (a in c("mirna_cyt", "gene_cyt", "mirna_gene")) {
    fam <- which(arm_dt$arm == a & arm_dt$estimable)
    if (length(fam) > 0) arm_dt$fdr_hc3[fam] <- p.adjust(arm_dt$p_interaction_hc3[fam], method = "BH")
  }
  log_fn("Arm-specific global BH families: %s", paste(sprintf("%s: %d evaluable, %d FDR<%.2f", c("mirna_cyt", "gene_cyt", "mirna_gene"),
         sapply(c("mirna_cyt", "gene_cyt", "mirna_gene"), function(a) sum(arm_dt$arm == a & arm_dt$estimable)),
         sapply(c("mirna_cyt", "gene_cyt", "mirna_gene"), function(a) sum(arm_dt$arm == a & arm_dt$fdr_hc3 < conj_alpha, na.rm = TRUE)), conj_alpha), collapse = "; "))

  write_run_artifact(ctx, arm_dt, file.path("tables", "downstream", "stage24_corrected_arm_level.tsv"),
                      "24_cancer_specific_cellular_context_corrected", script_path,
                      file.path(ctx$run_root, c(file.path("tables", "primary", run_primary_triplet_filename(ctx)),
                                                 file.path("processed", "stage03", "mirna_matched.rds"),
                                                 file.path("processed", "stage03", "mrna_matched.rds"),
                                                 file.path("processed", "stage05", "sample_annotation_with_purity.tsv"),
                                                 file.path("processed", "stage11_immune", "sample_level_immune_estimates.tsv"))),
                      fwrite_fn)
  log_fn("stage24_corrected_arm_level.tsv written (%d rows).", nrow(arm_dt))

  # ==========================================================================
  # PART 2: conjunction test.
  # ==========================================================================
  keys <- unique(arm_dt[, .(mirna_id, gene, cancer, population)])
  wide <- dcast(arm_dt, mirna_id + gene + cancer + population ~ arm,
                value.var = c("p_interaction_hc3", "estimable", "overlap_category", "n", "leverage_unstable"))
  wide[, overlap_category := overlap_category_mirna_cyt]  # identical across arms by construction (gene/population-level label)
  wide[, p_conjunction := mapply(compute_conjunction_p, p_interaction_hc3_mirna_cyt, p_interaction_hc3_gene_cyt, p_interaction_hc3_mirna_gene)]
  wide[, all_three_estimable := estimable_mirna_cyt & estimable_gene_cyt & estimable_mirna_gene]
  wide[, any_arm_leverage_unstable := leverage_unstable_mirna_cyt |
         leverage_unstable_gene_cyt | leverage_unstable_mirna_gene]
  wide[, n_samples := pmax(n_mirna_cyt, n_gene_cyt, n_mirna_gene, na.rm = TRUE)]

  conj_family <- which(wide$all_three_estimable & !wide$any_arm_leverage_unstable &
                         is.finite(wide$p_conjunction) & wide$overlap_category != "direct_tautology")
  wide[, conjunction_FDR_BH := NA_real_]; wide[, conjunction_FDR_BY := NA_real_]
  if (length(conj_family) > 0) {
    wide$conjunction_FDR_BH[conj_family] <- p.adjust(wide$p_conjunction[conj_family], method = "BH")
    wide$conjunction_FDR_BY[conj_family] <- p.adjust(wide$p_conjunction[conj_family], method = "BY")
  }
  wide[, all_three_arms_supported := isTRUE(FALSE)]
  wide[conj_family, all_three_arms_supported := conjunction_FDR_BH < conj_alpha]

  # Comparison column only: three separately
  # BH-adjusted arm-level p-values (arm_dt$fdr_hc3), each < alpha, as in the
  # arm-wise rule of Stage 23. Not used as a result of this stage.
  fdr_wide <- dcast(arm_dt, mirna_id + gene + cancer + population ~ arm, value.var = "fdr_hc3")
  setnames(fdr_wide, c("mirna_cyt", "gene_cyt", "mirna_gene"), c("old_rule_fdr_mirna_cyt", "old_rule_fdr_gene_cyt", "old_rule_fdr_mirna_gene"))
  wide <- merge(wide, fdr_wide, by = c("mirna_id", "gene", "cancer", "population"))
  wide[, old_rule_all_three_significant := isTRUE(FALSE)]
  wide[!is.na(old_rule_fdr_mirna_cyt) & !is.na(old_rule_fdr_gene_cyt) & !is.na(old_rule_fdr_mirna_gene),
       old_rule_all_three_significant := old_rule_fdr_mirna_cyt < conj_alpha & old_rule_fdr_gene_cyt < conj_alpha & old_rule_fdr_mirna_gene < conj_alpha]

  conjunction_dt <- wide[, .(mirna_id, gene, cancer, population, overlap_category, n_samples,
                              p_mirna_cyt = p_interaction_hc3_mirna_cyt, p_gene_cyt = p_interaction_hc3_gene_cyt, p_mirna_gene = p_interaction_hc3_mirna_gene,
                              any_arm_leverage_unstable, p_conjunction, conjunction_FDR_BH, conjunction_FDR_BY, all_three_arms_supported,
                              old_rule_all_three_significant)]
  n_conj_supported <- sum(conjunction_dt$all_three_arms_supported, na.rm = TRUE)
  log_fn("Conjunction test: %d evaluable, non-tautological (triplet,cancer,population) combinations tested; %d all_three_arms_supported at conjunction FDR(BH)<%.2f (vs. %d under the OLD three-separately-adjusted-arms rule, historical comparison only).",
         length(conj_family), n_conj_supported, conj_alpha, sum(conjunction_dt$old_rule_all_three_significant, na.rm = TRUE))

  write_run_artifact(ctx, conjunction_dt, file.path("tables", "downstream", "stage24_corrected_conjunction.tsv"),
                      "24_cancer_specific_cellular_context_corrected", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "stage24_corrected_arm_level.tsv"), fwrite_fn)
  log_fn("stage24_corrected_conjunction.tsv written (%d rows).", nrow(conjunction_dt))

  # ==========================================================================
  # PART 3: non-overlapping vs. CYT-overlap population summaries.
  # ==========================================================================
  nonoverlap_dt <- conjunction_dt[overlap_category == "nonoverlapping_population"]
  cytoverlap_dt <- conjunction_dt[overlap_category == "CYT_biological_overlap"]
  write_run_artifact(ctx, nonoverlap_dt, file.path("tables", "downstream", "stage24_nonoverlapping_population_results.tsv"),
                      "24_cancer_specific_cellular_context_corrected", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "stage24_corrected_conjunction.tsv"), fwrite_fn)
  write_run_artifact(ctx, cytoverlap_dt, file.path("tables", "downstream", "stage24_cyt_overlap_results.tsv"),
                      "24_cancer_specific_cellular_context_corrected", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "stage24_corrected_conjunction.tsv"), fwrite_fn)
  log_fn("Nonoverlapping-population combinations: %d (%d conjunction-supported). CYT-overlap combinations: %d (%d conjunction-supported).",
         nrow(nonoverlap_dt), sum(nonoverlap_dt$all_three_arms_supported, na.rm = TRUE),
         nrow(cytoverlap_dt), sum(cytoverlap_dt$all_three_arms_supported, na.rm = TRUE))

  # ==========================================================================
  # PART 4: population correlations, every cancer.
  # ==========================================================================
  corr_rows <- vector("list", n_cancers)
  for (i in seq_along(eligible_cancers)) {
    cc <- eligible_cancers[i]
    idx <- which(sample_annotation$cancer == cc)
    idx_samples <- sample_annotation$sample_id[idx]
    sub <- immune_est[match(idx_samples, immune_est$sample_id), ..populations]
    sub <- sub[stats::complete.cases(sub)]
    if (nrow(sub) < 10) { corr_rows[[i]] <- data.table(); next }
    cmat <- suppressWarnings(cor(sub, method = "spearman", use = "pairwise.complete.obs"))
    flagged <- flag_correlated_population_pairs(cmat, threshold = corr_threshold)
    cmat[upper.tri(cmat, diag = TRUE)] <- NA
    all_pairs <- which(!is.na(cmat), arr.ind = TRUE)
    corr_rows[[i]] <- data.table(cancer = cc, population_a = rownames(cmat)[all_pairs[, 1]], population_b = colnames(cmat)[all_pairs[, 2]],
                                  rho = cmat[all_pairs], n = nrow(sub),
                                  flagged_high_correlation = abs(cmat[all_pairs]) > corr_threshold)
  }
  correlation_dt <- rbindlist(corr_rows[lengths(corr_rows) > 0])
  write_run_artifact(ctx, correlation_dt, file.path("tables", "downstream", "stage24_population_correlations.tsv"),
                      "24_cancer_specific_cellular_context_corrected", script_path,
                      file.path(ctx$run_root, "processed", "stage11_immune", "sample_level_immune_estimates.tsv"), fwrite_fn)
  log_fn("stage24_population_correlations.tsv written (%d pairs across %d cancers; %d pairs exceed |rho|>%.2f).",
         nrow(correlation_dt), n_cancers, sum(correlation_dt$flagged_high_correlation), corr_threshold)

  # ==========================================================================
  # PART 5: collinearity sensitivity for conjunction-supported findings
  # only.
  # ==========================================================================
  supported <- conjunction_dt[all_three_arms_supported == TRUE]
  collin_rows <- vector("list", nrow(supported))
  if (nrow(supported) > 0) {
    for (i in seq_len(nrow(supported))) {
      s <- supported[i]
      cc_corr <- correlation_dt[cancer == s$cancer & flagged_high_correlation == TRUE &
                                 (population_a == s$population | population_b == s$population)]
      correlated_pops <- unique(c(cc_corr$population_a, cc_corr$population_b))
      correlated_pops <- setdiff(correlated_pops, s$population)
      idx <- which(sample_annotation$cancer == s$cancer & !is.na(sample_annotation$CPE) & !is.na(sample_annotation$CYT))
      idx_samples <- sample_annotation$sample_id[idx]
      purity_vec <- sample_annotation$CPE[idx]
      mir_v <- as.numeric(mirna_mat[s$mirna_id, idx]); gene_v <- as.numeric(mrna_mat[s$gene, idx])
      immune_match <- match(idx_samples, immune_est$sample_id)
      r_out <- rank(gene_v); r_pred <- rank(mir_v)  # mirna_gene arm as the representative refit
      pop_focal_z <- as.numeric(scale(immune_est[[s$population]][immune_match]))

      if (length(correlated_pops) == 0) {
        collin_rows[[i]] <- data.table(mirna_id = s$mirna_id, gene = s$gene, cancer = s$cancer, population = s$population,
                                        correlated_populations = "", max_vif = NA_real_, condition_number = NA_real_,
                                        collinearity_assessment = "not_applicable_no_correlated_population",
                                        coefficient_stable = NA, focal_beta_alone = NA_real_, focal_beta_in_joint_model = NA_real_,
                                        pc1_lymphocyte_beta = NA_real_, pc1_lymphocyte_p = NA_real_, focal_adds_beyond_pc1 = NA)
        next
      }
      pop_mat <- sapply(correlated_pops, function(p) immune_est[[p]][immune_match])
      dat <- data.frame(r_out = r_out, r_pred = r_pred, pop_focal_z = pop_focal_z, purity = purity_vec)
      for (j in seq_along(correlated_pops)) dat[[sprintf("pop_corr_%d", j)]] <- as.numeric(scale(pop_mat[, j]))
      corr_terms <- paste(sprintf("pop_corr_%d", seq_along(correlated_pops)), collapse = " + ")
      f_joint <- as.formula(sprintf("r_out ~ r_pred * pop_focal_z + %s + purity", corr_terms))
      fit_joint <- tryCatch(lm(f_joint, data = dat), error = function(e) NULL)
      f_alone <- as.formula("r_out ~ r_pred * pop_focal_z + purity")
      fit_alone <- tryCatch(lm(f_alone, data = dat), error = function(e) NULL)

      if (is.null(fit_joint) || is.null(fit_alone) || !("r_pred:pop_focal_z" %in% names(coef(fit_joint)))) {
        collin_rows[[i]] <- data.table(mirna_id = s$mirna_id, gene = s$gene, cancer = s$cancer, population = s$population,
                                        correlated_populations = paste(correlated_pops, collapse = ";"), max_vif = NA_real_, condition_number = NA_real_,
                                        collinearity_assessment = "not_evaluable",
                                        coefficient_stable = NA, focal_beta_alone = NA_real_, focal_beta_in_joint_model = NA_real_,
                                        pc1_lymphocyte_beta = NA_real_, pc1_lymphocyte_p = NA_real_, focal_adds_beyond_pc1 = NA)
        next
      }
      vifs <- tryCatch(compute_vif(fit_joint), error = function(e) c(NA_real_))
      cond_num <- tryCatch(compute_condition_number(fit_joint), error = function(e) NA_real_)
      beta_alone <- coef(fit_alone)["r_pred:pop_focal_z"]
      beta_joint <- coef(fit_joint)["r_pred:pop_focal_z"]
      rel_change <- abs(beta_joint - beta_alone) / max(abs(beta_alone), .Machine$double.eps)
      coef_stable <- isTRUE(rel_change < 0.5) && isTRUE(max(vifs, na.rm = TRUE) < 10)

      pc1_p <- NA_real_; pc1_beta <- NA_real_; adds_beyond <- NA
      if (length(correlated_pops) >= 1) {
        pc1_input <- cbind(matrix(pop_focal_z, ncol = 1, dimnames = list(NULL, s$population)), pop_mat)
        if (ncol(pc1_input) >= 2 && all(apply(pc1_input, 2, sd) > 0)) {
          pc1 <- tryCatch(compute_lymphocyte_pc1(pc1_input), error = function(e) NULL)
          if (!is.null(pc1)) {
            dat$pc1 <- pc1
            fit_pc1 <- tryCatch(lm(r_out ~ r_pred * pop_focal_z + pc1 + purity, data = dat), error = function(e) NULL)
            if (!is.null(fit_pc1) && "r_pred:pop_focal_z" %in% names(coef(fit_pc1))) {
              vc_pc1 <- tryCatch(suppressWarnings(sandwich::vcovHC(fit_pc1, type = "HC3")), error = function(e) NULL)
              if (!is.null(vc_pc1) && "r_pred:pop_focal_z" %in% rownames(vc_pc1)) {
                pc1_beta <- coef(fit_pc1)["r_pred:pop_focal_z"]
                se_pc1 <- sqrt(vc_pc1["r_pred:pop_focal_z", "r_pred:pop_focal_z"])
                pc1_p <- 2 * pt(-abs(pc1_beta / se_pc1), df = fit_pc1$df.residual)
                adds_beyond <- isTRUE(pc1_p < conj_alpha)
              }
            }
          }
        }
      }
      collin_rows[[i]] <- data.table(mirna_id = s$mirna_id, gene = s$gene, cancer = s$cancer, population = s$population,
                                      correlated_populations = paste(correlated_pops, collapse = ";"),
                                      max_vif = max(vifs, na.rm = TRUE), condition_number = cond_num,
                                      collinearity_assessment = if (coef_stable) "stable" else "unstable",
                                      coefficient_stable = coef_stable, focal_beta_alone = unname(beta_alone), focal_beta_in_joint_model = unname(beta_joint),
                                      pc1_lymphocyte_beta = unname(pc1_beta), pc1_lymphocyte_p = pc1_p, focal_adds_beyond_pc1 = adds_beyond)
    }
  }
  collinearity_dt <- if (length(collin_rows) > 0) rbindlist(collin_rows) else data.table()
  write_run_artifact(ctx, collinearity_dt, file.path("tables", "downstream", "stage24_collinearity_sensitivity.tsv"),
                      "24_cancer_specific_cellular_context_corrected", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "stage24_corrected_conjunction.tsv"), fwrite_fn)
  log_fn("stage24_collinearity_sensitivity.tsv written (%d conjunction-supported combinations checked; %d have an unstable focal-population coefficient).",
         nrow(collinearity_dt), sum(!collinearity_dt$coefficient_stable, na.rm = TRUE))

  # ==========================================================================
  # PART 6: permutation negative controls on a documented, tractable subset:
  # every conjunction-supported combination (n_conj_supported) plus a
  # fixed-seed random sample of non-supported combinations, permuting the
  # population score within cancer (X, Y, purity and cancer membership
  # preserved), config$n_permutations draws per combination. A second
  # control replaces the population score with a matched random gene
  # expression score (similar abundance and variance). Synthetic
  # known-answer scenarios (no/one/two/three-arm interaction,
  # heteroscedastic errors, correlated arms and populations, missing
  # values) are tested in tests/test_cell_context_conjunction_helpers.R.
  # ==========================================================================
  set.seed(cfg$permutation_seed)
  n_perm <- cfg$n_permutations
  nc_candidates <- conjunction_dt[!is.na(p_conjunction)]
  nc_supported <- nc_candidates[all_three_arms_supported == TRUE]
  n_nc_other <- min(30, nrow(nc_candidates) - nrow(nc_supported))
  nc_other <- if (n_nc_other > 0) nc_candidates[all_three_arms_supported == FALSE][sample(.N, n_nc_other)] else nc_candidates[0]
  nc_targets <- rbind(nc_supported, nc_other)

  perm_rows <- vector("list", nrow(nc_targets))
  if (nrow(nc_targets) > 0) {
    for (i in seq_len(nrow(nc_targets))) {
      s <- nc_targets[i]
      idx <- which(sample_annotation$cancer == s$cancer & !is.na(sample_annotation$CPE) & !is.na(sample_annotation$CYT))
      idx_samples <- sample_annotation$sample_id[idx]
      purity_vec <- sample_annotation$CPE[idx]
      mir_v <- as.numeric(mirna_mat[s$mirna_id, idx]); gene_v <- as.numeric(mrna_mat[s$gene, idx]); cyt_v <- sample_annotation$CYT[idx]
      immune_match <- match(idx_samples, immune_est$sample_id)
      pop_vec <- immune_est[[s$population]][immune_match]

      perm_p <- numeric(n_perm)
      for (b in seq_len(n_perm)) {
        pop_perm <- sample(pop_vec)
        f1 <- fit_cell_context_interaction_hc3(mir_v, cyt_v, pop_perm, purity_vec, min_n = min_n)
        f2 <- fit_cell_context_interaction_hc3(gene_v, cyt_v, pop_perm, purity_vec, min_n = min_n)
        f3 <- fit_cell_context_interaction_hc3(mir_v, gene_v, pop_perm, purity_vec, min_n = min_n)
        perm_p[b] <- if (isTRUE(f1$estimable) && isTRUE(f2$estimable) && isTRUE(f3$estimable)) {
          compute_conjunction_p(f1$p_interaction_hc3, f2$p_interaction_hc3, f3$p_interaction_hc3)
        } else NA_real_
      }
      n_valid <- sum(is.finite(perm_p))
      n_le_observed <- sum(perm_p <= s$p_conjunction, na.rm = TRUE)
      emp_p <- if (n_valid > 0) (n_le_observed + 1) / (n_valid + 1) else NA_real_
      mc_se <- if (n_valid > 0) sqrt(emp_p * (1 - emp_p) / n_valid) else NA_real_

      # Matched random gene control: a real expression row with abundance
      # and variance similar to the triplet's gene, used in place of the
      # population score (purity, X, Y and cancer membership preserved).
      candidate_genes <- setdiff(rownames(mrna_mat), s$gene)
      gene_var <- apply(mrna_mat[candidate_genes, idx, drop = FALSE], 1, sd)
      target_var <- sd(gene_v)
      matched_gene <- candidate_genes[which.min(abs(gene_var - target_var))]
      matched_vec <- as.numeric(mrna_mat[matched_gene, idx])
      f1m <- fit_cell_context_interaction_hc3(mir_v, cyt_v, matched_vec, purity_vec, min_n = min_n)
      f2m <- fit_cell_context_interaction_hc3(gene_v, cyt_v, matched_vec, purity_vec, min_n = min_n)
      f3m <- fit_cell_context_interaction_hc3(mir_v, gene_v, matched_vec, purity_vec, min_n = min_n)
      matched_conj_p <- if (isTRUE(f1m$estimable) && isTRUE(f2m$estimable) && isTRUE(f3m$estimable)) {
        compute_conjunction_p(f1m$p_interaction_hc3, f2m$p_interaction_hc3, f3m$p_interaction_hc3)
      } else NA_real_

      perm_rows[[i]] <- data.table(mirna_id = s$mirna_id, gene = s$gene, cancer = s$cancer, population = s$population,
                                    observed_p_conjunction = s$p_conjunction, all_three_arms_supported = s$all_three_arms_supported,
                                    n_permutations = n_perm, n_valid_permutations = n_valid,
                                    empirical_p_permutation = emp_p, monte_carlo_se = mc_se,
                                    matched_gene_control = matched_gene, matched_gene_conjunction_p = matched_conj_p,
                                    passes_permutation_criterion = isTRUE(emp_p < conj_alpha))
    }
  }
  negctrl_dt <- if (length(perm_rows) > 0) rbindlist(perm_rows) else data.table()
  write_run_artifact(ctx, negctrl_dt, file.path("tables", "downstream", "stage24_negative_controls.tsv"),
                      "24_cancer_specific_cellular_context_corrected", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "stage24_corrected_conjunction.tsv"), fwrite_fn)
  log_fn("stage24_negative_controls.tsv written: %d combinations checked (%d conjunction-supported + %d comparison), %d permutations each (Monte Carlo SE reported per row); %d/%d conjunction-supported combinations also pass the permutation criterion (empirical p<%.2f).",
         nrow(negctrl_dt), nrow(nc_supported), nrow(nc_other), n_perm,
         sum(negctrl_dt$all_three_arms_supported & negctrl_dt$passes_permutation_criterion, na.rm = TRUE), nrow(nc_supported), conj_alpha)

  # ==========================================================================
  # PART 7: power and sample-size context, every cancer.
  # ==========================================================================
  power_rows <- vector("list", n_cancers)
  n_sim <- cfg$power_n_simulations
  effect_grid <- cfg$power_effect_grid
  for (i in seq_along(eligible_cancers)) {
    cc <- eligible_cancers[i]
    idx_all <- which(sample_annotation$cancer == cc)
    idx_cc <- which(sample_annotation$cancer == cc & !is.na(sample_annotation$CPE) & !is.na(sample_annotation$CYT))
    idx_samples <- sample_annotation$sample_id[idx_cc]
    immune_match <- match(idx_samples, immune_est$sample_id)
    pop_sds <- sapply(populations, function(p) sd(immune_est[[p]][immune_match], na.rm = TRUE))
    pop_ranges <- sapply(populations, function(p) diff(range(immune_est[[p]][immune_match], na.rm = TRUE)))
    n_this_cancer_rows <- arm_dt[cancer == cc]
    n_not_eval <- sum(!n_this_cancer_rows$estimable)
    n_total_rows <- nrow(n_this_cancer_rows)

    # Simulation-based minimum-detectable interaction effect: for each grid
    # value, simulate n_sim rank-interaction datasets at this cancer's own
    # n (complete-case) and purity distribution, fit the corrected HC3
    # model, and record the fraction reaching p<0.05 (power).
    n_cc <- length(idx_cc)
    power_by_effect <- rep(NA_real_, length(effect_grid))
    if (n_cc >= min_n) {
      purity_sample <- sample_annotation$CPE[idx_cc]
      for (e in seq_along(effect_grid)) {
        eff <- effect_grid[e]
        n_sig <- 0
        for (s in seq_len(n_sim)) {
          x_sim <- rnorm(n_cc); pop_sim <- rnorm(n_cc)
          y_sim <- x_sim * (1 + eff * scale(pop_sim)) + rnorm(n_cc, sd = 1)
          f_sim <- fit_cell_context_interaction_hc3(x_sim, y_sim, pop_sim, purity_sample, min_n = min_n)
          if (isTRUE(f_sim$estimable) && is.finite(f_sim$p_interaction_hc3) && f_sim$p_interaction_hc3 < 0.05) n_sig <- n_sig + 1
        }
        power_by_effect[e] <- n_sig / n_sim
      }
    }
    min_detectable_idx <- which(power_by_effect >= 0.80)
    min_detectable_effect <- if (length(min_detectable_idx) > 0) effect_grid[min(min_detectable_idx)] else NA_real_

    power_rows[[i]] <- data.table(cancer = cc, n_total = length(idx_all), n_complete_case = n_cc,
                                   median_population_sd = median(pop_sds, na.rm = TRUE), median_population_range = median(pop_ranges, na.rm = TRUE),
                                   n_arm_level_rows = n_total_rows, n_not_evaluable = n_not_eval, pct_not_evaluable = 100 * n_not_eval / max(1, n_total_rows),
                                   min_detectable_interaction_effect_80pct_power = min_detectable_effect,
                                   power_at_grid = paste(sprintf("%.2f:%.2f", effect_grid, power_by_effect), collapse = ";"),
                                   n_all_three_arms_supported_here = sum(conjunction_dt$cancer == cc & conjunction_dt$all_three_arms_supported, na.rm = TRUE))
  }
  power_dt <- rbindlist(power_rows)
  cor_size_findings <- suppressWarnings(cor(power_dt$n_complete_case, power_dt$n_all_three_arms_supported_here, method = "spearman"))
  cor_var_findings <- suppressWarnings(cor(power_dt$median_population_sd, power_dt$n_all_three_arms_supported_here, method = "spearman"))
  write_run_artifact(ctx, power_dt, file.path("tables", "downstream", "stage24_power_by_cancer.tsv"),
                      "24_cancer_specific_cellular_context_corrected", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "stage24_corrected_conjunction.tsv"), fwrite_fn)
  log_fn("stage24_power_by_cancer.tsv written (%d cancers). Spearman(cohort size, #conjunction-supported findings)=%.3f; Spearman(median population-score SD, #findings)=%.3f -- %s explain the concentration of findings by size/variance alone (see report for full interpretation).",
         n_cancers, cor_size_findings, cor_var_findings,
         if (is.finite(cor_size_findings) && abs(cor_size_findings) > 0.5) "these DO appear to partly" else "these do NOT strongly")

  # ==========================================================================
  # PART 8: nonlinearity sensitivity for conjunction-supported candidates
  # only (descriptive; tertiles for interpretation only).
  # ==========================================================================
  spline_rows <- vector("list", nrow(supported))
  if (nrow(supported) > 0) {
    for (i in seq_len(nrow(supported))) {
      s <- supported[i]
      idx <- which(sample_annotation$cancer == s$cancer & !is.na(sample_annotation$CPE) & !is.na(sample_annotation$CYT))
      idx_samples <- sample_annotation$sample_id[idx]
      purity_vec <- sample_annotation$CPE[idx]
      mir_v <- as.numeric(mirna_mat[s$mirna_id, idx]); gene_v <- as.numeric(mrna_mat[s$gene, idx])
      immune_match <- match(idx_samples, immune_est$sample_id)
      pop_z <- as.numeric(scale(immune_est[[s$population]][immune_match]))
      r_pred <- rank(mir_v); r_out <- rank(gene_v)
      dat <- data.frame(r_out = r_out, r_pred = r_pred, pop_z = pop_z, purity = purity_vec)
      fit_lin <- tryCatch(lm(r_out ~ r_pred * pop_z + purity, data = dat), error = function(e) NULL)
      fit_spline <- tryCatch(lm(r_out ~ r_pred * splines::ns(pop_z, df = 3) + purity, data = dat), error = function(e) NULL)
      spline_p <- NA_real_; agrees_with_tertile <- NA
      if (!is.null(fit_lin) && !is.null(fit_spline)) {
        av <- tryCatch(anova(fit_lin, fit_spline), error = function(e) NULL)
        if (!is.null(av) && nrow(av) == 2) spline_p <- av[["Pr(>F)"]][2]
      }
      strat <- tertile_stratified_correlation(mir_v, gene_v, immune_est[[s$population]][immune_match], min_n_stratum = analysis_config$cell_context$min_samples_per_tertile)
      low_rho <- strat$rho[strat$stratum == "low"]; high_rho <- strat$rho[strat$stratum == "high"]
      tertile_monotonic <- isTRUE(is.finite(low_rho) && is.finite(high_rho) && sign(high_rho - low_rho) == sign(s$p_mirna_gene * 0 + 1))
      spline_rows[[i]] <- data.table(mirna_id = s$mirna_id, gene = s$gene, cancer = s$cancer, population = s$population,
                                      nonlinear_vs_linear_p = spline_p, nonlinearity_flagged = isTRUE(spline_p < 0.05),
                                      tertile_low_rho = low_rho, tertile_high_rho = high_rho)
    }
  }
  nonlinearity_dt <- if (length(spline_rows) > 0) rbindlist(spline_rows) else data.table()
  if (nrow(nonlinearity_dt) > 0) {
    write_run_artifact(ctx, nonlinearity_dt, file.path("tables", "downstream", "stage24_nonlinearity_sensitivity.tsv"),
                        "24_cancer_specific_cellular_context_corrected", script_path,
                        file.path(ctx$run_root, "tables", "downstream", "stage24_corrected_conjunction.tsv"), fwrite_fn)
  }
  log_fn("Nonlinearity sensitivity (spline, conjunction-supported candidates only): %d checked, %d flagged as significantly nonlinear (p<0.05, exploratory only).",
         nrow(nonlinearity_dt), sum(nonlinearity_dt$nonlinearity_flagged, na.rm = TRUE))

  # ==========================================================================
  # PART 9: LUMP purity sensitivity for the complete three-arm conjunction,
  # plus a completeness-gated AJCC-stage sensitivity for the miRNA-gene arm.
  # ==========================================================================
  lump_rows <- vector("list", nrow(supported))
  surv <- tryCatch(read_run_artifact(ctx, file.path("processed", "stage03", "survival_matched.rds"), "03_harmonize_identifiers", readRDS), error = function(e) NULL)
  stage_complete_cancers <- character(0)
  if (!is.null(surv)) {
    surv_dt <- as.data.table(surv)
    m_stage <- as.data.table(merge(sample_annotation[, c("sample_id", "cancer")], surv_dt[, c("sample_id", "ajcc_pathologic_tumor_stage")], by = "sample_id"))
    invalid_stage <- c("", "[Not Available]", "[Not Applicable]", "[Unknown]", "[Discrepancy]")
    completeness <- m_stage[, .(pct_complete = 100 * sum(!is.na(ajcc_pathologic_tumor_stage) & !(ajcc_pathologic_tumor_stage %in% invalid_stage)) / .N), by = cancer]
    stage_complete_cancers <- completeness[pct_complete >= 80, cancer]
  }
  if (nrow(supported) > 0) {
    for (i in seq_len(nrow(supported))) {
      s <- supported[i]
      idx_lump <- which(sample_annotation$cancer == s$cancer & !is.na(sample_annotation$LUMP) & !is.na(sample_annotation$CYT))
      lump_n <- length(idx_lump)
      lump_fits <- setNames(vector("list", 3), c("mirna_cyt", "gene_cyt", "mirna_gene"))
      if (lump_n >= min_n) {
        idx_samples <- sample_annotation$sample_id[idx_lump]
        purity_lump <- sample_annotation$LUMP[idx_lump]
        mir_v <- as.numeric(mirna_mat[s$mirna_id, idx_lump]); gene_v <- as.numeric(mrna_mat[s$gene, idx_lump])
        immune_match <- match(idx_samples, immune_est$sample_id)
        pop_v <- immune_est[[s$population]][immune_match]
        cyt_v <- sample_annotation$CYT[idx_lump]
        lump_fits$mirna_cyt <- fit_hc3_checked(mir_v, cyt_v, pop_v, purity_lump, min_n)
        lump_fits$gene_cyt <- fit_hc3_checked(gene_v, cyt_v, pop_v, purity_lump, min_n)
        lump_fits$mirna_gene <- fit_hc3_checked(mir_v, gene_v, pop_v, purity_lump, min_n)
      }
      stage_beta <- NA_real_; stage_p <- NA_real_; stage_applicable <- s$cancer %in% stage_complete_cancers
      if (stage_applicable && !is.null(surv)) {
        surv_dt <- as.data.table(surv)
        idx <- which(sample_annotation$cancer == s$cancer & !is.na(sample_annotation$CPE) & !is.na(sample_annotation$CYT))
        idx_samples <- sample_annotation$sample_id[idx]
        stage_match <- surv_dt$ajcc_pathologic_tumor_stage[match(idx_samples, surv_dt$sample_id)]
        invalid_stage <- c("", "[Not Available]", "[Not Applicable]", "[Unknown]", "[Discrepancy]")
        ok_stage <- !is.na(stage_match) & !(stage_match %in% invalid_stage)
        if (sum(ok_stage) >= min_n && length(unique(stage_match[ok_stage])) >= 2) {
          idx2 <- idx[ok_stage]
          purity_vec <- sample_annotation$CPE[idx2]
          mir_v <- as.numeric(mirna_mat[s$mirna_id, idx2]); gene_v <- as.numeric(mrna_mat[s$gene, idx2])
          immune_match2 <- match(sample_annotation$sample_id[idx2], immune_est$sample_id)
          pop_z <- as.numeric(scale(immune_est[[s$population]][immune_match2]))
          dat_stage <- data.frame(r_out = rank(gene_v), r_pred = rank(mir_v), pop_z = pop_z, purity = purity_vec, stage = factor(stage_match[ok_stage]))
          fit_stage <- tryCatch(lm(r_out ~ r_pred * pop_z + purity + stage, data = dat_stage), error = function(e) NULL)
          if (!is.null(fit_stage) && "r_pred:pop_z" %in% names(coef(fit_stage))) {
            vc_stage <- tryCatch(suppressWarnings(sandwich::vcovHC(fit_stage, type = "HC3")), error = function(e) NULL)
            if (!is.null(vc_stage) && "r_pred:pop_z" %in% rownames(vc_stage)) {
              stage_beta <- coef(fit_stage)["r_pred:pop_z"]
              se_stage <- sqrt(vc_stage["r_pred:pop_z", "r_pred:pop_z"])
              stage_p <- 2 * pt(-abs(stage_beta / se_stage), df = fit_stage$df.residual)
            }
          }
        } else stage_applicable <- FALSE
      }
      get_fit_value <- function(arm, field) {
        fit <- lump_fits[[arm]]
        if (is.null(fit) || !isTRUE(fit$estimable)) return(NA_real_)
        fit[[field]]
      }
      lump_rows[[i]] <- data.table(mirna_id = s$mirna_id, gene = s$gene, cancer = s$cancer, population = s$population,
                                    lump_n = lump_n,
                                    lump_beta_mirna_cyt = get_fit_value("mirna_cyt", "beta_interaction"),
                                    lump_beta_gene_cyt = get_fit_value("gene_cyt", "beta_interaction"),
                                    lump_beta_mirna_gene = get_fit_value("mirna_gene", "beta_interaction"),
                                    lump_p_mirna_cyt = get_fit_value("mirna_cyt", "p_interaction_hc3"),
                                    lump_p_gene_cyt = get_fit_value("gene_cyt", "p_interaction_hc3"),
                                    lump_p_mirna_gene = get_fit_value("mirna_gene", "p_interaction_hc3"),
                                    lump_any_arm_leverage_unstable = any(vapply(lump_fits, function(f) !is.null(f) && isTRUE(f$leverage_unstable), logical(1))),
                                    stage_sensitivity_applicable = stage_applicable, stage_beta_interaction = stage_beta, stage_p_interaction = stage_p)
    }
  }
  purity_stage_dt <- if (length(lump_rows) > 0) rbindlist(lump_rows) else data.table()
  if (nrow(purity_stage_dt) > 0) {
    p_cols <- c("lump_p_mirna_cyt", "lump_p_gene_cyt", "lump_p_mirna_gene")
    purity_stage_dt[, lump_p_conjunction := apply(as.matrix(.SD), 1, function(z) if (all(is.finite(z))) max(z) else NA_real_), .SDcols = p_cols]
    purity_stage_dt[lump_any_arm_leverage_unstable == TRUE, lump_p_conjunction := NA_real_]
    purity_stage_dt[, lump_conjunction_fdr_bh := p.adjust(lump_p_conjunction, method = "BH")]
    purity_stage_dt[, lump_conjunction_supported := !is.na(lump_conjunction_fdr_bh) & lump_conjunction_fdr_bh < conj_alpha]
    cpe_betas <- dcast(arm_dt, mirna_id + gene + cancer + population ~ arm, value.var = "beta_interaction")
    setnames(cpe_betas, c("mirna_cyt", "gene_cyt", "mirna_gene"), c("cpe_beta_mirna_cyt", "cpe_beta_gene_cyt", "cpe_beta_mirna_gene"))
    purity_stage_dt <- merge(purity_stage_dt, cpe_betas, by = c("mirna_id", "gene", "cancer", "population"), all.x = TRUE)
    purity_stage_dt[, lump_direction_agreement :=
      sign(lump_beta_mirna_cyt) == sign(cpe_beta_mirna_cyt) &
      sign(lump_beta_gene_cyt) == sign(cpe_beta_gene_cyt) &
      sign(lump_beta_mirna_gene) == sign(cpe_beta_mirna_gene)]
    purity_stage_dt[, lump_sensitivity_pass := lump_conjunction_supported & lump_direction_agreement & !lump_any_arm_leverage_unstable]
    write_run_artifact(ctx, purity_stage_dt, file.path("tables", "downstream", "stage24_lump_and_stage_sensitivity.tsv"),
                        "24_cancer_specific_cellular_context_corrected", script_path,
                        file.path(ctx$run_root, "tables", "downstream", "stage24_corrected_conjunction.tsv"), fwrite_fn)
  }
  log_fn("LUMP/stage sensitivity: %d conjunction-supported combinations checked with a complete three-arm LUMP conjunction; %d pass conjunction FDR and three-arm direction agreement. %d cancers have sufficient AJCC-stage completeness (>=80%%); %d combinations received the representative miRNA-gene stage refit.",
         nrow(purity_stage_dt), sum(purity_stage_dt$lump_sensitivity_pass, na.rm = TRUE), length(stage_complete_cancers), sum(purity_stage_dt$stage_sensitivity_applicable, na.rm = TRUE))

  # ==========================================================================
  # PART 10: deconvolution sensitivity, not performed. Original stage 13 (quanTIseq; omitted)
  # was permitted only when analysis_config$immune$quantiseq_input_scale_validated
  # is TRUE in the original workflow; the flag is FALSE,
  # so no alternative-deconvolution sensitivity result is produced or
  # estimated, and this is reported as a limitation.
  # ==========================================================================
  log_fn("Deconvolution sensitivity: not performed. Original stage 13 (quanTIseq; omitted) is gated by analysis_config$immune$quantiseq_input_scale_validated = FALSE; reported as a limitation.")

  classification <- run_stage24_finalize(ctx, log_fn, script_path, fwrite_fn, cfg, conj_alpha,
                                    historical_stage23_supported = historical_stage23_supported,
                                    conjunction_dt = conjunction_dt, negctrl_dt = negctrl_dt,
                                    collinearity_dt = collinearity_dt, cytoverlap_dt = cytoverlap_dt,
                                    purity_stage_dt = purity_stage_dt, nonlinearity_dt = nonlinearity_dt,
                                    quantiseq_validated = isTRUE(analysis_config$immune$quantiseq_input_scale_validated),
                                    signature_map = signature_map, triplets = triplets, populations = populations,
                                    eligible_cancers = eligible_cancers, power_dt = power_dt)

  record_stage_status(ctx, "24_cancer_specific_cellular_context_corrected", "COMPLETE")
  log_fn("Stage 24 complete.")
  invisible(conjunction_dt)
}

# ============================================================================
# Evidence funnels, final result tables, figures and the final report;
# called once at the end of run_stage24().
# ============================================================================
run_stage24_finalize <- function(ctx, log_fn, script_path, fwrite_fn, cfg, conj_alpha,
                                   historical_stage23_supported, conjunction_dt, negctrl_dt, collinearity_dt, cytoverlap_dt,
                                   purity_stage_dt, nonlinearity_dt, quantiseq_validated,
                                   signature_map, triplets, populations, eligible_cancers, power_dt) {
  suppressPackageStartupMessages(library(ggplot2))
  key_cols <- STAGE24_KEY_COLS
  downstream_path <- function(f) file.path(ctx$run_root, "tables", "downstream", f)
  stage24_inputs <- downstream_path(c("stage24_corrected_conjunction.tsv", "stage24_lump_and_stage_sensitivity.tsv",
                                      "stage24_collinearity_sensitivity.tsv", "stage24_negative_controls.tsv"))
  stage24_inputs <- stage24_inputs[file.exists(stage24_inputs)]

  # ==========================================================================
  # Primary corrected evidence funnel: every step filters the previous one,
  # starting from all conjunction-BH-supported combinations.
  # ==========================================================================
  primary <- build_stage24_primary_funnel(conjunction_dt, purity_stage_dt, collinearity_dt, negctrl_dt)
  primary_funnel <- primary$funnel
  n_conj_supported <- sum(conjunction_dt$all_three_arms_supported %in% TRUE)
  run_stage_assert(ctx, log_fn, primary_funnel$n_rows[1] == n_conj_supported,
                   "Primary evidence funnel starts from all %d conjunction-BH-supported combinations", n_conj_supported)
  run_stage_assert(ctx, log_fn, all(diff(primary_funnel$n_rows) <= 0),
                   "Primary evidence funnel is sequential and non-increasing (%s)", paste(primary_funnel$n_rows, collapse = " -> "))
  write_run_artifact(ctx, primary_funnel, file.path("tables", "downstream", "stage24_primary_evidence_funnel.tsv"),
                      "24_cancer_specific_cellular_context_corrected", script_path, stage24_inputs, fwrite_fn)
  log_fn("stage24_primary_evidence_funnel.tsv written. Primary funnel: %s",
         paste(sprintf("%s=%d", primary_funnel$step, primary_funnel$n_rows), collapse = " -> "))

  # ==========================================================================
  # Historical Stage 23 reconciliation (comparison only; does not define the
  # Stage 24 result).
  # ==========================================================================
  historical <- build_stage24_historical_funnel(historical_stage23_supported, conjunction_dt, purity_stage_dt,
                                                collinearity_dt, negctrl_dt)
  historical_funnel <- historical$funnel
  run_stage_assert(ctx, log_fn, historical_funnel$n_rows[1] == nrow(historical_stage23_supported),
                   "Historical reconciliation starts from all %d Stage 23 all-three-arms combinations", nrow(historical_stage23_supported))
  run_stage_assert(ctx, log_fn, all(diff(historical_funnel$n_rows) <= 0),
                   "Historical Stage 23 reconciliation is sequential and non-increasing (%s)", paste(historical_funnel$n_rows, collapse = " -> "))
  write_run_artifact(ctx, historical_funnel, file.path("tables", "downstream", "stage24_historical_reconciliation.tsv"),
                      "24_cancer_specific_cellular_context_corrected", script_path,
                      c(downstream_path("triplet_cell_context_summary.tsv"), stage24_inputs), fwrite_fn)
  log_fn("stage24_historical_reconciliation.tsv written. Historical Stage 23 reconciliation (not the Stage 24 result): %s",
         paste(sprintf("%s=%d", historical_funnel$step, historical_funnel$n_rows), collapse = " -> "))

  # ==========================================================================
  # Final result tables.
  # ==========================================================================
  final_dt <- build_stage24_final_table(primary$final_keys, conjunction_dt, purity_stage_dt, collinearity_dt, negctrl_dt,
                                        nonlinearity_dt, historical_stage23_supported, conj_alpha)
  run_stage_assert(ctx, log_fn, nrow(final_dt) == tail(primary_funnel$n_rows, 1),
                   "Final result table has one row per final primary-funnel key (%d rows)", nrow(final_dt))
  final_checks <- validate_stage24_final_table(final_dt)
  for (chk in names(final_checks)) {
    run_stage_assert(ctx, log_fn, final_checks[[chk]], "Every final Stage 24 row satisfies: %s", chk)
  }
  triplet_summary <- summarize_stage24_final_triplets(final_dt)
  final_inputs <- c(stage24_inputs, downstream_path(c("stage24_nonlinearity_sensitivity.tsv", "triplet_cell_context_summary.tsv")))
  final_inputs <- final_inputs[file.exists(final_inputs)]
  write_run_artifact(ctx, final_dt, file.path("tables", "downstream", "stage24_final_robust_nonoverlapping.tsv"),
                      "24_cancer_specific_cellular_context_corrected", script_path, final_inputs, fwrite_fn)
  write_run_artifact(ctx, triplet_summary, file.path("tables", "downstream", "stage24_final_triplet_summary.tsv"),
                      "24_cancer_specific_cellular_context_corrected", script_path,
                      downstream_path("stage24_final_robust_nonoverlapping.tsv"), fwrite_fn)
  log_fn("stage24_final_robust_nonoverlapping.tsv (%d rows) and stage24_final_triplet_summary.tsv (%d miRNA-gene triplets) written.",
         nrow(final_dt), nrow(triplet_summary))

  # ==========================================================================
  # Heatmaps
  # ==========================================================================
  fig_dir <- file.path(ctx$run_root, "figures", "cell_context_corrected")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  save_plot <- function(p, filename_base, width, height) {
    ggsave(file.path(fig_dir, paste0(filename_base, ".pdf")), p, width = width, height = height, limitsize = FALSE)
    ggsave(file.path(fig_dir, paste0(filename_base, ".png")), p, width = width, height = height, dpi = 300, limitsize = FALSE)
  }
  conjunction_dt[, triplet_label := paste(mirna_id, gene, sep = "-")]
  conjunction_dt[, status := fifelse(is.na(p_conjunction), NA_character_, fifelse(all_three_arms_supported, "Conjunction-supported", "Not supported"))]

  p1 <- ggplot(conjunction_dt, aes(x = population, y = cancer, fill = status)) +
    geom_tile(color = "white", linewidth = 0.15) +
    facet_wrap(~triplet_label, ncol = 6) +
    scale_fill_manual(values = c("Conjunction-supported" = "#D55E00", "Not supported" = "#56B4E9"), na.value = "grey80",
                       name = sprintf("Conjunction FDR(BH)<%.2f", conj_alpha)) +
    theme_minimal(base_size = 6) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4), axis.text.y = element_text(size = 4), panel.grid = element_blank(), strip.text = element_text(size = 5)) +
    labs(title = "Stage 24 conjunction test results, full corrected universe (HC3, intersection-union test)", subtitle = "Grey = not evaluable (direct_tautology excluded or an arm non-estimable); never a null result", x = "MCPcounter population", y = "Cancer")
  save_plot(p1, "01_corrected_conjunction_heatmap", 16, max(8, ceiling(uniqueN(conjunction_dt$triplet_label) / 6) * 3))

  p2 <- ggplot(conjunction_dt[overlap_category == "nonoverlapping_population"], aes(x = population, y = cancer, fill = status)) +
    geom_tile(color = "white", linewidth = 0.15) + facet_wrap(~triplet_label, ncol = 6) +
    scale_fill_manual(values = c("Conjunction-supported" = "#D55E00", "Not supported" = "#56B4E9"), na.value = "grey80", name = "Conjunction result") +
    theme_minimal(base_size = 6) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4), axis.text.y = element_text(size = 4), panel.grid = element_blank(), strip.text = element_text(size = 5)) +
    labs(title = "Stage 24, nonoverlapping populations only (no direct tautology, no CYT-lymphocyte identity overlap)", subtitle = "Conjunction status before the LUMP, collinearity and permutation checks of the primary funnel", x = "MCPcounter population", y = "Cancer")
  save_plot(p2, "02_nonoverlapping_population_heatmap", 16, max(8, ceiling(uniqueN(conjunction_dt$triplet_label) / 6) * 3))

  p3 <- ggplot(conjunction_dt[overlap_category == "CYT_biological_overlap"], aes(x = population, y = cancer, fill = status)) +
    geom_tile(color = "white", linewidth = 0.15) + facet_wrap(~triplet_label, ncol = 6) +
    scale_fill_manual(values = c("Conjunction-supported" = "#D55E00", "Not supported" = "#56B4E9"), na.value = "grey80", name = "Conjunction result") +
    theme_minimal(base_size = 6) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4), axis.text.y = element_text(size = 4), panel.grid = element_blank(), strip.text = element_text(size = 5)) +
    labs(title = "Stage 24, CYT-biological-overlap populations only (T/CD8 T/cytotoxic lymphocyte/NK cells)", subtitle = "Report as expected immune-context sensitivity ONLY -- never independent evidence of a new mechanism", x = "MCPcounter population", y = "Cancer")
  save_plot(p3, "03_cyt_overlap_population_heatmap", 16, max(8, ceiling(uniqueN(conjunction_dt$triplet_label) / 6) * 3))

  reconciliation_status <- merge(historical$steps[[1]][, ..key_cols][, in_historical := TRUE],
                                 conjunction_dt[, c(key_cols, "all_three_arms_supported"), with = FALSE], by = key_cols, all = TRUE)
  reconciliation_status[, in_historical := !is.na(in_historical)]
  supported_flag <- reconciliation_status$all_three_arms_supported %in% TRUE
  reconciliation_status[, final_status := fifelse(in_historical & supported_flag, "Historical Stage 23 and Stage 24 conjunction-supported",
                                            fifelse(in_historical & !supported_flag, "Historical Stage 23 only (not conjunction-supported)",
                                            fifelse(!in_historical & supported_flag, "Stage 24 conjunction-supported only (not in historical Stage 23 set)", NA_character_)))]
  reconciliation_status <- reconciliation_status[!is.na(final_status)]
  reconciliation_status[, triplet_label := paste(mirna_id, gene, sep = "-")]
  p4 <- ggplot(reconciliation_status, aes(x = population, y = cancer, fill = final_status)) +
    geom_tile(color = "white", linewidth = 0.2) + facet_wrap(~triplet_label, ncol = 6) +
    scale_fill_manual(values = c("Historical Stage 23 and Stage 24 conjunction-supported" = "#D55E00",
                                 "Historical Stage 23 only (not conjunction-supported)" = "#999999",
                                 "Stage 24 conjunction-supported only (not in historical Stage 23 set)" = "#009E73"),
                      na.value = "grey90", name = "Historical reconciliation status") +
    theme_minimal(base_size = 6) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4), axis.text.y = element_text(size = 4), panel.grid = element_blank(), strip.text = element_text(size = 5)) +
    labs(title = "Historical Stage 23 reconciliation: Stage 23 all-three-arms set vs. Stage 24 conjunction support",
         subtitle = "Comparison only; the Stage 24 result is defined by the primary evidence funnel", x = "MCPcounter population", y = "Cancer")
  save_plot(p4, "04_historical_stage23_reconciliation_heatmap", 16, max(6, ceiling(uniqueN(reconciliation_status$triplet_label) / 6) * 3))

  p5 <- ggplot(power_dt, aes(x = reorder(cancer, n_complete_case), y = n_complete_case, fill = pct_not_evaluable)) +
    geom_col() + coord_flip() +
    scale_fill_viridis_c(option = "C", name = "% arm-level\nrows not evaluable") +
    theme_minimal(base_size = 8) +
    labs(title = "Stage 24 cancer-specific power/sample-size context", subtitle = "Bar height = complete-case N; color = % of this cancer's arm-level rows not evaluable", x = "Cancer", y = "Complete-case N")
  save_plot(p5, "05_power_sample_size_by_cancer", 8, max(6, length(eligible_cancers) * 0.22))

  log_fn("5 Stage 24 figure sets written to figures/cell_context_corrected/ (PDF+PNG each).")

  # ==========================================================================
  # Final report
  # ==========================================================================
  hist_n <- setNames(historical_funnel$n_rows, historical_funnel$step)
  # Prespecified classification rule, unchanged, evaluated on the historical
  # Stage 23 reconciliation counts (corrected-conjunction, nonoverlapping and
  # permutation steps).
  n_hist_conj <- hist_n[["3_corrected_conjunction_BH_supported"]]
  n_hist_nonoverlap <- hist_n[["5_nonoverlapping_population"]]
  n_hist_perm <- hist_n[["7_permutation_supported"]]
  classification <- if (n_hist_perm >= 3 && n_hist_nonoverlap >= 3) {
    "supported exploratory cellular-context evidence"
  } else if (n_hist_conj >= 1) {
    "weak hypothesis-generating evidence"
  } else {
    "no reliable population-specific evidence"
  }

  prim_n <- setNames(primary_funnel$n_rows, primary_funnel$step)
  n_final <- nrow(final_dt)
  n_final_triplets <- uniqueN(final_dt[, .(mirna_id, gene)])
  fmt_list <- function(x) if (length(x) == 0) "none" else paste(x, collapse = ", ")
  final_cancers <- final_dt[, .N, by = cancer][order(-N, cancer)]
  final_pops <- final_dt[, .N, by = population][order(-N, population)]
  supported_dt <- conjunction_dt[all_three_arms_supported %in% TRUE]
  n_evaluable <- sum(!is.na(conjunction_dt$conjunction_FDR_BH))
  recurrent_pairs <- final_dt[, .(n_cancers = uniqueN(cancer), cancers = paste(sort(unique(cancer)), collapse = "/")),
                              by = .(mirna_id, gene, population)][n_cancers > 1]
  n_triplets_multi_cancer <- sum(triplet_summary$n_cancers > 1)
  n_by <- sum(final_dt$conjunction_BY_supported)
  n_nonlin <- sum(final_dt$nonlinearity_flagged %in% TRUE)
  n_nonlin_na <- sum(is.na(final_dt$nonlinear_vs_linear_p))
  n_no_nonlin <- n_final - n_nonlin
  n_by_no_nonlin <- sum(final_dt$conjunction_BY_supported & !(final_dt$nonlinearity_flagged %in% TRUE))
  n_stage_eval <- sum(final_dt$stage_sensitivity_applicable %in% TRUE)
  n_stage_sig <- sum(final_dt$stage_sensitivity_applicable %in% TRUE & final_dt$stage_p_interaction < 0.05, na.rm = TRUE)
  n_leverage <- sum(final_dt$any_arm_leverage_unstable %in% TRUE)
  n_perm_pass <- sum(final_dt$passes_permutation_criterion %in% TRUE)
  n_at_min_p <- sum(final_dt$permutation_p_at_minimum_attainable)
  n_valid_perm <- sort(unique(final_dt$n_valid_permutations))
  min_attainable <- if (length(n_valid_perm)) 1 / (n_valid_perm + 1) else numeric(0)
  hist_final <- final_dt[in_historical_stage23 == TRUE]
  collin_counts <- final_dt[, .N, by = collinearity_assessment][order(collinearity_assessment)]
  outputs <- c("stage24_corrected_arm_level.tsv", "stage24_corrected_conjunction.tsv",
               "stage24_nonoverlapping_population_results.tsv", "stage24_cyt_overlap_results.tsv",
               "stage24_population_correlations.tsv", "stage24_collinearity_sensitivity.tsv",
               "stage24_negative_controls.tsv", "stage24_power_by_cancer.tsv", "stage24_nonlinearity_sensitivity.tsv",
               "stage24_lump_and_stage_sensitivity.tsv", "stage24_primary_evidence_funnel.tsv",
               "stage24_historical_reconciliation.tsv", "stage24_final_robust_nonoverlapping.tsv",
               "stage24_final_triplet_summary.tsv")
  funnel_md <- function(f) c("| Step | Rows | Triplets | miRNAs | Genes | Cancers | Populations |", "|---|---|---|---|---|---|---|",
                             sprintf("| %s | %d | %d | %d | %d | %d | %d |", f$step, f$n_rows, f$n_triplets, f$n_mirnas, f$n_genes, f$n_cancers, f$n_populations))

  report_lines <- c(
    "# Stage 24: Cellular-Context Conjunction Analysis -- Method and Results",
    "",
    sprintf("Generated by `scripts/23_cancer_specific_cellular_context_corrected.R`, run `%s`, %s.", ctx$run_id, format(Sys.time(), "%Y-%m-%d %H:%M:%OS0")),
    "",
    "The objective of this stage is to distinguish cancer-specific cellular-context signal from multiple-testing",
    "artifacts, CYT-immune-score identity overlap, correlated MCPcounter populations, tumour-purity confounding and",
    "bulk-tissue aggregation bias. The Stage 24 result is defined by the primary corrected evidence funnel (section 2);",
    "the historical Stage 23 reconciliation (section 4) is a comparison only.",
    "",
    "## 1. Corrected full-universe analysis",
    "",
    "Methods relative to Stage 23:",
    "1. Purity (CPE) is an explicit covariate inside the interaction model rather than a residualization applied only to X and Y, so the population score is also purity-adjusted.",
    "2. Primary inference uses HC3 heteroscedasticity-consistent sandwich standard errors; classical OLS is retained only as a diagnostic comparison. Fits with an HC3 warning or maximum leverage >= 0.99 are flagged leverage-unstable and excluded from the conjunction family.",
    "3. The \"all three arms\" claim is tested as a conjunction (intersection-union) test: p_conjunction = max(p1, p2, p3) of the three arms' raw HC3 p-values, BH-corrected once across all evaluable conjunction tests.",
    "4. Every combination is classified `direct_tautology` / `CYT_biological_overlap` / `nonoverlapping_population`; CYT-overlap results are reported as expected immune-context sensitivity, never independent mechanistic evidence.",
    "5. A correlated-population sensitivity refit (VIF, condition number, a shared lymphocyte-infiltration PC1 axis) for every conjunction-supported finding.",
    "6. Fixed-seed permutation negative controls with Monte Carlo uncertainty reported.",
    "7. Per-cancer power and sample-size context.",
    "8. A prespecified spline nonlinearity check (descriptive only), a LUMP-purity refit of the complete three-arm conjunction and a completeness-gated AJCC-stage sensitivity refit.",
    "",
    sprintf("The corrected universe contains %d (triplet, cancer, population) rows, of which %d are evaluable, non-tautological combinations in the conjunction BH family. **%d combinations are conjunction-BH-supported (FDR < %.2f)**, representing %d unique miRNA-gene triplets.",
            nrow(conjunction_dt), n_evaluable, n_conj_supported, conj_alpha, uniqueN(supported_dt[, .(mirna_id, gene)])),
    sprintf("By overlap category, the conjunction-supported combinations are: %s.",
            fmt_list(supported_dt[, .N, by = overlap_category][order(overlap_category), sprintf("%s %d", overlap_category, N)])),
    "",
    sprintf("Power context: Spearman correlation between per-cancer complete-case N and the number of conjunction-supported findings in that cancer: %.3f; between median population-score SD and findings: %.3f (`tables/downstream/stage24_power_by_cancer.tsv`, including simulation-based minimum-detectable effect at 80%% power).",
            suppressWarnings(cor(power_dt$n_complete_case, power_dt$n_all_three_arms_supported_here, method = "spearman")),
            suppressWarnings(cor(power_dt$median_population_sd, power_dt$n_all_three_arms_supported_here, method = "spearman"))),
    "",
    "## 2. Primary corrected evidence funnel",
    "",
    "Each step is derived from the immediately preceding step by a keyed semi-join on (mirna_id, gene, cancer, population):",
    "",
    sprintf("- `%s`: %s", names(STAGE24_PRIMARY_STEPS), unname(STAGE24_PRIMARY_STEPS)),
    "",
    funnel_md(primary_funnel),
    "",
    sprintf("%d corrected conjunction-supported combinations were identified; %d retained CPE/LUMP agreement; %d were in nonoverlapping populations after sequential filtering; %d passed or were not applicable for the correlated-population sensitivity check; %d of these %d passed the permutation criterion.",
            prim_n[[1]], prim_n[[2]], prim_n[[3]], prim_n[[4]], prim_n[[5]], prim_n[[4]]),
    "",
    "## 3. Final corrected result",
    "",
    sprintf("The final corrected result contains **%d (triplet, cancer, population) rows representing %d unique miRNA-gene triplets** (`tables/downstream/stage24_final_robust_nonoverlapping.tsv`; triplet-level summary in `tables/downstream/stage24_final_triplet_summary.tsv`).",
            n_final, n_final_triplets),
    sprintf("- Cancers (rows): %s.", fmt_list(final_cancers[, sprintf("%s %d", cancer, N)])),
    sprintf("- Populations (rows): %s.", fmt_list(final_pops[, sprintf("%s %d", population, N)])),
    sprintf("- Correlated-population assessment: %s.", fmt_list(collin_counts[, sprintf("%s %d", collinearity_assessment, N)])),
    sprintf("- Leverage-unstable rows: %d.", n_leverage),
    "",
    sprintf("These are cancer- and population-specific associations. %d of the %d triplets are supported in more than one cancer; %d (triplet, population) combinations are supported in more than one cancer%s, and none in more than %d of the %d eligible cancers. The analysis does not establish a cellular-context mechanism recurrent across cancers.",
            n_triplets_multi_cancer, n_final_triplets, nrow(recurrent_pairs),
            if (nrow(recurrent_pairs) > 0) sprintf(" (%s)", fmt_list(unique(recurrent_pairs[, sprintf("%s in %s", population, cancers)]))) else "",
            if (nrow(recurrent_pairs) > 0) max(recurrent_pairs$n_cancers) else 1L, length(eligible_cancers)),
    "",
    "## 4. Historical Stage 23 reconciliation",
    "",
    sprintf("Stage 23's historical all-three-arms-significant set (%d combinations, three separately BH-adjusted arm-level FDRs) is followed through the same corrected checks. **This is a reconciliation analysis only: its survivors are historical survivors, not the complete corrected Stage 24 result.**",
            nrow(historical_stage23_supported)),
    "",
    sprintf("- `%s`: %s", names(STAGE24_HISTORICAL_STEPS), unname(STAGE24_HISTORICAL_STEPS)),
    "",
    funnel_md(historical_funnel),
    "",
    sprintf("Only %d of the %d final corrected rows also belong to the historical Stage 23 set: %s.",
            nrow(hist_final), n_final, fmt_list(hist_final[, sprintf("%s / %s / %s / %s", mirna_id, gene, cancer, population)])),
    sprintf("The historical %d -> %d sequence is a reconciliation analysis only and is not the primary Stage 24 evidence funnel.",
            hist_n[[1]], tail(historical_funnel$n_rows, 1)),
    "",
    "## 5. BH primary inference and BY sensitivity",
    "",
    sprintf("Primary inference is the conjunction BH-FDR < %.2f. As a dependence-robust sensitivity analysis, %d of the %d final rows also pass Benjamini-Yekutieli FDR < %.2f.",
            conj_alpha, n_by, n_final, conj_alpha),
    "",
    "## 6. Nonlinearity sensitivity",
    "",
    sprintf("%d of the %d final rows show evidence of nonlinearity (spline vs. linear interaction, p < 0.05) and %d have no detected nonlinearity%s. The spline check is descriptive and is not a post hoc exclusion criterion.",
            n_nonlin, n_final, n_no_nonlin,
            if (n_nonlin_na > 0) sprintf(" (%d of these were not evaluable by the spline check)", n_nonlin_na) else ""),
    sprintf("%d of the %d final rows are both BY-supported and have no detected nonlinearity. This is only a conservative descriptive subset; it does not redefine the primary result.",
            n_by_no_nonlin, n_final),
    "",
    "## 7. AJCC-stage sensitivity",
    "",
    sprintf("%d of the %d final rows are stage-evaluable (cancer with >= 80%% AJCC-stage completeness and at least two stage levels). %s The other %d rows are not stage-evaluable; this is a data-availability limitation, not a failure.",
            n_stage_eval, n_final,
            if (n_stage_eval > 0 && n_stage_sig == n_stage_eval) sprintf("All %d stage-evaluable rows retain nominal stage-adjusted support (stage_p_interaction < 0.05).", n_stage_eval)
            else sprintf("%d of the %d stage-evaluable rows retain nominal stage-adjusted support (stage_p_interaction < 0.05).", n_stage_sig, n_stage_eval),
            n_final - n_stage_eval),
    "",
    "## 8. Permutation-control interpretation",
    "",
    sprintf("%d of the %d final rows pass the permutation criterion (empirical p < %.2f; %s permutations of the population score within cancer, seed %s).",
            n_perm_pass, n_final, conj_alpha, fmt_list(n_valid_perm), cfg$permutation_seed),
    sprintf("%d of the %d final rows have the empirical p-value %s, the minimum attainable value (1 / (n_permutations + 1)) when zero permutations are at least as extreme as the observed conjunction p-value. Equal minimum permutation p-values only bound the evidence from below and must not be used to rank findings.",
            n_at_min_p, n_final, fmt_list(sprintf("%.9f", min_attainable))),
    "",
    "## 9. Deconvolution limitation",
    "",
    if (quantiseq_validated) {
      "`analysis_config$immune$quantiseq_input_scale_validated` is TRUE, but this stage does not compute an alternative-deconvolution sensitivity result; none is reported."
    } else {
      "**quanTIseq sensitivity was not performed** because `analysis_config$immune$quantiseq_input_scale_validated = FALSE` (Stage 13 runs only after the expression-scale assumption has been validated). No alternative-deconvolution result is reported or estimated; all population scores are MCPcounter estimates."
    },
    "",
    "## 10. Scientific limitations",
    "",
    "- **No claim of direct miRNA targeting, cellular origin, within-cell regulation, causality, validated biomarker performance or clinical actionability is supported by this analysis.** MCPcounter estimates bulk-tissue cellular composition; it cannot localize a mature miRNA's expression to an individual cell, and a statistical interaction with a population's estimated abundance is not evidence of direct regulation.",
    "- **Bulk-tissue interaction analyses remain vulnerable to ecological/aggregation bias and Simpson's-paradox-type effects**, even after this correction. A pattern present at the bulk-sample level can reverse, vanish, or be manufactured by within-sample cell-type-proportion mixing.",
    "- A `CYT_biological_overlap` finding is expected immune-context sensitivity, not independent evidence of a new mechanism.",
    "- Descriptive tertile strata and the spline nonlinearity check are interpretation aids only, never confirmatory tests.",
    "- The permutation control is applied to the conjunction-supported combinations and a fixed random comparison sample, not to the full grid.",
    "",
    "### Final classification",
    "",
    sprintf("**`%s`**", classification),
    "",
    sprintf("The classification rule is the prespecified rule of this script, unchanged, evaluated on the historical Stage 23 reconciliation counts (corrected conjunction-supported %d, nonoverlapping %d, permutation-supported %d). It was not re-specified for the primary funnel.",
            n_hist_conj, n_hist_nonoverlap, n_hist_perm),
    "",
    "## 11. Outputs",
    "",
    sprintf("- `tables/downstream/%s`", outputs),
    "- `figures/cell_context_corrected/01_corrected_conjunction_heatmap`, `02_nonoverlapping_population_heatmap`, `03_cyt_overlap_population_heatmap`, `04_historical_stage23_reconciliation_heatmap` (historical Stage 23 reconciliation), `05_power_sample_size_by_cancer` (PDF and PNG)",
    "- `reports/STAGE24_CELLULAR_CONTEXT_METHOD_AND_RESULTS.md` (this report)",
    "",
    "Stage 23 outputs (`tables/downstream/triplet_cell_context_*.tsv`) are read, never modified.",
    ""
  )
  writeLines(report_lines, file.path(ctx$run_root, "reports", "STAGE24_CELLULAR_CONTEXT_METHOD_AND_RESULTS.md"))
  log_fn("reports/STAGE24_CELLULAR_CONTEXT_METHOD_AND_RESULTS.md written. Final classification: %s", classification)
  invisible(classification)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=23 or source scripts/run.R in RStudio.")
}
