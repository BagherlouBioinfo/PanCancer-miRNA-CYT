# ============================================================================
# test_cell_context_helpers.R
#
# Known-answer, boundary, missing-value, non-estimability, tautology, and
# shuffled-order/label synthetic tests for R/cell_context_helpers.R --
# Stage 21's core cancer-specific cellular-composition-context statistical
# logic. Calls the REAL, unmodified production functions against
# constructed synthetic data with a KNOWN designed answer.
#
# SYNTHETIC DATA ONLY.
#
# Usage: Rscript tests/test_cell_context_helpers.R <project_root>
# ============================================================================

run_cell_context_helper_tests <- function(project_root) {
  source(file.path(project_root, "scripts", "lib", "cell_context_helpers.R"))

  results <- list()
  add <- function(id, desc, pass, detail) {
    results[[length(results) + 1]] <<- data.frame(test_id = id, description = desc,
      result = if (isTRUE(pass)) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE)
    cat(sprintf("[%s] %s -- %s\n", if (isTRUE(pass)) "PASS" else "FAIL", id, detail))
  }

  set.seed(2026)
  n <- 200

  # ==========================================================================
  # Scenario 1: an interaction IS truly present -- the predictor-outcome
  # association genuinely strengthens as the population score increases.
  # Constructed directly: outcome = predictor * (1 + population_z), so the
  # effective slope grows with population_z by design.
  # ==========================================================================
  predictor1 <- rnorm(n)
  population1 <- rnorm(n)
  pop_z1 <- as.numeric(scale(population1))
  outcome1 <- predictor1 * (1 + 1.5 * pop_z1) + rnorm(n, sd = 0.3)
  fit1 <- fit_cell_context_interaction(predictor1, outcome1, population1, min_n = 30)
  add("S1_interaction_truly_present_detected", "A designed, genuine interaction (association strengthens with population) is detected: estimable, interaction p-value small, |slope_high| > |slope_low|",
      fit1$estimable && fit1$p_interaction < 0.01 && abs(fit1$slope_high) > abs(fit1$slope_low),
      sprintf("estimable=%s p_int=%.4g slope_low=%.3f slope_high=%.3f", fit1$estimable, fit1$p_interaction, fit1$slope_low, fit1$slope_high))

  # ==========================================================================
  # Scenario 2: NO interaction is present -- a constant-strength
  # association regardless of population.
  # ==========================================================================
  predictor2 <- rnorm(n)
  population2 <- rnorm(n)
  outcome2 <- 0.6 * predictor2 + rnorm(n, sd = 0.6)  # no population term at all
  fit2 <- fit_cell_context_interaction(predictor2, outcome2, population2, min_n = 30)
  add("S2_no_interaction_present_not_significant", "No designed interaction: the interaction p-value is NOT small (fails to reject a true null), confirming the test does not spuriously find effect modification that was never constructed",
      fit2$estimable && fit2$p_interaction > 0.05,
      sprintf("estimable=%s p_int=%.4g", fit2$estimable, fit2$p_interaction))

  # ==========================================================================
  # Scenario 3: the relationship exists ONLY in high-composition samples
  # (predictor/outcome are pure noise, UNCORRELATED, when population is low;
  # a real, strong association only when population is high).
  # ==========================================================================
  population3 <- runif(n, -2, 2)
  predictor3 <- rnorm(n)
  is_high3 <- population3 > 0
  outcome3 <- ifelse(is_high3, 0.9 * predictor3, 0) + rnorm(n, sd = 0.3)
  fit3 <- fit_cell_context_interaction(predictor3, outcome3, population3, min_n = 30)
  strat3 <- tertile_stratified_correlation(predictor3, outcome3, population3, min_n_stratum = 10)
  add("S3_relationship_only_in_high_composition", "A relationship existing ONLY in high-population samples is detected as a significant interaction with |slope_high| > |slope_low|, and the descriptive tertile table shows a much stronger 'high' stratum rho than 'low'",
      fit3$estimable && fit3$p_interaction < 0.05 && abs(fit3$slope_high) > abs(fit3$slope_low) &&
        abs(strat3$rho[strat3$stratum == "high"]) > abs(strat3$rho[strat3$stratum == "low"]),
      sprintf("p_int=%.4g slope_low=%.3f slope_high=%.3f | tertile rho low=%.3f high=%.3f",
              fit3$p_interaction, fit3$slope_low, fit3$slope_high,
              strat3$rho[strat3$stratum == "low"], strat3$rho[strat3$stratum == "high"]))

  # ==========================================================================
  # Scenario 4: the relationship exists ONLY in low-composition samples
  # (mirror image of Scenario 3).
  # ==========================================================================
  population4 <- runif(n, -2, 2)
  predictor4 <- rnorm(n)
  is_low4 <- population4 <= 0
  outcome4 <- ifelse(is_low4, 0.9 * predictor4, 0) + rnorm(n, sd = 0.3)
  fit4 <- fit_cell_context_interaction(predictor4, outcome4, population4, min_n = 30)
  strat4 <- tertile_stratified_correlation(predictor4, outcome4, population4, min_n_stratum = 10)
  add("S4_relationship_only_in_low_composition", "A relationship existing ONLY in low-population samples is detected as a significant interaction with |slope_low| > |slope_high|, and the descriptive tertile table shows a much stronger 'low' stratum rho than 'high'",
      fit4$estimable && fit4$p_interaction < 0.05 && abs(fit4$slope_low) > abs(fit4$slope_high) &&
        abs(strat4$rho[strat4$stratum == "low"]) > abs(strat4$rho[strat4$stratum == "high"]),
      sprintf("p_int=%.4g slope_low=%.3f slope_high=%.3f | tertile rho low=%.3f high=%.3f",
              fit4$p_interaction, fit4$slope_low, fit4$slope_high,
              strat4$rho[strat4$stratum == "low"], strat4$rho[strat4$stratum == "high"]))

  # ==========================================================================
  # Scenario 5: missingness makes a model non-estimable.
  # ==========================================================================
  predictor5 <- rep(NA_real_, 40); outcome5 <- rnorm(40); population5 <- rnorm(40)
  fit5a <- fit_cell_context_interaction(predictor5, outcome5, population5, min_n = 30)
  add("S5a_all_missing_predictor_not_estimable", "A predictor vector that is entirely NA/missing produces estimable=FALSE, never a spurious fit or an error",
      !fit5a$estimable, sprintf("estimable=%s n=%d", fit5a$estimable, fit5a$n))

  predictor5b <- c(rnorm(25), rep(NA_real_, 15)); outcome5b <- rnorm(40); population5b <- rnorm(40)
  fit5b <- fit_cell_context_interaction(predictor5b, outcome5b, population5b, min_n = 30)
  add("S5b_below_min_n_after_missingness_removed_not_estimable", "After removing missing values, remaining n (25) is below min_n (30) -- correctly not estimable, not silently fit on too few samples",
      !fit5b$estimable, sprintf("estimable=%s n=%d (expected n=25, below min_n=30)", fit5b$estimable, fit5b$n))

  zero_var_pred <- rep(5, 40); outcome5c <- rnorm(40); population5c <- rnorm(40)
  fit5c <- fit_cell_context_interaction(zero_var_pred, outcome5c, population5c, min_n = 30)
  add("S5c_zero_variance_predictor_not_estimable", "A zero-variance predictor (constant value) is correctly flagged not estimable rather than producing a degenerate/undefined model",
      !fit5c$estimable, sprintf("estimable=%s", fit5c$estimable))

  # ==========================================================================
  # Scenario 6: a signature-gene tautology must be excluded.
  # ==========================================================================
  sig_map <- data.frame(gene = c("CD3D", "CD3G", "CD28", "CSF1R"), population = c("T cells", "T cells", "T cells", "Monocytic lineage"), stringsAsFactors = FALSE)
  add("S6a_signature_gene_flagged_tautological", "A gene that IS a marker of the queried population is correctly flagged as a tautology",
      is_gene_signature_member("CD3D", "T cells", sig_map), "checked CD3D against T cells (a real marker)")
  add("S6b_signature_gene_not_flagged_for_different_population", "The SAME gene is correctly NOT flagged when queried against a DIFFERENT population it is not a marker of",
      !is_gene_signature_member("CD3D", "Monocytic lineage", sig_map), "checked CD3D against Monocytic lineage (not a marker there)")
  add("S6c_non_signature_gene_never_flagged", "A gene that is not a marker of ANY population in the signature map is never flagged",
      !is_gene_signature_member("SOCS1", "T cells", sig_map), "checked SOCS1 (not in sig_map) against T cells")
  add("S6d_case_insensitive_gene_match", "The tautology check is case-insensitive (a real-world safety margin against symbol-casing mismatches between data sources)",
      is_gene_signature_member("cd3d", "T cells", sig_map), "checked lowercase 'cd3d' against T cells")

  # Full classify_cell_context() tautology-exclusion integration:
  fake_arm_ok <- list(estimable = TRUE, n = 100, beta_main = 1, beta_interaction = 0, se_interaction = 0.1,
                       p_interaction = 0.5, ci_lower = -0.1, ci_upper = 0.1, slope_low = 1, slope_high = 1, fdr = 0.5)
  verdict_tauto <- classify_cell_context(
    list(mirna_cyt = fake_arm_ok, gene_cyt = fake_arm_ok, mirna_gene = fake_arm_ok),
    tautological_arms = c("gene_cyt", "mirna_gene"),
    baseline_pass = c(mirna_cyt = TRUE, gene_cyt = TRUE, mirna_gene = TRUE))
  add("S6e_classify_cell_context_returns_tautological_exclusion", "classify_cell_context() returns 'tautological_exclusion' whenever any required arm is tautological, regardless of how strong the (biologically meaningless) interaction signal looks",
      identical(verdict_tauto, "tautological_exclusion"), sprintf("verdict=%s", verdict_tauto))

  # ==========================================================================
  # Scenario 7: cancer labels and sample orders are shuffled -- must not
  # change the result (order-independence / no positional-alignment bug).
  # ==========================================================================
  ord <- sample(seq_len(n))
  fit1_shuffled <- fit_cell_context_interaction(predictor1[ord], outcome1[ord], population1[ord], min_n = 30)
  add("S7a_shuffled_sample_order_gives_identical_result", "Shuffling predictor/outcome/population's sample order IN LOCKSTEP (a valid re-ordering, not a mismatch) gives numerically identical results to the unshuffled fit",
      isTRUE(all.equal(fit1$p_interaction, fit1_shuffled$p_interaction)) && isTRUE(all.equal(fit1$beta_interaction, fit1_shuffled$beta_interaction)),
      sprintf("unshuffled p_int=%.6g shuffled p_int=%.6g", fit1$p_interaction, fit1_shuffled$p_interaction))

  strat1 <- tertile_stratified_correlation(predictor1, outcome1, population1, min_n_stratum = 10)
  strat1_shuffled <- tertile_stratified_correlation(predictor1[ord], outcome1[ord], population1[ord], min_n_stratum = 10)
  add("S7b_shuffled_order_tertile_stratification_identical", "tertile_stratified_correlation() is also order-independent under a lockstep shuffle (re-derives strata from the values themselves, not positional assumptions)",
      isTRUE(all.equal(sort(strat1$rho), sort(strat1_shuffled$rho))),
      sprintf("unshuffled rhos=%s shuffled rhos=%s", paste(round(strat1$rho,4),collapse=","), paste(round(strat1_shuffled$rho,4),collapse=",")))

  # ==========================================================================
  # classify_cell_context(): remaining verdict categories, hand-worked.
  # ==========================================================================
  enriched_arm <- list(estimable = TRUE, n = 100, slope_low = 0.1, slope_high = 0.9, fdr = 0.001)
  flat_ns <- list(estimable = TRUE, n = 100, slope_low = 0.5, slope_high = 0.5, fdr = 0.6)
  v_enriched <- classify_cell_context(list(mirna_cyt = enriched_arm, gene_cyt = enriched_arm, mirna_gene = enriched_arm),
                                       tautological_arms = character(0), baseline_pass = c(mirna_cyt=TRUE,gene_cyt=TRUE,mirna_gene=TRUE))
  add("V1_composition_enriched_hand_worked", "All 3 arms show a significant interaction with |slope_high|>|slope_low| -> composition_enriched",
      identical(v_enriched, "composition_enriched"), sprintf("verdict=%s", v_enriched))

  depleted_arm <- list(estimable = TRUE, n = 100, slope_low = 0.9, slope_high = 0.1, fdr = 0.001)
  v_depleted <- classify_cell_context(list(mirna_cyt = depleted_arm, gene_cyt = depleted_arm, mirna_gene = depleted_arm),
                                       tautological_arms = character(0), baseline_pass = c(mirna_cyt=TRUE,gene_cyt=TRUE,mirna_gene=TRUE))
  add("V2_composition_depleted_hand_worked", "All 3 arms show a significant interaction with |slope_low|>|slope_high| -> composition_depleted",
      identical(v_depleted, "composition_depleted"), sprintf("verdict=%s", v_depleted))

  v_stable <- classify_cell_context(list(mirna_cyt = flat_ns, gene_cyt = flat_ns, mirna_gene = flat_ns),
                                     tautological_arms = character(0), baseline_pass = c(mirna_cyt=TRUE,gene_cyt=TRUE,mirna_gene=TRUE))
  add("V3_stable_across_composition_hand_worked", "No significant interaction in any arm, but the baseline (population-independent) association IS detected -> stable_across_composition",
      identical(v_stable, "stable_across_composition"), sprintf("verdict=%s", v_stable))

  v_not_detected <- classify_cell_context(list(mirna_cyt = flat_ns, gene_cyt = flat_ns, mirna_gene = flat_ns),
                                           tautological_arms = character(0), baseline_pass = c(mirna_cyt=FALSE,gene_cyt=TRUE,mirna_gene=TRUE))
  add("V4_not_detected_hand_worked", "No significant interaction, AND the baseline association is not detected in all three arms -> not_detected",
      identical(v_not_detected, "not_detected"), sprintf("verdict=%s", v_not_detected))

  not_est_arm <- list(estimable = FALSE)
  v_not_eval <- classify_cell_context(list(mirna_cyt = enriched_arm, gene_cyt = enriched_arm, mirna_gene = not_est_arm),
                                       tautological_arms = character(0), baseline_pass = c(mirna_cyt=TRUE,gene_cyt=TRUE,mirna_gene=TRUE))
  add("V5_not_evaluable_when_any_required_arm_non_estimable", "ANY of the 3 required arms being non-estimable (e.g. insufficient n) makes the whole triplet-cancer-population combination not_evaluable, regardless of the other two arms' results",
      identical(v_not_eval, "not_evaluable"), sprintf("verdict=%s", v_not_eval))

  df <- do.call(rbind, results)
  out_path <- file.path(project_root, "tests", "output", "test_cell_context_helpers_results.tsv")
  write.table(df, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
  n_fail <- sum(df$result == "FAIL")
  cat(sprintf("\n%d/%d PASS, %d FAIL\n", sum(df$result == "PASS"), nrow(df), n_fail))
  if (n_fail > 0) stop(sprintf("%d cell-context helper test(s) FAILED. See %s.", n_fail, out_path), call. = FALSE)
  invisible(df)
}

if (sys.nframe() == 0) {
  cli_args <- commandArgs(trailingOnly = TRUE)
  if (length(cli_args) < 1) stop("Usage: Rscript test_cell_context_helpers.R <project_root>")
  invisible(run_cell_context_helper_tests(normalizePath(cli_args[[1]], mustWork = TRUE)))
}
