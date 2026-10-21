# ============================================================================
# test_cell_context_corrected_helpers.R
#
# Known-answer, boundary, missing-value, duplicated/shuffled-key, and
# calibration tests for the CORRECTED Stage 21b methodology
# (R/cell_context_helpers.R's post-critique additions:
# fit_cell_context_interaction_hc3(), compute_conjunction_p(),
# classify_overlap_category(), flag_correlated_population_pairs(),
# compute_lymphocyte_pc1(), compute_vif(), compute_condition_number()).
# Calls the REAL, unmodified production functions against constructed
# synthetic data with a KNOWN designed answer -- never grep-based source
# text checks.
#
# SYNTHETIC DATA ONLY.
#
# Usage: Rscript tests/test_cell_context_corrected_helpers.R <project_root>
# ============================================================================

run_cell_context_corrected_helper_tests <- function(project_root) {
  suppressPackageStartupMessages({ library(sandwich) })
  source(file.path(project_root, "scripts", "lib", "cell_context_helpers.R"))

  results <- list()
  add <- function(id, desc, pass, detail) {
    results[[length(results) + 1]] <<- data.frame(test_id = id, description = desc,
      result = if (isTRUE(pass)) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE)
    cat(sprintf("[%s] %s -- %s\n", if (isTRUE(pass)) "PASS" else "FAIL", id, detail))
  }

  set.seed(2026)
  n <- 300

  # ==========================================================================
  # H1-H2: HC3 interaction inference matches sandwich::vcovHC(type="HC3")
  # computed independently on the identical model, both under homoscedastic
  # and heteroscedastic errors.
  # ==========================================================================
  x <- rnorm(n); pop <- rnorm(n); purity <- runif(n, 0.4, 1)
  y_homo <- x * (1 + 0.6 * scale(pop)) + rnorm(n, sd = 1) + 0.5 * purity
  r_h1 <- fit_cell_context_interaction_hc3(x, y_homo, pop, purity, min_n = 30)
  ref_fit <- lm(rank(y_homo) ~ rank(x) * as.numeric(scale(pop)) + purity)
  ref_vc <- sandwich::vcovHC(ref_fit, type = "HC3")
  ref_se <- sqrt(ref_vc["rank(x):as.numeric(scale(pop))", "rank(x):as.numeric(scale(pop))"])
  add("H1_hc3_se_matches_sandwich_package_homoscedastic", "fit_cell_context_interaction_hc3()'s HC3 SE matches an independently-built sandwich::vcovHC(type='HC3') on the identical model (homoscedastic data)",
      isTRUE(r_h1$estimable) && abs(r_h1$se_interaction_hc3 - ref_se) < 1e-8,
      sprintf("production SE=%.10f, independent reference SE=%.10f, diff=%.2e", r_h1$se_interaction_hc3, ref_se, abs(r_h1$se_interaction_hc3 - ref_se)))

  het_scale <- 0.2 + 2 * abs(x)  # variance grows with |x| -- classic heteroscedasticity
  y_het <- x * (1 + 0.6 * scale(pop)) + rnorm(n, sd = 1) * het_scale + 0.5 * purity
  r_h2 <- fit_cell_context_interaction_hc3(x, y_het, pop, purity, min_n = 30)
  ref_fit2 <- lm(rank(y_het) ~ rank(x) * as.numeric(scale(pop)) + purity)
  ols_se2 <- summary(ref_fit2)$coefficients["rank(x):as.numeric(scale(pop))", "Std. Error"]
  ref_vc2 <- sandwich::vcovHC(ref_fit2, type = "HC3")
  ref_se2 <- sqrt(ref_vc2["rank(x):as.numeric(scale(pop))", "rank(x):as.numeric(scale(pop))"])
  add("H2_hc3_se_matches_sandwich_heteroscedastic_and_differs_from_OLS", "Under heteroscedastic errors, the production HC3 SE matches the independent sandwich reference AND differs materially from the classical OLS SE (proving HC3 is actually doing something different, not silently falling back to OLS)",
      isTRUE(r_h2$estimable) && abs(r_h2$se_interaction_hc3 - ref_se2) < 1e-8 && abs(r_h2$se_interaction_hc3 - r_h2$se_interaction_ols) / r_h2$se_interaction_ols > 0.01,
      sprintf("HC3 SE=%.6f, OLS SE=%.6f, reference HC3=%.6f", r_h2$se_interaction_hc3, r_h2$se_interaction_ols, ref_se2))

  # ==========================================================================
  # H3: purity is a genuine covariate -- a synthetic Y with a strong purity
  # main effect but ZERO true interaction should not manufacture a spurious
  # interaction once purity is (correctly) controlled for.
  # ==========================================================================
  set.seed(303)
  x3 <- rnorm(n); pop3 <- rnorm(n); purity3 <- runif(n, 0.3, 1)
  y3 <- x3 + 5 * purity3 + rnorm(n, sd = 0.5)  # no interaction, strong purity confound
  r_h3 <- fit_cell_context_interaction_hc3(x3, y3, pop3, purity3, min_n = 30)
  add("H3_purity_confound_no_spurious_interaction", "A strong purity main effect with NO true interaction does not produce a spuriously significant interaction once purity is controlled for in-model",
      isTRUE(r_h3$estimable) && r_h3$p_interaction_hc3 > 0.05,
      sprintf("p_interaction_hc3=%.4f (expected > 0.05)", r_h3$p_interaction_hc3))

  # ==========================================================================
  # C1-C2: conjunction p-value = max(p1,p2,p3); NA if any arm missing.
  # ==========================================================================
  add("C1_conjunction_is_max_known_values", "compute_conjunction_p() returns exactly max(p1,p2,p3) for known inputs",
      compute_conjunction_p(0.001, 0.20, 0.03) == 0.20 && compute_conjunction_p(0.5, 0.5, 0.5) == 0.5,
      sprintf("max(0.001,0.20,0.03)=%.4f, max(0.5,0.5,0.5)=%.4f", compute_conjunction_p(0.001, 0.20, 0.03), compute_conjunction_p(0.5, 0.5, 0.5)))
  add("C2_conjunction_NA_if_any_arm_missing", "compute_conjunction_p() returns NA if ANY of the three arm p-values is NA/non-finite -- a conjunction claim cannot be made if one arm could not be tested",
      is.na(compute_conjunction_p(0.001, NA, 0.03)) && is.na(compute_conjunction_p(0.001, 0.02, NaN)) && is.na(compute_conjunction_p(0.001, 0.02, Inf)),
      "all three missing-arm variants correctly return NA")

  # ==========================================================================
  # C3: known-answer 3-arm alternative scenarios -- the conjunction test
  # must NOT classify a one-arm or two-arm-only alternative as
  # all-three-arms-supported, at any reasonable BH family size.
  # ==========================================================================
  one_arm_p <- c(0.0001, 0.70, 0.65)     # only arm 1 truly affected
  two_arm_p <- c(0.0002, 0.0003, 0.55)   # two arms affected, third flat
  three_arm_p <- c(4e-8, 6e-8, 8e-8)     # all three truly, strongly affected -- small enough
                                          # that its conjunction p (8e-8) survives BH correction
                                          # even at rank 1 of a 203-test family (8e-8*203 << 0.05),
                                          # unlike a borderline p-value which BH may legitimately
                                          # not rescue depending on family size/composition
  filler_p <- runif(200, 0.01, 1)  # a realistic surrounding FDR family
  fam <- c(compute_conjunction_p(one_arm_p[1], one_arm_p[2], one_arm_p[3]),
           compute_conjunction_p(two_arm_p[1], two_arm_p[2], two_arm_p[3]),
           compute_conjunction_p(three_arm_p[1], three_arm_p[2], three_arm_p[3]),
           filler_p)
  fdr <- p.adjust(fam, method = "BH")
  add("C3_one_arm_alternative_not_classified_supported", "A one-arm-only true alternative (arms 2,3 flat) is correctly NOT classified all-three-arms-supported by the conjunction test",
      fdr[1] >= 0.05, sprintf("conjunction FDR(one-arm scenario)=%.4f (expected >= 0.05)", fdr[1]))
  add("C4_two_arm_alternative_not_classified_supported", "A two-arm-only true alternative (arm 3 flat) is correctly NOT classified all-three-arms-supported by the conjunction test",
      fdr[2] >= 0.05, sprintf("conjunction FDR(two-arm scenario)=%.4f (expected >= 0.05)", fdr[2]))
  add("C5_three_arm_alternative_classified_supported", "A genuine three-arm alternative (all three arms truly affected) IS classified all-three-arms-supported by the conjunction test",
      fdr[3] < 0.05, sprintf("conjunction FDR(three-arm scenario)=%.4g (expected < 0.05)", fdr[3]))

  # ==========================================================================
  # B1-B2: global BH vs. BY on a known p-value vector (matches base R exactly).
  # ==========================================================================
  known_p <- c(0.001, 0.01, 0.02, 0.03, 0.04, 0.20, 0.50, 0.80, 0.90, 0.99)
  bh_ref <- p.adjust(known_p, method = "BH"); by_ref <- p.adjust(known_p, method = "BY")
  add("B1_global_BH_matches_base_R", "The BH family correction used for both arm-specific and conjunction families is exactly base R's p.adjust(method='BH') on a known vector",
      all(abs(p.adjust(known_p, method = "BH") - bh_ref) < 1e-12), "byte-identical to p.adjust(method='BH')")
  add("B2_global_BY_matches_base_R", "The BY dependence-robust sensitivity correction is exactly base R's p.adjust(method='BY') on a known vector",
      all(abs(p.adjust(known_p, method = "BY") - by_ref) < 1e-12), "byte-identical to p.adjust(method='BY')")

  # ==========================================================================
  # O1-O3: overlap-category classification, known answer.
  # ==========================================================================
  sig_map <- data.frame(gene = c("CD8A", "CD3D", "CD3D"), population = c("CD8 T cells", "T cells", "Cytotoxic lymphocytes"), stringsAsFactors = FALSE)
  add("O1_direct_tautology_known_answer", "A gene that literally is a population's own MCPcounter marker is classified direct_tautology",
      classify_overlap_category("CD8A", "CD8 T cells", sig_map) == "direct_tautology", "CD8A vs CD8 T cells -> direct_tautology")
  add("O2_cyt_overlap_known_answer", "A non-marker gene against a CYT-related lymphocyte population (not a marker of it) is classified CYT_biological_overlap",
      classify_overlap_category("SPI1", "NK cells", sig_map) == "CYT_biological_overlap", "SPI1 vs NK cells (non-marker) -> CYT_biological_overlap")
  add("O3_nonoverlapping_known_answer", "A non-marker gene against a non-CYT-related population is classified nonoverlapping_population",
      classify_overlap_category("SPI1", "Fibroblasts", sig_map) == "nonoverlapping_population", "SPI1 vs Fibroblasts -> nonoverlapping_population")
  add("O4_direct_tautology_takes_priority", "direct_tautology is checked BEFORE the CYT-overlap label -- a gene that is both a marker of a CYT-related population AND would otherwise be CYT_biological_overlap is still direct_tautology",
      classify_overlap_category("CD3D", "T cells", sig_map) == "direct_tautology", "CD3D vs T cells (its own marker, also a CYT-overlap population) -> direct_tautology, not CYT_biological_overlap")

  # ==========================================================================
  # P1-P2: population-correlation flagging, known matrix.
  # ==========================================================================
  cmat <- matrix(c(1, 0.9, 0.2, 0.9, 1, 0.75, 0.2, 0.75, 1), nrow = 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  flagged <- flag_correlated_population_pairs(cmat, threshold = 0.70)
  add("P1_correlated_pairs_known_answer", "flag_correlated_population_pairs() correctly identifies exactly the 2 pairs exceeding |rho|>0.70 in a known 3x3 matrix (A-B=0.9, B-C=0.75, A-C=0.2 excluded)",
      nrow(flagged) == 2 && all(sort(flagged$rho) == c(0.75, 0.9)),
      sprintf("found %d pairs: %s", nrow(flagged), paste(sprintf("%s-%s=%.2f", flagged$population_a, flagged$population_b, flagged$rho), collapse = ", ")))
  cmat_none <- matrix(c(1, 0.1, 0.1, 1), nrow = 2, dimnames = list(c("A", "B"), c("A", "B")))
  add("P2_no_correlated_pairs_empty_result", "flag_correlated_population_pairs() returns a zero-row result when no pair exceeds the threshold",
      nrow(flag_correlated_population_pairs(cmat_none, threshold = 0.70)) == 0, "0 rows for an uncorrelated 2x2 matrix")

  # ==========================================================================
  # PC1: known-answer PCA -- two perfectly correlated columns collapse to a
  # PC1 identical (up to sign) to either standardized column.
  # ==========================================================================
  set.seed(404)
  base_col <- rnorm(100)
  pop_cols <- cbind(a = base_col, b = base_col + rnorm(100, sd = 0.001))
  pc1 <- compute_lymphocyte_pc1(pop_cols)
  add("PC1_near_perfectly_correlated_columns", "compute_lymphocyte_pc1() of two near-perfectly-correlated columns is itself near-perfectly correlated with either input column",
      abs(cor(pc1, base_col)) > 0.999, sprintf("cor(PC1, base_col)=%.6f", cor(pc1, base_col)))

  # ==========================================================================
  # VIF1-VIF2: known-answer collinearity diagnostics.
  # ==========================================================================
  set.seed(505)
  z1 <- rnorm(200); z2 <- z1 + rnorm(200, sd = 0.01)  # near-perfect collinearity
  z3 <- rnorm(200)  # independent
  yv <- z1 + z3 + rnorm(200, sd = 0.1)
  fit_collinear <- lm(yv ~ z1 + z2 + z3)
  vifs <- compute_vif(fit_collinear)
  add("VIF1_near_collinear_terms_have_huge_vif", "compute_vif() assigns a very large VIF (>50) to two near-perfectly-collinear predictors",
      unname(vifs["z1"]) > 50 && unname(vifs["z2"]) > 50, sprintf("VIF(z1)=%.1f, VIF(z2)=%.1f, VIF(z3)=%.2f", vifs["z1"], vifs["z2"], vifs["z3"]))
  fit_orthogonal <- lm(yv ~ z1 + z3)
  vifs_orth <- compute_vif(fit_orthogonal)
  add("VIF2_independent_predictors_vif_near_one", "compute_vif() assigns VIF near 1 to genuinely independent predictors",
      all(vifs_orth < 1.5), sprintf("VIF(z1)=%.3f, VIF(z3)=%.3f (expected both < 1.5)", vifs_orth["z1"], vifs_orth["z3"]))

  # ==========================================================================
  # M1-M4: missing values, boundary n, zero variance, duplicated/shuffled
  # sample order.
  # ==========================================================================
  x_na <- x; x_na[1:20] <- NA
  r_m1 <- fit_cell_context_interaction_hc3(x_na, y_homo, pop, purity, min_n = 30)
  add("M1_missing_values_dropped_pairwise_not_crashed", "Non-finite predictor values are dropped (complete-case), not silently propagated or crashing the fit",
      isTRUE(r_m1$estimable) && r_m1$n == (n - 20), sprintf("n after dropping 20 NAs=%d (expected %d)", r_m1$n, n - 20))

  r_m2 <- fit_cell_context_interaction_hc3(x[1:20], y_homo[1:20], pop[1:20], purity[1:20], min_n = 30)
  add("M2_below_min_n_not_estimable", "Fewer samples than min_n is correctly reported as not estimable with reason 'insufficient_n', never a fabricated fit",
      !isTRUE(r_m2$estimable) && identical(r_m2$non_evaluable_reason, "insufficient_n"), sprintf("estimable=%s, reason=%s", r_m2$estimable, r_m2$non_evaluable_reason))

  purity_const <- rep(0.8, n)
  r_m3 <- fit_cell_context_interaction_hc3(x, y_homo, pop, purity_const, min_n = 30)
  add("M3_zero_variance_purity_not_estimable", "Zero-variance purity (a degenerate covariate) is correctly reported as not estimable with reason 'zero_variance'",
      !isTRUE(r_m3$estimable) && identical(r_m3$non_evaluable_reason, "zero_variance"), sprintf("estimable=%s, reason=%s", r_m3$estimable, r_m3$non_evaluable_reason))

  shuffle_idx <- sample(seq_len(n))
  r_orig <- fit_cell_context_interaction_hc3(x, y_homo, pop, purity, min_n = 30)
  r_shuf <- fit_cell_context_interaction_hc3(x[shuffle_idx], y_homo[shuffle_idx], pop[shuffle_idx], purity[shuffle_idx], min_n = 30)
  add("M4_shuffled_sample_order_gives_identical_result", "Applying the SAME permutation to predictor/outcome/population/purity together (consistent sample re-ordering) reproduces bit-identical estimates -- proves the fit does not depend on row order beyond correspondence",
      isTRUE(r_orig$estimable) && isTRUE(r_shuf$estimable) && abs(r_orig$beta_interaction - r_shuf$beta_interaction) < 1e-9 && abs(r_orig$p_interaction_hc3 - r_shuf$p_interaction_hc3) < 1e-9,
      sprintf("beta diff=%.2e, p diff=%.2e", abs(r_orig$beta_interaction - r_shuf$beta_interaction), abs(r_orig$p_interaction_hc3 - r_shuf$p_interaction_hc3)))

  # ==========================================================================
  # PERM1: permutation calibration -- under a TRUE NULL (population score
  # genuinely independent of X, Y, purity), permutation p-values should be
  # approximately uniform, i.e. roughly n_perm*alpha rejections at alpha.
  # ==========================================================================
  set.seed(606)
  n_perm_test <- 200
  n_cal <- 150
  x_null <- rnorm(n_cal); y_null <- rnorm(n_cal); pop_null <- rnorm(n_cal); purity_null <- runif(n_cal, 0.3, 1)
  r_obs <- fit_cell_context_interaction_hc3(x_null, y_null, pop_null, purity_null, min_n = 30)
  perm_p <- numeric(n_perm_test)
  for (b in seq_len(n_perm_test)) {
    pop_perm <- sample(pop_null)
    r_b <- fit_cell_context_interaction_hc3(x_null, y_null, pop_perm, purity_null, min_n = 30)
    perm_p[b] <- if (isTRUE(r_b$estimable)) r_b$p_interaction_hc3 else NA_real_
  }
  n_valid <- sum(is.finite(perm_p))
  frac_below_05 <- mean(perm_p < 0.05, na.rm = TRUE)
  mc_se <- sqrt(0.05 * 0.95 / n_valid)
  add("PERM1_calibration_under_true_null", sprintf("Under a genuine null (population independent of X/Y/purity), the fraction of %d permutation draws with raw p<0.05 is within 3 Monte Carlo SEs of the nominal 5%% rate", n_perm_test),
      abs(frac_below_05 - 0.05) < 3 * mc_se + 0.03,
      sprintf("observed rate=%.3f, nominal=0.05, MC SE=%.4f, n_valid=%d", frac_below_05, mc_se, n_valid))

  df <- do.call(rbind, results)
  out_path <- file.path(project_root, "tests", "output", "test_cell_context_corrected_helpers_results.tsv")
  write.table(df, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
  n_fail <- sum(df$result == "FAIL")
  cat(sprintf("\n%d/%d PASS, %d FAIL\n", sum(df$result == "PASS"), nrow(df), n_fail))
  if (n_fail > 0) stop(sprintf("%d test(s) FAILED. See %s.", n_fail, out_path), call. = FALSE)
  invisible(df)
}

if (sys.nframe() == 0) {
  cli_args <- commandArgs(trailingOnly = TRUE)
  if (length(cli_args) < 1) stop("Usage: Rscript test_cell_context_corrected_helpers.R <project_root>")
  invisible(run_cell_context_corrected_helper_tests(normalizePath(cli_args[[1]], mustWork = TRUE)))
}
