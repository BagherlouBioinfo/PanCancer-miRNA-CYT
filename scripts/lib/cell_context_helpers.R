# ============================================================================
# cell_context_helpers.R
#
# Helpers for the cancer-specific cellular-composition context analyses
# (Stages 23 and 24). Every function is a pure function of its arguments,
# with no file path, run context or project-specific literal, so each can
# be unit-tested with synthetic data (tests/test_cell_context_helpers.R,
# tests/test_cell_context_conjunction_helpers.R).
#
# Stage 23 design:
#
#   Primary test -- continuous effect modification. For one arm (a
#   predictor/outcome pair, e.g. miRNA residual vs. CYT residual) within
#   one cancer, with one MCPcounter population's standardized abundance as
#   the candidate effect modifier, fit
#       rank(outcome) ~ rank(predictor) * population_z
#   by ordinary least squares on rank-transformed residualized values. The
#   rank-based interaction model is consistent with the Spearman-based
#   discovery analysis (Spearman's rho is the Pearson correlation of ranks)
#   while allowing a continuous interaction term. `outcome` and `predictor`
#   are expected to be CPE-purity-residualized already (the
#   residualize_matrix() output used by Stage 10), computed once per
#   (triplet, cancer) and reused across populations; the function itself
#   performs no purity adjustment. The population is not dichotomized.
#
#   Descriptive summary -- tertile strata. Samples of the same cancer are
#   split into low/middle/high tertiles of the population's raw abundance,
#   and the Spearman correlation (with a Fisher-z confidence interval) of
#   the same residualized predictor and outcome is reported within each
#   stratum. This summary is interpretive only and is not a significance
#   test.
# ============================================================================

# Fits the Stage 23 continuous-interaction model for one arm x cancer x
# population combination. `predictor` and `outcome` must be CPE-residualized
# values of equal length and sample order with no NA (the caller restricts
# to complete cases); `population_raw` is the population's raw abundance for
# the same samples.
#
# Returns a list:
#   $estimable (logical) -- FALSE if n < min_n, any of predictor, outcome or
#     population has zero variance, or the fit is singular or yields a
#     non-finite interaction estimate; all other fields are then NA.
#   $n, $beta_main, $beta_interaction, $se_interaction, $p_interaction,
#   $ci_lower, $ci_upper (95% CI for the interaction coefficient),
#   $slope_low, $slope_high -- model-implied conditional slopes of
#     rank(outcome) on rank(predictor) at the 10th and 90th percentiles of
#     population_z (beta_main + beta_interaction * population_z). Callers
#     use them to classify enriched (|slope_high| > |slope_low|) versus
#     depleted (|slope_high| < |slope_low|) direction, because the sign of
#     beta_interaction alone is ambiguous without the sign of the baseline
#     association.
fit_cell_context_interaction <- function(predictor, outcome, population_raw, min_n = 30) {
  ok <- is.finite(predictor) & is.finite(outcome) & is.finite(population_raw)
  predictor <- predictor[ok]; outcome <- outcome[ok]; population_raw <- population_raw[ok]
  n <- length(predictor)
  not_estimable <- function() list(estimable = FALSE, n = n, beta_main = NA_real_, beta_interaction = NA_real_,
                                    se_interaction = NA_real_, p_interaction = NA_real_, ci_lower = NA_real_,
                                    ci_upper = NA_real_, slope_low = NA_real_, slope_high = NA_real_)
  if (n < min_n) return(not_estimable())
  if (sd(predictor) == 0 || sd(outcome) == 0 || sd(population_raw) == 0) return(not_estimable())

  r_pred <- rank(predictor); r_out <- rank(outcome)
  pop_z <- as.numeric(scale(population_raw))

  fit <- tryCatch(lm(r_out ~ r_pred * pop_z), error = function(e) NULL)
  if (is.null(fit)) return(not_estimable())
  co <- summary(fit)$coefficients
  int_term <- "r_pred:pop_z"
  if (!(int_term %in% rownames(co)) || !("r_pred" %in% rownames(co))) return(not_estimable())
  beta_main <- co["r_pred", "Estimate"]
  beta_int <- co[int_term, "Estimate"]
  se_int <- co[int_term, "Std. Error"]
  p_int <- co[int_term, "Pr(>|t|)"]
  if (!is.finite(beta_int) || !is.finite(se_int) || !is.finite(p_int)) return(not_estimable())

  ci <- tryCatch(confint(fit, int_term, level = 0.95), error = function(e) matrix(c(NA_real_, NA_real_), nrow = 1))

  pop_low <- as.numeric(quantile(pop_z, 0.10, na.rm = TRUE, type = 7))
  pop_high <- as.numeric(quantile(pop_z, 0.90, na.rm = TRUE, type = 7))
  slope_low <- beta_main + beta_int * pop_low
  slope_high <- beta_main + beta_int * pop_high

  list(estimable = TRUE, n = n, beta_main = beta_main, beta_interaction = beta_int,
       se_interaction = se_int, p_interaction = p_int,
       ci_lower = ci[1, 1], ci_upper = ci[1, 2],
       slope_low = slope_low, slope_high = slope_high)
}

# Descriptive tertile-stratified Spearman correlation. Splits samples into
# low/middle/high tertiles of `population_raw` within the cancer
# (probs = c(0, 1/3, 2/3, 1), quantile type 7) and computes the Spearman
# correlation with a Fisher-z 95% CI of `predictor` vs. `outcome` within
# each stratum. A stratum with fewer than `min_n_stratum` samples (default
# 10, smaller than the primary min_n = 30 because these subsets are
# descriptive only) is reported as not estimable.
tertile_stratified_correlation <- function(predictor, outcome, population_raw, min_n_stratum = 10) {
  ok <- is.finite(predictor) & is.finite(outcome) & is.finite(population_raw)
  predictor <- predictor[ok]; outcome <- outcome[ok]; population_raw <- population_raw[ok]
  n_total <- length(predictor)

  one_stratum <- function(idx, label) {
    n_s <- length(idx)
    if (n_s < min_n_stratum || sd(predictor[idx]) == 0 || sd(outcome[idx]) == 0) {
      return(data.frame(stratum = label, n = n_s, rho = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                         p_value = NA_real_, estimable = FALSE, stringsAsFactors = FALSE))
    }
    rho <- suppressWarnings(cor(predictor[idx], outcome[idx], method = "spearman"))
    test <- suppressWarnings(cor.test(predictor[idx], outcome[idx], method = "spearman", exact = FALSE))
    z <- atanh(pmin(pmax(rho, -0.999999), 0.999999))
    se <- 1 / sqrt(n_s - 3)
    ci_lo <- tanh(z - 1.96 * se); ci_hi <- tanh(z + 1.96 * se)
    data.frame(stratum = label, n = n_s, rho = rho, ci_lower = ci_lo, ci_upper = ci_hi,
               p_value = test$p.value, estimable = TRUE, stringsAsFactors = FALSE)
  }

  if (n_total < 3 * min_n_stratum) {
    return(rbind(one_stratum(integer(0), "low"), one_stratum(integer(0), "middle"), one_stratum(integer(0), "high")))
  }
  qs <- quantile(population_raw, probs = c(1/3, 2/3), na.rm = TRUE, type = 7)
  low_idx <- which(population_raw <= qs[1])
  high_idx <- which(population_raw > qs[2])
  mid_idx <- setdiff(seq_len(n_total), union(low_idx, high_idx))
  rbind(one_stratum(low_idx, "low"), one_stratum(mid_idx, "middle"), one_stratum(high_idx, "high"))
}

# Tautology safeguard: TRUE if `gene` is one of `population`'s MCPcounter
# marker genes. `signature_map` is a data.frame/data.table with columns
# `gene` and `population`, one row per marker gene, read by the caller from
# reference/mcpcounter_signatures/genes.txt. When TRUE, any arm involving
# `gene` against this population (gene-CYT, miRNA-gene) is circular, because
# the population's abundance estimate is partly constructed from this
# gene's expression; such arms are marked tautological and not evaluated.
is_gene_signature_member <- function(gene, population, signature_map) {
  any(toupper(signature_map$gene) == toupper(gene) & signature_map$population == population)
}

# Combines up to three per-arm interaction results for one
# (triplet, cancer, population) into one verdict:
#   "tautological_exclusion", "not_evaluable", "composition_enriched",
#   "composition_depleted", "stable_across_composition", "not_detected".
#
# `arm_results`: named list with any subset of "mirna_cyt", "gene_cyt",
#   "mirna_gene", each the return value of fit_cell_context_interaction()
#   plus an `fdr` field (the interaction p-value BH-corrected by the caller
#   across the run-level interaction-test family).
# `tautological_arms`: arm names (subset of "gene_cyt", "mirna_gene") that
#   are tautological for this population, from is_gene_signature_member().
# `baseline_pass`: named logical vector (one per arm) indicating whether the
#   arm's population-independent association meets the primary rule
#   (|rho| >= threshold and FDR < 0.05) in this cancer, taken from Stage 10
#   output. Used only to distinguish stable_across_composition from
#   not_detected when no interaction is significant.
# `fdr_alpha`: significance threshold for the interaction FDR (0.05).
classify_cell_context <- function(arm_results, tautological_arms, baseline_pass, fdr_alpha = 0.05) {
  required_arms <- c("mirna_cyt", "gene_cyt", "mirna_gene")

  if (length(tautological_arms) > 0) {
    return("tautological_exclusion")
  }

  evaluable <- vapply(required_arms, function(a) {
    r <- arm_results[[a]]
    !is.null(r) && isTRUE(r$estimable)
  }, logical(1))
  if (!all(evaluable)) return("not_evaluable")

  fdrs <- vapply(required_arms, function(a) arm_results[[a]]$fdr, numeric(1))
  directions <- vapply(required_arms, function(a) {
    r <- arm_results[[a]]
    if (!is.finite(r$slope_low) || !is.finite(r$slope_high)) return(NA_character_)
    if (abs(r$slope_high) > abs(r$slope_low)) "enriched" else if (abs(r$slope_high) < abs(r$slope_low)) "depleted" else "flat"
  }, character(1))

  sig_arms <- fdrs < fdr_alpha
  if (any(sig_arms, na.rm = TRUE)) {
    sig_directions <- directions[sig_arms & !is.na(directions)]
    n_enriched <- sum(sig_directions == "enriched")
    n_depleted <- sum(sig_directions == "depleted")
    if (n_enriched > n_depleted) return("composition_enriched")
    if (n_depleted > n_enriched) return("composition_depleted")
    # Tied directions (e.g. one significant arm in each direction, or all
    # flat): fall through to the baseline check rather than assigning a
    # direction without majority support.
  }

  baseline_detected <- all(vapply(required_arms, function(a) isTRUE(baseline_pass[[a]]), logical(1)))
  if (baseline_detected) return("stable_across_composition")
  "not_detected"
}

# ============================================================================
# Conjunction methodology used by Stage 24
# (scripts/23_cancer_specific_cellular_context_corrected.R). Stage 23 uses
# only the functions above.
# ============================================================================

# Interaction model for one arm x cancer x population combination:
#   rank(outcome) ~ rank(predictor) + population_z + rank(predictor):population_z + purity
# fit on the winsorized, matched values without prior purity
# residualization; purity enters as an explicit covariate, so the
# population score is also purity-adjusted. Inference for the interaction
# coefficient uses HC3 heteroscedasticity-consistent standard errors
# (sandwich::vcovHC(type = "HC3")); classical OLS standard errors and
# p-values are reported only as a diagnostic comparison.
#
# `predictor_raw`, `outcome_raw`, `population_raw` and `purity` must have
# equal length and sample order; non-finite entries are dropped jointly
# (complete-case analysis) before fitting. Returns estimable = FALSE with a
# `non_evaluable_reason` if n < min_n, any input has zero within-cancer
# variance, the model is singular, or the HC3 covariance cannot be
# computed.
fit_cell_context_interaction_hc3 <- function(predictor_raw, outcome_raw, population_raw, purity, min_n = 30) {
  ok <- is.finite(predictor_raw) & is.finite(outcome_raw) & is.finite(population_raw) & is.finite(purity)
  predictor_raw <- predictor_raw[ok]; outcome_raw <- outcome_raw[ok]
  population_raw <- population_raw[ok]; purity <- purity[ok]
  n <- length(predictor_raw)
  n_complete_case <- n

  not_estimable <- function(reason) list(
    estimable = FALSE, non_evaluable_reason = reason, n = n, n_complete_case = n_complete_case,
    beta_main = NA_real_, beta_interaction = NA_real_,
    se_interaction_hc3 = NA_real_, p_interaction_hc3 = NA_real_,
    ci_lower_hc3 = NA_real_, ci_upper_hc3 = NA_real_,
    beta_interaction_ols = NA_real_, se_interaction_ols = NA_real_, p_interaction_ols = NA_real_,
    slope_low = NA_real_, slope_high = NA_real_, df_residual = NA_integer_
  )

  if (n < min_n) return(not_estimable("insufficient_n"))
  if (sd(predictor_raw) == 0 || sd(outcome_raw) == 0 || sd(population_raw) == 0 || sd(purity) == 0) {
    return(not_estimable("zero_variance"))
  }

  r_pred <- rank(predictor_raw); r_out <- rank(outcome_raw)
  pop_z <- as.numeric(scale(population_raw))
  dat <- data.frame(r_out = r_out, r_pred = r_pred, pop_z = pop_z, purity = purity)

  fit <- tryCatch(lm(r_out ~ r_pred * pop_z + purity, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(not_estimable("singular_fit"))
  co <- summary(fit)$coefficients
  int_term <- "r_pred:pop_z"
  if (!(int_term %in% rownames(co)) || !("r_pred" %in% rownames(co)) || fit$rank < ncol(model.matrix(fit))) {
    return(not_estimable("singular_fit"))
  }

  vc <- tryCatch(sandwich::vcovHC(fit, type = "HC3"), error = function(e) NULL)
  if (is.null(vc) || !(int_term %in% rownames(vc)) || !is.finite(vc[int_term, int_term]) || vc[int_term, int_term] <= 0) {
    return(not_estimable("hc3_not_computable"))
  }

  beta_int <- unname(coef(fit)[int_term])
  beta_main <- unname(coef(fit)["r_pred"])
  se_hc3 <- sqrt(vc[int_term, int_term])
  df_resid <- fit$df.residual
  if (!is.finite(beta_int) || !is.finite(se_hc3) || df_resid < 1) return(not_estimable("nonfinite_hc3"))
  t_stat <- beta_int / se_hc3
  p_hc3 <- 2 * pt(-abs(t_stat), df = df_resid)
  tcrit <- qt(0.975, df = df_resid)
  ci_lo <- beta_int - tcrit * se_hc3
  ci_hi <- beta_int + tcrit * se_hc3

  beta_int_ols <- co[int_term, "Estimate"]; se_ols <- co[int_term, "Std. Error"]; p_ols <- co[int_term, "Pr(>|t|)"]

  pop_low <- as.numeric(quantile(pop_z, 0.10, na.rm = TRUE, type = 7))
  pop_high <- as.numeric(quantile(pop_z, 0.90, na.rm = TRUE, type = 7))
  slope_low <- beta_main + beta_int * pop_low
  slope_high <- beta_main + beta_int * pop_high

  list(estimable = TRUE, non_evaluable_reason = NA_character_, n = n, n_complete_case = n_complete_case,
       beta_main = beta_main, beta_interaction = beta_int,
       se_interaction_hc3 = se_hc3, p_interaction_hc3 = p_hc3,
       ci_lower_hc3 = ci_lo, ci_upper_hc3 = ci_hi,
       beta_interaction_ols = unname(beta_int_ols), se_interaction_ols = unname(se_ols), p_interaction_ols = unname(p_ols),
       slope_low = slope_low, slope_high = slope_high, df_residual = df_resid)
}

# Intersection-union (conjunction) test for "all three arms show effect
# modification": p_conjunction = max(p1, p2, p3) of the three arms' raw HC3
# interaction p-values. NA if any of the three is missing or non-finite,
# because a conjunction claim requires all three arms to be tested. Simes'
# test addresses a different hypothesis (at least one non-null arm) and is
# not used for this claim.
compute_conjunction_p <- function(p_mirna_cyt, p_gene_cyt, p_mirna_gene) {
  ps <- c(p_mirna_cyt, p_gene_cyt, p_mirna_gene)
  if (any(!is.finite(ps))) return(NA_real_)
  max(ps)
}

# Overlap classification for one (gene, population) pair:
#   "direct_tautology"          -- gene is an MCPcounter marker of this
#                                  population (is_gene_signature_member());
#                                  these models remain excluded.
#   "CYT_biological_overlap"    -- population is one of the prespecified
#                                  CYT-related lymphocyte populations. The
#                                  MCPcounter reference has no literal marker
#                                  overlap with GZMA/PRF1, but strong
#                                  biological overlap with cytotoxic effector
#                                  activity is expected; results are reported
#                                  as immune-context sensitivity, not as
#                                  independent mechanistic evidence.
#   "nonoverlapping_population" -- neither of the above.
classify_overlap_category <- function(gene, population, signature_map,
                                       cyt_overlap_populations = c("T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells")) {
  if (is_gene_signature_member(gene, population, signature_map)) return("direct_tautology")
  if (population %in% cyt_overlap_populations) return("CYT_biological_overlap")
  "nonoverlapping_population"
}

# Returns a long-format data.frame of population pairs whose absolute
# Spearman correlation exceeds `threshold` (default 0.70, prespecified) in a
# square, named population x population correlation matrix (as produced by
# cor(..., method = "spearman")). Each unordered pair is reported once.
flag_correlated_population_pairs <- function(cor_mat, threshold = 0.70) {
  cor_mat[upper.tri(cor_mat, diag = TRUE)] <- NA
  idx <- which(abs(cor_mat) > threshold, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(data.frame(population_a = character(0), population_b = character(0), rho = numeric(0), stringsAsFactors = FALSE))
  }
  data.frame(population_a = rownames(cor_mat)[idx[, 1]], population_b = colnames(cor_mat)[idx[, 2]],
             rho = cor_mat[idx], stringsAsFactors = FALSE)
}

# First principal component of a set of mutually correlated population
# columns, centred and scaled before PCA: a general lymphocytic-infiltration
# axis used to test whether a focal population's interaction effect adds
# information beyond this shared axis. `pop_cols` is a numeric matrix or
# data.frame, samples x populations (at least two columns).
compute_lymphocyte_pc1 <- function(pop_cols) {
  pop_cols <- as.matrix(pop_cols)
  if (ncol(pop_cols) < 2) stop("compute_lymphocyte_pc1(): need at least 2 population columns to form a shared axis.", call. = FALSE)
  const_cols <- apply(pop_cols, 2, function(x) sd(x) == 0)
  if (any(const_cols)) stop(sprintf("compute_lymphocyte_pc1(): column(s) %s have zero variance.", paste(colnames(pop_cols)[const_cols], collapse = ", ")), call. = FALSE)
  pr <- prcomp(pop_cols, center = TRUE, scale. = TRUE)
  as.numeric(pr$x[, 1])
}

# Variance inflation factors for every non-intercept term of a fitted lm,
# VIF_j = 1 / (1 - R_j^2) from regressing column j of the model matrix on
# the other columns. Returns a named numeric vector.
compute_vif <- function(fit) {
  mm <- model.matrix(fit)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  if (ncol(mm) < 2) return(setNames(rep(1, ncol(mm)), colnames(mm)))
  vifs <- vapply(colnames(mm), function(j) {
    other <- setdiff(colnames(mm), j)
    r2 <- summary(lm(mm[, j] ~ mm[, other, drop = FALSE]))$r.squared
    if (!is.finite(r2) || r2 >= 1) return(Inf)
    1 / (1 - r2)
  }, numeric(1))
  vifs
}

# Condition number of the design matrix excluding the intercept (kappa() on
# the non-intercept columns), a scale-sensitive collinearity diagnostic
# complementary to VIF.
compute_condition_number <- function(fit) {
  mm <- model.matrix(fit)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  kappa(mm, exact = TRUE)
}
