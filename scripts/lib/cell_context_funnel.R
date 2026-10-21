# ============================================================================
# cell_context_funnel.R
#
# Reporting helpers for Stage 24
# (scripts/23_cancer_specific_cellular_context_corrected.R). These functions
# perform no statistical test: they combine the tables Stage 24 has already
# computed into
#   - the primary corrected evidence funnel, which starts from every
#     conjunction-BH-supported (triplet, cancer, population) combination;
#   - the historical Stage 23 reconciliation, which starts from Stage 23's
#     all-three-arms-significant combinations and is a comparison only;
#   - the row-level final result table and its triplet-level summary.
# Every funnel step is derived from the immediately preceding step by a
# keyed semi-join on (mirna_id, gene, cancer, population).
# ============================================================================

suppressPackageStartupMessages(library(data.table))

STAGE24_KEY_COLS <- c("mirna_id", "gene", "cancer", "population")

STAGE24_PRIMARY_STEPS <- c(
  "1_corrected_conjunction_BH_supported" = "conjunction FDR(BH) below the cutoff (all_three_arms_supported)",
  "2_LUMP_supported_with_direction_agreement" = "previous step AND lump_sensitivity_pass (LUMP conjunction FDR(BH) supported, three-arm direction agreement with CPE, no leverage-unstable arm)",
  "3_nonoverlapping_population" = "previous step AND overlap_category == nonoverlapping_population",
  "4_collinearity_stable_or_not_applicable" = "previous step AND (coefficient_stable OR collinearity_assessment == not_applicable_no_correlated_population)",
  "5_permutation_supported" = "previous step AND passes_permutation_criterion"
)

STAGE24_HISTORICAL_STEPS <- c(
  "1_original_stage23_all_three_arms" = "Stage 23 all_three_arms_fdr_significant (historical rule, three separately adjusted arm FDRs)",
  "2_HC3_evaluable" = "previous step AND a finite Stage 24 conjunction p-value",
  "3_corrected_conjunction_BH_supported" = "previous step AND all_three_arms_supported",
  "4_LUMP_supported_with_direction_agreement" = "previous step AND lump_sensitivity_pass",
  "5_nonoverlapping_population" = "previous step AND overlap_category == nonoverlapping_population",
  "6_collinearity_stable_or_not_applicable" = "previous step AND (coefficient_stable OR not_applicable_no_correlated_population)",
  "7_permutation_supported" = "previous step AND passes_permutation_criterion"
)

# Checks that a source table has the key and required columns and unique
# keys. A zero-row table without columns (written by Stage 24 when nothing
# was checked) is given the missing columns so that it filters to zero rows.
stage24_prepare_source <- function(dt, name, required_cols, key_cols = STAGE24_KEY_COLS) {
  dt <- as.data.table(dt)
  needed <- c(key_cols, required_cols)
  if (nrow(dt) == 0) {
    for (col in setdiff(needed, names(dt))) dt[, (col) := logical(0)]
  }
  missing <- setdiff(needed, names(dt))
  if (length(missing) > 0) {
    stop(sprintf("%s lacks required column(s): %s", name, paste(missing, collapse = ", ")))
  }
  n_dup <- sum(duplicated(dt[, ..key_cols]))
  if (n_dup > 0) stop(sprintf("%s has %d duplicated (%s) keys", name, n_dup, paste(key_cols, collapse = ", ")))
  dt
}

# Rows of `prev` whose key appears in `keep`. Fails if the result duplicates
# a key, gains rows, or contains a key absent from `prev`.
stage24_semi_join <- function(prev, keep, step, key_cols = STAGE24_KEY_COLS) {
  keep <- unique(as.data.table(keep)[, ..key_cols])
  if (nrow(keep) == 0) return(prev[0])
  out <- prev[keep, on = key_cols, nomatch = NULL]
  if (anyDuplicated(out[, ..key_cols])) stop(sprintf("Funnel step %s introduced duplicated keys", step))
  if (nrow(out) > nrow(prev)) stop(sprintf("Funnel step %s has more rows (%d) than its previous step (%d)", step, nrow(out), nrow(prev)))
  absent <- out[!prev, on = key_cols]
  if (nrow(absent) > 0) stop(sprintf("Funnel step %s introduced %d key(s) absent from the previous step", step, nrow(absent)))
  out
}

stage24_funnel_stats <- function(dt) {
  if (nrow(dt) == 0) {
    return(data.table(n_rows = 0L, n_triplets = 0L, n_mirnas = 0L, n_genes = 0L, n_cancers = 0L, n_populations = 0L))
  }
  data.table(n_rows = nrow(dt), n_triplets = uniqueN(dt[, .(mirna_id, gene)]), n_mirnas = uniqueN(dt$mirna_id),
             n_genes = uniqueN(dt$gene), n_cancers = uniqueN(dt$cancer), n_populations = uniqueN(dt$population))
}

stage24_funnel_table <- function(steps, definitions) {
  stopifnot(identical(names(steps), names(definitions)))
  funnel <- rbindlist(lapply(names(steps), function(s) {
    cbind(data.table(step = s, step_definition = unname(definitions[[s]])), stage24_funnel_stats(steps[[s]]))
  }))
  if (any(diff(funnel$n_rows) > 0)) stop("Funnel counts are not non-increasing")
  funnel
}

# Robustness steps shared by both funnels, applied in order to `start`
# (which must carry overlap_category). Each step filters the previous one.
stage24_robustness_steps <- function(start, sources, step_names) {
  k <- STAGE24_KEY_COLS
  s_lump <- stage24_semi_join(start, sources$purity_stage_dt[lump_sensitivity_pass %in% TRUE, ..k], step_names[1])
  s_nonoverlap <- s_lump[overlap_category == "nonoverlapping_population"]
  s_collin <- stage24_semi_join(s_nonoverlap,
    sources$collinearity_dt[coefficient_stable %in% TRUE |
                              collinearity_assessment %in% "not_applicable_no_correlated_population", ..k],
    step_names[3])
  s_perm <- stage24_semi_join(s_collin, sources$negctrl_dt[passes_permutation_criterion %in% TRUE, ..k], step_names[4])
  setNames(list(s_lump, s_nonoverlap, s_collin, s_perm), step_names)
}

stage24_prepare_funnel_sources <- function(conjunction_dt, purity_stage_dt, collinearity_dt, negctrl_dt) {
  list(
    conjunction_dt = stage24_prepare_source(conjunction_dt, "Stage 24 conjunction table",
                                            c("all_three_arms_supported", "overlap_category", "p_conjunction")),
    purity_stage_dt = stage24_prepare_source(purity_stage_dt, "Stage 24 LUMP/stage sensitivity table", "lump_sensitivity_pass"),
    collinearity_dt = stage24_prepare_source(collinearity_dt, "Stage 24 collinearity table",
                                             c("coefficient_stable", "collinearity_assessment")),
    negctrl_dt = stage24_prepare_source(negctrl_dt, "Stage 24 negative-control table", "passes_permutation_criterion")
  )
}

# Primary corrected evidence funnel: starts from every
# conjunction-BH-supported combination in the full corrected universe.
build_stage24_primary_funnel <- function(conjunction_dt, purity_stage_dt, collinearity_dt, negctrl_dt) {
  src <- stage24_prepare_funnel_sources(conjunction_dt, purity_stage_dt, collinearity_dt, negctrl_dt)
  step_names <- names(STAGE24_PRIMARY_STEPS)
  s1 <- src$conjunction_dt[all_three_arms_supported %in% TRUE, c(STAGE24_KEY_COLS, "overlap_category"), with = FALSE]
  steps <- c(setNames(list(s1), step_names[1]), stage24_robustness_steps(s1, src, step_names[2:5]))
  list(steps = steps, funnel = stage24_funnel_table(steps, STAGE24_PRIMARY_STEPS), final_keys = steps[[step_names[5]]])
}

# Historical Stage 23 reconciliation: follows Stage 23's own
# all-three-arms-significant set through the same corrected checks. It does
# not define the Stage 24 result.
build_stage24_historical_funnel <- function(historical_stage23_supported, conjunction_dt, purity_stage_dt,
                                            collinearity_dt, negctrl_dt) {
  src <- stage24_prepare_funnel_sources(conjunction_dt, purity_stage_dt, collinearity_dt, negctrl_dt)
  hist <- stage24_prepare_source(historical_stage23_supported, "Historical Stage 23 all-three-arms set", character(0))
  step_names <- names(STAGE24_HISTORICAL_STEPS)
  h1 <- hist[, ..STAGE24_KEY_COLS]
  h2 <- stage24_semi_join(h1, src$conjunction_dt[!is.na(p_conjunction), ..STAGE24_KEY_COLS], step_names[2])
  h3 <- stage24_semi_join(h2, src$conjunction_dt[all_three_arms_supported %in% TRUE, ..STAGE24_KEY_COLS], step_names[3])
  h3 <- merge(h3, src$conjunction_dt[, c(STAGE24_KEY_COLS, "overlap_category"), with = FALSE], by = STAGE24_KEY_COLS, sort = FALSE)
  steps <- c(setNames(list(h1, h2, h3), step_names[1:3]), stage24_robustness_steps(h3, src, step_names[4:7]))
  list(steps = steps, funnel = stage24_funnel_table(steps, STAGE24_HISTORICAL_STEPS), final_keys = steps[[step_names[7]]])
}

# Left-joins the per-combination columns of one source table onto `base`,
# requiring a one-to-one key match.
stage24_attach <- function(base, src, cols, name) {
  src <- stage24_prepare_source(src, name, cols)
  if (nrow(src) == 0) {
    out <- copy(base)
    for (col in cols) out[, (col) := NA]
    return(out)
  }
  out <- merge(base, src[, c(STAGE24_KEY_COLS, cols), with = FALSE], by = STAGE24_KEY_COLS, all.x = TRUE, sort = FALSE)
  if (nrow(out) != nrow(base)) stop(sprintf("Attaching %s changed the final row count (%d -> %d)", name, nrow(base), nrow(out)))
  out
}

# Row-level final result: one row per final (mirna_id, gene, cancer,
# population) key with the evidence and sensitivity columns of every
# Stage 24 table.
build_stage24_final_table <- function(final_keys, conjunction_dt, purity_stage_dt, collinearity_dt, negctrl_dt,
                                      nonlinearity_dt, historical_stage23_supported, conj_alpha) {
  base <- unique(as.data.table(final_keys)[, ..STAGE24_KEY_COLS])
  if (nrow(base) != nrow(final_keys)) stop("Final keys are not unique")
  out <- stage24_attach(base, conjunction_dt,
                        c("overlap_category", "n_samples", "any_arm_leverage_unstable", "p_conjunction",
                          "conjunction_FDR_BH", "conjunction_FDR_BY", "all_three_arms_supported"),
                        "Stage 24 conjunction table")
  out <- stage24_attach(out, purity_stage_dt,
                        c("lump_sensitivity_pass", "lump_conjunction_fdr_bh", "lump_direction_agreement",
                          "stage_sensitivity_applicable", "stage_beta_interaction", "stage_p_interaction"),
                        "Stage 24 LUMP/stage sensitivity table")
  out <- stage24_attach(out, collinearity_dt,
                        c("coefficient_stable", "collinearity_assessment", "focal_adds_beyond_pc1"),
                        "Stage 24 collinearity table")
  out <- stage24_attach(out, negctrl_dt,
                        c("n_permutations", "n_valid_permutations", "empirical_p_permutation", "monte_carlo_se",
                          "passes_permutation_criterion"),
                        "Stage 24 negative-control table")
  out <- stage24_attach(out, nonlinearity_dt, c("nonlinearity_flagged", "nonlinear_vs_linear_p"),
                        "Stage 24 nonlinearity table")
  hist <- stage24_prepare_source(historical_stage23_supported, "Historical Stage 23 all-three-arms set", character(0))
  out[, in_historical_stage23 := FALSE]
  if (nrow(hist) > 0 && nrow(out) > 0) out[hist, on = STAGE24_KEY_COLS, in_historical_stage23 := TRUE]
  out[, conjunction_BY_supported := !is.na(conjunction_FDR_BY) & conjunction_FDR_BY < conj_alpha]
  # With zero permutations at least as extreme, (0 + 1) / (n_valid + 1) is
  # the smallest attainable empirical p-value.
  out[, permutation_p_at_minimum_attainable :=
        is.finite(empirical_p_permutation) & n_valid_permutations > 0 &
        abs(empirical_p_permutation - 1 / (n_valid_permutations + 1)) < 1e-12]
  setcolorder(out, c(STAGE24_KEY_COLS, "overlap_category", "n_samples", "p_conjunction", "conjunction_FDR_BH",
                     "conjunction_FDR_BY", "conjunction_BY_supported", "lump_sensitivity_pass", "coefficient_stable",
                     "collinearity_assessment", "focal_adds_beyond_pc1", "empirical_p_permutation",
                     "passes_permutation_criterion", "permutation_p_at_minimum_attainable", "nonlinearity_flagged",
                     "nonlinear_vs_linear_p", "stage_sensitivity_applicable", "stage_beta_interaction",
                     "stage_p_interaction", "in_historical_stage23"))
  setorderv(out, c("cancer", "population", "mirna_id", "gene"))
  out[]
}

# Named logical vector: each element is one property every final row must
# satisfy.
validate_stage24_final_table <- function(final_dt) {
  c(
    keys_unique = !anyDuplicated(final_dt[, ..STAGE24_KEY_COLS]),
    conjunction_supported = all(final_dt$all_three_arms_supported %in% TRUE),
    lump_sensitivity_pass = all(final_dt$lump_sensitivity_pass %in% TRUE),
    nonoverlapping_population = all(final_dt$overlap_category %in% "nonoverlapping_population"),
    collinearity_stable_or_not_applicable = all(final_dt$coefficient_stable %in% TRUE |
      final_dt$collinearity_assessment %in% "not_applicable_no_correlated_population"),
    permutation_supported = all(final_dt$passes_permutation_criterion %in% TRUE)
  )
}

# Triplet-level summary of the final result, aggregated by (mirna_id, gene).
summarize_stage24_final_triplets <- function(final_dt, stage_alpha = 0.05) {
  final_dt[, .(
    n_cancer_population_rows = .N,
    n_cancers = uniqueN(cancer),
    cancers = paste(sort(unique(cancer)), collapse = ";"),
    n_populations = uniqueN(population),
    populations = paste(sort(unique(population)), collapse = ";"),
    min_conjunction_FDR_BH = min(conjunction_FDR_BH),
    max_conjunction_FDR_BH = max(conjunction_FDR_BH),
    n_conjunction_BY_supported = sum(conjunction_BY_supported %in% TRUE),
    n_nonlinearity_flagged = sum(nonlinearity_flagged %in% TRUE),
    n_stage_evaluable = sum(stage_sensitivity_applicable %in% TRUE),
    n_stage_nominal_supported = sum(stage_sensitivity_applicable %in% TRUE & stage_p_interaction < stage_alpha, na.rm = TRUE),
    n_in_historical_stage23 = sum(in_historical_stage23 %in% TRUE)
  ), by = .(mirna_id, gene)][order(-n_cancer_population_rows, mirna_id, gene)]
}
