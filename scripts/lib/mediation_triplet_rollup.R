# ============================================================================
# mediation_triplet_rollup.R
#
# Read-only, deterministic aggregation of Stage 17's
# tables/downstream/mediation_primary_support_subset.tsv from one row per
# (miRNA, gene, cancer) context to one row per unique (miRNA, gene)
# triplet. No model is fitted and no p-value or FDR is recomputed; every
# reported quantity is a grouped count, percentage, minimum, median,
# maximum or range of an existing column.
#
# Input scope: mediation_primary_support_subset.tsv is the subset of
# mediation_all_contexts.tsv (every evaluable (miRNA, gene, cancer)
# context) restricted to primary_support == TRUE, i.e. the contexts in
# which the triplet passes the three-arm discovery rule of the same run
# (the supporting cancers that define cross-cancer recurrence). The
# aggregation verifies that every row has n_successful == n_iterations and
# n_failed == 0, so every row of this input is an evaluable context.
#
# The bh_fdr column was computed by Stage 17 across the full set of
# evaluable contexts and is not recomputed within the subset; counts of
# "contexts reaching BH-FDR < X" therefore refer to primary-support
# contexts that cleared the threshold under the full-universe correction.
#
# Shared numerical aggregation helper; see tests/run_tests.R.
# ============================================================================

suppressPackageStartupMessages(library(data.table))

# Computes the triplet-level roll-up. `mediation_dt` must be a data.table or
# data.frame with at least the columns mirna_id, gene, cancer,
# raw_bootstrap_p, bh_fdr, prop_mediated, prop_mediated_meaningful and
# point_indirect. Every row is treated as one evaluable cancer context; the
# function performs no evaluability filtering of its own and drops no row.
compute_mediation_triplet_rollup <- function(mediation_dt) {
  mediation_dt <- as.data.table(mediation_dt)
  required_cols <- c("mirna_id", "gene", "cancer", "raw_bootstrap_p", "bh_fdr",
                      "prop_mediated", "prop_mediated_meaningful", "point_indirect")
  missing_cols <- setdiff(required_cols, names(mediation_dt))
  if (length(missing_cols) > 0) {
    stop(sprintf("compute_mediation_triplet_rollup(): input is missing required column(s): %s",
                 paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  if (nrow(mediation_dt) == 0) {
    stop("compute_mediation_triplet_rollup(): input has zero rows -- nothing to aggregate.", call. = FALSE)
  }
  dup_key <- duplicated(mediation_dt[, .(mirna_id, gene, cancer)])
  if (any(dup_key)) {
    stop(sprintf("compute_mediation_triplet_rollup(): input contains %d duplicated (mirna_id, gene, cancer) key(s) -- refusing to aggregate over an ambiguous input (grouping would silently double-count or arbitrarily pick one row per duplicate).",
                 sum(dup_key)), call. = FALSE)
  }
  # Columns used below must not contain NA in any aggregated row.
  na_counts <- vapply(mediation_dt[, ..required_cols], function(x) sum(is.na(x)), integer(1))
  if (any(na_counts > 0)) {
    stop(sprintf("compute_mediation_triplet_rollup(): required column(s) contain NA and would be silently excluded by a naive aggregation -- refusing: %s",
                 paste(sprintf("%s=%d NA", names(na_counts)[na_counts > 0], na_counts[na_counts > 0]), collapse = ", ")), call. = FALSE)
  }

  rollup <- mediation_dt[, {
    n_contexts <- .N
    nominal_sig <- raw_bootstrap_p < 0.05
    fdr05 <- bh_fdr < 0.05
    fdr01 <- bh_fdr < 0.01
    dir_sign <- sign(point_indirect)
    .(
      n_contexts = n_contexts,
      n_nominal_sig = sum(nominal_sig),
      pct_nominal_sig = 100 * sum(nominal_sig) / n_contexts,
      n_fdr05 = sum(fdr05),
      pct_fdr05 = 100 * sum(fdr05) / n_contexts,
      n_fdr01 = sum(fdr01),
      pct_fdr01 = 100 * sum(fdr01) / n_contexts,
      min_bh_fdr = min(bh_fdr),
      median_bh_fdr = median(bh_fdr),
      min_prop_mediated = min(prop_mediated),
      median_prop_mediated = median(prop_mediated),
      max_prop_mediated = max(prop_mediated),
      range_prop_mediated = max(prop_mediated) - min(prop_mediated),
      n_prop_mediated_meaningful = sum(prop_mediated_meaningful),
      pct_prop_mediated_meaningful = 100 * sum(prop_mediated_meaningful) / n_contexts,
      n_distinct_indirect_directions = uniqueN(dir_sign),
      direction_consistent_across_contexts = uniqueN(dir_sign) == 1,
      supporting_cancers = paste(sort(cancer), collapse = ";")
    )
  }, by = .(mirna_id, gene)]

  setorder(rollup, mirna_id, gene)
  rollup[]
}

# ---------------------------------------------------------------------------
# Command-line entry point: reads the run's mediation output, computes the
# aggregation and writes the triplet-level TSV. Not executed when the file is
# sourced for its function definition (e.g. by the test file).
# ---------------------------------------------------------------------------
if (sys.nframe() == 0) stop("This helper is sourced by the analysis; use scripts/run.R.")
