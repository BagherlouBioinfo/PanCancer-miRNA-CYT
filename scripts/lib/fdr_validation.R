# ============================================================================
# fdr_validation.R
#
# Generic Benjamini-Hochberg hypothesis-family recomputation, used by
# tests/test_fdr_hypothesis_families.R to validate every family documented
# in reports/FDR_FAMILY_SPECIFICATION.tsv with a single implementation.
# ============================================================================

# Validates that a stored FDR column equals the Benjamini-Hochberg
# correction of a stored raw p-value column, computed once per unique
# hypothesis key, separately within each level of `group_vars` (the family
# boundary).
#
# Arguments:
#   dt: a data.table (or data.frame) with one row per hypothesis instance.
#       A family may contain several rows per hypothesis key only as a
#       documented broadcast (e.g. a gene-CYT result attached to every
#       miRNA that targets that gene); all rows sharing a key must then
#       carry identical p and FDR values, otherwise the key is reported as
#       inconsistent.
#   group_vars: column names defining the FDR family boundary (e.g.
#       c("cancer") or c("arm", "adjustment")). BH is computed
#       independently within each unique combination.
#   key_vars: column names defining one unique hypothesis (e.g.
#       c("mirna_id", "cancer"), or c("gene", "cancer") for a family
#       deduplicated by gene). Must include every column in group_vars.
#   p_col, fdr_col: names of the raw-p and stored-FDR columns.
#   tolerance: maximum allowed |stored - recomputed| difference (numerical
#       tolerance only).
#
# Returns a list: $consistent (TRUE/FALSE), $max_diff, $n_duplicate_keys
# (keys whose repeated rows carry inconsistent p or FDR values),
# $n_hypotheses (unique keys checked) and $detail (human-readable).
validate_fdr_family <- function(dt, group_vars, key_vars, p_col = "p", fdr_col = "fdr", tolerance = 1e-8) {
  dt <- as.data.table(dt)
  stopifnot(all(group_vars %in% key_vars))
  cols <- unique(c(key_vars, group_vars, p_col, fdr_col))
  stopifnot(all(cols %in% names(dt)))

  key_dt <- unique(dt[, cols, with = FALSE])
  n_dup_keys <- sum(duplicated(key_dt[, key_vars, with = FALSE]) | duplicated(key_dt[, key_vars, with = FALSE], fromLast = TRUE))
  if (n_dup_keys > 0) {
    return(list(consistent = FALSE, max_diff = NA_real_, n_duplicate_keys = n_dup_keys, n_hypotheses = NA_integer_,
                detail = sprintf("%d row(s) share a hypothesis key (%s) but carry DIFFERING %s/%s values -- a hypothesis key must map to exactly one raw p-value and one stored FDR value.",
                                  n_dup_keys, paste(key_vars, collapse = ","), p_col, fdr_col)))
  }

  key_dt[, .recomputed_fdr := p.adjust(get(p_col), method = "BH"), by = group_vars]
  diffs <- abs(key_dt[[fdr_col]] - key_dt$.recomputed_fdr)
  max_diff <- suppressWarnings(max(diffs, na.rm = TRUE))
  consistent <- is.finite(max_diff) && max_diff <= tolerance
  list(consistent = consistent, max_diff = max_diff, n_duplicate_keys = 0L, n_hypotheses = nrow(key_dt),
       detail = sprintf("max |stored - recomputed| = %.2e across %d unique hypotheses, grouped by (%s)",
                         if (is.finite(max_diff)) max_diff else NA_real_, nrow(key_dt), paste(group_vars, collapse = ",")))
}
