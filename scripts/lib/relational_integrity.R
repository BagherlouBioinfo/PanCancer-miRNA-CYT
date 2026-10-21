# ============================================================================
# relational_integrity.R
#
# Explicitly keyed join and cardinality-assertion helpers. Every table built
# from a full pair-universe table must be restricted to the run's primary
# triplet keys before it is stacked into a triplet-scoped table; a missing
# restriction expands row counts without ever duplicating a key, so
# key-uniqueness checks alone cannot detect it. Each function returns
# $ok = FALSE with an explicit $problem message instead of tolerating an
# unexpected key expansion; merge(..., allow.cartesian = TRUE) is never used
# to silence a mismatch.
# ============================================================================

# Restricts `dt` (a data.table containing every column in `by`) to the rows
# whose `by` key combination appears in `keys` (a data.table with one row
# per allowed key combination, e.g. the primary triplets' (mirna_id, gene)
# pairs). This semi-join must be applied to every full pair-universe table
# before it is stacked into a triplet-scoped table.
#
# Before restricting, asserts that `keys` has unique `by` combinations, so a
# duplicated allowed-key list is never treated as several distinct keys.
# Row-level uniqueness of `dt` on `by` is not required (e.g. one row per
# cancer is expected).
#
# Returns list($data, $ok, $n_before, $n_after, $n_keys_expected,
# $n_keys_matched, $problem). The operation only narrows `dt`: a key in
# `keys` can match many existing rows of `dt` but never creates new rows.
restrict_to_keys <- function(dt, keys, by) {
  key_dup <- keys[, .N, by = by][N > 1]
  if (nrow(key_dup) > 0) {
    return(list(ok = FALSE, data = NULL, n_before = nrow(dt), n_after = NA_integer_,
                n_keys_expected = NA_integer_, n_keys_matched = NA_integer_,
                problem = sprintf(
                  "restrict_to_keys(): the allowed-key table itself has %d duplicated %s combination(s); an allowed-key list must be unique (e.g. use unique(triplets[, .(mirna_id, gene)]), never the raw triplets table, which may repeat a pair across cancers).",
                  nrow(key_dup), paste(by, collapse = ","))))
  }
  keys_unique <- unique(keys[, ..by])
  n_before <- nrow(dt)
  out <- merge(dt, keys_unique, by = by)
  n_after <- nrow(out)
  n_keys_matched <- uniqueN(out[, ..by])
  n_keys_expected <- nrow(keys_unique)
  list(ok = TRUE, data = out, n_before = n_before, n_after = n_after,
       n_keys_expected = n_keys_expected, n_keys_matched = n_keys_matched, problem = NA_character_)
}

# Post-assembly cardinality assertion: `dt` must have a unique combination
# of `key_cols` (the intended grain), and its row count must not exceed
# `max_expected`, a theoretical maximum derived by the caller from the
# run's own counts (e.g. n_triplets * n_arms * n_adjustments * n_cancers).
# Returns list($ok, $n_rows, $n_unique_keys, $n_duplicate_keys,
# $duplicate_examples, $problem).
assert_cardinality <- function(dt, key_cols, max_expected = Inf, label = "table") {
  n_rows <- nrow(dt)
  key_counts <- dt[, .N, by = key_cols]
  n_unique <- nrow(key_counts)
  dup_keys <- key_counts[N > 1]
  n_dup <- nrow(dup_keys)
  over_max <- n_rows > max_expected
  ok <- n_dup == 0 && !over_max
  problem <- if (n_dup > 0) {
    sprintf("%s: %d duplicated key combination(s) found at the intended grain (%s) -- e.g. %s",
            label, n_dup, paste(key_cols, collapse = ","),
            paste(utils::head(do.call(paste, c(as.list(dup_keys[, ..key_cols]), sep = "|")), 3), collapse = "; "))
  } else if (over_max) {
    sprintf("%s: %d rows exceeds the dynamically derived theoretical maximum of %d for this run's actual inputs",
            label, n_rows, max_expected)
  } else NA_character_
  list(ok = ok, n_rows = n_rows, n_unique_keys = n_unique, n_duplicate_keys = n_dup,
       duplicate_examples = dup_keys, problem = problem)
}

# Confirms that an intended one-to-many broadcast join (e.g. one miRNA-CYT
# estimate replicated across every triplet sharing that miRNA) multiplied
# each `by` group by exactly its expected fan-out, so that no unintended
# many-to-many join is present. `dt_before` is the pre-broadcast table,
# `dt_after` the broadcast result, and `fanout_keys` the table joined
# against (e.g. unique(triplets[, .(mirna_id, gene)])), used to compute the
# expected number of output rows per `by` group independently of the join.
# Returns list($ok, $problem) plus per-group expected and observed counts.
assert_broadcast_cardinality <- function(dt_before, dt_after, by, fanout_keys, fanout_by, label = "broadcast") {
  expected_fanout <- fanout_keys[, .N, by = fanout_by]
  before_groups <- dt_before[, .N, by = by]
  # Left join (all.x = TRUE): a `by` group in `dt_before` that is absent
  # from `fanout_keys` has an expected fan-out of 0 and is still checked.
  chk <- merge(before_groups, expected_fanout, by = intersect(by, fanout_by), all.x = TRUE)
  chk[is.na(N.y), N.y := 0L]
  chk[, expected_rows := N.x * N.y]
  observed <- dt_after[, .N, by = by]
  chk <- merge(chk, observed, by = by, all.x = TRUE)
  chk[is.na(N), N := 0L]
  bad <- chk[N != expected_rows]
  ok <- nrow(bad) == 0
  list(ok = ok, n_groups_checked = nrow(chk), n_groups_mismatched = nrow(bad),
       problem = if (!ok) sprintf("%s: %d group(s) had an observed row count different from the expected fan-out (e.g. %s)",
                                   label, nrow(bad), paste(utils::head(do.call(paste, c(as.list(bad[, ..by]), sep = "|")), 3), collapse = "; ")) else NA_character_)
}
