# ============================================================================
# expression_preprocessing.R
#
# Shared expression-matrix preprocessing: TCGA barcode standardization,
# restriction to primary tumor samples, validation of the declared input
# scale, duplicate-column and duplicate-row classification and resolution,
# and zero-variance feature removal. The miRNA and mRNA matrices receive
# the same structural treatment, while each source's scale is declared
# independently (config/analysis_config.R, $expression).
#
# Winsorization is not part of preprocess_expression_matrix(). It is
# unnecessary before rank-based Spearman correlation, where per-feature
# clipping can introduce ties, but is appropriate immediately before
# value-sensitive methods (OLS residualization, regression, Cox models,
# clustering). It is therefore exposed as a separate step
# (winsorize_matrix_rows(), winsorize_matrix_rows_by_group()) that each
# consuming stage applies explicitly.
# ============================================================================

# --- Scale declaration and validation ---------------------------------------

# Empirical helper: TRUE if every finite value in `mat` is <= threshold. Not
# used on its own to decide whether to transform; validate_expression_scale()
# treats the declared scale (config/analysis_config.R, $expression$mrna_scale
# and $expression$mirna_scale) as the decision and uses this function only
# to detect a contradiction.
is_log_scale_matrix <- function(mat, threshold = 50) {
  vals <- as.numeric(mat)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(NA)
  max(vals, na.rm = TRUE) <= threshold
}

# Validates a declared source scale ("log2" or "linear") against the
# matrix's empirical range as a consistency check. Returns a list with
# $consistent (TRUE/FALSE/NA), $empirical_max and $detail. The calling
# stage must run_stage_assert() on $consistent before applying any
# transform, so a mismatch (a misdeclared source or a corrupted input file)
# stops the stage.
validate_expression_scale <- function(mat, declared_scale, source_name, sanity_max = 50) {
  stopifnot(declared_scale %in% c("log2", "linear"))
  vals <- as.numeric(mat)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) {
    return(list(consistent = NA, empirical_max = NA_real_, empirical_min = NA_real_,
                detail = sprintf("%s: matrix has no finite values to validate against the declared scale.", source_name)))
  }
  mx <- max(vals); mn <- min(vals)
  if (declared_scale == "log2") {
    consistent <- mx <= sanity_max
    detail <- if (consistent)
      sprintf("%s: declared scale 'log2' is consistent with the observed range [%.4g, %.4g] (max <= sanity bound %.0f).", source_name, mn, mx, sanity_max)
    else
      sprintf("%s: declared scale is 'log2' but the observed maximum (%.4g) exceeds the sanity bound (%.0f) for already-log2-transformed expression -- this contradicts the declared source-dataset identity in config/analysis_config.R$expression and must be resolved by a human (wrong file supplied? scale declaration wrong for this source?) before proceeding.", source_name, mx, sanity_max)
  } else {
    consistent <- !(mx <= 1 && mn >= 0)
    detail <- if (consistent)
      sprintf("%s: declared scale 'linear' is consistent with the observed range [%.4g, %.4g] (not confined to [0,1], which would suggest an already-transformed or proportion-scale input).", source_name, mn, mx)
    else
      sprintf("%s: declared scale is 'linear' but the observed range ([%.4g, %.4g]) is entirely within [0,1] -- this looks like an already-transformed or proportion-scale input, contradicting the declared source-dataset identity, and must be resolved by a human before proceeding.", source_name, mn, mx)
  }
  list(consistent = consistent, empirical_max = mx, empirical_min = mn, detail = detail)
}

# --- Duplicate-column classification and resolution -------------------------

# Classifies exact-string duplicate standardized column identifiers in a
# feature-by-sample matrix (two or more columns whose standardize_tcga_id()
# output is identical) as:
#   - "technical_duplicate": all columns in the group are numerically
#     identical (NA matches NA), i.e. the same measurement represented
#     twice. Resolved deterministically by retaining the leftmost column.
#   - "conflicting_record": the columns differ, i.e. one identifier maps to
#     different data. Never resolved automatically; the caller stops and
#     requires a documented decision (config/sample_exclusions.tsv).
# Returns a data.frame with one row per duplicate group (zero rows if none).
classify_duplicate_columns <- function(mat) {
  cn <- colnames(mat)
  dup_ids <- unique(cn[duplicated(cn)])
  if (length(dup_ids) == 0) {
    return(data.frame(id = character(0), classification = character(0), n_columns = integer(0), detail = character(0), stringsAsFactors = FALSE))
  }
  rows <- lapply(dup_ids, function(id) {
    idx <- which(cn == id)
    cols <- mat[, idx, drop = FALSE]
    ref <- cols[, 1]
    identical_to_ref <- vapply(seq_len(ncol(cols))[-1], function(j) {
      b <- cols[, j]
      both_na <- is.na(ref) & is.na(b)
      eq <- both_na | (!is.na(ref) & !is.na(b) & abs(ref - b) < 1e-9)
      all(eq)
    }, logical(1))
    is_technical <- length(identical_to_ref) == 0 || all(identical_to_ref)
    data.frame(id = id,
               classification = if (is_technical) "technical_duplicate" else "conflicting_record",
               n_columns = length(idx),
               detail = if (is_technical) sprintf("%d columns share this identifier; values are numerically identical (safe to resolve to one column).", length(idx))
                        else sprintf("%d columns share this identifier; values DIFFER between columns -- requires a user decision, not auto-resolved.", length(idx)),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Resolution for technical duplicates only (conflicting groups must already
# have stopped the caller). Retains the leftmost column of each group;
# because the retained and dropped columns are numerically identical, the
# choice does not affect any downstream value.
resolve_technical_duplicate_columns <- function(mat, dup_report) {
  keep <- !logical(ncol(mat))
  cn <- colnames(mat)
  for (id in dup_report$id[dup_report$classification == "technical_duplicate"]) {
    idx <- which(cn == id)
    keep[idx[-1]] <- FALSE
  }
  mat[, keep, drop = FALSE]
}

# Detects aliquot duplicates: two or more distinct full standardized column
# identifiers sharing the same 15-character Project-TSS-Participant-
# SampleType prefix (tcga_sample_barcode_15()), i.e. several physical
# aliquots, vials or portions of one participant and sample type. TCGA
# legitimately resequences some samples, but keeping both columns as
# independent samples would double-count one tumor in every downstream
# correlation, regression and survival model. Returns a data.frame with one
# row per affected 15-character group (zero rows if none).
detect_aliquot_duplicates <- function(colnames_full) {
  id15 <- tcga_sample_barcode_15(colnames_full)
  dup15 <- unique(id15[duplicated(id15)])
  if (length(dup15) == 0) {
    return(data.frame(sample_id_15 = character(0), full_ids = character(0), n_aliquots = integer(0), stringsAsFactors = FALSE))
  }
  rows <- lapply(dup15, function(s15) {
    full <- colnames_full[id15 == s15]
    data.frame(sample_id_15 = s15, full_ids = paste(full, collapse = ";"), n_aliquots = length(full), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Deterministic resolution for aliquot duplicates: retain the aliquot column
# with the fewest missing or non-finite values, breaking ties by the
# alphabetically first full identifier. Returns the full identifiers to
# drop. This rule differs from the purity duplicate handling in
# scripts/05_define_eligible_features.R because an expression column carries
# its own completeness information whereas a purity record does not.
resolve_aliquot_duplicates <- function(mat, aliquot_report) {
  drop_ids <- character(0)
  for (i in seq_len(nrow(aliquot_report))) {
    full_ids <- strsplit(aliquot_report$full_ids[i], ";", fixed = TRUE)[[1]]
    n_missing <- vapply(full_ids, function(id) sum(!is.finite(mat[, id])), integer(1))
    best <- order(n_missing, full_ids)[1]
    drop_ids <- c(drop_ids, setdiff(full_ids, full_ids[best]))
  }
  drop_ids
}

# Same technical_duplicate / conflicting_record classification as
# classify_duplicate_columns(), applied to rows of a data.frame keyed by an
# identifier column (the clinical/survival table, one row per sample).
classify_duplicate_rows <- function(df, id_col) {
  ids <- df[[id_col]]
  dup_ids <- unique(ids[duplicated(ids)])
  if (length(dup_ids) == 0) {
    return(data.frame(id = character(0), classification = character(0), n_rows = integer(0), detail = character(0), stringsAsFactors = FALSE))
  }
  value_cols <- setdiff(colnames(df), id_col)
  rows <- lapply(dup_ids, function(id) {
    idx <- which(ids == id)
    sub <- df[idx, value_cols, drop = FALSE]
    ref <- sub[1, , drop = FALSE]
    identical_to_ref <- vapply(seq_len(nrow(sub))[-1], function(j) {
      isTRUE(all.equal(ref, sub[j, , drop = FALSE], check.attributes = FALSE))
    }, logical(1))
    is_technical <- length(identical_to_ref) == 0 || all(identical_to_ref)
    data.frame(id = id,
               classification = if (is_technical) "technical_duplicate" else "conflicting_record",
               n_rows = length(idx),
               detail = if (is_technical) sprintf("%d rows share this identifier; all other fields are identical (safe to resolve to one row).", length(idx))
                        else sprintf("%d rows share this identifier; at least one other field DIFFERS between rows -- requires a user decision, not auto-resolved.", length(idx)),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

resolve_technical_duplicate_rows <- function(df, id_col, dup_report) {
  keep <- !logical(nrow(df))
  ids <- df[[id_col]]
  for (id in dup_report$id[dup_report$classification == "technical_duplicate"]) {
    idx <- which(ids == id)
    keep[idx[-1]] <- FALSE
  }
  df[keep, , drop = FALSE]
}

# Non-numeric and infinite values are detected inside
# preprocess_expression_matrix() on the unclass()+list representation of
# the raw columns, because a loop over setdiff(colnames(df), id_col) would
# collapse duplicate raw column names before they are examined.

# --- Winsorization (explicit, opt-in; see the file header) -------------------

winsorize_vector <- function(x, probs = c(0.001, 0.999)) {
  qs <- suppressWarnings(quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
  x[x < qs[1]] <- qs[1]
  x[x > qs[2]] <- qs[2]
  x
}

# Per-feature winsorization (MARGIN = 1: each row, one miRNA or gene, is
# clipped across its values over the sample cohort). Intended for
# value-sensitive methods to apply to their own input immediately before
# the computation: OLS residualization (Stage 10), MCPcounter estimation and
# immune-adjustment regressions (Stage 12), Cox regression (Stages 15-16),
# mediation path models (Stage 17) and PCA/k-means clustering (Stage 18).
# It is not applied by default because the rank-based Spearman correlations
# of Stages 06 and 08 do not require it and clipping to a shared bound can
# introduce ties.
winsorize_matrix_rows <- function(mat, probs = c(0.001, 0.999)) {
  out <- t(apply(mat, 1, winsorize_vector, probs = probs))
  dimnames(out) <- dimnames(mat)
  out
}

# Per-cancer-group winsorization: applies winsorize_matrix_rows() to each
# group's column subset independently, so that no cancer type's clipping
# threshold is informed by another cancer's distribution. Used wherever
# winsorization precedes a per-cancer analysis (e.g. the CPE/LUMP-adjusted
# OLS residualization of Stage 10). `group` must be aligned to the columns
# of `mat` (length == ncol(mat)); columns with NA group are returned
# unchanged. Clipping probabilities and NA handling match
# winsorize_matrix_rows() and winsorize_vector().
winsorize_matrix_rows_by_group <- function(mat, group, probs = c(0.001, 0.999)) {
  stopifnot(length(group) == ncol(mat))
  out <- mat
  for (g in unique(group[!is.na(group)])) {
    cols <- which(group == g)
    if (length(cols) == 0) next
    out[, cols] <- winsorize_matrix_rows(mat[, cols, drop = FALSE], probs = probs)
  }
  out
}

# --- Main preprocessing entry point ------------------------------------------

# Structural preprocessing shared by both expression layers: barcode
# standardization, restriction to primary-tumor columns, exact-duplicate
# column classification and resolution, aliquot-duplicate detection and
# resolution, numeric coercion with an explicit count of unparseable values,
# declared-scale validation, and zero-variance feature removal. Does not
# winsorize. Returns a list exposing every diagnostic; $ok is FALSE whenever
# an unresolved conflict exists, and the caller must then stop
# (run_stage_assert() on $conflicting_columns and $scale_check).
preprocess_expression_matrix <- function(df, id_col, declared_scale, pseudo_count,
                                          scale_sanity_max = 50, remove_zero_variance = TRUE) {
  # All column selection below operates on a matrix built with
  # unclass()+cbind(), not on `df[, <selector>, drop = FALSE]`.
  # `[.data.frame` renames repeated column names in the selected result
  # (a second "X" becomes "X.1") for character, integer and logical
  # selectors alike, which would hide the duplicate standardized barcodes
  # that classify_duplicate_columns() must detect. Matrix subsetting
  # preserves duplicate dimnames, so the data stay in matrix form until the
  # deduplicated return value is constructed.
  lst <- unclass(df)
  feature_ids <- as.character(lst[[1]])
  raw_colnames <- names(lst)[-1]
  sample_ids <- standardize_tcga_id(raw_colnames)

  missing_tokens <- c("", "NA", "N/A", "NaN", "null", "NULL", "nan")
  n_newly_na <- 0L
  mat_cols <- vector("list", length(raw_colnames))
  for (j in seq_along(raw_colnames)) {
    orig <- lst[[j + 1]]
    orig_chr <- trimws(as.character(orig))
    was_missing <- is.na(orig) | orig_chr %in% missing_tokens
    coerced <- suppressWarnings(as.numeric(orig))
    n_newly_na <- n_newly_na + sum(is.na(coerced) & !was_missing)
    mat_cols[[j]] <- coerced
  }
  mat_all <- do.call(cbind, mat_cols)
  colnames(mat_all) <- sample_ids

  primary_keep <- is_primary_tumor_barcode(sample_ids)
  mat_all <- mat_all[, primary_keep, drop = FALSE]

  dup_report <- classify_duplicate_columns(mat_all)
  conflicting <- dup_report[dup_report$classification == "conflicting_record", , drop = FALSE]
  if (nrow(conflicting) > 0) {
    return(list(ok = FALSE, stage = "duplicate_column_conflict", duplicate_column_report = dup_report,
                message = sprintf("%d standardized sample identifier(s) map to conflicting (non-identical) values across duplicate columns: %s. This requires a human decision (which record is correct, or whether both should be excluded) and is not auto-resolved.",
                                   nrow(conflicting), paste(conflicting$id, collapse = ", "))))
  }
  if (nrow(dup_report) > 0) {
    mat_all <- resolve_technical_duplicate_columns(mat_all, dup_report)
  }

  aliquot_report <- detect_aliquot_duplicates(colnames(mat_all))
  aliquot_drop_ids <- character(0)
  if (nrow(aliquot_report) > 0) {
    aliquot_drop_ids <- resolve_aliquot_duplicates(mat_all, aliquot_report)
    mat_all <- mat_all[, !(colnames(mat_all) %in% aliquot_drop_ids), drop = FALSE]
  }

  scale_check <- validate_expression_scale(mat_all, declared_scale, id_col, sanity_max = scale_sanity_max)

  transformed_by_pipeline <- FALSE
  mat <- mat_all
  if (isTRUE(scale_check$consistent) && declared_scale == "linear") {
    mat <- log2(mat + pseudo_count)
    transformed_by_pipeline <- TRUE
  }

  n_inf <- sum(is.infinite(mat))
  n_nan <- sum(is.nan(mat))

  rownames(mat) <- feature_ids

  if (remove_zero_variance) {
    row_sd <- apply(mat, 1, function(z) sd(z, na.rm = TRUE))
    zero_var_removed <- sum(!(is.finite(row_sd) & row_sd > 0))
    mat <- mat[is.finite(row_sd) & row_sd > 0, , drop = FALSE]
  } else {
    zero_var_removed <- NA_integer_
  }

  out <- data.frame(feature_id = rownames(mat), mat, check.names = FALSE)
  colnames(out)[1] <- id_col

  list(ok = TRUE, data = out, declared_scale = declared_scale, scale_check = scale_check,
       transformed_by_pipeline = transformed_by_pipeline,
       duplicate_column_report = dup_report, aliquot_report = aliquot_report,
       aliquot_dropped = aliquot_drop_ids, n_newly_na_from_coercion = n_newly_na,
       n_infinite = n_inf, n_nan = n_nan, n_zero_variance_removed = zero_var_removed)
}

highest_variance_row <- function(mat, idx) {
  if (length(idx) == 1) return(idx)
  row_var <- apply(mat[idx, , drop = FALSE], 1, var, na.rm = TRUE)
  best <- which.max(row_var)
  # When every candidate row's variance is NA (all-NA rows, or too few
  # non-NA values), which.max() returns integer(0); in that case the first
  # candidate row is retained so that the duplicate group is never dropped
  # entirely.
  if (length(best) == 0) return(idx[1])
  idx[best]
}

extract_gene_row <- function(mat, gene_symbol) {
  clean_id <- function(x) toupper(sub("\\..*$", "", sub("\\|.*$", "", trimws(as.character(x)))))
  rn_clean <- clean_id(rownames(mat))
  idx <- which(rn_clean == toupper(gene_symbol) | toupper(rownames(mat)) == toupper(gene_symbol))
  if (length(idx) == 0) stop(sprintf("Gene %s not found in expression matrix.", gene_symbol))
  idx <- highest_variance_row(mat, idx)
  as.numeric(mat[idx, ])
}

# Row-level deduplication for a feature-by-sample matrix whose rownames may
# contain the same symbol more than once (e.g. SLC35E2 in the matched mRNA
# matrix). Row-indexing a matrix by a duplicated dimname returns only the
# first match, which would make lookups depend on row order; this function
# applies one shared rule (retain the duplicate with the highest variance)
# wherever a gene-keyed matrix is built.
deduplicate_matrix_rows <- function(mat) {
  rn <- rownames(mat)
  if (!any(duplicated(rn))) return(mat)
  keep <- !logical(nrow(mat))
  for (s in unique(rn[duplicated(rn)])) {
    idx <- which(rn == s)
    keep[setdiff(idx, highest_variance_row(mat, idx))] <- FALSE
  }
  mat[keep, , drop = FALSE]
}
