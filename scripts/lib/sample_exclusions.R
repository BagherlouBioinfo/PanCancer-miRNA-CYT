# ============================================================================
# sample_exclusions.R
#
# Version-controlled, prespecified sample-level exclusions. Loads and
# validates config/sample_exclusions.tsv, a reviewed list of samples to
# exclude from some or all analytical layers, each with a documented reason
# and evidence reference. No sample identity is hard-coded in this file;
# every exclusion is driven by the configuration table passed in at call
# time.
#
# Matching uses the 15-character TCGA barcode
# (Project-TSS-Participant-SampleType; R/tcga_identifiers.R::
# tcga_sample_barcode_15()), the grain used for cross-file sample identity
# throughout the pipeline.
#
# Column removal from feature-by-sample data.frames uses unclass() and list
# reconstruction rather than `df[, <selector>, drop = FALSE]`, because
# `[.data.frame` applies make.unique()-style renaming to repeated column
# names in the selected result, which would interfere with the
# duplicate-barcode conflict detection in scripts/lib/expression_preprocessing.R.
# ============================================================================

# Loads and validates config/sample_exclusions.tsv. Stops on:
#   - a missing file
#   - a missing required column
#   - a malformed sample_id (not a valid TCGA barcode once standardized to
#     its 15-character form)
#   - an unknown `scope` value
#   - a duplicate entry (same 15-character sample_id listed more than once)
#   - an empty `reason` or `evidence_report`
# Returns a data.frame with columns sample_id_raw, sample_id_15, scope,
# reason and evidence_report, one row per validated entry.
load_sample_exclusions <- function(path, valid_scopes = c("all_layers", "expression_only", "clinical_only", "purity_only")) {
  if (!file.exists(path)) {
    stop(sprintf("Sample exclusion config not found: %s", path), call. = FALSE)
  }
  raw <- read.delim(path, stringsAsFactors = FALSE, colClasses = "character")
  required_cols <- c("sample_id", "scope", "reason", "evidence_report")
  missing_cols <- setdiff(required_cols, colnames(raw))
  if (length(missing_cols) > 0) {
    stop(sprintf("%s is missing required column(s): %s", path, paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  if (nrow(raw) == 0) {
    return(data.frame(sample_id_raw = character(0), sample_id_15 = character(0), scope = character(0),
                       reason = character(0), evidence_report = character(0), stringsAsFactors = FALSE))
  }

  raw$sample_id <- trimws(raw$sample_id)
  raw$scope <- trimws(raw$scope)

  if (any(raw$sample_id == "" | is.na(raw$sample_id))) {
    stop(sprintf("%s contains a row with an empty sample_id.", path), call. = FALSE)
  }

  standardized <- tcga_sample_barcode_15(standardize_tcga_id(raw$sample_id))
  is_valid_pattern <- grepl("^TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}-[0-9]{2}$", standardized)
  if (any(!is_valid_pattern)) {
    bad <- unique(raw$sample_id[!is_valid_pattern])
    stop(sprintf(
      "%s contains malformed sample_id value(s) that do not standardize to a valid 15-character TCGA barcode (Project-TSS-Participant-SampleType, e.g. 'TCGA-23-1023-01'): %s",
      path, paste(bad, collapse = ", ")), call. = FALSE)
  }

  if (any(!raw$scope %in% valid_scopes)) {
    bad_scope <- unique(raw$scope[!raw$scope %in% valid_scopes])
    stop(sprintf("%s contains unknown scope value(s): %s (valid: %s)",
                 path, paste(bad_scope, collapse = ", "), paste(valid_scopes, collapse = ", ")), call. = FALSE)
  }

  if (any(duplicated(standardized))) {
    dup_ids <- unique(standardized[duplicated(standardized)])
    stop(sprintf(
      "%s contains duplicate exclusion entries for: %s -- each sample may be listed exactly once, with one authoritative reason, not multiple conflicting or redundant rows.",
      path, paste(dup_ids, collapse = ", ")), call. = FALSE)
  }

  if (any(raw$reason == "" | is.na(raw$reason))) {
    stop(sprintf("%s has a row with an empty reason -- every exclusion must be justified.", path), call. = FALSE)
  }
  if (any(raw$evidence_report == "" | is.na(raw$evidence_report))) {
    stop(sprintf("%s has a row with an empty evidence_report -- every exclusion must cite its supporting evidence.", path), call. = FALSE)
  }

  data.frame(sample_id_raw = raw$sample_id, sample_id_15 = standardized, scope = raw$scope,
             reason = raw$reason, evidence_report = raw$evidence_report, stringsAsFactors = FALSE)
}

# Matches sample identifiers (full, possibly unstandardized barcodes, e.g.
# raw matrix column names or a clinical/purity identifier column) against a
# validated exclusion table from load_sample_exclusions(), by 15-character
# barcode. `applicable_scopes` restricts which scope values apply to the
# calling layer (e.g. a clinical table honours "all_layers" and
# "clinical_only" but not "expression_only").
#
# Returns a list:
#   $keep  -- logical vector in the order of `ids`; TRUE = retain.
#   $audit -- data.frame with one row per excluded identifier (`layer` is
#             the caller-supplied label, e.g. "miRNA", "mRNA", "clinical",
#             "purity"): original_id, sample_id_15, reason, evidence_report.
apply_sample_exclusions <- function(ids, exclusions, layer, applicable_scopes = c("all_layers")) {
  applicable <- exclusions[exclusions$scope %in% applicable_scopes, , drop = FALSE]
  ids_15 <- tcga_sample_barcode_15(standardize_tcga_id(ids))
  match_idx <- match(ids_15, applicable$sample_id_15)
  excluded <- !is.na(match_idx)
  audit <- if (any(excluded)) {
    data.frame(layer = layer, original_id = ids[excluded], sample_id_15 = ids_15[excluded],
               reason = applicable$reason[match_idx[excluded]],
               evidence_report = applicable$evidence_report[match_idx[excluded]],
               stringsAsFactors = FALSE)
  } else {
    data.frame(layer = character(0), original_id = character(0), sample_id_15 = character(0),
               reason = character(0), evidence_report = character(0), stringsAsFactors = FALSE)
  }
  list(keep = !excluded, audit = audit)
}

# Drops excluded sample columns from a feature-by-sample data.frame whose
# first column is the feature identifier (the shape of the raw miRNA and
# mRNA matrices). Uses unclass() and list reconstruction so that unrelated
# duplicate column names are left intact for the duplicate-conflict
# detection in scripts/lib/expression_preprocessing.R. Returns a list: $df (reduced
# data.frame), $audit, $n_removed, $removed_names.
drop_excluded_sample_columns <- function(df, exclusions, layer, applicable_scopes = c("all_layers", "expression_only")) {
  lst <- unclass(df)
  id_col_name <- names(lst)[1]
  raw_names <- names(lst)[-1]
  res <- apply_sample_exclusions(raw_names, exclusions, layer, applicable_scopes)
  kept_lst <- c(lst[1], lst[-1][res$keep])
  df_out <- as.data.frame(kept_lst, check.names = FALSE, stringsAsFactors = FALSE)
  names(df_out)[1] <- id_col_name
  list(df = df_out, audit = res$audit, n_removed = sum(!res$keep), removed_names = raw_names[!res$keep])
}

# Drops excluded sample rows from a sample-by-field data.frame (the shape of
# the clinical/survival and purity tables), matching on `id_col`. Logical
# row indexing is not affected by the duplicate column-name behaviour noted
# above.
drop_excluded_sample_rows <- function(df, id_col, exclusions, layer, applicable_scopes = c("all_layers")) {
  res <- apply_sample_exclusions(df[[id_col]], exclusions, layer, applicable_scopes)
  df_out <- df[res$keep, , drop = FALSE]
  list(df = df_out, audit = res$audit, n_removed = sum(!res$keep), removed_ids = df[[id_col]][!res$keep])
}
