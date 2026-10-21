# ============================================================================
# tcga_identifiers.R
#
# TCGA barcode standardization and sample-type classification. A TCGA
# barcode has the form TCGA-XX-XXXX-01A-...; the two-digit sample-type
# code (positions 14-15) distinguishes primary tumor (01) from normal,
# metastatic, and other sample types.
# ============================================================================

standardize_tcga_id <- function(x) {
  x <- trimws(as.character(x))
  gsub("\\.", "-", x)
}

tcga_sample_type_code <- function(x) {
  out <- rep(NA_character_, length(x))
  ok <- !is.na(x) & grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}", x)
  out[ok] <- sub("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-([0-9]{2}).*$", "\\1", x[ok])
  out
}

# %in% rather than == because tcga_sample_type_code() returns NA for an
# identifier that does not match the standard barcode pattern, and
# `== "01"` would propagate that NA. The result is used directly as a row
# or column selector (scripts/03_harmonize_identifiers.R,
# scripts/lib/expression_preprocessing.R); `[.data.frame` inserts an all-NA row for
# an NA row selector and errors for an NA column selector, so a malformed
# identifier must resolve to FALSE.
is_primary_tumor_barcode <- function(x) tcga_sample_type_code(x) %in% "01"

# The first 15 characters of a standardized TCGA barcode
# (Project-TSS-Participant-SampleType, e.g. "TCGA-BL-A0C8-01") identify the
# project, tissue source site, participant and two-digit sample type
# (01 = primary solid tumor). Truncation removes the vial letter at
# position 16 and every later field (portion, analyte, plate, center).
# Different vials of the same participant and sample type therefore share
# one 15-character identifier, which is why TCGAbiolinks::Tumor.purity can
# contain several rows per identifier; see the duplicate-barcode handling
# in scripts/05_define_eligible_features.R.
tcga_sample_barcode_15 <- function(x) substr(standardize_tcga_id(x), 1, 15)

# Parses a tumor-purity value (CPE or LUMP) from any of the representations
# a source may provide: numeric, character with a decimal point, character
# with a comma decimal separator (the format of the TCGAbiolinks::Tumor.purity
# factor labels), or a factor with either separator. Values always pass
# through as.character() and decimal-separator normalization before
# as.numeric(); calling as.numeric() directly on a factor would return the
# internal level code instead of the labelled value. Covered by
# tests/test_purity_value_parsing.R.
parse_purity_value <- function(x) as.numeric(gsub(",", ".", as.character(x), fixed = TRUE))

# Integrity check applied to every parsed purity vector before use (see
# scripts/05_define_eligible_features.R). Returns a zero-row result when no
# problem is found, or a one-row data.frame describing the problem, so that
# callers can run_stage_assert() on nrow(result) == 0 with a specific
# message. Three independent checks:
#   1. integer_level_codes: every non-missing value is a whole number and
#      the range extends far beyond [0, 1], the pattern produced when
#      factor level codes are read as values.
#   2. out_of_range: any non-missing value lies outside [0, 1].
#   3. rank_equivalent_to_factor_order: if the raw factor is supplied,
#      checks whether parsed values are numerically equal to the factor's
#      internal level codes for a large share of samples. Exact equality is
#      used rather than rank correlation because CPE/LUMP labels are
#      fixed-width "0,dddd" strings whose alphabetical order coincides with
#      numeric order, so a correctly parsed vector is also rank-concordant
#      with the level codes; a correctly parsed fraction (e.g. 0.9246) is
#      never equal to its level code (e.g. 3368).
check_purity_value_corruption <- function(parsed, raw_factor = NULL, rank_cor_threshold = 0.95) {
  ok <- !is.na(parsed) & !is.nan(parsed)
  problems <- list()

  looks_like_level_code <- any(ok) && all(parsed[ok] == round(parsed[ok])) && any(parsed[ok] > 1 | parsed[ok] < 0)
  if (looks_like_level_code) {
    problems[["integer_level_codes"]] <- sprintf(
      "All %d non-missing parsed values are whole numbers and at least one lies outside [0,1] (range [%.4g, %.4g]) -- this is exactly the shape produced by calling as.numeric() directly on a factor instead of parsing its character labels.",
      sum(ok), min(parsed[ok]), max(parsed[ok])
    )
  }

  out_of_range <- any(ok) && any(parsed[ok] < 0 | parsed[ok] > 1)
  if (out_of_range && !looks_like_level_code) {
    problems[["out_of_range"]] <- sprintf(
      "%d of %d non-missing parsed values fall outside the valid [0,1] purity-fraction range (range [%.4g, %.4g]).",
      sum(parsed[ok] < 0 | parsed[ok] > 1), sum(ok), min(parsed[ok]), max(parsed[ok])
    )
  }

  if (!is.null(raw_factor) && is.factor(raw_factor) && sum(ok) >= 10 && !looks_like_level_code) {
    level_code <- as.integer(raw_factor)[ok]
    frac_equal_to_level_code <- mean(abs(parsed[ok] - level_code) < 1e-9)
    if (frac_equal_to_level_code >= rank_cor_threshold) {
      problems[["rank_equivalent_to_factor_order"]] <- sprintf(
        "%.1f%% of non-missing parsed values are numerically IDENTICAL to the raw factor's own internal level-code integer (not merely rank-correlated) -- this is the direct signature of a factor-level-code leak, even in the (here, already excluded by the checks above) case where those integers happen to fall in [0,1].",
        100 * frac_equal_to_level_code
      )
    }
  }

  if (length(problems) == 0) return(data.frame(check = character(0), detail = character(0), stringsAsFactors = FALSE))
  data.frame(check = names(problems), detail = unlist(problems, use.names = FALSE), stringsAsFactors = FALSE)
}

# Resolves one purity metric (CPE or LUMP) for a sample_id_15 group that may
# contain several vial records, without a vial-preference tie-break. Used by
# scripts/05_define_eligible_features.R for each metric independently (a
# sample can have a resolvable CPE and an unresolvable LUMP, or vice versa).
#
# `vals`: every raw (possibly NA) value for this metric across the vial
#   records in the group (length 1 for a non-duplicated sample_id_15).
# `is_conflict`: TRUE if the group's non-missing values disagree (decided
#   by the caller via exact inequality among unique non-NA values).
#
# Returns NA_real_ if is_conflict is TRUE (no vial is preferred, and the
# sample is excluded only from analyses that require this metric as a
# covariate); otherwise the single available non-NA value, or NA_real_ if
# every record is NA.
resolve_purity_duplicate_value <- function(vals, is_conflict) {
  if (isTRUE(is_conflict)) return(NA_real_)
  nn <- unique(vals[!is.na(vals)])
  if (length(nn) == 0) return(NA_real_)
  nn[1]
}
