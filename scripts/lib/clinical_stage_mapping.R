# ============================================================================
# clinical_stage_mapping.R
#
# Single authoritative mapping from a raw clinical tumor-stage string (from
# either ajcc_pathologic_tumor_stage or clinical_stage, whichever a given
# cancer type populates) to an ordinal Stage I-IV value. Accepts:
#   - the "Stage X" form and a bare Roman-numeral form (observed for THYM,
#     which records stage as "I", "IIa", "IIb", "III", "IVa", "IVb" with no
#     "Stage " prefix);
#   - FIGO-style trailing sub-sub-stage digits (observed for the
#     gynecological cancers, e.g. CESC's "Stage IB1"/"Stage IIA2" and UCS's
#     "Stage IIIC1"/"Stage IIIC2").
# Non-ordinal placeholder values ("[Discrepancy]", "[Unknown]", "IS" for
# in situ, "I/II NOS", "Stage X") return NA rather than being forced onto
# the ordinal scale.
# ============================================================================

stage_to_numeric <- function(x) {
  x <- toupper(trimws(x))
  x <- sub("^STAGE\\s*", "", x)
  out <- rep(NA_real_, length(x))
  out[grepl("^IV[ABC]?[0-9]*$", x)] <- 4
  out[grepl("^III[ABC]?[0-9]*$", x)] <- 3
  out[grepl("^II[ABC]?[0-9]*$", x)] <- 2
  out[grepl("^I[ABC]?[0-9]*$", x)] <- 1
  out
}
