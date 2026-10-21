# ============================================================================
# targetscan_mapping.R
#
# Maps TCGA miRNA identifiers to human TargetScan prediction identifiers
# while preserving the original TCGA identifier used in the expression
# matrix.
#
# Mapping priority:
#   1. direct
#   2. canonical_accession_alias
#   3. seed_family_representative
#   4. unresolved
#
# Ambiguous mappings stop execution instead of being resolved arbitrarily.
# ============================================================================

normalize_targetscan_id <- function(x) {
  sub("[.][0-9]+$", "", trimws(x))
}


build_targetscan_mapping <- function(candidate_ids, prediction_ids, family_table) {

  candidate_ids <- sort(unique(trimws(candidate_ids)))
  prediction_ids <- sort(unique(trimws(prediction_ids)))

  required_columns <- c(
    "miR family",
    "Seed+m8",
    "Species ID",
    "MiRBase ID",
    "MiRBase Accession"
  )

  missing_columns <- setdiff(required_columns, names(family_table))

  if (length(missing_columns) > 0) {
    stop(
      "TargetScan family table is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  fam <- data.table::as.data.table(family_table)

  fam <- fam[
    `Species ID` == 9606,
    .(
      family = trimws(`miR family`),
      seed = trimws(`Seed+m8`),
      targetscan_id = trimws(`MiRBase ID`),
      accession = trimws(`MiRBase Accession`)
    )
  ]

  fam[, normalized_id := normalize_targetscan_id(targetscan_id)]
  fam <- unique(fam)

  mapping_rows <- vector("list", length(candidate_ids))

  for (i in seq_along(candidate_ids)) {

    candidate <- candidate_ids[i]

    # ----------------------------------------------------------------------
    # 1. Exact identifier match
    # ----------------------------------------------------------------------
    if (candidate %in% prediction_ids) {

      mapping_rows[[i]] <- data.table::data.table(
        mirna_id = candidate,
        targetscan_mirna_id = candidate,
        mapping_type = "direct",
        family = NA_character_,
        seed = NA_character_,
        accession = NA_character_
      )

      next
    }

    candidate_normalized <- normalize_targetscan_id(candidate)

    candidate_family_rows <- fam[
      targetscan_id == candidate |
      normalized_id == candidate_normalized
    ]

    # ----------------------------------------------------------------------
    # 2. Canonical accession-bearing alias
    #
    # TargetScan may split one miRNA identifier into variants such as .1
    # and .2. When exactly one prediction-bearing variant carries the
    # canonical miRBase accession, that variant is used.
    # ----------------------------------------------------------------------
    canonical_aliases <- unique(
      candidate_family_rows[
        accession != "" &
        !is.na(accession) &
        targetscan_id %in% prediction_ids,
        targetscan_id
      ]
    )

    if (length(canonical_aliases) == 1) {

      selected <- candidate_family_rows[
        targetscan_id == canonical_aliases,
      ][1]

      mapping_rows[[i]] <- data.table::data.table(
        mirna_id = candidate,
        targetscan_mirna_id = canonical_aliases,
        mapping_type = "canonical_accession_alias",
        family = selected$family,
        seed = selected$seed,
        accession = selected$accession
      )

      next
    }

    if (length(canonical_aliases) > 1) {
      stop(
        "Ambiguous canonical TargetScan alias for ", candidate, ": ",
        paste(sort(canonical_aliases), collapse = ", ")
      )
    }

    # ----------------------------------------------------------------------
    # 3. TargetScan seed-family representative
    #
    # Some mature miRNAs have no prediction row under their own name because
    # TargetScan stores predictions under one representative of the same
    # seed family.
    # ----------------------------------------------------------------------
    candidate_families <- unique(candidate_family_rows$family)
    candidate_families <- candidate_families[
      !is.na(candidate_families) & candidate_families != ""
    ]

    family_prediction_ids <- unique(
      fam[
        family %in% candidate_families &
        targetscan_id %in% prediction_ids,
        targetscan_id
      ]
    )

    if (length(family_prediction_ids) == 1) {

      selected <- fam[
        targetscan_id == family_prediction_ids &
        family %in% candidate_families,
      ][1]

      mapping_rows[[i]] <- data.table::data.table(
        mirna_id = candidate,
        targetscan_mirna_id = family_prediction_ids,
        mapping_type = "seed_family_representative",
        family = selected$family,
        seed = selected$seed,
        accession = selected$accession
      )

      next
    }

    if (length(family_prediction_ids) > 1) {
      stop(
        "Ambiguous TargetScan family mapping for ", candidate, ": ",
        paste(sort(family_prediction_ids), collapse = ", ")
      )
    }

    # ----------------------------------------------------------------------
    # 4. No valid mapping
    # ----------------------------------------------------------------------
    mapping_rows[[i]] <- data.table::data.table(
      mirna_id = candidate,
      targetscan_mirna_id = NA_character_,
      mapping_type = "unresolved",
      family = if (length(candidate_families) > 0) {
        paste(sort(candidate_families), collapse = ";")
      } else {
        NA_character_
      },
      seed = NA_character_,
      accession = NA_character_
    )
  }

  mapping <- data.table::rbindlist(mapping_rows, fill = TRUE)

  if (anyDuplicated(mapping$mirna_id)) {
    stop("Internal error: duplicate candidate miRNAs in TargetScan mapping.")
  }

  allowed_mapping_types <- c(
    "direct",
    "canonical_accession_alias",
    "seed_family_representative",
    "unresolved"
  )

  if (!all(mapping$mapping_type %in% allowed_mapping_types)) {
    stop("Internal error: unknown TargetScan mapping type.")
  }

  mapping[]
}
