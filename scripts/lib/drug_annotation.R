# ============================================================================
# drug_annotation.R
#
# Helpers for Stage 22 (drug-gene annotation through the DGIdb public
# GraphQL API, https://dgidb.org/api/graphql). Response parsing and
# deduplication are pure functions of an already-fetched, parsed JSON list,
# so they can be unit-tested against synthetic responses without network
# access.
#
# The DGIdb `genes(names: [...])` query omits unrecognized names from the
# returned `nodes` array; there is no explicit "not found" entry. A queried
# gene absent from `nodes` is therefore classified as unresolved, which is
# distinct from a recognized gene with an empty `interactions` array (zero
# interactions, a reportable negative result).
# ============================================================================

# Builds the GraphQL query body (an R list ready for jsonlite::toJSON) for a
# batch of gene symbols. Symbols should be de-duplicated by the caller.
build_dgidb_query <- function(gene_symbols) {
  query <- "query($names: [String!]) { genes(names: $names) { nodes { name conceptId longName interactions { drug { name conceptId approved } interactionScore interactionTypes { type directionality } interactionAttributes { name value } sources { sourceDbName } } } } }"
  list(query = query, variables = list(names = as.list(gene_symbols)))
}

# Parses one DGIdb GraphQL response (parsed with jsonlite::fromJSON(...,
# simplifyVector = FALSE)) for a batch of `queried_symbols` (the exact
# strings sent: official symbols or, in a retry batch, aliases). Returns a
# list with:
#   $interactions -- one row per (queried_symbol, drug, interaction)
#     record; deduplication happens separately in dedup_interactions().
#   $zero_hit -- queried symbols recognized by DGIdb (present in `nodes`)
#     with an empty `interactions` array.
#   $unresolved -- queried symbols absent from `nodes`.
parse_dgidb_response <- function(resp, queried_symbols) {
  nodes <- resp$data$genes$nodes
  if (is.null(nodes)) nodes <- list()
  found_names <- vapply(nodes, function(n) n$name, character(1))
  unresolved <- setdiff(queried_symbols, found_names)

  interaction_rows <- list()
  zero_hit <- character(0)
  for (node in nodes) {
    gene_name <- node$name
    inter <- node$interactions
    if (length(inter) == 0) {
      zero_hit <- c(zero_hit, gene_name)
      next
    }
    for (ix in inter) {
      drug_name <- ix$drug$name
      drug_concept_id <- ix$drug$conceptId
      approved <- ix$drug$approved
      types <- ix$interactionTypes
      if (length(types) == 0) {
        type_str <- NA_character_
        dir_str <- NA_character_
      } else {
        type_vec <- vapply(types, function(t) t$type, character(1))
        dir_vec <- vapply(types, function(t) if (is.null(t$directionality)) "NA" else t$directionality, character(1))
        # Sort the (type, directionality) pairs jointly, preserving pairing,
        # so that records describing the same interaction merge in
        # dedup_interactions() regardless of the API's return order.
        ord <- order(type_vec, dir_vec)
        type_str <- paste(type_vec[ord], collapse = ";")
        dir_str <- paste(dir_vec[ord], collapse = ";")
      }
      srcs <- ix$sources
      src_str <- if (length(srcs) == 0) NA_character_ else paste(vapply(srcs, function(s) s$sourceDbName, character(1)), collapse = ";")
      interaction_rows[[length(interaction_rows) + 1]] <- data.frame(
        queried_symbol = gene_name, gene_concept_id = node$conceptId, gene_long_name = ifelse(is.null(node$longName), NA_character_, node$longName),
        drug_name = drug_name, drug_concept_id = drug_concept_id, drug_approved = isTRUE(approved),
        interaction_score = ifelse(is.null(ix$interactionScore), NA_real_, ix$interactionScore),
        interaction_types = type_str, interaction_directionality = dir_str, sources = src_str,
        stringsAsFactors = FALSE
      )
    }
  }
  interactions <- if (length(interaction_rows)) do.call(rbind, interaction_rows) else
    data.frame(queried_symbol = character(0), gene_concept_id = character(0), gene_long_name = character(0),
               drug_name = character(0), drug_concept_id = character(0), drug_approved = logical(0),
               interaction_score = numeric(0), interaction_types = character(0), interaction_directionality = character(0),
               sources = character(0), stringsAsFactors = FALSE)
  list(interactions = interactions, zero_hit = zero_hit, unresolved = unresolved)
}

# Deduplicates identical (queried_symbol, drug_name, drug_concept_id,
# interaction_types, interaction_directionality) records while preserving
# every distinct source database (merged into one sorted, semicolon-joined
# list per record) and keeping each distinct interaction type and
# directionality combination as a separate row.
dedup_interactions <- function(dt) {
  if (nrow(dt) == 0) return(dt)
  key_cols <- c("queried_symbol", "drug_name", "drug_concept_id", "interaction_types", "interaction_directionality")
  key <- do.call(paste, c(dt[key_cols], sep = ""))
  split_idx <- split(seq_len(nrow(dt)), key)
  rows <- lapply(split_idx, function(idx) {
    sub <- dt[idx, , drop = FALSE]
    all_sources <- unique(unlist(strsplit(sub$sources[!is.na(sub$sources)], ";", fixed = TRUE)))
    out <- sub[1, , drop = FALSE]
    out$sources <- if (length(all_sources)) paste(sort(all_sources), collapse = ";") else NA_character_
    out$n_source_records_collapsed <- nrow(sub)
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Re-keys an alias string returned by DGIdb for a gene queried in the alias
# retry step (DGIdb echoes the name that resolved the query, not
# necessarily the official symbol) onto the official gene symbol. Used by
# both re-keying call sites in scripts/21_drug_gene_annotation.R (the
# interactions loop and the zero-hit classification).
#
# `alias_map` must have columns `gene_symbol` (official symbol),
# `alias_used` (queried alias) and `resolved` (whether the alias resolved
# through the API). Only rows with `resolved %in% TRUE` are considered, so
# NA aliases of unresolved genes do not propagate.
#
# An alias may belong to more than one official symbol. Returns:
#   - the single official gene_symbol, if exactly one match;
#   - `alias_string` unchanged, if there is no match (handled by the
#     unresolved-gene path elsewhere);
#   - the sentinel "UNRESOLVED_AMBIGUOUS_ALIAS", if more than one match;
#     such records are retained downstream.
resolve_alias_rekey <- function(alias_string, alias_map) {
  orig <- alias_map$gene_symbol[alias_map$alias_used == alias_string & alias_map$resolved %in% TRUE]
  if (length(orig) == 1) return(orig[1])
  if (length(orig) == 0) return(alias_string)
  "UNRESOLVED_AMBIGUOUS_ALIAS"
}
