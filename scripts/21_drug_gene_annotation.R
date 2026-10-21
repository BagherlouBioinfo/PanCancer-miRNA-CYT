# ============================================================================
# 22_drug_gene_annotation.R
#
# Purpose: drug-gene annotation for the unique target genes of the current
#          run's primary triplets through the DGIdb public GraphQL API
#          (https://dgidb.org/api/graphql). This is database annotation
#          only; no claim of drug repurposing, treatment-response
#          prediction, clinical actionability, or drug action through the
#          miRNA relationship is made.
# Inputs:  <run_root>/tables/primary/triplets_rho030_CPE.tsv (primary
#          triplets, via load_run_primary_triplets(ctx)).
# Outputs: <run_root>/tables/downstream/drug_gene_id_resolution.tsv,
#          drug_gene_interactions.tsv, drug_gene_coverage_summary.tsv;
#          <run_root>/checkpoints/dgidb_raw_response_<timestamp>.json
#          (raw API responses, cached for reproducibility);
#          reports/26_drug_gene_annotation_report.md.
# Method:  all unique target genes are queried in one batched request by
#          official gene symbol; genes not recognized by DGIdb are retried
#          once with their org.Hs.eg.db alias symbols, recording which query
#          term resolved each gene (R/drug_annotation.R). Each raw JSON
#          response is hashed (SHA-256) and cached under the run's
#          checkpoints/ directory. Every gene is represented in the final
#          table as resolved with interactions, resolved with zero
#          interactions (a retained negative result) or unresolved (neither
#          the symbol nor any alias was recognized). Records are
#          deduplicated on (gene, drug, interaction_types,
#          interaction_directionality) while preserving every supporting
#          source database (R/drug_annotation.R::dedup_interactions());
#          interaction types are kept distinct rather than collapsed into a
#          single indicator. Because DGIdb content changes over time, the
#          access date and response hashes are recorded in the report.
# Depends on: 09_recurrence_analysis.R, 10_purity_adjusted_analysis.R
#             (via load_run_primary_triplets(ctx)).
# ============================================================================

run_stage22 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "22_drug_gene_annotation")
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "run_discovery_source.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "drug_annotation.R"), local = TRUE)
  suppressPackageStartupMessages({
    library(data.table)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(httr2)
    library(jsonlite)
    library(digest)
  })
  script_path <- file.path(ctx$data_root, "scripts", "21_drug_gene_annotation.R")
  fwrite_fn <- function(x, p) fwrite(x, p, sep = "\t")

  # ==========================================================================
  # 1. Genes: the current run's primary triplets.
  # ==========================================================================
  triplets <- load_run_primary_triplets(ctx)
  triplets_path <- file.path(ctx$run_root, "tables", "primary", run_primary_triplet_filename(ctx))
  genes <- sort(unique(triplets$gene))
  run_stage_assert(ctx, log_fn, length(genes) > 0, "This run's primary triplets contain at least one target gene (%d found)", length(genes))
  log_fn("Querying DGIdb for this run's %d unique primary-triplet target gene(s).", length(genes))

  api_url <- analysis_config$drug_annotation$api_url
  package_versions <- sprintf("httr2 %s, jsonlite %s", as.character(packageVersion("httr2")), as.character(packageVersion("jsonlite")))
  access_date <- format(Sys.time(), "%Y-%m-%d %H:%M:%OS0 UTC", tz = "UTC")

  # ==========================================================================
  # 2-4. Query DGIdb: official symbols first, then a retry batch for
  # anything unresolved using org.Hs.eg.db ALIAS symbols.
  # ==========================================================================
  do_query <- function(symbols) {
    body <- build_dgidb_query(symbols)
    req <- httr2::request(api_url) |> httr2::req_headers(`Content-Type` = "application/json") |>
      httr2::req_body_json(body) |> httr2::req_retry(max_tries = 3)
    resp <- httr2::req_perform(req)
    raw_text <- httr2::resp_body_string(resp)
    list(parsed = jsonlite::fromJSON(raw_text, simplifyVector = FALSE), raw_text = raw_text, status = httr2::resp_status(resp))
  }

  primary <- do_query(genes)
  run_stage_assert(ctx, log_fn, primary$status == 200, "DGIdb API request succeeded (HTTP status %d)", primary$status)
  primary_parsed <- parse_dgidb_response(primary$parsed, genes)
  log_fn("DGIdb primary query (official symbols): %d resolved with interactions, %d resolved with zero interactions, %d unresolved.",
         uniqueN(primary_parsed$interactions$queried_symbol), length(primary_parsed$zero_hit), length(primary_parsed$unresolved))

  raw_responses <- list(primary = primary$raw_text)
  query_log <- data.table(batch = "primary", query_terms = paste(genes, collapse = ";"), n_terms = length(genes), http_status = primary$status)

  alias_map <- data.table(gene_symbol = character(0), alias_used = character(0), resolved = logical(0))
  retry_parsed <- list(interactions = primary_parsed$interactions[0, ], zero_hit = character(0), unresolved = character(0))
  if (length(primary_parsed$unresolved) > 0) {
    alias_lookup <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = primary_parsed$unresolved, keytype = "SYMBOL", columns = "ALIAS"))
    alias_lookup <- alias_lookup[!is.na(alias_lookup$ALIAS) & alias_lookup$ALIAS != alias_lookup$SYMBOL, ]
    if (nrow(alias_lookup) > 0) {
      alias_terms <- unique(alias_lookup$ALIAS)
      log_fn("Retrying %d unresolved gene(s) using %d org.Hs.eg.db ALIAS symbol(s).", length(primary_parsed$unresolved), length(alias_terms))
      retry <- do_query(alias_terms)
      run_stage_assert(ctx, log_fn, retry$status == 200, "DGIdb alias-retry API request succeeded (HTTP status %d)", retry$status)
      retry_parsed_raw <- parse_dgidb_response(retry$parsed, alias_terms)
      raw_responses$alias_retry <- retry$raw_text
      query_log <- rbind(query_log, data.table(batch = "alias_retry", query_terms = paste(alias_terms, collapse = ";"), n_terms = length(alias_terms), http_status = retry$status))
      # Map each resolved alias back to its original unresolved official symbol.
      resolved_aliases <- unique(c(retry_parsed_raw$interactions$queried_symbol, retry_parsed_raw$zero_hit))
      for (a in resolved_aliases) {
        orig <- unique(alias_lookup$SYMBOL[alias_lookup$ALIAS == a])
        for (o in orig) alias_map <- rbind(alias_map, data.table(gene_symbol = o, alias_used = a, resolved = TRUE))
      }
      still_unresolved <- setdiff(primary_parsed$unresolved, alias_map$gene_symbol)
      for (o in still_unresolved) alias_map <- rbind(alias_map, data.table(gene_symbol = o, alias_used = NA_character_, resolved = FALSE))
      retry_parsed <- retry_parsed_raw
    } else {
      alias_map <- data.table(gene_symbol = primary_parsed$unresolved, alias_used = NA_character_, resolved = FALSE)
    }
  }

  # Raw response provenance: each batch's raw JSON is hashed and cached
  # under the run's checkpoints/ directory.
  checkpoints_dir <- file.path(ctx$run_root, "checkpoints")
  if (!dir.exists(checkpoints_dir)) dir.create(checkpoints_dir, recursive = TRUE)
  response_hashes <- list()
  for (nm in names(raw_responses)) {
    fpath <- file.path(checkpoints_dir, sprintf("dgidb_raw_response_%s.json", nm))
    writeLines(raw_responses[[nm]], fpath)
    response_hashes[[nm]] <- digest::digest(fpath, algo = "sha256", file = TRUE)
    log_fn("Cached immutable raw DGIdb response (%s batch): %s (SHA-256 %s).", nm, fpath, response_hashes[[nm]])
  }

  # ==========================================================================
  # 6. Gene-ID / alias resolution audit table.
  # ==========================================================================
  resolution_dt <- data.table(gene_symbol = genes)
  resolution_dt[, resolved_directly := gene_symbol %in% c(primary_parsed$interactions$queried_symbol, primary_parsed$zero_hit)]
  resolution_dt <- merge(resolution_dt, alias_map, by = "gene_symbol", all.x = TRUE)
  resolution_dt[, status := fifelse(resolved_directly, "resolved_official_symbol",
                              fifelse(resolved %in% TRUE, "resolved_via_alias",
                              "unresolved"))]
  write_run_artifact(ctx, resolution_dt, file.path("tables", "downstream", "drug_gene_id_resolution.tsv"),
                      "22_drug_gene_annotation", script_path, triplets_path, fwrite_fn)

  # ==========================================================================
  # 7-9. Combine, deduplicate, and ensure that every queried gene has a row
  # (interactions, explicit zero-hit, or unresolved).
  # ==========================================================================
  all_interactions_raw <- rbind(primary_parsed$interactions, retry_parsed$interactions)
  # Re-key alias-resolved rows onto the original official gene symbol. An
  # alias shared by more than one official symbol has no single correct
  # target and is flagged with a retained sentinel value
  # (R/drug_annotation.R::resolve_alias_rekey(), tested in
  # tests/test_drug_alias_rekeying.R).
  if (nrow(retry_parsed$interactions) > 0) {
    for (i in seq_len(nrow(all_interactions_raw))) {
      all_interactions_raw$queried_symbol[i] <- resolve_alias_rekey(all_interactions_raw$queried_symbol[i], alias_map)
    }
  }
  deduped <- as.data.table(dedup_interactions(all_interactions_raw))
  if (nrow(deduped) > 0) {
    deduped[, record_status := fifelse(queried_symbol == "UNRESOLVED_AMBIGUOUS_ALIAS",
                                        "unresolved_ambiguous_alias", "database_listed_interaction")]
  }

  resolved_zero_hit <- unique(c(primary_parsed$zero_hit, {
    if (nrow(retry_parsed$interactions) >= 0 && length(retry_parsed$zero_hit) > 0) {
      vapply(retry_parsed$zero_hit, resolve_alias_rekey, character(1), alias_map = alias_map)
    } else character(0)
  }))
  final_unresolved <- resolution_dt$gene_symbol[resolution_dt$status == "unresolved"]

  zero_hit_rows <- if (length(resolved_zero_hit) > 0) data.table(
    queried_symbol = resolved_zero_hit, gene_concept_id = NA_character_, gene_long_name = NA_character_,
    drug_name = NA_character_, drug_concept_id = NA_character_, drug_approved = NA,
    interaction_score = NA_real_, interaction_types = NA_character_, interaction_directionality = NA_character_,
    sources = NA_character_, n_source_records_collapsed = 0L,
    record_status = fifelse(resolved_zero_hit == "UNRESOLVED_AMBIGUOUS_ALIAS",
                             "unresolved_ambiguous_alias", "resolved_zero_interactions")
  ) else data.table()
  unresolved_rows <- if (length(final_unresolved) > 0) data.table(
    queried_symbol = final_unresolved, gene_concept_id = NA_character_, gene_long_name = NA_character_,
    drug_name = NA_character_, drug_concept_id = NA_character_, drug_approved = NA,
    interaction_score = NA_real_, interaction_types = NA_character_, interaction_directionality = NA_character_,
    sources = NA_character_, n_source_records_collapsed = 0L, record_status = "unresolved_by_dgidb"
  ) else data.table()

  final_dt <- rbindlist(list(deduped, zero_hit_rows, unresolved_rows), use.names = TRUE, fill = TRUE)
  setnames(final_dt, "queried_symbol", "gene_symbol")
  setcolorder(final_dt, c("gene_symbol", "record_status", "drug_name", "drug_concept_id", "drug_approved",
                           "interaction_types", "interaction_directionality", "sources", "interaction_score",
                           "n_source_records_collapsed", "gene_concept_id", "gene_long_name"))
  # Every queried gene must have at least one row in the final table
  # (interaction, zero-hit or unresolved). The only permitted gene value
  # outside the queried set is the "UNRESOLVED_AMBIGUOUS_ALIAS" sentinel;
  # any other value fails the assertion.
  final_gene_symbols <- unique(final_dt$gene_symbol)
  every_gene_present <- all(genes %in% final_gene_symbols)
  extra_symbols <- setdiff(final_gene_symbols, genes)
  extras_are_only_sentinel <- all(extra_symbols %in% "UNRESOLVED_AMBIGUOUS_ALIAS")
  run_stage_assert(ctx, log_fn, every_gene_present && extras_are_only_sentinel,
                    "Every one of this run's %d genes has at least one row in the final drug-gene table (interaction, zero-hit, or unresolved), with only the 'UNRESOLVED_AMBIGUOUS_ALIAS' sentinel allowed as an additional gene_symbol value", length(genes))

  write_run_artifact(ctx, final_dt, file.path("tables", "downstream", "drug_gene_interactions.tsv"),
                      "22_drug_gene_annotation", script_path, triplets_path, fwrite_fn)

  # ==========================================================================
  # 14. Exact coverage report.
  # ==========================================================================
  n_with_interactions <- uniqueN(deduped$queried_symbol)
  n_zero <- length(resolved_zero_hit)
  n_unresolved <- length(final_unresolved)
  coverage_dt <- data.table(
    metric = c("n_genes_queried", "n_genes_mapped_resolved", "n_genes_with_ge1_interaction", "n_genes_zero_interactions", "n_genes_unresolved", "n_total_interaction_records_after_dedup"),
    value = c(length(genes), length(genes) - n_unresolved, n_with_interactions, n_zero, n_unresolved, nrow(deduped))
  )
  write_run_artifact(ctx, coverage_dt, file.path("tables", "downstream", "drug_gene_coverage_summary.tsv"),
                      "22_drug_gene_annotation", script_path, triplets_path, fwrite_fn)
  log_fn("Coverage: %d genes queried, %d mapped/resolved, %d with >=1 interaction, %d zero-interaction, %d unresolved; %d interaction records after dedup.",
         length(genes), length(genes) - n_unresolved, n_with_interactions, n_zero, n_unresolved, nrow(deduped))

  # ==========================================================================
  # Report
  # ==========================================================================
  report <- c(
    "# Stage 22: Drug-Gene Annotation (DGIdb)",
    "",
    sprintf("Generated by `scripts/21_drug_gene_annotation.R`, run `%s`, %s.", ctx$run_id, format(Sys.time(), "%Y-%m-%d %H:%M:%OS0")),
    "",
    "**This is database annotation only.** No claim of drug repurposing,",
    "treatment-response prediction, or clinical actionability is made, and no",
    "drug is claimed to act through the discovered miRNA relationship. The",
    "table describes the target genes themselves; evidence on miRNA-target",
    "pair specificity is reported separately by the Stage 19 matched-null",
    "target-specificity test (tables/audit/matched_null_target_specificity_results.tsv).",
    "",
    "## Provenance",
    sprintf("- Database: DGIdb, GraphQL API (`%s`)", api_url),
    sprintf("- Access date: %s", access_date),
    sprintf("- Package versions: %s", package_versions),
    sprintf("- Raw response SHA-256: %s", paste(sprintf("%s=%s", names(response_hashes), unlist(response_hashes)), collapse = "; ")),
    "- Raw responses cached under this run's own `checkpoints/` directory (immutable, not committed to source control or the public GitHub copy).",
    "",
    "## Coverage",
    sprintf("- %d genes queried (this run's own unique primary-triplet target genes).", length(genes)),
    sprintf("- %d mapped/resolved by DGIdb (official symbol or alias); %d unresolved.", length(genes) - n_unresolved, n_unresolved),
    sprintf("- %d genes with >=1 reported interaction; %d genes resolved with zero reported interactions (retained, not dropped).", n_with_interactions, n_zero),
    sprintf("- %d distinct interaction records after deduplication (sources preserved, interaction types never collapsed).", nrow(deduped)),
    "",
    "## What this table does NOT claim",
    "- Database presence does not mean a gene is \"druggable\" or clinically actionable.",
    "- `drug_approved` reflects DGIdb's own field for that drug in general use, not approval for any indication related to this project's cancer types or this specific gene-drug pair.",
    "- Interaction type/directionality/source are reported exactly as DGIdb provides them; fields DGIdb does not provide are left as `NA` (unknown), never inferred.",
    "- No mechanism, direction of action, or effect on cytolytic immunity is inferred for any gene-drug pair."
  )
  write_run_artifact(ctx, report, file.path("reports", "26_drug_gene_annotation_report.md"),
                      "22_drug_gene_annotation", script_path, character(0), function(x, p) writeLines(x, p))
  log_fn("reports/26_drug_gene_annotation_report.md written.")

  record_stage_status(ctx, "22_drug_gene_annotation", "COMPLETE")
  log_fn("Stage complete.")
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=21 or source scripts/run.R in RStudio.")
}
