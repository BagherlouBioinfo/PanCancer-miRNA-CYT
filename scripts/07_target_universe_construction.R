# ============================================================================
# 07_target_universe_construction.R
#
# Purpose: for each threshold's independently-derived candidate-miRNA list,
#          retrieve the complete set of eligible human TargetScan predicted
#          target genes and intersect it with the matched mRNA feature space.
# Inputs:  <run_root>/processed/stage06/candidates_threshold_*.tsv (same-run);
#          <run_root>/processed/stage03/mrna_matched.rds (same-run);
#          reference/targets/Predicted_Targets_Context_Scores.default_predictions.txt.zip;
#          reference/targets/miR_Family_Info.txt.zip
#          (immutable external references, allowlist-verified)
# Outputs: <run_root>/processed/stage07/pair_universe_threshold_0.30.tsv,
#          pair_universe_threshold_0.25.tsv, pair_universe_threshold_0.20.tsv;
#          targetscan_mapping_threshold_0.30.tsv,
#          targetscan_mapping_threshold_0.25.tsv,
#          targetscan_mapping_threshold_0.20.tsv;
#          target_universe_summary.tsv
# Method:  TargetScan records and family metadata are restricted to human
#          entries (Tax ID 9606). Candidate identifiers are mapped in a fixed
#          priority order: direct prediction identifier, canonical
#          accession-bearing alias, and unique seed-family representative.
#          Ambiguous mappings stop execution. Pair tables retain the original
#          expression-matrix miRNA identifier; the TargetScan identifier is
#          recorded separately in mapping audit tables. Target genes are
#          restricted to symbols present in the matched mRNA feature space.
#          No expression correlation is computed in this stage.
# Depends on: 03_harmonize_identifiers.R, 06_mirna_cyt_discovery.R
# ============================================================================

run_stage07 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "07_target_universe_construction")
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "targetscan_mapping.R"), local = TRUE)
  suppressPackageStartupMessages(library(data.table))
  script_path <- file.path(ctx$data_root, "scripts", "07_target_universe_construction.R")
  mapping_helper_path <- file.path(ctx$data_root, "scripts", "lib", "targetscan_mapping.R")

  mrna_matched <- read_run_artifact(
    ctx,
    file.path("processed", "stage03", "mrna_matched.rds"),
    "03_harmonize_identifiers",
    readRDS
  )
  matched_genes <- unique(toupper(trimws(mrna_matched$gene_id_raw)))
  mrna_input_path <- file.path(ctx$run_root, "processed", "stage03", "mrna_matched.rds")

  ts_path <- file.path(ctx$data_root, analysis_config$paths$targetscan)
  family_path <- file.path(ctx$data_root, analysis_config$paths$mir_family_info)
  verify_allowed_raw_input(ctx$data_root, ts_path)
  verify_allowed_raw_input(ctx$data_root, family_path)

  log_fn("Reading TargetScan prediction reference (%s) ...", ts_path)
  ts_raw <- fread(ts_path)
  ts_raw <- ts_raw[
    `Gene Tax ID` == 9606,
    .(
      targetscan_mirna_id = trimws(miRNA),
      gene_symbol = toupper(trimws(`Gene Symbol`))
    )
  ]

  log_fn("Reading TargetScan miRNA family reference (%s) ...", family_path)
  family_raw <- fread(family_path)

  run_stage_assert(
    ctx,
    log_fn,
    nrow(ts_raw) > 0,
    "Human TargetScan prediction reference contains at least one record"
  )
  run_stage_assert(
    ctx,
    log_fn,
    any(family_raw$`Species ID` == 9606),
    "TargetScan family reference contains human records"
  )

  log_fn(
    "Human TargetScan records: %d; unique prediction identifiers: %d; unique target genes: %d",
    nrow(ts_raw),
    uniqueN(ts_raw$targetscan_mirna_id),
    uniqueN(ts_raw$gene_symbol)
  )

  write_fn <- function(x, p) {
    write.table(x, p, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  }

  summary_rows <- list()
  pair_lists <- list()

  for (T in analysis_config$thresholds$all) {
    cand_rel <- file.path(
      "processed",
      "stage06",
      sprintf("candidates_threshold_%.2f.tsv", T)
    )
    candidates <- read_run_artifact(
      ctx,
      cand_rel,
      "06_mirna_cyt_discovery",
      function(p) read.delim(p, stringsAsFactors = FALSE)
    )$mirna_id
    candidates <- sort(unique(trimws(candidates)))
    cand_input_path <- file.path(ctx$run_root, cand_rel)

    mapping <- build_targetscan_mapping(
      candidate_ids = candidates,
      prediction_ids = unique(ts_raw$targetscan_mirna_id),
      family_table = family_raw
    )

    run_stage_assert(
      ctx,
      log_fn,
      nrow(mapping) == length(candidates),
      "Threshold %.2f mapping contains exactly one row for every candidate miRNA (%d)",
      T,
      length(candidates)
    )
    run_stage_assert(
      ctx,
      log_fn,
      !anyDuplicated(mapping$mirna_id),
      "Threshold %.2f mapping contains no duplicated candidate identifiers",
      T
    )

    mapping_fname <- sprintf("targetscan_mapping_threshold_%.2f.tsv", T)
    stage_inputs <- c(
      mrna_input_path,
      ts_path,
      family_path,
      cand_input_path,
      mapping_helper_path
    )
    write_run_artifact(
      ctx,
      mapping,
      file.path("processed", "stage07", mapping_fname),
      "07_target_universe_construction",
      script_path,
      stage_inputs,
      write_fn
    )

    mapped <- mapping[
      mapping_type != "unresolved" & !is.na(targetscan_mirna_id)
    ]

    ts_records <- merge(
      mapped[, .(mirna_id, targetscan_mirna_id)],
      ts_raw,
      by = "targetscan_mirna_id",
      allow.cartesian = TRUE,
      sort = FALSE
    )

    pairs <- unique(ts_records[, .(mirna_id, gene = gene_symbol)])
    pairs <- pairs[gene %in% matched_genes]
    setorder(pairs, mirna_id, gene)

    run_stage_assert(
      ctx,
      log_fn,
      all(pairs$mirna_id %in% candidates),
      "Threshold %.2f pair universe retains only original candidate-miRNA identifiers",
      T
    )
    run_stage_assert(
      ctx,
      log_fn,
      all(pairs$gene %in% matched_genes),
      "Threshold %.2f pair universe contains only genes in the matched mRNA feature space",
      T
    )
    run_stage_assert(
      ctx,
      log_fn,
      !anyDuplicated(pairs[, .(mirna_id, gene)]),
      "Threshold %.2f pair universe contains no duplicated miRNA-gene pairs",
      T
    )

    pair_fname <- sprintf("pair_universe_threshold_%.2f.tsv", T)
    write_run_artifact(
      ctx,
      pairs,
      file.path("processed", "stage07", pair_fname),
      "07_target_universe_construction",
      script_path,
      stage_inputs,
      write_fn
    )
    pair_lists[[as.character(T)]] <- pairs

    mapping_counts <- mapping[, .N, by = mapping_type]
    count_type <- function(type) {
      value <- mapping_counts[mapping_type == type, N]
      if (length(value) == 0) 0L else as.integer(value)
    }

    n_direct <- count_type("direct")
    n_alias <- count_type("canonical_accession_alias")
    n_family <- count_type("seed_family_representative")
    n_unresolved <- count_type("unresolved")
    n_mapped <- nrow(mapped)

    log_fn(
      paste0(
        "Threshold %.2f: %d candidate miRNAs -> %d mapped ",
        "(direct=%d, canonical alias=%d, seed-family representative=%d; unresolved=%d) ",
        "-> %d TargetScan records -> %d unique pairs -> %d unique target genes"
      ),
      T,
      length(candidates),
      n_mapped,
      n_direct,
      n_alias,
      n_family,
      n_unresolved,
      nrow(ts_records),
      nrow(pairs),
      uniqueN(pairs$gene)
    )

    summary_rows[[as.character(T)]] <- data.frame(
      threshold = T,
      n_candidate_mirnas = length(candidates),
      n_mapped_mirnas = n_mapped,
      n_direct_mappings = n_direct,
      n_canonical_alias_mappings = n_alias,
      n_seed_family_mappings = n_family,
      n_unresolved_mirnas = n_unresolved,
      n_targetscan_records = nrow(ts_records),
      n_unique_pairs = nrow(pairs),
      n_unique_genes = uniqueN(pairs$gene),
      stringsAsFactors = FALSE
    )
  }

  sorted_thresholds <- sort(analysis_config$thresholds$all, decreasing = TRUE)
  for (i in seq_len(length(sorted_thresholds) - 1)) {
    tighter <- pair_lists[[as.character(sorted_thresholds[i])]]
    looser <- pair_lists[[as.character(sorted_thresholds[i + 1])]]
    not_subset <- fsetdiff(tighter, looser)
    run_stage_assert(
      ctx,
      log_fn,
      nrow(not_subset) == 0,
      "Threshold %.2f pair universe is a subset of threshold %.2f's pair universe (0 exceptions)",
      sorted_thresholds[i],
      sorted_thresholds[i + 1]
    )
  }

  summary_table <- rbindlist(summary_rows, use.names = TRUE)
  write_run_artifact(
    ctx,
    summary_table,
    file.path("processed", "stage07", "target_universe_summary.tsv"),
    "07_target_universe_construction",
    script_path,
    c(mrna_input_path, ts_path, family_path, mapping_helper_path),
    write_fn
  )

  record_stage_status(ctx, "07_target_universe_construction", "COMPLETE")
  log_fn("Stage complete.")
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=7 or source scripts/run.R in RStudio.")
}
