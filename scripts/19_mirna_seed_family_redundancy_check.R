# ============================================================================
# 20_mirna_seed_family_redundancy_check.R  (Stage 20)
#
# Purpose: determine whether any primary miRNAs of the current run are
#          biologically redundant (same TargetScan seed family or seed
#          sequence, implying largely overlapping predicted target sets),
#          and if so whether they carry independent, related or clearly
#          redundant information. Several human miRNA paralogs share a seed
#          family (e.g. miR-146a-5p/miR-146b-5p, miR-199a-3p/miR-199b-3p,
#          miR-221-3p/miR-222-3p); whether such pairs occur in the candidate
#          set, and how redundant they are, is determined empirically from
#          the run's data against the TargetScan seed-family reference
#          (reference/targets/miR_Family_Info.txt.zip). The candidate set is
#          read from the current run's primary table through
#          load_run_primary_triplets(ctx); no miRNA identity or count is
#          hard-coded.
#
# Prespecification:
#   - Candidate set: every unique miRNA in the current run's primary
#     triplet table (qualifies == TRUE), read without modification.
#   - Redundancy criterion: two miRNAs are flagged as same-family if the
#     TargetScan family assignment (miR_Family_Info.txt, Species ID 9606)
#     places them in the same "miR family"; seed-sequence identity
#     ("Seed+m8") is reported in a separate column. This is TargetScan's
#     published family definition.
#   - Target-overlap statistic: Jaccard index of each flagged pair's
#     predicted-target gene sets (Gene Tax ID 9606), computed from the same
#     TargetScan file with the human filter and upper-case symbol convention
#     of scripts/07_target_universe_construction.R.
#   - Expression-similarity statistics (two, reported separately):
#       (a) Within-cancer: for each eligible cancer in which both miRNAs
#           have valid expression, the Spearman correlation of their
#           expression across that cancer's samples. Reported per cancer
#           (detail table) and as the median across shared eligible
#           cancers (summary). This statistic addresses whether the two
#           miRNAs carry independent information within one tissue
#           context.
#       (b) Pan-cohort: the same correlation computed once across the full
#           matched cohort. Pooling cancers can inflate the correlation of
#           two features whose typical expression differs between cancers
#           (a Simpson's-paradox-type confound; see
#           tests/test_synthetic_scenarios.R, scenario S07), so this value
#           is not adjusted for cancer of origin and is not interpreted as
#           within-cancer co-expression.
#   - Redundancy verdict (fixed rule): a pair is labelled "clearly
#     redundant" if Seed+m8 is identical, the Jaccard index is >= 0.5 and
#     the median within-cancer |rho| is >= 0.5; "merely related" if it shares
#     a family but any of these conditions fails. Every pair sharing a
#     TargetScan family is reported regardless of verdict.
#   - The analysis characterizes the primary set only; it does not remove,
#     merge or re-select any miRNA or triplet and writes nothing under
#     tables/primary/.
# ============================================================================

run_stage20 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "20_mirna_seed_family_redundancy_check")
  source(file.path(ctx$data_root, "scripts", "lib", "run_discovery_source.R"), local = TRUE)
  suppressPackageStartupMessages(library(data.table))
  script_path <- file.path(ctx$data_root, "scripts", "19_mirna_seed_family_redundancy_check.R")
  fwrite_fn <- function(x, p) fwrite(x, p, sep = "\t")

  log_fn("miRNA seed-family/target-overlap redundancy check starting.")

  # --- 1. This run's own primary candidate set (read-only) -----------------
  triplets <- load_run_primary_triplets(ctx)
  run_primary_mirnas <- sort(unique(triplets$mirna_id))
  log_fn("This run's primary miRNAs (n=%d): %s", length(run_primary_mirnas), paste(run_primary_mirnas, collapse = ", "))

  # --- 2. Seed-family assignment (TargetScan's own definition) -------------
  fam_zip <- file.path(ctx$data_root, "reference", "targets", "miR_Family_Info.txt.zip")
  verify_allowed_raw_input(ctx$data_root, fam_zip)
  fam_dir <- file.path(tempdir(), sprintf("mir_family_info_%s", ctx$run_id))
  dir.create(fam_dir, showWarnings = FALSE)
  unzip(fam_zip, exdir = fam_dir, overwrite = TRUE)
  fam_raw <- fread(file.path(fam_dir, "miR_Family_Info.txt"), sep = "\t", quote = "")
  fam_hs <- fam_raw[`Species ID` == 9606, .(
    mirna_id = trimws(`MiRBase ID`),
    family = trimws(`miR family`),
    seed = trimws(`Seed+m8`)
  )]
  fam_hs <- unique(fam_hs)

  fam_run <- fam_hs[mirna_id %in% run_primary_mirnas]
  missing_from_family_ref <- setdiff(run_primary_mirnas, fam_run$mirna_id)
  if (length(missing_from_family_ref) > 0) {
    log_fn("WARNING: %d primary miRNA(s) not found in miR_Family_Info.txt: %s",
           length(missing_from_family_ref), paste(missing_from_family_ref, collapse = ", "))
  }
  write_run_artifact(ctx, fam_run[order(family, mirna_id)], file.path("tables", "audit", "mirna_redundancy_seed_family_assignment.tsv"),
                      "20_mirna_seed_family_redundancy_check", script_path, fam_zip, fwrite_fn)
  log_fn("Seed-family assignment for %d/%d primary miRNAs written.", nrow(fam_run), length(run_primary_mirnas))

  # Pairs sharing a seed family
  fam_counts <- fam_run[, .N, by = family]
  redundant_families <- fam_counts[N > 1, family]
  pair_rows <- list()
  for (fm in redundant_families) {
    members <- fam_run[family == fm, mirna_id]
    combos <- combn(sort(members), 2, simplify = FALSE)
    for (cb in combos) {
      pair_rows[[length(pair_rows) + 1]] <- data.table(
        mirna_a = cb[1], mirna_b = cb[2], shared_family = fm,
        seed_a = fam_run[mirna_id == cb[1], seed][1],
        seed_b = fam_run[mirna_id == cb[2], seed][1]
      )
    }
  }
  same_family_pairs <- if (length(pair_rows)) rbindlist(pair_rows) else
    data.table(mirna_a = character(), mirna_b = character(), shared_family = character(), seed_a = character(), seed_b = character())
  same_family_pairs[, identical_seed := seed_a == seed_b]
  log_fn("Same-seed-family pairs found within the %d primary miRNAs: %d.", length(run_primary_mirnas), nrow(same_family_pairs))
  for (i in seq_len(nrow(same_family_pairs))) {
    log_fn("  %s / %s -- family '%s', identical Seed+m8: %s",
           same_family_pairs$mirna_a[i], same_family_pairs$mirna_b[i],
           same_family_pairs$shared_family[i], same_family_pairs$identical_seed[i])
  }

  # --- 3. Predicted-target-set Jaccard overlap for each flagged pair -------
  ts_zip <- file.path(ctx$data_root, "reference", "targets", "Predicted_Targets_Context_Scores.default_predictions.txt.zip")
  overlap_rows <- list()
  if (nrow(same_family_pairs) > 0) {
    verify_allowed_raw_input(ctx$data_root, ts_zip)
    ts_dir <- file.path(tempdir(), sprintf("targetscan_predictions_%s", ctx$run_id))
    dir.create(ts_dir, showWarnings = FALSE)
    unzip(ts_zip, exdir = ts_dir, overwrite = TRUE)
    ts_raw <- fread(file.path(ts_dir, "Predicted_Targets_Context_Scores.default_predictions.txt"),
                     sep = "\t", quote = "", select = c("Gene Symbol", "Gene Tax ID", "miRNA"))
    ts_hs <- ts_raw[`Gene Tax ID` == 9606, .(gene_symbol = toupper(trimws(`Gene Symbol`)), mirna_id = trimws(miRNA))]
    ts_hs <- unique(ts_hs)
    log_fn("TargetScan human predicted-target rows (unique gene-mirna pairs): %d.", nrow(ts_hs))

    for (i in seq_len(nrow(same_family_pairs))) {
      ma <- same_family_pairs$mirna_a[i]; mb <- same_family_pairs$mirna_b[i]
      targets_a <- unique(ts_hs[mirna_id == ma, gene_symbol])
      targets_b <- unique(ts_hs[mirna_id == mb, gene_symbol])
      inter <- intersect(targets_a, targets_b)
      uni <- union(targets_a, targets_b)
      jaccard <- if (length(uni) > 0) length(inter) / length(uni) else NA_real_

      run_genes_a <- sort(unique(triplets[mirna_id == ma, gene]))
      run_genes_b <- sort(unique(triplets[mirna_id == mb, gene]))
      run_shared <- intersect(run_genes_a, run_genes_b)

      overlap_rows[[length(overlap_rows) + 1]] <- data.table(
        mirna_a = ma, mirna_b = mb,
        n_targets_a = length(targets_a), n_targets_b = length(targets_b),
        n_targets_intersection = length(inter), n_targets_union = length(uni),
        jaccard_index = jaccard,
        run_primary_genes_a = paste(run_genes_a, collapse = ";"),
        run_primary_genes_b = paste(run_genes_b, collapse = ";"),
        run_primary_genes_shared = paste(run_shared, collapse = ";"),
        n_run_primary_genes_shared = length(run_shared)
      )
      log_fn("  %s vs %s: %d predicted targets each side overlap, Jaccard=%.4f; this run's own primary-triplet genes shared: %s",
             ma, mb, length(inter), jaccard, if (length(run_shared)) paste(run_shared, collapse = ",") else "(none)")
    }
  }
  overlap_dt <- if (length(overlap_rows)) rbindlist(overlap_rows) else
    data.table(mirna_a = character(), mirna_b = character())
  write_run_artifact(ctx, overlap_dt, file.path("tables", "audit", "mirna_redundancy_target_overlap.tsv"),
                      "20_mirna_seed_family_redundancy_check", script_path, ts_zip, fwrite_fn)

  # --- 4. Expression similarity: within-cancer (primary statistic) and
  #        pan-cohort (subject to between-cancer confounding), reported separately --
  per_cancer_rows <- list()
  pan_cohort_rows <- list()
  mirna_matched_path <- file.path(ctx$run_root, "processed", "stage03", "mirna_matched.rds")
  if (nrow(same_family_pairs) > 0) {
    mirna_expr <- read_run_artifact(ctx, file.path("processed", "stage03", "mirna_matched.rds"), "03_harmonize_identifiers", readRDS)
    setDT(mirna_expr)
    sample_annot <- read_run_artifact(ctx, file.path("processed", "stage03", "sample_annotation.tsv"), "03_harmonize_identifiers",
                                       function(p) read.delim(p, stringsAsFactors = FALSE))
    sample_cols <- setdiff(names(mirna_expr), "mirna_id")
    cancer_of_sample <- sample_annot$cancer[match(sample_cols, sample_annot$sample_id)]

    for (i in seq_len(nrow(same_family_pairs))) {
      ma <- same_family_pairs$mirna_a[i]; mb <- same_family_pairs$mirna_b[i]
      row_a <- mirna_expr[mirna_id == ma]; row_b <- mirna_expr[mirna_id == mb]
      if (nrow(row_a) != 1 || nrow(row_b) != 1) {
        log_fn("  WARNING: could not locate a unique expression row for %s and/or %s in mirna_matched.rds.", ma, mb)
        next
      }
      vec_a <- as.numeric(row_a[, ..sample_cols]); vec_b <- as.numeric(row_b[, ..sample_cols])
      complete <- is.finite(vec_a) & is.finite(vec_b)

      # Within-cancer: one Spearman rho per cancer where both miRNAs have
      # expression, computed on that cancer's own samples only.
      cancers_here <- sort(unique(cancer_of_sample[complete]))
      within_rhos <- numeric(0)
      for (cc in cancers_here) {
        idx <- complete & !is.na(cancer_of_sample) & cancer_of_sample == cc
        if (sum(idx) < 10) next  # too few samples in this cancer for a stable within-cancer estimate
        rho_cc <- suppressWarnings(cor(vec_a[idx], vec_b[idx], method = "spearman"))
        within_rhos <- c(within_rhos, rho_cc)
        per_cancer_rows[[length(per_cancer_rows) + 1]] <- data.table(
          mirna_a = ma, mirna_b = mb, cancer = cc, n_samples = sum(idx), within_cancer_spearman_rho = rho_cc
        )
      }
      median_within <- if (length(within_rhos)) median(within_rhos, na.rm = TRUE) else NA_real_

      # Pan-cohort: pooled across all cancers (see the caveat in the header).
      rho_pooled <- suppressWarnings(cor(vec_a[complete], vec_b[complete], method = "spearman"))
      pan_cohort_rows[[length(pan_cohort_rows) + 1]] <- data.table(
        mirna_a = ma, mirna_b = mb, n_samples_pooled = sum(complete),
        pan_cohort_spearman_rho_UNADJUSTED_FOR_CANCER = rho_pooled,
        n_cancers_with_within_cancer_estimate = length(within_rhos),
        median_within_cancer_spearman_rho = median_within
      )
      log_fn("  %s vs %s: median within-cancer rho=%.4f (%d cancers), pan-cohort (unadjusted) rho=%.4f (n=%d).",
             ma, mb, median_within, length(within_rhos), rho_pooled, sum(complete))
    }
  }
  per_cancer_dt <- if (length(per_cancer_rows)) rbindlist(per_cancer_rows) else
    data.table(mirna_a = character(), mirna_b = character(), cancer = character(), n_samples = integer(), within_cancer_spearman_rho = numeric())
  pan_cohort_dt <- if (length(pan_cohort_rows)) rbindlist(pan_cohort_rows) else
    data.table(mirna_a = character(), mirna_b = character())
  write_run_artifact(ctx, per_cancer_dt, file.path("tables", "audit", "mirna_redundancy_expression_correlation_within_cancer.tsv"),
                      "20_mirna_seed_family_redundancy_check", script_path, mirna_matched_path, fwrite_fn)
  write_run_artifact(ctx, pan_cohort_dt, file.path("tables", "audit", "mirna_redundancy_expression_correlation_pan_cohort.tsv"),
                      "20_mirna_seed_family_redundancy_check", script_path, mirna_matched_path, fwrite_fn)

  # --- 5. Combined summary table with the prespecified redundancy verdict --
  summary_dt <- same_family_pairs
  if (nrow(summary_dt) > 0) {
    summary_dt <- merge(summary_dt, overlap_dt, by = c("mirna_a", "mirna_b"), all.x = TRUE)
    summary_dt <- merge(summary_dt, pan_cohort_dt, by = c("mirna_a", "mirna_b"), all.x = TRUE)
    summary_dt[, redundancy_verdict := ifelse(
      identical_seed & !is.na(jaccard_index) & jaccard_index >= 0.5 &
        !is.na(median_within_cancer_spearman_rho) & abs(median_within_cancer_spearman_rho) >= 0.5,
      "clearly redundant", "merely related (same TargetScan family, but does not meet all three redundancy conditions)"
    )]
  } else {
    summary_dt[, redundancy_verdict := character(0)]
  }
  write_run_artifact(ctx, summary_dt, file.path("tables", "audit", "mirna_redundancy_summary.tsv"),
                      "20_mirna_seed_family_redundancy_check", script_path, mirna_matched_path, fwrite_fn)

  run_stage_assert(ctx, log_fn, !any(duplicated(paste(summary_dt$mirna_a, summary_dt$mirna_b))),
                    "No duplicate (mirna_a, mirna_b) pair in the redundancy summary")

  log_fn("DONE. %d same-seed-family pair(s) identified and characterized among this run's %d primary miRNAs (%d 'clearly redundant').",
         nrow(same_family_pairs), length(run_primary_mirnas),
         if (nrow(summary_dt) > 0) sum(summary_dt$redundancy_verdict == "clearly redundant") else 0L)

  record_stage_status(ctx, "20_mirna_seed_family_redundancy_check", "COMPLETE")
  invisible(summary_dt)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=19 or source scripts/run.R in RStudio.")
}
