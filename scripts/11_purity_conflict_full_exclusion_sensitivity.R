# ============================================================================
# 11_purity_conflict_full_exclusion_sensitivity.R
#
# Purpose: sensitivity analysis for the handling of purity (CPE/LUMP)
#          conflicts between duplicate vials in Stage 05. The primary
#          analysis (Stage 05) sets CPE and/or LUMP to missing for the
#          affected samples and excludes them only from analyses that
#          require that purity metric, keeping their expression, clinical
#          and survival data available elsewhere ("Option 1"). This stage
#          evaluates the alternative of excluding the affected samples from
#          every layer ("Option 2") for comparison only; it never replaces
#          the primary analysis. All inputs derive from the current run's
#          Stage 05 conflict census; no sample count or identifier is
#          hard-coded.
# Inputs:  <run_root>/processed/stage03/mirna_matched.rds, mrna_matched.rds;
#          stage05/sample_annotation_with_purity.tsv;
#          stage07/pair_universe_threshold_*.tsv (all same-run)
# Outputs: <run_root>/tables/audit/purity_sensitivity_opt2_triplets_unadjusted.tsv,
#          purity_sensitivity_opt2_triplets_CPE.tsv,
#          purity_sensitivity_opt2_triplets_LUMP.tsv (written under
#          tables/audit/ to keep them distinct from the prespecified
#          correlation-threshold sensitivity tables);
#          tables/audit/purity_sensitivity_option1_vs_option2_comparison.tsv,
#          purity_sensitivity_option1_vs_option2_per_triplet_changes.tsv;
#          reports/27_purity_conflict_sensitivity_report.md
# Method:  re-derives the eligible-cancer sets, the full miRNA-CYT family
#          (as in Stage 06), the unadjusted three-arm test (as in Stages 08
#          and 09) and the CPE/LUMP-adjusted three-arm test (as in Stage 10,
#          including per-cancer winsorization) at the primary threshold,
#          using the production functions of scripts/lib/stats_engine.R and
#          scripts/lib/expression_preprocessing.R on the cohort with affected samples
#          excluded from every layer. Results are compared with the Stage 09
#          and 10 primary outputs of the same run.
# Depends on: 03_harmonize_identifiers.R, 05_define_eligible_features.R,
#             07_target_universe_construction.R, 09_recurrence_analysis.R,
#             10_purity_adjusted_analysis.R
# ============================================================================

run_stage11 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "11_purity_conflict_full_exclusion_sensitivity")
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "stats_engine.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "expression_preprocessing.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "run_discovery_source.R"), local = TRUE)
  suppressPackageStartupMessages(library(data.table))
  script_path <- file.path(ctx$data_root, "scripts", "11_purity_conflict_full_exclusion_sensitivity.R")
  fwrite_fn <- function(x, p) fwrite(x, p, sep = "\t")
  min_n <- analysis_config$discovery$min_samples_per_cancer
  T <- analysis_config$thresholds$primary

  sample_annotation <- read_run_artifact(ctx, file.path("processed", "stage05", "sample_annotation_with_purity.tsv"), "05_define_eligible_features",
                                          function(p) read.delim(p, stringsAsFactors = FALSE))
  run_stage_assert(ctx, log_fn,
                    all(c("purity_cpe_forced_missing_conflict", "purity_lump_forced_missing_conflict") %in% names(sample_annotation)),
                    "sample_annotation_with_purity.tsv carries the Stage 05 purity-conflict flag columns this sensitivity analysis depends on")

  excluded_ids <- sample_annotation$sample_id[sample_annotation$purity_cpe_forced_missing_conflict |
                                                 sample_annotation$purity_lump_forced_missing_conflict]
  n_excluded <- length(excluded_ids)
  log_fn("Option-2 sensitivity: %d sample(s) with an unresolvable CPE/LUMP vial conflict (from THIS run's own Stage 05 census) excluded from EVERY layer for this comparison only -- not from the primary analysis.", n_excluded)

  mirna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mirna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mirna_mat_full <- as.matrix(mirna_matched[, -1, drop = FALSE]); rownames(mirna_mat_full) <- mirna_matched$mirna_id
  mrna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mrna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mrna_mat_full <- as.matrix(mrna_matched[, -1, drop = FALSE]); rownames(mrna_mat_full) <- toupper(trimws(mrna_matched$gene_id_raw))

  base_input_paths <- file.path(ctx$run_root, c(
    file.path("processed", "stage03", "mirna_matched.rds"), file.path("processed", "stage03", "mrna_matched.rds"),
    file.path("processed", "stage05", "sample_annotation_with_purity.tsv")
  ))

  if (n_excluded == 0) {
    log_fn("No samples have an unresolvable purity conflict in this run -- Option 1 and Option 2 are identical by construction. Writing a trivial, explicit comparison record rather than skipping this stage.")
    comparison <- data.table(comparison = "option1_vs_option2", n_excluded_samples = 0L,
                              note = "No purity duplicate-vial conflicts exist in this run's matched cohort; Option 1 and Option 2 are identical.")
    write_run_artifact(ctx, comparison, file.path("tables", "audit", "purity_sensitivity_option1_vs_option2_comparison.tsv"),
                        "11_purity_conflict_full_exclusion_sensitivity", script_path, base_input_paths, fwrite_fn)
    record_stage_status(ctx, "11_purity_conflict_full_exclusion_sensitivity", "COMPLETE")
    log_fn("Stage complete.")
    return(invisible(ctx))
  }

  keep_mirna <- !(colnames(mirna_mat_full) %in% excluded_ids)
  keep_mrna <- !(colnames(mrna_mat_full) %in% excluded_ids)
  mirna_mat <- mirna_mat_full[, keep_mirna, drop = FALSE]
  mrna_mat <- mrna_mat_full[, keep_mrna, drop = FALSE]
  sa2 <- sample_annotation[!(sample_annotation$sample_id %in% excluded_ids), , drop = FALSE]

  # --- Recompute eligible-cancer sets on the Option-2 (fully excluded) cohort ---
  cancer_counts2 <- as.data.table(sa2)[, .N, by = cancer]
  eligible_unadj2 <- sort(cancer_counts2[N >= min_n, cancer])
  cpe_counts2 <- as.data.table(sa2)[!is.na(CPE) & !is.na(CYT), .N, by = cancer]
  eligible_cpe2 <- sort(cpe_counts2[N >= min_n, cancer])
  lump_counts2 <- as.data.table(sa2)[!is.na(LUMP) & !is.na(CYT), .N, by = cancer]
  eligible_lump2 <- sort(lump_counts2[N >= min_n, cancer])

  eligible_unadj1 <- read_run_artifact(ctx, file.path("processed", "stage05", "eligible_cancers_unadjusted.tsv"), "05_define_eligible_features",
                                        function(p) read.delim(p, stringsAsFactors = FALSE))$cancer
  eligible_cpe1 <- read_run_artifact(ctx, file.path("processed", "stage05", "eligible_cancers_cpe.tsv"), "05_define_eligible_features",
                                      function(p) read.delim(p, stringsAsFactors = FALSE))$cancer
  eligible_lump1 <- read_run_artifact(ctx, file.path("processed", "stage05", "eligible_cancers_lump.tsv"), "05_define_eligible_features",
                                       function(p) read.delim(p, stringsAsFactors = FALSE))$cancer

  log_fn("Eligible-cancer-type counts: unadjusted Option1=%d Option2=%d; CPE Option1=%d Option2=%d; LUMP Option1=%d Option2=%d.",
         length(eligible_unadj1), length(eligible_unadj2), length(eligible_cpe1), length(eligible_cpe2), length(eligible_lump1), length(eligible_lump2))

  union_rel <- file.path("processed", "stage07", sprintf("pair_universe_threshold_%.2f.tsv", min(analysis_config$thresholds$all)))
  union_pairs <- read_run_artifact(ctx, union_rel, "07_target_universe_construction", fread)
  mirna_to_genes <- split(union_pairs$gene, union_pairs$mirna_id)
  needed_genes <- sort(unique(union_pairs$gene))
  pairs_T_rel <- file.path("processed", "stage07", sprintf("pair_universe_threshold_%.2f.tsv", T))
  pairs_T <- read_run_artifact(ctx, pairs_T_rel, "07_target_universe_construction", fread)
  genes_T <- unique(pairs_T$gene); mirnas_T <- unique(pairs_T$mirna_id)

  # --- Unadjusted three-arm test (Stage 06+08+09 pattern) on Option-2 cohort, primary threshold only ---
  mircyt_rows <- vector("list", length(eligible_unadj2))
  mg_rows <- vector("list", length(eligible_unadj2))
  gc_rows <- vector("list", length(eligible_unadj2))
  for (i in seq_along(eligible_unadj2)) {
    cc <- eligible_unadj2[i]
    idx <- which(sa2$cancer == cc)
    mc_res <- one_vs_many_spearman_pairwise(sa2$CYT[idx], mirna_mat[, idx, drop = FALSE], min_n = min_n)
    mircyt_rows[[i]] <- data.table(cancer = cc, mirna_id = rownames(mirna_mat), n_samples = mc_res$n, cor_mirna_cyt = mc_res$cor, p_mirna_cyt = mc_res$p)
    gc_res <- one_vs_many_spearman_pairwise(sa2$CYT[idx], mrna_mat[needed_genes, idx, drop = FALSE], min_n = min_n)
    gc_rows[[i]] <- data.table(gene = needed_genes, cancer = cc, n_samples = gc_res$n, cor_gene_cyt = gc_res$cor, p_gene_cyt = gc_res$p)
    mg_list <- vector("list", length(mirna_to_genes))
    for (j in seq_along(mirna_to_genes)) {
      mir <- names(mirna_to_genes)[j]; targ <- mirna_to_genes[[j]]
      mg_res <- one_vs_many_spearman_pairwise(mirna_mat[mir, idx], mrna_mat[targ, idx, drop = FALSE], min_n = min_n)
      mg_list[[j]] <- data.table(mirna_id = mir, gene = targ, cancer = cc, n_samples = mg_res$n, cor_mir_gene = mg_res$cor, p_mir_gene = mg_res$p)
    }
    mg_rows[[i]] <- rbindlist(mg_list)
  }
  master_mircyt2 <- rbindlist(mircyt_rows); master_mircyt2[, fdr_mirna_cyt := p.adjust(p_mirna_cyt, "BH"), by = cancer]
  master_mirgene2 <- rbindlist(mg_rows)
  master_genecyt2 <- rbindlist(gc_rows)

  mg2 <- master_mirgene2[paste(mirna_id, gene) %in% paste(pairs_T$mirna_id, pairs_T$gene)]
  mg2[, fdr_mir_gene := p.adjust(p_mir_gene, "BH"), by = cancer]
  gc2 <- master_genecyt2[gene %in% genes_T]
  gc2[, fdr_gene_cyt := p.adjust(p_gene_cyt, "BH"), by = cancer]
  mc2 <- master_mircyt2[mirna_id %in% mirnas_T]
  long_unadj2 <- merge(mg2, gc2, by = c("gene", "cancer"), all.x = TRUE)
  long_unadj2 <- merge(long_unadj2, mc2[, .(mirna_id, cancer, cor_mirna_cyt, fdr_mirna_cyt)], by = c("mirna_id", "cancer"), all.x = TRUE)
  long_unadj2[, three_arm_pass := (abs(cor_mirna_cyt) >= T & fdr_mirna_cyt < analysis_config$discovery$fdr_cutoff) &
                (abs(cor_mir_gene) >= T & fdr_mir_gene < analysis_config$discovery$fdr_cutoff) &
                (abs(cor_gene_cyt) >= T & fdr_gene_cyt < analysis_config$discovery$fdr_cutoff)]
  summary_unadj2 <- long_unadj2[, .(strict_recurrence = sum(three_arm_pass, na.rm = TRUE),
                                     supporting_cancers = paste(sort(unique(cancer[three_arm_pass %in% TRUE])), collapse = "; ")),
                                 by = .(mirna_id, gene)]
  summary_unadj2[, qualifies := strict_recurrence >= analysis_config$discovery$min_recurrence]
  write_run_artifact(ctx, summary_unadj2, file.path("tables", "audit", "purity_sensitivity_opt2_triplets_unadjusted.tsv"),
                      "11_purity_conflict_full_exclusion_sensitivity", script_path, base_input_paths, fwrite_fn)

  # --- CPE/LUMP-adjusted three-arm test (as in Stage 10, with per-cancer
  #     winsorization) on the Option-2 cohort, primary threshold only ---------
  run_metric2 <- function(metric_name, eligible_cancers2) {
    cancer_group <- sa2$cancer[match(colnames(mirna_mat), sa2$sample_id)]
    mirna_mat_w <- winsorize_matrix_rows_by_group(mirna_mat, cancer_group)
    mrna_cancer_group <- sa2$cancer[match(colnames(mrna_mat), sa2$sample_id)]
    mrna_mat_w <- winsorize_matrix_rows_by_group(mrna_mat, mrna_cancer_group)
    mc_full_rows <- vector("list", length(eligible_cancers2))
    mg_rows2 <- vector("list", length(eligible_cancers2))
    gc_rows2 <- vector("list", length(eligible_cancers2))
    all_mirnas <- rownames(mirna_mat_w)
    union_mirnas_here <- sort(unique(union_pairs$mirna_id))
    for (i in seq_along(eligible_cancers2)) {
      cc <- eligible_cancers2[i]
      idx <- which(sa2$cancer == cc & !is.na(sa2[[metric_name]]) & !is.na(sa2$CYT))
      purity_vec <- sa2[[metric_name]][idx]
      cyt_resid <- as.numeric(residualize_matrix(matrix(sa2$CYT[idx], nrow = 1), purity_vec))
      mir_resid_full <- residualize_matrix(mirna_mat_w[all_mirnas, idx, drop = FALSE], purity_vec)
      mc_res <- one_vs_many_spearman_pairwise(cyt_resid, mir_resid_full, min_n = min_n)
      mc_full_rows[[i]] <- data.table(mirna_id = all_mirnas, cancer = cc, n_samples = mc_res$n, cor_mirna_cyt_adj = mc_res$cor, p_mirna_cyt_adj = mc_res$p)
      mir_resid_union <- mir_resid_full[match(union_mirnas_here, all_mirnas), , drop = FALSE]
      gene_resid <- residualize_matrix(mrna_mat_w[needed_genes, idx, drop = FALSE], purity_vec)
      gc_res <- one_vs_many_spearman_pairwise(cyt_resid, gene_resid, min_n = min_n)
      gc_rows2[[i]] <- data.table(gene = needed_genes, cancer = cc, n_samples = gc_res$n, cor_gene_cyt_adj = gc_res$cor, p_gene_cyt_adj = gc_res$p)
      mg_list <- vector("list", length(mirna_to_genes))
      for (j in seq_along(mirna_to_genes)) {
        mir <- names(mirna_to_genes)[j]; targ <- mirna_to_genes[[j]]
        ridx <- match(targ, needed_genes); mridx <- match(mir, union_mirnas_here)
        mg_res <- one_vs_many_spearman_pairwise(mir_resid_union[mridx, ], gene_resid[ridx, , drop = FALSE], min_n = min_n)
        mg_list[[j]] <- data.table(mirna_id = mir, gene = targ, cancer = cc, n_samples = mg_res$n, cor_mir_gene_adj = mg_res$cor, p_mir_gene_adj = mg_res$p)
      }
      mg_rows2[[i]] <- rbindlist(mg_list)
    }
    mircyt_full2 <- rbindlist(mc_full_rows); mircyt_full2[, fdr_mirna_cyt_adj := p.adjust(p_mirna_cyt_adj, "BH"), by = cancer]
    list(mircyt_full = mircyt_full2, mirgene = rbindlist(mg_rows2), genecyt = rbindlist(gc_rows2))
  }

  build_adj_summary2 <- function(res, metric_name) {
    mg <- res$mirgene[paste(mirna_id, gene) %in% paste(pairs_T$mirna_id, pairs_T$gene)]
    mg[, fdr_mir_gene_adj := p.adjust(p_mir_gene_adj, "BH"), by = cancer]
    gc <- res$genecyt[gene %in% genes_T]
    gc[, fdr_gene_cyt_adj := p.adjust(p_gene_cyt_adj, "BH"), by = cancer]
    mc <- res$mircyt_full[mirna_id %in% mirnas_T]
    long <- merge(mg, gc, by = c("gene", "cancer"), all.x = TRUE)
    long <- merge(long, mc[, .(mirna_id, cancer, cor_mirna_cyt_adj, fdr_mirna_cyt_adj)], by = c("mirna_id", "cancer"), all.x = TRUE)
    long[, three_arm_pass := (abs(cor_mirna_cyt_adj) >= T & fdr_mirna_cyt_adj < analysis_config$discovery$fdr_cutoff) &
           (abs(cor_mir_gene_adj) >= T & fdr_mir_gene_adj < analysis_config$discovery$fdr_cutoff) &
           (abs(cor_gene_cyt_adj) >= T & fdr_gene_cyt_adj < analysis_config$discovery$fdr_cutoff)]
    summary_dt <- long[, .(strict_recurrence = sum(three_arm_pass, na.rm = TRUE),
                            supporting_cancers = paste(sort(unique(cancer[three_arm_pass %in% TRUE])), collapse = "; ")),
                        by = .(mirna_id, gene)]
    summary_dt[, qualifies := strict_recurrence >= analysis_config$discovery$min_recurrence]
    summary_dt
  }

  res_cpe2 <- run_metric2("CPE", eligible_cpe2)
  summary_cpe2 <- build_adj_summary2(res_cpe2, "CPE")
  write_run_artifact(ctx, summary_cpe2, file.path("tables", "audit", "purity_sensitivity_opt2_triplets_CPE.tsv"),
                      "11_purity_conflict_full_exclusion_sensitivity", script_path, base_input_paths, fwrite_fn)

  res_lump2 <- run_metric2("LUMP", eligible_lump2)
  summary_lump2 <- build_adj_summary2(res_lump2, "LUMP")
  write_run_artifact(ctx, summary_lump2, file.path("tables", "audit", "purity_sensitivity_opt2_triplets_LUMP.tsv"),
                      "11_purity_conflict_full_exclusion_sensitivity", script_path, base_input_paths, fwrite_fn)

  # --- Compare with the Option-1 (primary) output of the same run ---------
  opt1_cpe <- load_run_primary_triplets(ctx)  # qualifies == TRUE rows of the primary CPE table
  opt1_unadj_full <- read_run_artifact(ctx, file.path("tables", "primary", sprintf("triplets_rho%03d_unadjusted.tsv", round(T * 100))),
                                        "09_recurrence_analysis", fread)
  opt1_unadj <- opt1_unadj_full[qualifies == TRUE]
  opt2_cpe <- summary_cpe2[qualifies == TRUE]
  opt2_unadj <- summary_unadj2[qualifies == TRUE]

  key <- function(dt) if (nrow(dt) == 0) character(0) else paste(dt$mirna_id, dt$gene)
  gained_cpe <- setdiff(key(opt2_cpe), key(opt1_cpe))
  lost_cpe <- setdiff(key(opt1_cpe), key(opt2_cpe))
  gained_unadj <- setdiff(key(opt2_unadj), key(opt1_unadj))
  lost_unadj <- setdiff(key(opt1_unadj), key(opt2_unadj))

  comparison <- data.table(
    comparison = c("CPE_adjusted_primary_table", "unadjusted_table"),
    n_excluded_samples = n_excluded,
    n_option1_qualifying = c(nrow(opt1_cpe), nrow(opt1_unadj)),
    n_option2_qualifying = c(nrow(opt2_cpe), nrow(opt2_unadj)),
    n_gained_under_option2 = c(length(gained_cpe), length(gained_unadj)),
    n_lost_under_option2 = c(length(lost_cpe), length(lost_unadj)),
    gained_triplets = c(paste(gained_cpe, collapse = "; "), paste(gained_unadj, collapse = "; ")),
    lost_triplets = c(paste(lost_cpe, collapse = "; "), paste(lost_unadj, collapse = "; "))
  )
  write_run_artifact(ctx, comparison, file.path("tables", "audit", "purity_sensitivity_option1_vs_option2_comparison.tsv"),
                      "11_purity_conflict_full_exclusion_sensitivity", script_path, base_input_paths, fwrite_fn)

  # Per-triplet recurrence comparison for every Option-1 CPE-adjusted qualifying triplet (the reported headline set)
  opt1_cpe_full <- read_run_artifact(ctx, file.path("tables", "primary", run_primary_triplet_filename(ctx)), "10_purity_adjusted_analysis", fread)
  per_triplet <- merge(
    opt1_cpe_full[, .(mirna_id, gene, recurrence_option1 = strict_recurrence, qualifies_option1 = qualifies)],
    summary_cpe2[, .(mirna_id, gene, recurrence_option2 = strict_recurrence, qualifies_option2 = qualifies)],
    by = c("mirna_id", "gene"), all = TRUE
  )
  per_triplet[is.na(recurrence_option1), recurrence_option1 := 0L]
  per_triplet[is.na(recurrence_option2), recurrence_option2 := 0L]
  per_triplet[is.na(qualifies_option1), qualifies_option1 := FALSE]
  per_triplet[is.na(qualifies_option2), qualifies_option2 := FALSE]
  per_triplet[, recurrence_delta := recurrence_option2 - recurrence_option1]
  per_triplet[, qualification_changed := qualifies_option1 != qualifies_option2]
  per_triplet_changed <- per_triplet[recurrence_delta != 0 | qualification_changed == TRUE]
  write_run_artifact(ctx, per_triplet_changed, file.path("tables", "audit", "purity_sensitivity_option1_vs_option2_per_triplet_changes.tsv"),
                      "11_purity_conflict_full_exclusion_sensitivity", script_path, base_input_paths, fwrite_fn)

  log_fn("Option 1 vs Option 2 (CPE-adjusted primary table, threshold %.2f): Option1=%d qualifying, Option2=%d qualifying, %d gained, %d lost under Option2. %d triplet(s) show any recurrence-count or qualification change. Full detail: tables/audit/purity_sensitivity_option1_vs_option2_comparison.tsv, purity_sensitivity_option1_vs_option2_per_triplet_changes.tsv.",
         T, nrow(opt1_cpe), nrow(opt2_cpe), length(gained_cpe), length(lost_cpe), nrow(per_triplet_changed))

  report_lines <- c(
    "# Stage 11: Purity-Conflict Handling Sensitivity Report (Option 1 vs Option 2)",
    "",
    "**This is a sensitivity analysis only. It never replaces the primary (Option 1)",
    "analysis and is never used to select whichever result looks more favorable.**",
    "Option 1 (primary): CPE/LUMP values that conflict between duplicate vials are set",
    "to missing. Option 2 (this comparison): affected samples are excluded from every layer.",
    "",
    sprintf("- Samples with an unresolvable CPE/LUMP vial conflict in this run's matched cohort: %d", n_excluded),
    sprintf("- Eligible cancer types (unadjusted): Option 1 = %d, Option 2 (full exclusion) = %d", length(eligible_unadj1), length(eligible_unadj2)),
    sprintf("- Eligible cancer types (CPE-adjusted): Option 1 = %d, Option 2 = %d", length(eligible_cpe1), length(eligible_cpe2)),
    sprintf("- Eligible cancer types (LUMP-adjusted): Option 1 = %d, Option 2 = %d", length(eligible_lump1), length(eligible_lump2)),
    "",
    "## CPE-adjusted primary table (the table consumed by Stages 12 onward via load_run_primary_triplets())",
    sprintf("- Option 1 (primary) qualifying triplets: %d", nrow(opt1_cpe)),
    sprintf("- Option 2 (full-exclusion sensitivity) qualifying triplets: %d", nrow(opt2_cpe)),
    sprintf("- Gained under Option 2: %d (%s)", length(gained_cpe), if (length(gained_cpe) > 0) paste(gained_cpe, collapse = ", ") else "none"),
    sprintf("- Lost under Option 2: %d (%s)", length(lost_cpe), if (length(lost_cpe) > 0) paste(lost_cpe, collapse = ", ") else "none"),
    sprintf("- Does the CPE-adjusted qualifying triplet list change between Option 1 and Option 2? %s",
            if (length(gained_cpe) == 0 && length(lost_cpe) == 0) "NO -- identical qualifying set under both options." else "YES -- see gained/lost lists above."),
    "",
    "## Unadjusted table (diagnostic only; not the table downstream stages consume)",
    sprintf("- Option 1 qualifying triplets: %d", nrow(opt1_unadj)),
    sprintf("- Option 2 qualifying triplets: %d", nrow(opt2_unadj)),
    sprintf("- Gained under Option 2: %d", length(gained_unadj)),
    sprintf("- Lost under Option 2: %d", length(lost_unadj)),
    "",
    "## Downstream reach (Stages 12 onward)",
    "Stages 12 onward read the run's tables/primary/ CPE-adjusted output via",
    "load_run_primary_triplets(ctx).",
    if (length(gained_cpe) == 0 && length(lost_cpe) == 0)
      "Because the CPE-adjusted qualifying set is identical under Option 1 and Option 2, downstream stages would analyse the same (miRNA, gene) set under either option; any residual difference is confined to the affected samples within each stage's own purity-adjusted computations, not to which miRNAs and genes are analysed."
      else "The CPE-adjusted qualifying set differs between Option 1 and Option 2 (see gained/lost above); downstream stages would analyse a different (miRNA, gene) set depending on which option is primary. Only Option 1 is used as primary in this run.",
    "",
    sprintf("Per-triplet recurrence/qualification changes: %d of %d Option-1 CPE-adjusted qualifying triplets show any change under Option 2 (full detail: tables/audit/purity_sensitivity_option1_vs_option2_per_triplet_changes.tsv).",
            nrow(per_triplet_changed), nrow(opt1_cpe))
  )
  writeLines(report_lines, file.path(ctx$run_root, "reports", "27_purity_conflict_sensitivity_report.md"))
  log_fn("reports/27_purity_conflict_sensitivity_report.md written.")

  record_stage_status(ctx, "11_purity_conflict_full_exclusion_sensitivity", "COMPLETE")
  log_fn("Stage complete.")
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=11 or source scripts/run.R in RStudio.")
}
