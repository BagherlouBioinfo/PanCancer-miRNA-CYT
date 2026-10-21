# ============================================================================
# 10_purity_adjusted_analysis.R
#
# Purpose: repeat the same-cancer three-arm test on tumor-purity-
#          residualized expression and CYT values, using the primary
#          purity metric (CPE) consistently for all three arms, and
#          separately using the sensitivity metric (LUMP) with strict
#          labeling so the two are never combined.
# Inputs:  <run_root>/processed/stage03/mirna_matched.rds, mrna_matched.rds;
#          stage05/sample_annotation_with_purity.tsv;
#          stage07/pair_universe_threshold_*.tsv (all same-run)
# Outputs: <run_root>/tables/primary/triplets_rho030_CPE.tsv;
#          tables/sensitivity/triplets_rho0{25,20}_CPE.tsv,
#          triplets_rho0{30,25,20}_LUMP.tsv;
#          tables/audit/purity_concordance_summary.tsv
# Method:  for each CPE-eligible (or LUMP-eligible) cancer type, miRNA
#          expression, target-gene expression and CYT are each residualized
#          on the purity value by closed-form OLS
#          (scripts/lib/stats_engine.R::residualize_matrix(), validated against
#          per-row lm() fits at run time; pairwise-complete for rows with
#          missing values), and pairwise-complete Spearman correlations
#          (scripts/lib/stats_engine.R::one_vs_many_spearman_pairwise()) are computed
#          on the residuals. All three arms use the same purity metric
#          within an adjusted analysis; CPE- and LUMP-adjusted results are
#          stored and reported in separate files and never merged.
#
#          The adjusted miRNA-CYT correlation and FDR are computed once per
#          purity metric across every miRNA eligible for that metric (not
#          only the TargetScan-candidate pool), mirroring the full-family
#          design of Stage 06. This table (res$mircyt_full) is reused by
#          every threshold's three-arm test and never recomputed on a
#          threshold-restricted subset. The miRNA-gene and gene-CYT families
#          are scoped to each threshold's pair universe.
# Depends on: 03_harmonize_identifiers.R, 05_define_eligible_features.R,
#             07_target_universe_construction.R
# ============================================================================

run_stage10 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "10_purity_adjusted_analysis")
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "stats_engine.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "expression_preprocessing.R"), local = TRUE)
  suppressPackageStartupMessages(library(data.table))
  script_path <- file.path(ctx$data_root, "scripts", "10_purity_adjusted_analysis.R")
  fwrite_fn <- function(x, p) fwrite(x, p, sep = "\t")

  mirna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mirna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mirna_mat <- as.matrix(mirna_matched[, -1, drop = FALSE]); rownames(mirna_mat) <- mirna_matched$mirna_id
  mrna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mrna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mrna_mat <- as.matrix(mrna_matched[, -1, drop = FALSE]); rownames(mrna_mat) <- toupper(trimws(mrna_matched$gene_id_raw))
  sample_annotation <- read_run_artifact(ctx, file.path("processed", "stage05", "sample_annotation_with_purity.tsv"), "05_define_eligible_features",
                                          function(p) read.delim(p, stringsAsFactors = FALSE))

  # Per-feature winsorization within each cancer type, immediately before
  # OLS purity residualization. The rank-based Spearman arms of Stages 06
  # and 08 use the unwinsorized matrices, because clipping is unnecessary
  # for rank-based methods and can introduce ties, whereas the OLS fit in
  # residualize_matrix() is sensitive to the leverage of extreme samples.
  # Clipping points are computed independently per cancer type
  # (winsorize_matrix_rows_by_group(), scripts/lib/expression_preprocessing.R) with
  # probabilities 0.001/0.999 and winsorize_vector()'s NA handling, so no
  # cancer's clip point is informed by another cancer's distribution.
  cancer_group <- sample_annotation$cancer[match(colnames(mirna_mat), sample_annotation$sample_id)]
  run_stage_assert(ctx, log_fn, !anyNA(cancer_group),
                    "Every column of the matched miRNA matrix has a known cancer-type assignment in sample_annotation_with_purity.tsv (required for per-cancer winsorization)")
  mrna_cancer_group <- sample_annotation$cancer[match(colnames(mrna_mat), sample_annotation$sample_id)]
  run_stage_assert(ctx, log_fn, !anyNA(mrna_cancer_group),
                    "Every column of the matched mRNA matrix has a known cancer-type assignment in sample_annotation_with_purity.tsv (required for per-cancer winsorization)")
  mirna_mat <- winsorize_matrix_rows_by_group(mirna_mat, cancer_group)
  mrna_mat <- winsorize_matrix_rows_by_group(mrna_mat, mrna_cancer_group)

  loosest_threshold <- min(analysis_config$thresholds$all)
  union_rel <- file.path("processed", "stage07", sprintf("pair_universe_threshold_%.2f.tsv", loosest_threshold))
  union_pairs <- read_run_artifact(ctx, union_rel, "07_target_universe_construction", fread)
  mirna_to_genes <- split(union_pairs$gene, union_pairs$mirna_id)
  needed_genes <- sort(unique(union_pairs$gene))
  all_mirnas <- rownames(mirna_mat)  # full eligible-miRNA universe, as in Stage 06

  # Eligible-cancer sets per metric are read from Stage 05's
  # provenance-verified output, so that the eligibility rule is defined in
  # one place.
  eligible_cancers_by_metric <- list(
    CPE = read_run_artifact(ctx, file.path("processed", "stage05", "eligible_cancers_cpe.tsv"), "05_define_eligible_features",
                             function(p) read.delim(p, stringsAsFactors = FALSE)$cancer),
    LUMP = read_run_artifact(ctx, file.path("processed", "stage05", "eligible_cancers_lump.tsv"), "05_define_eligible_features",
                              function(p) read.delim(p, stringsAsFactors = FALSE)$cancer)
  )

  base_input_paths <- file.path(ctx$run_root, c(
    file.path("processed", "stage03", "mirna_matched.rds"), file.path("processed", "stage03", "mrna_matched.rds"),
    file.path("processed", "stage05", "sample_annotation_with_purity.tsv"),
    file.path("processed", "stage05", "eligible_cancers_cpe.tsv"),
    file.path("processed", "stage05", "eligible_cancers_lump.tsv"), union_rel
  ))

  # Validate the Spearman-correlation engine and the OLS-residualization
  # engine against their reference implementations (cor.test() and lm()) on
  # real data before any purity-adjusted computation, as Stage 06 does for
  # the correlation engine.
  validation_cancer <- eligible_cancers_by_metric$CPE[1]
  validation_idx <- which(sample_annotation$cancer == validation_cancer &
                             !is.na(sample_annotation$CPE) & !is.na(sample_annotation$CYT))
  validation <- validate_stats_engine(mirna_mat[, validation_idx, drop = FALSE],
                                       sample_annotation$CYT[validation_idx],
                                       covariate = sample_annotation$CPE[validation_idx],
                                       seed = analysis_config$seeds$master_seed)
  run_stage_assert(ctx, log_fn,
                    validation$max_cor_diff < 1e-9 && validation$max_p_diff < 1e-6 && validation$max_resid_diff < 1e-8,
                    "Vectorized Spearman correlation AND OLS residualization match their reference (cor.test()/lm()) implementations on real data (max cor diff=%.2e, max p diff=%.2e, max resid diff=%.2e)",
                    validation$max_cor_diff, validation$max_p_diff, validation$max_resid_diff)

  run_metric <- function(metric_name) {
    eligible_cancers <- sort(eligible_cancers_by_metric[[metric_name]])
    # Independent recomputation of eligibility from the sample annotation,
    # asserted to agree with Stage 05's output.
    recomputed <- sort(as.data.table(sample_annotation)[!is.na(get(metric_name)) & !is.na(CYT), .N, by = cancer][
      N >= analysis_config$discovery$min_samples_per_cancer, cancer])
    run_stage_assert(ctx, log_fn, setequal(eligible_cancers, recomputed),
                      "[%s] Eligible-cancer set read from Stage 05's output matches an independent recomputation from sample_annotation_with_purity.tsv", metric_name)
    log_fn("[%s] %d eligible cancer types.", metric_name, length(eligible_cancers))

    mc_full_rows <- vector("list", length(eligible_cancers))
    mg_rows <- vector("list", length(eligible_cancers))
    gc_rows <- vector("list", length(eligible_cancers))
    for (i in seq_along(eligible_cancers)) {
      cc <- eligible_cancers[i]
      idx <- which(sample_annotation$cancer == cc & !is.na(sample_annotation[[metric_name]]) & !is.na(sample_annotation$CYT))
      purity_vec <- sample_annotation[[metric_name]][idx]

      # residualize_matrix() leaves a per-row pairwise-complete NA pattern,
      # and every ranking and correlation step below uses
      # one_vs_many_spearman_pairwise() (vectorized for fully finite rows,
      # pairwise-complete otherwise) rather than the strict
      # rank_rows()/safe_rank()/one_vs_many_spearman() path.
      cyt_resid <- as.numeric(residualize_matrix(matrix(sample_annotation$CYT[idx], nrow = 1), purity_vec))

      # Full-family miRNA-CYT arm: every eligible miRNA, not only the TargetScan-candidate pool.
      mir_resid_full <- residualize_matrix(mirna_mat[all_mirnas, idx, drop = FALSE], purity_vec)
      mc_full_res <- one_vs_many_spearman_pairwise(cyt_resid, mir_resid_full,
                                                    min_n = analysis_config$discovery$min_samples_per_cancer)
      mc_full_rows[[i]] <- data.table(mirna_id = all_mirnas, cancer = cc, n_samples = mc_full_res$n,
                                       cor_mirna_cyt_adj = mc_full_res$cor, p_mirna_cyt_adj = mc_full_res$p)

      # miRNA-gene and gene-CYT arms are restricted to the union TargetScan
      # candidate/target pool, their prespecified hypothesis family.
      mir_resid_union <- mir_resid_full[match(sort(unique(union_pairs$mirna_id)), all_mirnas), , drop = FALSE]
      union_mirnas_here <- sort(unique(union_pairs$mirna_id))

      gene_resid <- residualize_matrix(mrna_mat[needed_genes, idx, drop = FALSE], purity_vec)
      gc_res <- one_vs_many_spearman_pairwise(cyt_resid, gene_resid,
                                               min_n = analysis_config$discovery$min_samples_per_cancer)
      gc_rows[[i]] <- data.table(gene = needed_genes, cancer = cc, n_samples = gc_res$n,
                                  cor_gene_cyt_adj = gc_res$cor, p_gene_cyt_adj = gc_res$p)

      mg_list <- vector("list", length(mirna_to_genes))
      for (j in seq_along(mirna_to_genes)) {
        mir <- names(mirna_to_genes)[j]; targ <- mirna_to_genes[[j]]
        ridx <- match(targ, needed_genes); mridx <- match(mir, union_mirnas_here)
        mg_res <- one_vs_many_spearman_pairwise(mir_resid_union[mridx, ], gene_resid[ridx, , drop = FALSE],
                                                 min_n = analysis_config$discovery$min_samples_per_cancer)
        mg_list[[j]] <- data.table(mirna_id = mir, gene = targ, cancer = cc, n_samples = mg_res$n,
                                    cor_mir_gene_adj = mg_res$cor, p_mir_gene_adj = mg_res$p)
      }
      mg_rows[[i]] <- rbindlist(mg_list)
      log_fn("  [%s] Cancer %s done (%d/%d).", metric_name, cc, i, length(eligible_cancers))
    }
    mircyt_full <- rbindlist(mc_full_rows)
    mircyt_full[, fdr_mirna_cyt_adj := p.adjust(p_mirna_cyt_adj, "BH"), by = cancer]
    list(mircyt_full = mircyt_full, mirgene = rbindlist(mg_rows), genecyt = rbindlist(gc_rows), eligible = eligible_cancers)
  }

  metric_results <- list()
  # Qualifying sets are captured in memory for the concordance check below.
  # They are not re-read from disk because read_run_artifact() requires the
  # generating stage (this one) to be recorded COMPLETE first.
  primary_cpe_qualifying <- NULL
  primary_lump_qualifying <- NULL
  for (metric_name in c(analysis_config$purity$primary_metric, analysis_config$purity$sensitivity_metric)) {
    log_fn("=== Running purity-adjusted analysis for metric: %s ===", metric_name)
    res <- run_metric(metric_name)
    metric_results[[metric_name]] <- res
    write_run_artifact(ctx, res, file.path("processed", "stage10", sprintf("master_%s.rds", metric_name)),
                        "10_purity_adjusted_analysis", script_path, base_input_paths, saveRDS)

    run_stage_assert(ctx, log_fn, !any(grepl(setdiff(c("CPE", "LUMP"), metric_name), names(res$mircyt_full), ignore.case = TRUE)),
                      "[%s] No column from the other purity metric is present in this metric's tables", metric_name)

    for (T in analysis_config$thresholds$all) {
      pairs_T_rel <- file.path("processed", "stage07", sprintf("pair_universe_threshold_%.2f.tsv", T))
      pairs_T <- read_run_artifact(ctx, pairs_T_rel, "07_target_universe_construction", fread)
      genes_T <- unique(pairs_T$gene); mirnas_T <- unique(pairs_T$mirna_id)

      mg <- res$mirgene[paste(mirna_id, gene) %in% paste(pairs_T$mirna_id, pairs_T$gene)]
      mg[, fdr_mir_gene_adj := p.adjust(p_mir_gene_adj, "BH"), by = cancer]
      gc <- res$genecyt[gene %in% genes_T]
      gc[, fdr_gene_cyt_adj := p.adjust(p_gene_cyt_adj, "BH"), by = cancer]
      # Restrict rows for the join, reusing the full-family fdr_mirna_cyt_adj
      # computed above rather than recomputing BH on this subset.
      mc <- res$mircyt_full[mirna_id %in% mirnas_T]

      long <- merge(mg, gc, by = c("gene", "cancer"), all.x = TRUE)
      long <- merge(long, mc[, .(mirna_id, cancer, cor_mirna_cyt_adj, fdr_mirna_cyt_adj)], by = c("mirna_id", "cancer"), all.x = TRUE)
      long[, pass_mirna_cyt := abs(cor_mirna_cyt_adj) >= T & fdr_mirna_cyt_adj < analysis_config$discovery$fdr_cutoff]
      long[, pass_mirna_gene := abs(cor_mir_gene_adj) >= T & fdr_mir_gene_adj < analysis_config$discovery$fdr_cutoff]
      long[, pass_gene_cyt := abs(cor_gene_cyt_adj) >= T & fdr_gene_cyt_adj < analysis_config$discovery$fdr_cutoff]
      long[, three_arm_pass := pass_mirna_cyt & pass_mirna_gene & pass_gene_cyt]
      long[, purity_setting := metric_name]; long[, threshold := T]; long[, run_id := ctx$run_id]

      summary_dt <- long[, .(
        strict_recurrence = sum(three_arm_pass, na.rm = TRUE),
        supporting_cancers = paste(sort(unique(cancer[three_arm_pass %in% TRUE])), collapse = "; ")
      ), by = .(mirna_id, gene)]
      summary_dt[, qualifies := strict_recurrence >= analysis_config$discovery$min_recurrence]
      summary_dt$purity_setting <- metric_name; summary_dt$threshold <- T; summary_dt$run_id <- ctx$run_id
      summary_dt$generating_script <- "10_purity_adjusted_analysis.R"

      thr_tag <- sprintf("rho%03d", round(T * 100))
      out_table <- sprintf("triplets_%s_%s.tsv", thr_tag, metric_name)
      out_subdir <- if (T == analysis_config$thresholds$primary && metric_name == analysis_config$purity$primary_metric) "primary" else "sensitivity"
      stage_inputs <- c(base_input_paths, file.path(ctx$run_root, pairs_T_rel))
      write_run_artifact(ctx, summary_dt, file.path("tables", out_subdir, out_table),
                          "10_purity_adjusted_analysis", script_path, stage_inputs, fwrite_fn)
      write_run_artifact(ctx, long, file.path("tables", "audit", sprintf("cancer_level_evidence_%s_%s.tsv", metric_name, thr_tag)),
                          "10_purity_adjusted_analysis", script_path, stage_inputs, fwrite_fn)
      log_fn("[%s] Threshold %.2f: %d qualifying triplets written to tables/%s/%s",
             metric_name, T, sum(summary_dt$qualifies), out_subdir, out_table)

      if (T == analysis_config$thresholds$primary && metric_name == analysis_config$purity$primary_metric) {
        primary_cpe_qualifying <- summary_dt[qualifies == TRUE]
      }
      if (T == analysis_config$thresholds$primary && metric_name == analysis_config$purity$sensitivity_metric) {
        primary_lump_qualifying <- summary_dt[qualifies == TRUE]
      }
    }
  }

  concordance <- data.table(
    metric_pair = "CPE_vs_LUMP", threshold = analysis_config$thresholds$primary,
    n_cpe_qualifying = nrow(primary_cpe_qualifying), n_lump_qualifying = nrow(primary_lump_qualifying),
    n_overlap = sum(paste(primary_cpe_qualifying$mirna_id, primary_cpe_qualifying$gene) %in%
                      paste(primary_lump_qualifying$mirna_id, primary_lump_qualifying$gene))
  )
  # Derived from analysis_config, so that this declared input path always
  # matches the file written by the primary-threshold iteration above.
  primary_thr_tag <- sprintf("rho%03d", round(analysis_config$thresholds$primary * 100))
  cpe_subdir <- if (analysis_config$purity$primary_metric == "CPE") "primary" else "sensitivity"
  lump_subdir <- if (analysis_config$purity$primary_metric == "LUMP") "primary" else "sensitivity"
  write_run_artifact(ctx, concordance, file.path("tables", "audit", "purity_concordance_summary.tsv"),
                      "10_purity_adjusted_analysis", script_path,
                      file.path(ctx$run_root, c(file.path("tables", cpe_subdir, sprintf("triplets_%s_CPE.tsv", primary_thr_tag)),
                                                 file.path("tables", lump_subdir, sprintf("triplets_%s_LUMP.tsv", primary_thr_tag)))),
                      fwrite_fn)

  record_stage_status(ctx, "10_purity_adjusted_analysis", "COMPLETE")
  log_fn("Stage complete.")
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=10 or source scripts/run.R in RStudio.")
}
