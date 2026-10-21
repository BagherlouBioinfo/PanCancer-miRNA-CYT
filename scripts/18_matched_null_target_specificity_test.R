# ============================================================================
# 19_matched_null_target_specificity_test.R  (Stage 19)
#
# Purpose: prespecified matched-null permutation test of miRNA-target pair
#          specificity: whether the real target gene of each primary triplet
#          of the current run shows more same-cancer three-arm recurrence
#          than comparably matched genes that are not predicted targets of
#          the same miRNA. This distinguishes pair-specific evidence from a
#          general immune-associated co-expression signal. The complete
#          design (matched null genes, number of draws, empirical-p formula,
#          FDR family, exclusion of predicted targets and random seed) is
#          specified in this header.
#
#          The test evaluates exactly the current run's primary set (read
#          through load_run_primary_triplets(ctx) and read_run_artifact())
#          and never adds to, removes from or re-ranks it; the result is
#          reported as a descriptive specificity signal alongside primary
#          membership. No triplet, miRNA or gene count is hard-coded.
#
#
# Prespecification:
#   - Permutation unit: one primary triplet (real miRNA + real target gene)
#     -> replaced by (same real miRNA + one matched null gene), independently
#     for each primary triplet. The miRNA is never permuted, only the gene
#     -- this isolates gene-side target specificity.
#   - Null-candidate exclusions: (a) the real target gene; (b) every gene
#     predicted as a TargetScan target of that miRNA at the most permissive
#     threshold (processed/stage07/pair_universe_threshold_0.20.tsv), so no
#     predicted target enters a miRNA's null pool.
#   - Matching variables (5, capturing expression level, variability,
#     CYT association, and immune-population association):
#       1. mean log2 expression (pan-cohort, over the matched mRNA matrix)
#       2. variance of log2 expression (pan-cohort)
#       3. gene-CYT association magnitude: median |cor_gene_cyt_adj| (CPE-
#          adjusted, from processed/stage10/master_CPE.rds$genecyt) across
#          the gene's own >=30-sample cancers
#       4. immune-cell association 1: pan-cohort Spearman |rho| between gene
#          expression and the MCPcounter "T cells" score (Stage 12 output)
#       5. immune-cell association 2: the same against the MCPcounter
#          "Monocytic lineage" score (a lymphoid and a myeloid population,
#          used as separate matching variables and never combined)
#   - Cancer-availability hard filter: a null gene is eligible for a given
#     real triplet only if it has a valid (n_samples>=30) CPE-adjusted
#     gene-CYT estimate in every cancer where the real target gene does.
#   - Matching procedure: for each triplet, the five matching variables are
#     standardized (z-scores) across the eligible-gene universe; the
#     Euclidean distance from the real gene to every eligible,
#     non-predicted-target candidate is computed; a caliper of the 200
#     nearest candidates (or all candidates if fewer pass the filters) is
#     formed; and up to 99 genes are drawn from it without replacement.
#     Caliper-based matched sampling follows Rosenbaum & Rubin (1985, Am
#     Stat 39(1):33-38), applied directly to the five matching covariates.
#   - Number of draws: up to 99 matched null genes per triplet (fewer only
#     if the caliper contains fewer candidates, which is reported).
#   - Random seed: analysis_config$seeds$permutation, the project's single
#     standing permutation seed, applied once per triplet in triplet order
#     so the draw is reproducible run-to-run.
#   - Per-cancer pass rule (identical for real and null pairs, using raw
#     nominal significance rather than a full-universe FDR family per null
#     draw): |rho| >= 0.30 and raw p < 0.05 in all three arms (miRNA-CYT,
#     miRNA-gene, gene-CYT). miRNA-CYT and gene-CYT reuse the CPE-adjusted
#     correlations in processed/stage10/master_CPE.rds (gene-CYT covers the
#     full eligible-gene universe independently of any miRNA). miRNA-gene
#     was precomputed by Stage 10 only for the TargetScan-predicted pair
#     universe, from which null genes are excluded by design, so the
#     (real miRNA, null gene) arm is computed here with the same method as
#     Stage 10: single-covariate CPE residualization
#     (scripts/lib/stats_engine.R::residualize_matrix()) of both expression vectors
#     followed by pairwise-complete Spearman correlation of the residuals
#     (scripts/lib/stats_engine.R::one_vs_many_spearman_pairwise()).
#     nominal_recurrence = number of the real gene's eligible cancers
#     meeting this rule.
#   - Empirical p-value (fixed formula, North et al. 2002 continuity
#     correction): p = (1 + #{null nominal_recurrence >= real nominal_recurrence}) /
#     (1 + #{null genes drawn}).
#   - Two recurrence statistics are distinguished. "nominal_recurrence" is
#     the raw-p statistic defined above and the only statistic used in the
#     empirical-p test, for real and null genes alike. The output also
#     reports "fdr_qualified_recurrence", based on Stage 10's per-cancer
#     BH-FDR < 0.05 definition (the definition behind strict_recurrence and
#     qualifies in the primary triplet table), for cross-reference only; it
#     never enters the empirical-p formula or null-gene scoring.
#   - FDR family: Benjamini-Hochberg across all primary-triplet empirical
#     p-values of this test, a single family distinct from the discovery
#     families and the Stage 12 per-population sensitivity families.
#   - Decision rule: a triplet outperforms its matched null pool if its
#     empirical FDR < 0.05. The result is reported as the number and
#     identity of triplets meeting this rule.
#
# Inputs:  <run_root>/tables/primary/<primary triplet file> (via
#          load_run_primary_triplets(), same-run); processed/stage10/
#          master_CPE.rds; processed/stage07/pair_universe_threshold_0.20.tsv;
#          processed/stage03/mrna_matched.rds, mirna_matched.rds;
#          processed/stage05/sample_annotation_with_purity.tsv;
#          processed/stage11_immune/sample_level_immune_estimates.tsv
#          (all same-run, provenance-verified).
# Outputs: <run_root>/tables/audit/matched_null_target_specificity_results.tsv
#          (one row per triplet: observed nominal_recurrence, the
#          reporting-only fdr_qualified_recurrence, null-distribution
#          summary, empirical p, FDR and verdict);
#          tables/audit/matched_null_target_specificity_null_draws.tsv
#          (one row per triplet x null gene drawn); the stage log
#          logs/19_matched_null_target_specificity.log.
# Depends on (same-run, enforced by read_run_artifact()): 03, 05, 07, 09,
#          10, 12.
# ============================================================================

run_stage19 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "19_matched_null_target_specificity")
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "run_discovery_source.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "expression_preprocessing.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "stats_engine.R"), local = TRUE)
  suppressPackageStartupMessages(library(data.table))
  script_path <- file.path(ctx$data_root, "scripts", "18_matched_null_target_specificity_test.R")
  fwrite_fn <- function(x, p) fwrite(x, p, sep = "\t")

  min_n <- analysis_config$discovery$min_samples_per_cancer      # 30
  rho_thr <- analysis_config$thresholds$primary                  # 0.30
  p_thr <- 0.05
  n_null_ceiling <- 99
  caliper_size <- 200
  seed_base <- analysis_config$seeds$permutation

  triplets <- load_run_primary_triplets(ctx)
  run_stage_assert(ctx, log_fn, all(c("mirna_id", "gene", "strict_recurrence") %in% names(triplets)),
                    "This run's primary triplet table has mirna_id, gene, and strict_recurrence columns (%d triplets loaded)", nrow(triplets))

  master_cpe <- read_run_artifact(ctx, file.path("processed", "stage10", "master_CPE.rds"), "10_purity_adjusted_analysis", readRDS)
  mirgene <- copy(master_cpe$mirgene); setkey(mirgene, mirna_id, gene, cancer)
  genecyt <- copy(master_cpe$genecyt); setkey(genecyt, gene, cancer)
  mircyt <- copy(master_cpe$mircyt); setkey(mircyt, mirna_id, cancer)

  eligible_genes_universe <- sort(unique(genecyt$gene))

  # --- matching variable 1-2: pan-cohort mean/variance of log2 expression ---
  mrna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mrna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mrna_mat_full <- as.matrix(mrna_matched[, -1, drop = FALSE])
  rownames(mrna_mat_full) <- toupper(trimws(mrna_matched$gene_id_raw))
  mrna_mat_full <- deduplicate_matrix_rows(mrna_mat_full)
  mrna_mat <- mrna_mat_full[rownames(mrna_mat_full) %in% eligible_genes_universe, , drop = FALSE]

  mirna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mirna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mirna_mat_full <- as.matrix(mirna_matched[, -1, drop = FALSE])
  rownames(mirna_mat_full) <- mirna_matched$mirna_id

  # Per-feature winsorization before the OLS purity residualization in
  # evaluate_null_genes_batch() below, as in Stage 10. The matching
  # variables above use the unwinsorized matrices so that mean and variance
  # matching reflects each gene's actual expression distribution; the
  # winsorized copy is therefore created only after that computation.
  mirna_mat_full <- winsorize_matrix_rows(mirna_mat_full)
  mrna_mat_full <- winsorize_matrix_rows(mrna_mat_full)

  sample_annotation <- read_run_artifact(ctx, file.path("processed", "stage05", "sample_annotation_with_purity.tsv"), "05_define_eligible_features",
                                          function(p) read.delim(p, stringsAsFactors = FALSE))

  gene_stats <- data.table(
    gene = rownames(mrna_mat),
    mean_expr = rowMeans(mrna_mat, na.rm = TRUE),
    var_expr = apply(mrna_mat, 1, var, na.rm = TRUE)
  )

  # --- matching variable 3: median |gene-CYT rho| across the gene's own eligible cancers ---
  gene_cyt_summary <- genecyt[n_samples >= min_n, .(gene_cyt_medianabs = median(abs(cor_gene_cyt_adj), na.rm = TRUE)), by = gene]
  gene_stats <- merge(gene_stats, gene_cyt_summary, by = "gene", all.x = TRUE)

  # --- matching variables 4-5: pan-cohort |rho| against two individual, standard,
  #     already-computed MCPcounter populations (never combined into one index) ---
  immune_est <- read_run_artifact(ctx, file.path("processed", "stage11_immune", "sample_level_immune_estimates.tsv"), "12_immune_composition_analysis", fread)
  mrna_samples <- colnames(mrna_mat)
  tcell_v <- immune_est$`T cells`[match(mrna_samples, immune_est$sample_id)]
  mono_v <- immune_est$`Monocytic lineage`[match(mrna_samples, immune_est$sample_id)]
  ok_immune <- is.finite(tcell_v) & is.finite(mono_v)
  tcell_rank <- rank(tcell_v[ok_immune]); mono_rank <- rank(mono_v[ok_immune])
  gene_ranks <- t(apply(mrna_mat[, ok_immune, drop = FALSE], 1, rank))   # computed once, reused for both immune matching variables
  gene_immune <- data.table(
    gene = rownames(mrna_mat),
    tcell_assoc_abs = abs(as.numeric(cor(t(gene_ranks), tcell_rank))),
    mono_assoc_abs = abs(as.numeric(cor(t(gene_ranks), mono_rank)))
  )
  gene_stats <- merge(gene_stats, gene_immune, by = "gene", all.x = TRUE)

  gene_stats <- gene_stats[complete.cases(gene_stats)]
  match_vars <- c("mean_expr", "var_expr", "gene_cyt_medianabs", "tcell_assoc_abs", "mono_assoc_abs")
  gene_stats_z <- copy(gene_stats)
  for (v in match_vars) gene_stats_z[[v]] <- as.numeric(scale(gene_stats_z[[v]]))
  setkey(gene_stats_z, gene)

  log_fn("Matched-null target-specificity test: %d primary triplets (this run's own), up to %d null genes each, caliper=%d, seed_base=%d.",
         nrow(triplets), n_null_ceiling, caliper_size, seed_base)

  # --- predicted-target exclusion universe (most permissive threshold) ---
  pred_targets <- read_run_artifact(ctx, file.path("processed", "stage07", "pair_universe_threshold_0.20.tsv"), "07_target_universe_construction", fread)
  setkey(pred_targets, mirna_id)

  # Real pairs belong to the TargetScan predicted-target universe evaluated
  # by Stage 10, so their miRNA-gene arm is read from master_CPE.rds.
  # pass_rule_real() wraps scripts/lib/stats_engine.R::three_arm_raw_p_pass_count()
  # (tested in tests/test_matched_null_helpers.R).
  pass_rule_real <- function(mir, gene_x, cancers) {
    three_arm_raw_p_pass_count(mir, gene_x, cancers, mirgene, genecyt, mircyt, rho_thr, p_thr)
  }

  # Null genes lie outside the predicted-target universe, so master_CPE.rds
  # has no miRNA-gene value for them; this function computes that arm per
  # cancer for a batch of null genes with the same closed-form CPE
  # residualization and Spearman method as
  # scripts/10_purity_adjusted_analysis.R. gene-CYT and miRNA-CYT are read
  # from master_CPE.rds, which covers the full eligible gene and miRNA
  # universes.
  evaluate_null_genes_batch <- function(mir, null_genes, cancers) {
    mc <- mircyt[.(mir, cancers), on = c("mirna_id", "cancer")]
    mc_pass <- !is.na(mc$p_mirna_cyt_adj) & abs(mc$cor_mirna_cyt_adj) >= rho_thr & mc$p_mirna_cyt_adj < p_thr

    all_pass <- matrix(FALSE, nrow = length(null_genes), ncol = length(cancers))
    for (ci in seq_along(cancers)) {
      cc <- cancers[ci]
      gc <- genecyt[.(null_genes, cc), on = c("gene", "cancer")]
      pass_gc <- !is.na(gc$p_gene_cyt_adj) & abs(gc$cor_gene_cyt_adj) >= rho_thr & gc$p_gene_cyt_adj < p_thr

      idx <- which(sample_annotation$cancer == cc & !is.na(sample_annotation$CPE) & !is.na(sample_annotation$CYT))
      purity_vec <- sample_annotation$CPE[idx]
      # Pairwise-complete computation: mrna_mat_full may contain native
      # missing values for some genes and cancers. residualize_matrix()
      # leaves a per-row pairwise-complete NA pattern, and
      # one_vs_many_spearman_pairwise() uses the vectorized path for fully
      # finite rows and a per-row computation otherwise.
      mir_resid <- residualize_matrix(mirna_mat_full[mir, idx, drop = FALSE], purity_vec)
      gene_resid <- residualize_matrix(mrna_mat_full[null_genes, idx, drop = FALSE], purity_vec)
      mg_res <- one_vs_many_spearman_pairwise(mir_resid[1, ], gene_resid, min_n = min_n)
      pass_mg <- !is.na(mg_res$p) & abs(mg_res$cor) >= rho_thr & mg_res$p < p_thr

      all_pass[, ci] <- mc_pass[ci] & pass_gc & pass_mg
    }
    rowSums(all_pass, na.rm = TRUE)
  }

  summary_rows <- list()
  draw_rows <- list()

  for (i in seq_len(nrow(triplets))) {
    mir <- triplets$mirna_id[i]; real_gene <- triplets$gene[i]
    set.seed(seed_base + i)   # per-triplet seed, deterministic in triplet order

    real_cancers <- sort(genecyt[gene == real_gene & n_samples >= min_n]$cancer)
    if (length(real_cancers) == 0) {
      log_fn("  [%d/%d] %s-%s: real gene has no eligible cancer at n>=%d; skipped.", i, nrow(triplets), mir, real_gene, min_n)
      next
    }
    nominal_recurrence <- pass_rule_real(mir, real_gene, real_cancers)

    # fdr_qualified_recurrence is reported for cross-reference only and is
    # not used in the empirical-p test. It is taken from the current run's
    # Stage 10 primary triplet table (triplets$strict_recurrence, row i;
    # `triplets` is iterated in the order loaded), computed under Stage 10's
    # per-cancer BH-FDR < 0.05 rule. Null genes are scored only with the
    # raw-p rule, so a comparison based on FDR-qualified recurrence would
    # require a full BH-FDR family for every null draw; comparing an
    # FDR-qualified real statistic with a raw-p null statistic would be
    # inconsistent.
    fdr_qualified_recurrence <- triplets$strict_recurrence[i]

    mir_predicted_targets <- pred_targets[.(mir), on = "mirna_id"]$gene
    cancer_ok_genes <- genecyt[n_samples >= min_n, .(n_covered = sum(cancer %in% real_cancers)), by = gene][n_covered == length(real_cancers)]$gene

    candidate_pool <- setdiff(cancer_ok_genes, c(real_gene, mir_predicted_targets))
    candidate_pool <- intersect(candidate_pool, gene_stats_z$gene)

    if (length(candidate_pool) == 0) {
      log_fn("  [%d/%d] %s-%s: empty matched-null candidate pool after exclusions; skipped.", i, nrow(triplets), mir, real_gene)
      next
    }

    # scripts/lib/stats_engine.R::build_caliper() (tested in
    # tests/test_matched_null_helpers.R).
    caliper_genes <- build_caliper(real_gene, candidate_pool, gene_stats_z, match_vars, caliper_size)

    n_draw <- min(n_null_ceiling, length(caliper_genes))
    null_genes <- sample(caliper_genes, n_draw, replace = FALSE)

    null_recurrence <- evaluate_null_genes_batch(mir, null_genes, real_cancers)

    # Empirical p-value, (1 + sum(null_dist >= observed)) / (1 + length(null_dist)),
    # computed with the shared scripts/lib/stats_engine.R::empirical_p_upper();
    # length(null_recurrence) equals n_draw by construction.
    emp_p <- empirical_p_upper(nominal_recurrence, null_recurrence)

    summary_rows[[length(summary_rows) + 1]] <- data.table(
      mirna_id = mir, gene = real_gene,
      n_eligible_cancers = length(real_cancers), nominal_recurrence = nominal_recurrence,
      fdr_qualified_recurrence = fdr_qualified_recurrence,
      n_null_drawn = n_draw, caliper_pool_size = length(candidate_pool),
      null_recurrence_mean = mean(null_recurrence), null_recurrence_median = median(null_recurrence),
      null_recurrence_max = max(null_recurrence), null_recurrence_min = min(null_recurrence),
      empirical_p = emp_p
    )
    draw_rows[[length(draw_rows) + 1]] <- data.table(
      mirna_id = mir, real_gene = real_gene, null_gene = null_genes, null_recurrence = null_recurrence
    )

    if (i %% 5 == 0) log_fn("  Processed %d/%d triplets.", i, nrow(triplets))
  }

  run_stage_assert(ctx, log_fn, length(summary_rows) > 0,
    "At least one primary triplet was evaluable by the matched-null test (0 evaluable would indicate an eligibility/data problem upstream, not a legitimate empty result)")

  summary_dt <- rbindlist(summary_rows)
  summary_dt[, fdr := p.adjust(empirical_p, "BH")]
  summary_dt[, outperforms_matched_null := fdr < 0.05]
  draws_dt <- rbindlist(draw_rows)

  input_paths <- file.path(ctx$run_root, c(
    file.path("processed", "stage10", "master_CPE.rds"),
    file.path("processed", "stage07", "pair_universe_threshold_0.20.tsv"),
    file.path("processed", "stage03", "mrna_matched.rds"),
    file.path("processed", "stage03", "mirna_matched.rds"),
    file.path("processed", "stage05", "sample_annotation_with_purity.tsv"),
    file.path("processed", "stage11_immune", "sample_level_immune_estimates.tsv")
  ))
  write_run_artifact(ctx, summary_dt, file.path("tables", "audit", "matched_null_target_specificity_results.tsv"),
                      "19_matched_null_target_specificity", script_path, input_paths, fwrite_fn)
  write_run_artifact(ctx, draws_dt, file.path("tables", "audit", "matched_null_target_specificity_null_draws.tsv"),
                      "19_matched_null_target_specificity", script_path, input_paths, fwrite_fn)

  log_fn("DONE. %d/%d evaluable triplets outperform their matched null pool at FDR<0.05.",
         sum(summary_dt$outperforms_matched_null), nrow(summary_dt))

  record_stage_status(ctx, "19_matched_null_target_specificity", "COMPLETE")
  invisible(summary_dt)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=18 or source scripts/run.R in RStudio.")
}
