# ============================================================================
# 08_same_cancer_three_arm_testing.R
#
# Purpose: compute miRNA-gene and gene-CYT Spearman correlations within each
#          eligible cancer type for the threshold-0.20 pair universe (the
#          loosest threshold, a superset of the 0.25 and 0.30 universes), and
#          combine them with the miRNA-CYT arm from Stage 06 to evaluate the
#          same-cancer three-arm condition independently for every
#          threshold.
# Inputs:  <run_root>/processed/stage03/mirna_matched.rds, mrna_matched.rds;
#          stage05/sample_annotation_with_purity.tsv, eligible_cancers_unadjusted.tsv;
#          stage06/master_mirna_cyt.tsv;
#          stage07/pair_universe_threshold_*.tsv (all same-run)
# Outputs: <run_root>/processed/stage08/master_mirna_gene.tsv,
#          master_gene_cyt.tsv (threshold-invariant, union universe);
#          three_arm_threshold_0.30.tsv, 0.25.tsv, 0.20.tsv
#          (one row per miRNA x gene x cancer, unadjusted)
# Method:  correlations are computed with scripts/lib/stats_engine.R::
#          one_vs_many_spearman_pairwise(). The matched mRNA matrix contains
#          native missing values for a minority of gene/cancer combinations,
#          so both variables of every miRNA-gene and gene-CYT test are
#          re-ranked within that test's complete-case subset (vectorized for
#          fully observed genes, per-gene otherwise). One call per
#          (miRNA, cancer) correlates the miRNA with all of its target genes,
#          and one call per cancer correlates CYT with every unique union
#          target gene, so a gene shared by several miRNAs is tested once per
#          cancer. Benjamini-Hochberg FDR for the miRNA-gene and gene-CYT
#          arms is computed independently for each threshold within each
#          cancer type, across that threshold's tested pairs (miRNA-gene) or
#          unique tested genes (gene-CYT); these are the prespecified
#          hypothesis families for the two arms. The miRNA-CYT FDR is taken
#          from Stage 06's master_mirna_cyt.tsv, computed across all eligible
#          miRNAs, and is not recomputed on a threshold-restricted subset.
#          The three-arm pass indicator is the logical AND of all three arms
#          evaluated in the same cancer.
# Depends on: 03_harmonize_identifiers.R, 05_define_eligible_features.R,
#             06_mirna_cyt_discovery.R, 07_target_universe_construction.R
# ============================================================================

run_stage08 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "08_same_cancer_three_arm_testing")
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "stats_engine.R"), local = TRUE)
  suppressPackageStartupMessages(library(data.table))
  script_path <- file.path(ctx$data_root, "scripts", "08_same_cancer_three_arm_testing.R")

  mirna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mirna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mirna_mat <- as.matrix(mirna_matched[, -1, drop = FALSE]); rownames(mirna_mat) <- mirna_matched$mirna_id
  mrna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mrna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mrna_mat <- as.matrix(mrna_matched[, -1, drop = FALSE]); rownames(mrna_mat) <- toupper(trimws(mrna_matched$gene_id_raw))

  sample_annotation <- read_run_artifact(ctx, file.path("processed", "stage05", "sample_annotation_with_purity.tsv"), "05_define_eligible_features",
                                          function(p) read.delim(p, stringsAsFactors = FALSE))
  eligible_cancers <- read_run_artifact(ctx, file.path("processed", "stage05", "eligible_cancers_unadjusted.tsv"), "05_define_eligible_features",
                                         function(p) read.delim(p, stringsAsFactors = FALSE))$cancer
  master_mircyt <- read_run_artifact(ctx, file.path("processed", "stage06", "master_mirna_cyt.tsv"), "06_mirna_cyt_discovery", fread)

  loosest_threshold <- min(analysis_config$thresholds$all)
  union_rel <- file.path("processed", "stage07", sprintf("pair_universe_threshold_%.2f.tsv", loosest_threshold))
  union_pairs <- read_run_artifact(ctx, union_rel, "07_target_universe_construction", fread)
  log_fn("Union pair universe (threshold %.2f, used once for the shared master computation): %d pairs, %d unique miRNAs, %d unique genes",
         loosest_threshold, nrow(union_pairs), uniqueN(union_pairs$mirna_id), uniqueN(union_pairs$gene))

  base_input_paths <- file.path(ctx$run_root, c(
    file.path("processed", "stage03", "mirna_matched.rds"), file.path("processed", "stage03", "mrna_matched.rds"),
    file.path("processed", "stage05", "sample_annotation_with_purity.tsv"), file.path("processed", "stage05", "eligible_cancers_unadjusted.tsv"),
    file.path("processed", "stage06", "master_mirna_cyt.tsv"), union_rel
  ))

  mirna_to_genes <- split(union_pairs$gene, union_pairs$mirna_id)
  needed_genes <- sort(unique(union_pairs$gene))

  log_fn("Computing master miRNA-gene and gene-CYT correlation tables (unadjusted) across %d cancers ...", length(eligible_cancers))
  mg_rows <- vector("list", length(eligible_cancers))
  gc_rows <- vector("list", length(eligible_cancers))
  for (i in seq_along(eligible_cancers)) {
    cc <- eligible_cancers[i]
    idx <- which(sample_annotation$cancer == cc)
    # Pairwise-complete correlation: the matched mRNA matrix contains
    # native TCGA missing values, so the strict
    # rank_rows()/safe_rank()/one_vs_many_spearman() path (which requires
    # fully finite input) cannot be applied to `mrna_mat` directly.
    # one_vs_many_spearman_pairwise() uses the vectorized path for genes
    # that are fully finite in this cancer and a per-gene pairwise-complete
    # computation for the remainder.
    cyt_here <- sample_annotation$CYT[idx]

    gc_res <- one_vs_many_spearman_pairwise(cyt_here, mrna_mat[needed_genes, idx, drop = FALSE],
                                             min_n = analysis_config$discovery$min_samples_per_cancer)
    gc_rows[[i]] <- data.table(gene = needed_genes, cancer = cc, n_samples = gc_res$n,
                                cor_gene_cyt = gc_res$cor, p_gene_cyt = gc_res$p)

    mg_list <- vector("list", length(mirna_to_genes))
    for (j in seq_along(mirna_to_genes)) {
      mir <- names(mirna_to_genes)[j]; targ <- mirna_to_genes[[j]]
      mg_res <- one_vs_many_spearman_pairwise(mirna_mat[mir, idx], mrna_mat[targ, idx, drop = FALSE],
                                               min_n = analysis_config$discovery$min_samples_per_cancer)
      mg_list[[j]] <- data.table(mirna_id = mir, gene = targ, cancer = cc, n_samples = mg_res$n,
                                  cor_mir_gene = mg_res$cor, p_mir_gene = mg_res$p)
    }
    mg_rows[[i]] <- rbindlist(mg_list)
    log_fn("  Cancer %s done (%d/%d).", cc, i, length(eligible_cancers))
  }
  master_mirgene <- rbindlist(mg_rows)
  master_genecyt <- rbindlist(gc_rows)

  run_stage_assert(ctx, log_fn, nrow(master_mirgene) == nrow(union_pairs) * length(eligible_cancers),
                    "Master miRNA-gene table has the expected row count (%s)",
                    format(nrow(union_pairs) * length(eligible_cancers), big.mark = ","))
  run_stage_assert(ctx, log_fn, nrow(master_genecyt) == length(needed_genes) * length(eligible_cancers),
                    "Master gene-CYT table has the expected row count (%s)",
                    format(length(needed_genes) * length(eligible_cancers), big.mark = ","))

  fwrite_fn <- function(x, p) fwrite(x, p, sep = "\t")
  write_run_artifact(ctx, master_mirgene, file.path("processed", "stage08", "master_mirna_gene.tsv"),
                      "08_same_cancer_three_arm_testing", script_path, base_input_paths, fwrite_fn)
  write_run_artifact(ctx, master_genecyt, file.path("processed", "stage08", "master_gene_cyt.tsv"),
                      "08_same_cancer_three_arm_testing", script_path, base_input_paths, fwrite_fn)

  for (T in analysis_config$thresholds$all) {
    pairs_T_rel <- file.path("processed", "stage07", sprintf("pair_universe_threshold_%.2f.tsv", T))
    pairs_T <- read_run_artifact(ctx, pairs_T_rel, "07_target_universe_construction", fread)
    genes_T <- unique(pairs_T$gene); mirnas_T <- unique(pairs_T$mirna_id)

    mg <- master_mirgene[paste(mirna_id, gene) %in% paste(pairs_T$mirna_id, pairs_T$gene)]
    mg[, fdr_mir_gene := p.adjust(p_mir_gene, "BH"), by = cancer]
    gc <- master_genecyt[gene %in% genes_T]
    gc[, fdr_gene_cyt := p.adjust(p_gene_cyt, "BH"), by = cancer]
    # miRNA-CYT: restrict to the threshold-T rows for the join, reusing Stage 06's
    # fdr_mirna_cyt (full eligible-miRNA family) rather than recomputing BH on
    # the restricted subset.
    mc <- master_mircyt[mirna_id %in% mirnas_T]

    long <- merge(mg, gc, by = c("gene", "cancer"), all.x = TRUE)
    long <- merge(long, mc[, .(mirna_id, cancer, cor_mirna_cyt, fdr_mirna_cyt)], by = c("mirna_id", "cancer"), all.x = TRUE)

    long[, pass_mirna_cyt := abs(cor_mirna_cyt) >= T & fdr_mirna_cyt < analysis_config$discovery$fdr_cutoff]
    long[, pass_mirna_gene := abs(cor_mir_gene) >= T & fdr_mir_gene < analysis_config$discovery$fdr_cutoff]
    long[, pass_gene_cyt := abs(cor_gene_cyt) >= T & fdr_gene_cyt < analysis_config$discovery$fdr_cutoff]
    long[, three_arm_pass := pass_mirna_cyt & pass_mirna_gene & pass_gene_cyt]
    long[, sign_mirna_cyt := sign(cor_mirna_cyt)]
    long[, sign_mirna_gene := sign(cor_mir_gene)]
    long[, sign_gene_cyt := sign(cor_gene_cyt)]
    long[, threshold := T]
    long[, adjustment := "unadjusted"]

    dup <- long[, .N, by = .(mirna_id, gene, cancer)]
    run_stage_assert(ctx, log_fn, all(dup$N == 1), "[threshold %.2f] No duplicate miRNA-gene-cancer rows", T)
    run_stage_assert(ctx, log_fn, all(long$cor_mir_gene[!is.na(long$cor_mir_gene)] >= -1 & long$cor_mir_gene[!is.na(long$cor_mir_gene)] <= 1) &&
                        all(long$cor_gene_cyt[!is.na(long$cor_gene_cyt)] >= -1 & long$cor_gene_cyt[!is.na(long$cor_gene_cyt)] <= 1),
                      "[threshold %.2f] All correlations lie in [-1,1]", T)

    fname <- sprintf("three_arm_threshold_%.2f.tsv", T)
    write_run_artifact(ctx, long, file.path("processed", "stage08", fname),
                        "08_same_cancer_three_arm_testing", script_path,
                        c(base_input_paths, file.path(ctx$run_root, pairs_T_rel)), fwrite_fn)
    log_fn("Threshold %.2f: three-arm long table written (%s rows, %d three-arm-passing rows).",
           T, format(nrow(long), big.mark = ","), sum(long$three_arm_pass, na.rm = TRUE))
  }

  record_stage_status(ctx, "08_same_cancer_three_arm_testing", "COMPLETE")
  log_fn("Stage complete.")
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=8 or source scripts/run.R in RStudio.")
}
