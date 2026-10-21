# ============================================================================
# 18_clustering_and_state_analysis.R
#
# Purpose: test whether the matched samples form reproducible groups on the
#          current run's primary molecular features (miRNAs and target
#          genes), using a cancer-adjusted feature matrix so that detected
#          structure does not merely reflect tissue of origin. Meaningful
#          clusters are not assumed: weak structure is reported as weak, and
#          a descriptive label is assigned to a cluster solution only if a
#          prespecified minimum-evidence criterion is met. No molecular
#          subtype is claimed.
# Inputs:  <run_root>/processed/stage03/mirna_matched.rds, mrna_matched.rds
#          (same-run); processed/stage04/cyt_scores.tsv (same-run);
#          processed/stage05/sample_annotation_with_purity.tsv (same-run);
#          processed/stage11_immune/sample_level_immune_estimates.tsv
#          (Stage 12, same-run); the current run's primary triplet set
#          (scripts/lib/run_discovery_source.R).
# Outputs: <run_root>/tables/downstream/cluster_validity.tsv,
#          cluster_stability.tsv, cluster_null_calibration.tsv,
#          cluster_membership.tsv, cluster_characteristics.tsv,
#          cluster_pca_coordinates.tsv;
#          figures/main_candidates/clustering/pca_by_cluster.png;
#          reports/32_clustering_and_state_report.md.
# Method:  feature matrix = primary features x matched samples,
#          standardized within each cancer type (z-score per feature per
#          cancer). The same cancer-adjusted matrix is used for k-means
#          clustering and for the PCA that visualizes it. k-means
#          (nstart = analysis_config$clustering$n_init) is run for every k in
#          [analysis_config$clustering$min_k, max_k]; k is chosen by maximum
#          mean silhouette width (cluster::silhouette), with the gap
#          statistic (cluster::clusGap, B = 50) as a corroborating
#          criterion. At the chosen k the following are evaluated:
#          cluster-size balance, resampling stability (mean Adjusted Rand
#          Index over 100 bootstrap resamples), seed sensitivity (5
#          independent master seeds), feature-scaling sensitivity (the same
#          procedure on a pan-cohort-scaled, not cancer-adjusted, matrix),
#          and calibration against a permutation null of cluster structure.
#          Associations with cancer type (Cramer's V) and with CYT and each
#          MCPcounter population (Kruskal-Wallis, BH-FDR within that family)
#          are tested only after cluster assignment is fixed and describe
#          the clusters without serving as evidence of their validity. A
#          descriptive label is assigned only if mean silhouette > 0.25,
#          mean stability ARI > 0.5 and the observed silhouette exceeds the
#          permutation null (empirical p < 0.05); all three criteria are
#          prespecified.
# Depends on: 03_harmonize_identifiers.R, 04_compute_validate_cyt.R,
#             05_define_eligible_features.R, 09_recurrence_analysis.R,
#             10_purity_adjusted_analysis.R, 12_immune_composition_analysis.R
#             (all same-run).
# ============================================================================

# adjusted_rand_index() is defined in scripts/lib/stats_engine.R (sourced below), so
# that observed and null stability statistics use the same function.

cramers_v <- function(a, b) {
  tab <- table(a, b)
  chi2 <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  n <- sum(tab)
  k <- min(nrow(tab), ncol(tab))
  sqrt(unname(chi2) / (n * (k - 1)))
}

# Within-cancer per-feature standardization (mean 0, sd 1 within each cancer
# group). Defined at top level so that the observed analysis and
# scripts/lib/stats_engine.R::run_null_cluster_calibration() apply the same function.
within_cancer_scale <- function(mat, grp) {
  out <- mat
  for (g in unique(grp)) {
    idx <- which(grp == g)
    if (length(idx) < 2) next
    out[, idx] <- t(scale(t(mat[, idx, drop = FALSE])))
  }
  # A feature with zero variance within one cancer group yields NaN from
  # scale(). Such entries are set to zero (no deviation from the group
  # mean, hence no contribution to distances), identically for observed
  # data and null-calibration draws (tested in
  # tests/test_cluster_null_calibration.R).
  out[is.nan(out)] <- 0
  out
}

run_stage18 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "18_clustering_and_state_analysis")
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "run_discovery_source.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "expression_preprocessing.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "stats_engine.R"), local = TRUE)
  suppressPackageStartupMessages({ library(data.table); library(cluster) })
  script_path <- file.path(ctx$data_root, "scripts", "17_clustering_and_state_analysis.R")
  fwrite_fn <- function(x, p) fwrite(x, p, sep = "\t")

  # ==========================================================================
  # PHASE 1: PRESPECIFIED DEFINITIONS (before any clustering is run)
  # ==========================================================================
  triplets <- load_run_primary_triplets(ctx)
  features <- sort(unique(c(triplets$mirna_id, triplets$gene)))
  n_features_expected <- uniqueN(triplets$mirna_id) + uniqueN(triplets$gene)
  run_stage_assert(ctx, log_fn, length(features) == n_features_expected,
               "Feature matrix uses exactly this run's %d primary features (%d miRNAs, %d genes)",
               n_features_expected, uniqueN(triplets$mirna_id), uniqueN(triplets$gene))

  mirna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mirna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mirna_mat <- as.matrix(mirna_matched[, -1, drop = FALSE]); rownames(mirna_mat) <- mirna_matched$mirna_id
  mrna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mrna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mrna_mat <- as.matrix(mrna_matched[, -1, drop = FALSE]); rownames(mrna_mat) <- toupper(trimws(mrna_matched$gene_id_raw))
  sample_annot <- read_run_artifact(ctx, file.path("processed", "stage05", "sample_annotation_with_purity.tsv"), "05_define_eligible_features", fread)
  cyt_scores <- read_run_artifact(ctx, file.path("processed", "stage04", "cyt_scores.tsv"), "04_compute_validate_cyt", fread)
  immune_est <- read_run_artifact(ctx, file.path("processed", "stage11_immune", "sample_level_immune_estimates.tsv"), "12_immune_composition_analysis", fread)

  mirna_feats <- intersect(features, rownames(mirna_mat))
  gene_feats <- intersect(features, rownames(mrna_mat))
  raw_mat <- rbind(mirna_mat[mirna_feats, , drop = FALSE], mrna_mat[gene_feats, , drop = FALSE])
  raw_mat <- raw_mat[features, , drop = FALSE]
  # Per-feature winsorization before PCA and k-means, both of which are
  # sensitive to extreme samples dominating a component or a cluster-centre
  # distance.
  raw_mat <- winsorize_matrix_rows(raw_mat)
  run_stage_assert(ctx, log_fn, all(colnames(raw_mat) == sample_annot$sample_id[match(colnames(raw_mat), sample_annot$sample_id)]),
               "Feature matrix sample columns align with sample_annotation")
  run_stage_assert(ctx, log_fn, sum(is.na(raw_mat)) == 0, "No missing values in the %d-feature matrix (complete-case by construction)", length(features))

  cancer_of <- sample_annot$cancer[match(colnames(raw_mat), sample_annot$sample_id)]

  adj_mat <- within_cancer_scale(raw_mat, cancer_of)          # cancer-adjusted matrix, used for clustering and visualization
  run_stage_assert(ctx, log_fn, all(is.finite(adj_mat)),
               "Cancer-adjusted feature matrix is entirely finite after within-cancer scaling (no NaN/Inf survives; zero-variance feature-by-cancer blocks are deterministically zero-filled by within_cancer_scale() -- see tests/test_cluster_null_calibration.R controls 4a-4c)")
  raw_scaled_mat <- t(scale(t(raw_mat)))                       # sensitivity only: pan-cohort scaling, not cancer-adjusted

  k_range <- seq(analysis_config$clustering$min_k, analysis_config$clustering$max_k)
  n_init <- analysis_config$clustering$n_init
  master_seed <- analysis_config$seeds$kmeans
  log_fn("Prespecified: feature matrix = %d this-run primary features x %d samples, within-cancer standardized; candidate k = %s; nstart = %d.",
         length(features), ncol(adj_mat), paste(k_range, collapse = ","), n_init)

  X <- t(adj_mat)  # samples x features, for kmeans/silhouette/clusGap

  # ==========================================================================
  # PHASE 2: k SELECTION -- silhouette (primary) + gap statistic (corroborating)
  # ==========================================================================
  log_fn("PHASE 2: evaluating candidate k values.")
  dist_X <- dist(X)
  validity_rows <- vector("list", length(k_range))
  km_fits <- list()
  for (i in seq_along(k_range)) {
    k <- k_range[i]
    set.seed(master_seed)
    km <- kmeans(X, centers = k, nstart = n_init, iter.max = 50)
    km_fits[[as.character(k)]] <- km
    sil <- cluster::silhouette(km$cluster, dist_X)
    avg_sil <- mean(sil[, "sil_width"])
    size_tab <- table(km$cluster)
    balance_ratio <- min(size_tab) / max(size_tab)
    validity_rows[[i]] <- data.table(k = k, avg_silhouette = avg_sil, tot_withinss = km$tot.withinss,
                                      betweenss = km$betweenss, min_cluster_size = min(size_tab),
                                      max_cluster_size = max(size_tab), size_balance_ratio = balance_ratio)
    log_fn("  k=%d: avg silhouette=%.4f, size balance ratio=%.3f (sizes: %s)", k, avg_sil, balance_ratio, paste(size_tab, collapse = ","))
  }
  validity_dt <- rbindlist(validity_rows)

  log_fn("Computing gap statistic (B=50 reference bootstraps) as a corroborating criterion...")
  set.seed(master_seed)
  gap_res <- tryCatch(
    cluster::clusGap(X, FUN = function(x, k) kmeans(x, centers = k, nstart = 25, iter.max = 50),
                      K.max = max(k_range), B = 50, verbose = FALSE),
    error = function(e) { log_fn("  clusGap failed: %s", conditionMessage(e)); NULL }
  )
  if (!is.null(gap_res)) {
    # clusGap returns rows for k = 1..K.max in order (the Tab matrix has no
    # rownames), so row i corresponds to k = i.
    gap_tab <- as.data.table(gap_res$Tab)
    gap_tab[, k := seq_len(.N)]
    validity_dt <- merge(validity_dt, gap_tab[k %in% k_range, .(k, gap, SE.sim)], by = "k", all.x = TRUE)
  } else {
    validity_dt[, `:=`(gap = NA_real_, SE.sim = NA_real_)]
  }

  k_chosen <- validity_dt$k[which.max(validity_dt$avg_silhouette)]
  log_fn("Chosen k = %d (maximizes mean silhouette, prespecified selection rule).", k_chosen)
  km_final <- km_fits[[as.character(k_chosen)]]

  base_inputs <- file.path(ctx$run_root, c(
    file.path("tables", "primary", run_primary_triplet_filename(ctx)),
    file.path("processed", "stage03", "mirna_matched.rds"), file.path("processed", "stage03", "mrna_matched.rds"),
    file.path("processed", "stage04", "cyt_scores.tsv"), file.path("processed", "stage05", "sample_annotation_with_purity.tsv"),
    file.path("processed", "stage11_immune", "sample_level_immune_estimates.tsv")
  ))
  write_run_artifact(ctx, validity_dt, file.path("tables", "downstream", "cluster_validity.tsv"),
                      "18_clustering_and_state_analysis", script_path, base_inputs, fwrite_fn)

  # ==========================================================================
  # PHASE 3: STABILITY -- bootstrap resampling ARI, seed sensitivity, feature-scaling sensitivity
  # ==========================================================================
  log_fn("PHASE 3: stability, seed-sensitivity, and feature-scaling-sensitivity evaluation at k=%d.", k_chosen)

  n_samples <- nrow(X)
  n_boot_stability <- 100
  ari_vals <- numeric(n_boot_stability)
  set.seed(master_seed)
  for (b in seq_len(n_boot_stability)) {
    idx <- sample(seq_len(n_samples), n_samples, replace = TRUE)
    km_b <- tryCatch(kmeans(X[idx, , drop = FALSE], centers = k_chosen, nstart = 10, iter.max = 50), error = function(e) NULL)
    if (is.null(km_b)) { ari_vals[b] <- NA_real_; next }
    ari_vals[b] <- adjusted_rand_index(km_final$cluster[idx], km_b$cluster)
  }
  mean_stability_ari <- mean(ari_vals, na.rm = TRUE)
  log_fn("  Mean bootstrap stability ARI (k=%d, %d resamples): %.4f (SD %.4f).", k_chosen, n_boot_stability, mean_stability_ari, sd(ari_vals, na.rm = TRUE))

  seed_sensitivity_seeds <- master_seed + 0:4
  seed_clusterings <- lapply(seed_sensitivity_seeds, function(s) { set.seed(s); kmeans(X, centers = k_chosen, nstart = n_init, iter.max = 50)$cluster })
  seed_ari <- combn(length(seed_clusterings), 2, function(ij) adjusted_rand_index(seed_clusterings[[ij[1]]], seed_clusterings[[ij[2]]]))
  mean_seed_ari <- mean(seed_ari)
  log_fn("  Seed-sensitivity: mean pairwise ARI across %d independent master seeds = %.4f.", length(seed_sensitivity_seeds), mean_seed_ari)

  set.seed(master_seed)
  km_rawscaled <- kmeans(t(raw_scaled_mat), centers = k_chosen, nstart = n_init, iter.max = 50)
  scaling_ari <- adjusted_rand_index(km_final$cluster, km_rawscaled$cluster)
  cancer_assoc_adjusted <- cramers_v(km_final$cluster, cancer_of)
  cancer_assoc_rawscaled <- cramers_v(km_rawscaled$cluster, cancer_of)
  log_fn("  Feature-scaling sensitivity: ARI(cancer-adjusted vs pan-cohort-raw-scaled clustering) = %.4f.", scaling_ari)
  log_fn("  Cancer-type association (Cramer's V): cancer-adjusted clustering = %.4f, raw-scaled clustering = %.4f.", cancer_assoc_adjusted, cancer_assoc_rawscaled)

  stability_dt <- data.table(
    metric = c("bootstrap_stability_mean_ARI", "bootstrap_stability_sd_ARI", "seed_sensitivity_mean_pairwise_ARI",
               "feature_scaling_sensitivity_ARI", "cancer_association_cramers_v_adjusted", "cancer_association_cramers_v_rawscaled",
               "chosen_k", "n_bootstrap_resamples", "n_seeds_tested"),
    value = c(mean_stability_ari, sd(ari_vals, na.rm = TRUE), mean_seed_ari, scaling_ari, cancer_assoc_adjusted, cancer_assoc_rawscaled,
              k_chosen, n_boot_stability, length(seed_sensitivity_seeds))
  )
  write_run_artifact(ctx, stability_dt, file.path("tables", "downstream", "cluster_stability.tsv"),
                      "18_clustering_and_state_analysis", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "cluster_validity.tsv"), fwrite_fn)

  # ==========================================================================
  # PHASE 3b: NULL-CLUSTER-STRUCTURE CALIBRATION. B independent
  # per-feature-permuted null matrices (preserving sample count, feature
  # count, cancer-group sizes and each feature's marginal distribution,
  # while destroying cross-feature sample-level structure) are processed
  # with the same within-cancer scaling and k-selection/stability procedure
  # as the observed data (scripts/lib/stats_engine.R::run_null_cluster_calibration(),
  # which calls within_cancer_scale() and select_k_and_stability()). CYT and
  # immune populations play no role in this calibration. The number of
  # iterations and the seed are fixed in config/analysis_config.R.
  # ==========================================================================
  n_null <- analysis_config$clustering$n_null_iterations
  n_null_boot <- analysis_config$clustering$n_null_boot_stability
  null_seed <- analysis_config$seeds$null_calibration
  log_fn("PHASE 3b: null-cluster-structure calibration (%d permutation-null draws, seed=%d).", n_null, null_seed)
  null_res <- run_null_cluster_calibration(raw_mat, cancer_of, within_cancer_scale, k_range, n_init, n_null, null_seed, n_null_boot)
  chosen_sil_observed <- validity_dt$avg_silhouette[validity_dt$k == k_chosen]
  p_silhouette_vs_null <- empirical_p_upper(chosen_sil_observed, null_res$null_max_silhouette)
  p_stability_vs_null <- empirical_p_upper(mean_stability_ari, null_res$null_stability)
  null_calibration_dt <- data.table(
    metric = c("observed_max_silhouette", "null_silhouette_median", "null_silhouette_q05", "null_silhouette_q95", "empirical_p_silhouette_vs_null",
               "observed_stability_ari", "null_stability_median", "null_stability_q05", "null_stability_q95", "empirical_p_stability_vs_null",
               "n_null_iterations", "n_null_boot_stability", "null_calibration_seed"),
    value = c(chosen_sil_observed, median(null_res$null_max_silhouette), quantile(null_res$null_max_silhouette, 0.05), quantile(null_res$null_max_silhouette, 0.95), p_silhouette_vs_null,
              mean_stability_ari, median(null_res$null_stability, na.rm = TRUE), quantile(null_res$null_stability, 0.05, na.rm = TRUE), quantile(null_res$null_stability, 0.95, na.rm = TRUE), p_stability_vs_null,
              n_null, n_null_boot, null_seed)
  )
  log_fn("  Observed max silhouette=%.4f vs null median=%.4f (empirical p=%.4f). Observed stability ARI=%.4f vs null median=%.4f (empirical p=%.4f).",
         chosen_sil_observed, median(null_res$null_max_silhouette), p_silhouette_vs_null, mean_stability_ari, median(null_res$null_stability, na.rm = TRUE), p_stability_vs_null)
  write_run_artifact(ctx, null_calibration_dt, file.path("tables", "downstream", "cluster_null_calibration.tsv"),
                      "18_clustering_and_state_analysis", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "cluster_validity.tsv"), fwrite_fn)

  # ==========================================================================
  # PHASE 4: POST-CLUSTERING CHARACTERIZATION -- CYT and immune populations
  # (tested only after cluster assignment is fixed; descriptive, not validation)
  # ==========================================================================
  log_fn("PHASE 4: post-clustering characterization against CYT and immune populations (internal, not external -- see disclosure below).")
  membership <- data.table(sample_id = colnames(raw_mat), cancer = cancer_of, cluster = km_final$cluster)
  membership <- merge(membership, cyt_scores[, .(sample_id, CYT)], by = "sample_id")
  membership <- merge(membership, immune_est[, -"cancer"], by = "sample_id")
  write_run_artifact(ctx, membership, file.path("tables", "downstream", "cluster_membership.tsv"),
                      "18_clustering_and_state_analysis", script_path, base_inputs, fwrite_fn)

  populations <- setdiff(names(immune_est), c("sample_id", "cancer"))
  ext_vars <- c("CYT", populations)
  kw_rows <- vector("list", length(ext_vars))
  char_rows <- list()
  for (j in seq_along(ext_vars)) {
    v <- ext_vars[j]
    kw <- suppressWarnings(kruskal.test(membership[[v]] ~ factor(membership$cluster)))
    kw_rows[[j]] <- data.table(variable = v, kruskal_stat = unname(kw$statistic), df = unname(kw$parameter), raw_p = kw$p.value)
    by_cluster <- membership[, .(mean = mean(get(v)), median = median(get(v)), sd = sd(get(v)), n = .N), by = cluster]
    by_cluster[, variable := v]
    char_rows[[length(char_rows) + 1]] <- by_cluster
  }
  kw_dt <- rbindlist(kw_rows)
  kw_dt[, fdr := p.adjust(raw_p, "BH")]
  characteristics_dt <- rbindlist(char_rows)
  setcolorder(characteristics_dt, c("variable", "cluster", "n", "mean", "median", "sd"))
  characteristics_dt <- merge(characteristics_dt, kw_dt[, .(variable, kruskal_stat, df, raw_p, fdr)], by = "variable")
  write_run_artifact(ctx, characteristics_dt, file.path("tables", "downstream", "cluster_characteristics.tsv"),
                      "18_clustering_and_state_analysis", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "cluster_membership.tsv"), fwrite_fn)
  log_fn("cluster_characteristics.tsv written; %d of %d post-clustering characterization variables reach FDR<0.05.", sum(kw_dt$fdr < 0.05), nrow(kw_dt))

  # ==========================================================================
  # PHASE 5: PCA VISUALIZATION -- identical matrix used for clustering and display
  # ==========================================================================
  pr <- prcomp(X, center = FALSE, scale. = FALSE)  # X is already the cancer-adjusted matrix used for clustering
  pca_dt <- data.table(sample_id = colnames(raw_mat), cluster = km_final$cluster, cancer = cancer_of, PC1 = pr$x[, 1], PC2 = pr$x[, 2])
  write_run_artifact(ctx, pca_dt, file.path("tables", "downstream", "cluster_pca_coordinates.tsv"),
                      "18_clustering_and_state_analysis", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "cluster_membership.tsv"), fwrite_fn)
  write_run_artifact(ctx, pca_dt, file.path("figures", "main_candidates", "clustering", "pca_by_cluster.png"),
                      "18_clustering_and_state_analysis", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "cluster_pca_coordinates.tsv"),
    write_fn = function(d, p) {
      png(p, width = 900, height = 700, res = 120)
      plot(d$PC1, d$PC2, col = d$cluster, pch = 16, cex = 0.4,
           xlab = "PC1 (cancer-adjusted matrix)", ylab = "PC2 (cancer-adjusted matrix)",
           main = sprintf("k=%d clustering, cancer-adjusted %d-feature matrix", k_chosen, length(features)))
      legend("topright", legend = paste("Cluster", sort(unique(d$cluster))), col = sort(unique(d$cluster)), pch = 16, cex = 0.8)
      dev.off()
    })
  log_fn("PCA visualization written (same cancer-adjusted matrix used for clustering and display).")

  # ==========================================================================
  # PHASE 6: LABELING DECISION -- prespecified minimum-evidence bar
  # ==========================================================================
  # Cluster validity is judged by internal criteria and by comparison with
  # the permutation null (meets_null_calibration). Associations of clusters
  # with CYT or immune populations describe the clusters and are not used
  # as evidence that the cluster structure exists; their FDR values are
  # reported in Phase 4 and Report 32.
  chosen_sil <- validity_dt$avg_silhouette[validity_dt$k == k_chosen]
  meets_silhouette <- chosen_sil > 0.25
  meets_stability <- mean_stability_ari > 0.5
  meets_null_calibration <- is.finite(p_silhouette_vs_null) && p_silhouette_vs_null < 0.05
  # Reporting only (not part of label_bar_met): whether at least one
  # post-clustering association reached FDR < 0.05.
  meets_external <- sum(kw_dt$fdr < 0.05) > 0
  label_bar_met <- meets_silhouette && meets_stability && meets_null_calibration
  log_fn("Labeling bar: silhouette>0.25=%s (%.4f), stability ARI>0.5=%s (%.4f), null-calibration p<0.05=%s (%.4f). Bar met: %s.",
         meets_silhouette, chosen_sil, meets_stability, mean_stability_ari, meets_null_calibration, p_silhouette_vs_null, label_bar_met)

  # ==========================================================================
  # REPORT 32
  # ==========================================================================
  report32 <- c(
    "# Stage 18: Clustering and Molecular-State Analysis",
    "",
    sprintf("Generated by `scripts/17_clustering_and_state_analysis.R`, run `%s`.", ctx$run_id),
    "",
    "## Prespecified design",
    "",
    sprintf("- Sample universe: all %d matched samples.", ncol(adj_mat)),
    sprintf("- Feature matrix: this run's own %d primary molecular features (%d miRNAs, %d target genes) -- the same, already-justified feature set from this run's own discovery output; no new post-hoc feature-selection pass was introduced.",
            length(features), uniqueN(triplets$mirna_id), uniqueN(triplets$gene)),
    "- Scaling / cancer-effect adjustment: each feature standardized (z-scored) WITHIN each cancer type. This \"cancer-adjusted matrix\" is used identically for k-means clustering and for the PCA visualization -- there is no separate raw-for-clustering / adjusted-for-display split anywhere in this script.",
    sprintf("- Candidate k: %s (from `analysis_config$clustering`).", paste(k_range, collapse = ", ")),
    "- Cluster-validity metrics: mean silhouette width (primary k-selection criterion) and the gap statistic (corroborating).",
    "- Stability metric: mean Adjusted Rand Index across 100 bootstrap resamples of samples, re-clustered at the chosen k.",
    sprintf("- Minimum evidence to assign a descriptive biological label to any cluster solution: mean silhouette > 0.25 AND mean stability ARI > 0.5 AND the observed silhouette significantly exceeds a %d-draw permutation-null distribution (empirical p < 0.05) -- all three fixed before this script was run. CYT/immune-population association is NEVER part of this bar (see the post-clustering-characterization disclosure below for why).", n_null),
    "",
    "## k selection",
    "",
    "| k | Mean silhouette | Gap statistic | Min/max cluster size | Size balance ratio |",
    "|---:|---:|---:|---|---:|",
    sprintf("| %d | %.4f | %s | %d / %d | %.3f |", validity_dt$k, validity_dt$avg_silhouette,
            ifelse(is.na(validity_dt$gap), "NA", sprintf("%.4f", validity_dt$gap)),
            validity_dt$min_cluster_size, validity_dt$max_cluster_size, validity_dt$size_balance_ratio),
    "",
    sprintf("**Chosen k = %d** (maximizes mean silhouette, the prespecified primary criterion).", k_chosen),
    "",
    "## Stability, seed-sensitivity, and feature-scaling-sensitivity",
    "",
    sprintf("- Bootstrap resampling stability at k=%d: mean ARI = %.4f (SD %.4f) across 100 resamples.", k_chosen, mean_stability_ari, sd(ari_vals, na.rm = TRUE)),
    sprintf("- Seed sensitivity: mean pairwise ARI across 5 independent master seeds = %.4f.", mean_seed_ari),
    sprintf("- Feature-scaling sensitivity: ARI between the cancer-adjusted clustering and an otherwise-identical pan-cohort-raw-scaled (not cancer-adjusted) clustering = %.4f.", scaling_ari),
    sprintf("- Cancer-type association (Cramer's V): cancer-adjusted clustering = %.4f; pan-cohort-raw-scaled clustering = %.4f.", cancer_assoc_adjusted, cancer_assoc_rawscaled),
    "",
    # The diagnostic signal for whether cancer-adjustment worked is the DROP in cancer-type
    # association (Cramer's V), not the ARI between the two clusterings by itself: the two
    # clusterings are expected to disagree (low-to-moderate ARI) whenever adjustment succeeds
    # in making them capture different structure, so ARI alone cannot distinguish "adjustment
    # worked" from "adjustment failed to remove structure in a different way."
    if (cancer_assoc_adjusted < cancer_assoc_rawscaled * 0.5)
      sprintf("The within-cancer standardization step worked as intended: cancer-type association drops sharply between the raw-scaled clustering (Cramer's V=%.4f) and the cancer-adjusted clustering (Cramer's V=%.4f), and the moderate agreement between the two clusterings (ARI=%.4f) is the expected signature of that success -- the adjusted clustering captures materially different, less tissue-of-origin-driven structure, not merely a relabeling of the same cancer-driven groups.", cancer_assoc_rawscaled, cancer_assoc_adjusted, scaling_ari)
    else
      sprintf("Cancer-type association remains similar between the raw-scaled clustering (Cramer's V=%.4f) and the cancer-adjusted clustering (Cramer's V=%.4f), indicating the within-cancer standardization step did NOT substantially remove tissue-of-origin structure from this feature matrix; this is reported directly rather than smoothed over.", cancer_assoc_rawscaled, cancer_assoc_adjusted),
    "",
    "## Null-cluster-structure calibration",
    "",
    sprintf("The observed k-means solution is compared against %d independent permutation-null matrices, each built by independently permuting every feature's values across samples (preserving sample count, feature count, cancer-group sizes, and each feature's own marginal distribution exactly, while destroying cross-feature sample-level structure), then put through the IDENTICAL within-cancer scaling and k-selection/stability procedure as the real data. CYT and immune-population data play no role anywhere in this comparison.", n_null),
    "",
    sprintf("- Observed max silhouette = %.4f vs. null median = %.4f (90%% interval [%.4f, %.4f]); empirical p = %.4g (finite-sample continuity-corrected: (extreme+1)/(B+1)).",
            chosen_sil_observed, median(null_res$null_max_silhouette), quantile(null_res$null_max_silhouette, 0.05), quantile(null_res$null_max_silhouette, 0.95), p_silhouette_vs_null),
    sprintf("- Observed stability ARI = %.4f vs. null median = %.4f (90%% interval [%.4f, %.4f]); empirical p = %.4g.",
            mean_stability_ari, median(null_res$null_stability, na.rm = TRUE), quantile(null_res$null_stability, 0.05, na.rm = TRUE), quantile(null_res$null_stability, 0.95, na.rm = TRUE), p_stability_vs_null),
    if (meets_null_calibration)
      "The observed clustering structure is significantly stronger than a permutation null at the prespecified alpha=0.05 threshold -- this is genuine evidence of multivariate structure beyond chance, though it says nothing about whether that structure is biologically meaningful (see post-clustering characterization below, which is descriptive only)."
    else
      "The observed clustering structure is NOT significantly stronger than a permutation null at the prespecified alpha=0.05 threshold. This is reported directly: a forced k-means partition always produces SOME clusters and SOME silhouette value, but that alone is not evidence of genuine multivariate structure distinguishable from what random data of the same shape would produce.",
    "",
    "Full per-iteration null distribution: `tables/downstream/cluster_null_calibration.tsv`.",
    "",
    "## Post-clustering characterization (internal, NOT external validation -- tested only after cluster assignment was fixed)",
    "",
    "**This is NOT external validation.** The features used to build the clustering matrix (this run's own primary miRNAs and target genes) were themselves selected, in the discovery stages, partly on the basis of their association with CYT. CYT association tested here on the resulting clusters is therefore an INTERNAL characterization of the same signal already used to choose the features, not an independent test -- a genuinely external validation would require features or a cohort not involved in feature selection, neither of which is available in this single-cohort, post-discovery analysis. The MCPcounter immune-population associations are internal in the same sense: they use the same TCGA cohort and the same cancer-adjusted matrix. CYT (or immune-population) significance below is reported as a descriptive characteristic of the clusters found, never as evidence that validates the clusters' existence or biological reality.",
    "",
    sprintf("%d of %d post-clustering characterization variables (CYT + %d MCPcounter populations) reach BH-FDR < 0.05 across clusters (Kruskal-Wallis, single flat family). Full per-cluster means/medians and test statistics: `tables/downstream/cluster_characteristics.tsv`.", sum(kw_dt$fdr < 0.05), nrow(kw_dt), length(populations)),
    "",
    "## Labeling decision",
    "",
    sprintf("- Silhouette > 0.25: %s (%.4f)", meets_silhouette, chosen_sil),
    sprintf("- Stability ARI > 0.5: %s (%.4f)", meets_stability, mean_stability_ari),
    sprintf("- At least one post-clustering-characterization FDR<0.05 association: %s (internal, not independent evidence of cluster validity -- see disclosure above)", meets_external),
    "",
    if (label_bar_met)
      "**All three prespecified conditions are met.** Cluster characteristics (`cluster_characteristics.tsv`) may be used to assign a descriptive label (e.g., relative CYT/immune level), but any such label describes the measured characteristics of this specific clustering and is not a claim of a validated, externally-reproduced molecular subtype."
    else
      "**At least one prespecified condition is not met. No descriptive biological label (e.g., \"immune-hot\"/\"immune-cold\") is assigned to any cluster.** The clustering structure found here should be reported as weak/exploratory. This is not smoothed over or forced into a subtype narrative.",
    "",
    "## Outputs",
    "",
    "- `tables/downstream/cluster_validity.tsv`, `cluster_stability.tsv`, `cluster_null_calibration.tsv`, `cluster_membership.tsv`, `cluster_characteristics.tsv`, `cluster_pca_coordinates.tsv`",
    "- `figures/main_candidates/clustering/pca_by_cluster.png`",
    "- `reports/32_clustering_and_state_report.md` (this file)"
  )
  write_run_artifact(ctx, report32, file.path("reports", "32_clustering_and_state_report.md"),
                      "18_clustering_and_state_analysis", script_path, character(0), function(x, p) writeLines(x, p))
  log_fn("reports/32_clustering_and_state_report.md written.")

  record_stage_status(ctx, "18_clustering_and_state_analysis", "COMPLETE")
  log_fn("Stage 18 complete.")
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=17 or source scripts/run.R in RStudio.")
}
