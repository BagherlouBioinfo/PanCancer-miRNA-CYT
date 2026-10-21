# ============================================================================
# 23_cancer_specific_cellular_context.R
#
# Purpose: for every primary triplet of the current run, determine in which
#          cancer types and MCPcounter cellular-composition contexts the
#          triplet's three-arm relationship varies with composition. The
#          analysis is cancer-specific and uses the run's sample-level
#          MCPcounter estimates (Stage 12) rather than the recurrence
#          summary.
#
# Statistical design (formulas in R/cell_context_helpers.R):
#   For every eligible (triplet x cancer x MCPcounter population), the three
#   arms (miRNA-CYT, gene-CYT, miRNA-gene) are evaluated on the
#   CPE-residualized values computed as in Stage 10 (same preprocessing,
#   eligibility rule analysis_config$discovery$min_samples_per_cancer, and
#   covariate CPE).
#
#   Primary test: a continuous rank-based interaction model,
#   rank(outcome) ~ rank(predictor) * standardize(population_abundance),
#   fitted once per (triplet, cancer, population, arm) with
#   R/cell_context_helpers.R::fit_cell_context_interaction(). The population
#   is not dichotomized. FDR family (prespecified): Benjamini-Hochberg across
#   every evaluable, non-tautological (triplet, cancer, population, arm)
#   interaction p-value of the run, a single family distinct from those of
#   other stages.
#
#   Descriptive summary: low/middle/high tertiles of the population's raw
#   within-cancer abundance, with the Spearman correlation and Fisher-z CI
#   of the same residualized predictor and outcome in each stratum
#   (R/cell_context_helpers.R::tertile_stratified_correlation()). This
#   summary is interpretive only and is not FDR-corrected.
#
#   Each MCPcounter population is analysed separately (its own interaction
#   model and tertile split); populations are not combined into a composite
#   score or placed together in one multivariable model.
#
#   Tautology safeguard: when the target gene is one of a population's
#   MCPcounter marker genes (R/cell_context_helpers.R::
#   is_gene_signature_member(), read from
#   reference/mcpcounter_signatures/genes.txt), the gene-CYT and miRNA-gene
#   arms for that (gene, population) pair are circular, because the
#   population estimate is partly constructed from the gene's expression.
#   They are marked tautological and excluded from the interaction-FDR
#   family and from evidence-based classification. The miRNA-CYT arm does
#   not involve the gene and is evaluated normally.
#
#   Six-category classification per (triplet, cancer, population)
#   (R/cell_context_helpers.R::classify_cell_context()):
#   tautological_exclusion, not_evaluable, composition_enriched,
#   composition_depleted, stable_across_composition, not_detected.
#   "Detected" (for the stable/not_detected distinction) uses Stage 10's
#   three-arm pass/fail per (triplet, cancer) from
#   tables/audit/cancer_level_evidence_CPE_rho030.tsv.
#
# Interpretation: MCPcounter estimates cellular composition of bulk tissue
# and cannot localize mature miRNA expression to individual cells. No
# result of this stage establishes the cellular origin of a miRNA or direct
# regulation.
#
# Inputs:  <run_root>/tables/primary/<primary triplet file> (via
#          load_run_primary_triplets(ctx), same-run);
#          tables/audit/cancer_level_evidence_CPE_rho030.tsv (via
#          load_run_evidence(ctx, "CPE"), same-run);
#          processed/stage03/mirna_matched.rds, mrna_matched.rds;
#          processed/stage05/sample_annotation_with_purity.tsv,
#          eligible_cancers_cpe.tsv;
#          processed/stage11_immune/sample_level_immune_estimates.tsv
#          (all same-run, provenance-verified);
#          reference/mcpcounter_signatures/genes.txt (allowlist-verified).
# Outputs: <run_root>/tables/downstream/triplet_cell_context_long.tsv,
#          triplet_cell_context_summary.tsv,
#          triplet_cell_context_all_three_arms_significant.tsv,
#          triplet_cell_context_by_cancer.tsv;
#          figures/cell_context/*.pdf, *.png;
#          reports/CANCER_SPECIFIC_CELLULAR_CONTEXT_REPORT.md.
# Depends on (same-run, enforced by read_run_artifact()): 03, 05, 09, 10, 12.
# ============================================================================

run_stage23 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "23_cancer_specific_cellular_context")
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "run_discovery_source.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "expression_preprocessing.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "stats_engine.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "cell_context_helpers.R"), local = TRUE)
  suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
  })
  script_path <- file.path(ctx$data_root, "scripts", "22_cancer_specific_cellular_context.R")
  fwrite_fn <- function(x, p) fwrite(x, p, sep = "\t")

  min_n <- analysis_config$discovery$min_samples_per_cancer
  min_n_stratum <- analysis_config$cell_context$min_samples_per_tertile
  fdr_alpha <- analysis_config$cell_context$interaction_fdr_cutoff
  rho_thr <- analysis_config$thresholds$primary

  triplets <- load_run_primary_triplets(ctx)
  run_stage_assert(ctx, log_fn, all(c("mirna_id", "gene") %in% names(triplets)) && nrow(triplets) > 0,
                    "This run's primary triplet table has mirna_id/gene columns and at least 1 row (%d triplets loaded)", nrow(triplets))

  baseline_evidence <- load_run_evidence(ctx, "CPE")
  setkey(baseline_evidence, mirna_id, gene, cancer)

  # --- MCPcounter signature map: the reference file used by Stage 12
  #     (allowlist-verified). ---
  sig_path <- file.path(ctx$data_root, "reference", "mcpcounter_signatures", "genes.txt")
  verify_allowed_raw_input(ctx$data_root, sig_path)
  sig_raw <- read.table(sig_path, sep = "\t", stringsAsFactors = FALSE, header = TRUE, colClasses = "character", check.names = FALSE)
  signature_map <- data.frame(gene = toupper(trimws(sig_raw[["HUGO symbols"]])), population = sig_raw[["Cell population"]], stringsAsFactors = FALSE)
  populations <- sort(unique(signature_map$population))
  log_fn("MCPcounter populations (this run's own reference file, %d marker genes): %s", nrow(signature_map), paste(populations, collapse = ", "))

  # --- matched matrices, winsorized within each cancer as in Stage 10
  #     (scripts/lib/expression_preprocessing.R::winsorize_matrix_rows_by_group()). ---
  mirna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mirna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mirna_mat <- as.matrix(mirna_matched[, -1, drop = FALSE]); rownames(mirna_mat) <- mirna_matched$mirna_id
  mrna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mrna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mrna_mat <- as.matrix(mrna_matched[, -1, drop = FALSE]); rownames(mrna_mat) <- toupper(trimws(mrna_matched$gene_id_raw))
  mrna_mat <- deduplicate_matrix_rows(mrna_mat)

  sample_annotation <- read_run_artifact(ctx, file.path("processed", "stage05", "sample_annotation_with_purity.tsv"), "05_define_eligible_features",
                                          function(p) read.delim(p, stringsAsFactors = FALSE))
  eligible_cancers <- sort(read_run_artifact(ctx, file.path("processed", "stage05", "eligible_cancers_cpe.tsv"), "05_define_eligible_features",
                                              function(p) read.delim(p, stringsAsFactors = FALSE))$cancer)

  cancer_group_mirna <- sample_annotation$cancer[match(colnames(mirna_mat), sample_annotation$sample_id)]
  run_stage_assert(ctx, log_fn, !anyNA(cancer_group_mirna), "Every column of the matched miRNA matrix has a known cancer-type assignment (required for per-cancer winsorization)")
  cancer_group_mrna <- sample_annotation$cancer[match(colnames(mrna_mat), sample_annotation$sample_id)]
  run_stage_assert(ctx, log_fn, !anyNA(cancer_group_mrna), "Every column of the matched mRNA matrix has a known cancer-type assignment (required for per-cancer winsorization)")
  mirna_mat <- winsorize_matrix_rows_by_group(mirna_mat, cancer_group_mirna)
  mrna_mat <- winsorize_matrix_rows_by_group(mrna_mat, cancer_group_mrna)

  immune_est <- read_run_artifact(ctx, file.path("processed", "stage11_immune", "sample_level_immune_estimates.tsv"), "12_immune_composition_analysis", fread)
  run_stage_assert(ctx, log_fn, all(populations %in% names(immune_est)),
                    "Every MCPcounter population in this run's own signature file has a matching column in Stage 12's sample-level immune estimates")

  n_triplets <- nrow(triplets)
  n_cancers <- length(eligible_cancers)
  n_populations <- length(populations)
  log_fn("Cancer-specific cellular-context analysis: %d primary triplets x up to %d eligible cancers x %d MCPcounter populations x 3 arms.",
         n_triplets, n_cancers, n_populations)

  long_rows <- vector("list", n_triplets * n_cancers * n_populations * 3)
  ridx <- 0L

  for (ti in seq_len(n_triplets)) {
    mir <- triplets$mirna_id[ti]; gene_x <- triplets$gene[ti]
    if (!(mir %in% rownames(mirna_mat)) || !(gene_x %in% rownames(mrna_mat))) {
      log_fn("  [%d/%d] %s-%s: miRNA or gene not present in the matched matrices; skipped entirely.", ti, n_triplets, mir, gene_x)
      next
    }

    for (cc in eligible_cancers) {
      idx <- which(sample_annotation$cancer == cc & !is.na(sample_annotation$CPE) & !is.na(sample_annotation$CYT))
      idx_samples <- sample_annotation$sample_id[idx]

      base_row <- baseline_evidence[.(mir, gene_x, cc), on = c("mirna_id", "gene", "cancer")]
      baseline_pass <- c(
        mirna_cyt = if (nrow(base_row) == 1) isTRUE(base_row$pass_mirna_cyt) else FALSE,
        gene_cyt  = if (nrow(base_row) == 1) isTRUE(base_row$pass_gene_cyt) else FALSE,
        mirna_gene = if (nrow(base_row) == 1) isTRUE(base_row$pass_mirna_gene) else FALSE
      )

      if (length(idx) < min_n) {
        for (pop in populations) {
          for (arm in c("mirna_cyt", "gene_cyt", "mirna_gene")) {
            ridx <- ridx + 1L
            long_rows[[ridx]] <- data.table(mirna_id = mir, gene = gene_x, cancer = cc, population = pop, arm = arm,
                                             n = length(idx), tautological = FALSE, estimable = FALSE,
                                             beta_main = NA_real_, beta_interaction = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                                             p_interaction = NA_real_, slope_low = NA_real_, slope_high = NA_real_,
                                             tertile_low_rho = NA_real_, tertile_low_n = NA_integer_,
                                             tertile_mid_rho = NA_real_, tertile_mid_n = NA_integer_,
                                             tertile_high_rho = NA_real_, tertile_high_n = NA_integer_,
                                             baseline_pass_arm = baseline_pass[[arm]])
          }
        }
        next
      }

      purity_vec <- sample_annotation$CPE[idx]
      mir_resid <- as.numeric(residualize_matrix(mirna_mat[mir, idx, drop = FALSE], purity_vec))
      gene_resid <- as.numeric(residualize_matrix(mrna_mat[gene_x, idx, drop = FALSE], purity_vec))
      cyt_resid <- as.numeric(residualize_matrix(matrix(sample_annotation$CYT[idx], nrow = 1), purity_vec))

      immune_match <- match(idx_samples, immune_est$sample_id)

      for (pop in populations) {
        pop_vec <- immune_est[[pop]][immune_match]
        gene_tautological <- is_gene_signature_member(gene_x, pop, signature_map)
        tautological_arm <- c(mirna_cyt = FALSE, gene_cyt = gene_tautological, mirna_gene = gene_tautological)

        arm_pairs <- list(
          mirna_cyt = list(predictor = mir_resid, outcome = cyt_resid),
          gene_cyt = list(predictor = gene_resid, outcome = cyt_resid),
          mirna_gene = list(predictor = mir_resid, outcome = gene_resid)
        )

        for (arm in names(arm_pairs)) {
          ridx <- ridx + 1L
          if (isTRUE(tautological_arm[[arm]])) {
            long_rows[[ridx]] <- data.table(mirna_id = mir, gene = gene_x, cancer = cc, population = pop, arm = arm,
                                             n = length(idx), tautological = TRUE, estimable = FALSE,
                                             beta_main = NA_real_, beta_interaction = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                                             p_interaction = NA_real_, slope_low = NA_real_, slope_high = NA_real_,
                                             tertile_low_rho = NA_real_, tertile_low_n = NA_integer_,
                                             tertile_mid_rho = NA_real_, tertile_mid_n = NA_integer_,
                                             tertile_high_rho = NA_real_, tertile_high_n = NA_integer_,
                                             baseline_pass_arm = baseline_pass[[arm]])
            next
          }
          pred <- arm_pairs[[arm]]$predictor; out <- arm_pairs[[arm]]$outcome
          fit <- fit_cell_context_interaction(pred, out, pop_vec, min_n = min_n)
          strat <- tertile_stratified_correlation(pred, out, pop_vec, min_n_stratum = min_n_stratum)
          long_rows[[ridx]] <- data.table(
            mirna_id = mir, gene = gene_x, cancer = cc, population = pop, arm = arm,
            n = fit$n, tautological = FALSE, estimable = fit$estimable,
            beta_main = fit$beta_main, beta_interaction = fit$beta_interaction,
            ci_lower = fit$ci_lower, ci_upper = fit$ci_upper, p_interaction = fit$p_interaction,
            slope_low = fit$slope_low, slope_high = fit$slope_high,
            tertile_low_rho = strat$rho[strat$stratum == "low"], tertile_low_n = strat$n[strat$stratum == "low"],
            tertile_mid_rho = strat$rho[strat$stratum == "middle"], tertile_mid_n = strat$n[strat$stratum == "middle"],
            tertile_high_rho = strat$rho[strat$stratum == "high"], tertile_high_n = strat$n[strat$stratum == "high"],
            baseline_pass_arm = baseline_pass[[arm]]
          )
        }
      }
    }
    if (ti %% 5 == 0) log_fn("  Processed %d/%d triplets.", ti, n_triplets)
  }

  long_dt <- rbindlist(long_rows[seq_len(ridx)])
  run_stage_assert(ctx, log_fn, nrow(long_dt) > 0, "At least one (triplet, cancer, population, arm) row was produced")

  # --- Single prespecified BH-FDR family: every evaluable, non-tautological
  #     interaction p-value of the run, corrected together. ---
  fam_idx <- which(long_dt$estimable & !long_dt$tautological)
  long_dt[, fdr_interaction := NA_real_]
  if (length(fam_idx) > 0) {
    long_dt$fdr_interaction[fam_idx] <- p.adjust(long_dt$p_interaction[fam_idx], method = "BH")
  }
  # The count below is arm-level (up to three per (triplet, cancer,
  # population) combination) and is a different unit from the
  # combination-level composition_enriched/depleted counts.
  n_sig_arms <- sum(long_dt$fdr_interaction < fdr_alpha, na.rm = TRUE)
  log_fn("Interaction-test FDR family: %d evaluable, non-tautological ARM-LEVEL tests (up to 3 per triplet x cancer x population combination); %d of those individual arm-level tests are significant at FDR<%.2f. This is NOT the same unit as the composition_enriched/composition_depleted combination counts reported below (a combination requires only ONE significant arm, by majority direction) -- see the all-three-arms-significant count for the stricter triplet-level result.",
         length(fam_idx), n_sig_arms, fdr_alpha)

  write_run_artifact(ctx, long_dt, file.path("tables", "downstream", "triplet_cell_context_long.tsv"),
                      "23_cancer_specific_cellular_context", script_path,
                      file.path(ctx$run_root, c(file.path("tables", "primary", run_primary_triplet_filename(ctx)),
                                                 file.path("processed", "stage03", "mirna_matched.rds"),
                                                 file.path("processed", "stage03", "mrna_matched.rds"),
                                                 file.path("processed", "stage05", "sample_annotation_with_purity.tsv"),
                                                 file.path("processed", "stage11_immune", "sample_level_immune_estimates.tsv"))),
                      fwrite_fn)

  # ==========================================================================
  # Summary table: one row per (triplet, cancer, population), combining the
  # 3 arms via classify_cell_context().
  # ==========================================================================
  summary_keys <- unique(long_dt[, .(mirna_id, gene, cancer, population)])
  summary_rows <- vector("list", nrow(summary_keys))
  for (i in seq_len(nrow(summary_keys))) {
    k <- summary_keys[i]
    rows_i <- long_dt[mirna_id == k$mirna_id & gene == k$gene & cancer == k$cancer & population == k$population]
    arm_results <- list()
    for (arm_name in c("mirna_cyt", "gene_cyt", "mirna_gene")) {
      r <- rows_i[arm == arm_name]
      arm_results[[arm_name]] <- list(estimable = isTRUE(r$estimable), n = r$n, beta_main = r$beta_main,
                                       beta_interaction = r$beta_interaction, slope_low = r$slope_low, slope_high = r$slope_high,
                                       fdr = r$fdr_interaction)
    }
    tautological_arms <- rows_i$arm[rows_i$tautological]
    baseline_pass_vec <- setNames(rows_i$baseline_pass_arm, rows_i$arm)
    verdict <- classify_cell_context(arm_results, tautological_arms, baseline_pass_vec, fdr_alpha = fdr_alpha)
    # n_arms_fdr_significant: the number of evaluable, non-tautological arms
    # (up to three) that individually reach FDR < fdr_alpha. This
    # descriptive column distinguishes combinations classified on one
    # significant arm from those classified on all three; `category` is not
    # affected.
    n_arms_sig <- sum(vapply(c("mirna_cyt", "gene_cyt", "mirna_gene"), function(a) {
      r <- arm_results[[a]]
      isTRUE(r$estimable) && !is.na(r$fdr) && r$fdr < fdr_alpha
    }, logical(1)))
    summary_rows[[i]] <- data.table(
      mirna_id = k$mirna_id, gene = k$gene, cancer = k$cancer, population = k$population,
      n_samples = max(rows_i$n, na.rm = TRUE),
      all_three_arms_estimable = all(vapply(c("mirna_cyt", "gene_cyt", "mirna_gene"), function(a) isTRUE(arm_results[[a]]$estimable), logical(1))),
      min_interaction_fdr = suppressWarnings(min(rows_i$fdr_interaction, na.rm = TRUE)),
      n_arms_fdr_significant = n_arms_sig,
      all_three_arms_fdr_significant = isTRUE(n_arms_sig == 3),
      category = verdict
    )
  }
  summary_dt <- rbindlist(summary_rows)
  summary_dt[!is.finite(min_interaction_fdr), min_interaction_fdr := NA_real_]
  n_all_three_sig <- sum(summary_dt$all_three_arms_fdr_significant)
  log_fn("Summary classification counts: %s", paste(sprintf("%s=%d", names(table(summary_dt$category)), table(summary_dt$category)), collapse = ", "))
  log_fn("Of %d total (triplet,cancer,population) combinations, %d (%.2f%%) have ALL THREE arms simultaneously FDR-significant -- this is the stricter, scientifically meaningful triplet-level result; a single significant arm is not, by itself, evidence that the complete miRNA-gene-CYT triplet is composition-dependent (%d of the %d composition_enriched/composition_depleted combinations, %.1f%%, are significant in fewer than all three arms).",
         nrow(summary_dt), n_all_three_sig, 100 * n_all_three_sig / nrow(summary_dt),
         sum(summary_dt$category %in% c("composition_enriched", "composition_depleted") & !summary_dt$all_three_arms_fdr_significant),
         sum(summary_dt$category %in% c("composition_enriched", "composition_depleted")),
         100 * sum(summary_dt$category %in% c("composition_enriched", "composition_depleted") & !summary_dt$all_three_arms_fdr_significant) / max(1, sum(summary_dt$category %in% c("composition_enriched", "composition_depleted"))))

  write_run_artifact(ctx, summary_dt, file.path("tables", "downstream", "triplet_cell_context_summary.tsv"),
                      "23_cancer_specific_cellular_context", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "triplet_cell_context_long.tsv"), fwrite_fn)

  # All-three-arms-significant table: the subset of
  # composition_enriched/composition_depleted combinations in which all three
  # arms are individually FDR-significant. A single significant arm is not
  # evidence that the complete miRNA-gene-CYT triplet is
  # composition-dependent; this table is the stricter summary. Stage 24
  # replaces arm-wise significance with a formal conjunction test.
  all_three_dt <- summary_dt[all_three_arms_fdr_significant == TRUE,
                              .(mirna_id, gene, cancer, population, category, n_samples, min_interaction_fdr)]
  setorder(all_three_dt, min_interaction_fdr)
  write_run_artifact(ctx, all_three_dt, file.path("tables", "downstream", "triplet_cell_context_all_three_arms_significant.tsv"),
                      "23_cancer_specific_cellular_context", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "triplet_cell_context_long.tsv"), fwrite_fn)

  # ==========================================================================
  # By-cancer readable table: one row per (triplet, cancer), listing which
  # populations fall into each category.
  # ==========================================================================
  by_cancer_rows <- unique(summary_dt[, .(mirna_id, gene, cancer)])
  by_cancer_list <- vector("list", nrow(by_cancer_rows))
  for (i in seq_len(nrow(by_cancer_rows))) {
    k <- by_cancer_rows[i]
    sub <- summary_dt[mirna_id == k$mirna_id & gene == k$gene & cancer == k$cancer]
    list_pops <- function(cat) paste(sort(sub$population[sub$category == cat]), collapse = "; ")
    by_cancer_list[[i]] <- data.table(
      mirna_id = k$mirna_id, gene = k$gene, cancer = k$cancer,
      composition_enriched_populations = list_pops("composition_enriched"),
      composition_depleted_populations = list_pops("composition_depleted"),
      stable_across_composition_populations = list_pops("stable_across_composition"),
      not_detected_populations = list_pops("not_detected"),
      not_evaluable_populations = list_pops("not_evaluable"),
      tautological_exclusion_populations = list_pops("tautological_exclusion")
    )
  }
  by_cancer_dt <- rbindlist(by_cancer_list)
  write_run_artifact(ctx, by_cancer_dt, file.path("tables", "downstream", "triplet_cell_context_by_cancer.tsv"),
                      "23_cancer_specific_cellular_context", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "triplet_cell_context_summary.tsv"), fwrite_fn)

  # ==========================================================================
  # Heatmaps
  # ==========================================================================
  fig_dir <- file.path(ctx$run_root, "figures", "cell_context")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  cat_colors <- c(
    composition_enriched = "#D55E00", composition_depleted = "#0072B2",
    stable_across_composition = "#009E73", not_detected = "#F0E442",
    not_evaluable = "grey80", tautological_exclusion = "grey50"
  )  # colour-blind-safe (Okabe-Ito); grey is used only for not_evaluable/tautological, not for not_detected

  summary_dt[, triplet_label := paste(mirna_id, gene, sep = "-")]
  summary_dt[, category := factor(category, levels = names(cat_colors))]

  save_plot <- function(p, filename_base, width, height) {
    ggsave(file.path(fig_dir, paste0(filename_base, ".pdf")), p, width = width, height = height, limitsize = FALSE)
    ggsave(file.path(fig_dir, paste0(filename_base, ".png")), p, width = width, height = height, dpi = 300, limitsize = FALSE)
  }

  p_main <- ggplot(summary_dt, aes(x = population, y = cancer, fill = category)) +
    geom_tile(color = "white", linewidth = 0.2) +
    facet_wrap(~triplet_label, ncol = 6) +
    scale_fill_manual(values = cat_colors, drop = FALSE, name = "Cellular-context category") +
    theme_minimal(base_size = 7) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
          axis.text.y = element_text(size = 5), strip.text = element_text(size = 6),
          panel.grid = element_blank()) +
    labs(title = "Cancer-specific cellular-composition context, all primary triplets",
         subtitle = "Grey = not evaluable / tautological exclusion (never a null result); descriptive over-representation only",
         x = "MCPcounter population", y = "Cancer type")
  save_plot(p_main, "01_all_triplets_faceted_heatmap", width = 16, height = max(8, ceiling(n_triplets / 6) * 3))
  log_fn("Figure written: figures/cell_context/01_all_triplets_faceted_heatmap.{pdf,png}")

  interaction_effect_dt <- long_dt[arm == "mirna_gene" & estimable == TRUE & tautological == FALSE,
                                     .(effect = mean(abs(slope_high) - abs(slope_low), na.rm = TRUE)), by = .(mirna_id, gene, cancer, population)]
  interaction_effect_dt[, triplet_label := paste(mirna_id, gene, sep = "-")]
  p_effect <- ggplot(interaction_effect_dt, aes(x = population, y = cancer, fill = effect)) +
    geom_tile(color = "white", linewidth = 0.2) +
    facet_wrap(~triplet_label, ncol = 6) +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, name = "|slope_high|-|slope_low|\n(miRNA-gene arm)", na.value = "grey85") +
    theme_minimal(base_size = 7) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5), axis.text.y = element_text(size = 5),
          strip.text = element_text(size = 6), panel.grid = element_blank()) +
    labs(title = "Interaction effect size (miRNA-gene arm): |slope at high population| - |slope at low population|",
         subtitle = "Positive (orange) = composition-enriched direction; negative (blue) = composition-depleted direction. Grey = not evaluable.",
         x = "MCPcounter population", y = "Cancer type")
  save_plot(p_effect, "02_interaction_effect_size_heatmap", width = 16, height = max(8, ceiling(n_triplets / 6) * 3))
  log_fn("Figure written: figures/cell_context/02_interaction_effect_size_heatmap.{pdf,png}")

  fdr_dt <- long_dt[arm == "mirna_gene", .(min_fdr = suppressWarnings(min(fdr_interaction, na.rm = TRUE)),
                                            any_tautological = any(tautological)), by = .(mirna_id, gene, cancer, population)]
  fdr_dt[!is.finite(min_fdr), min_fdr := NA_real_]
  fdr_dt[, fdr_supported := fifelse(any_tautological, NA, min_fdr < fdr_alpha)]
  fdr_dt[, triplet_label := paste(mirna_id, gene, sep = "-")]
  p_fdr <- ggplot(fdr_dt, aes(x = population, y = cancer, fill = fdr_supported)) +
    geom_tile(color = "white", linewidth = 0.2) +
    facet_wrap(~triplet_label, ncol = 6) +
    scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#56B4E9"), na.value = "grey80",
                       name = sprintf("Interaction\nFDR<%.2f", fdr_alpha), breaks = c(TRUE, FALSE),
                       labels = c("Supported", "Not supported")) +
    theme_minimal(base_size = 7) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5), axis.text.y = element_text(size = 5),
          strip.text = element_text(size = 6), panel.grid = element_blank()) +
    labs(title = "FDR-supported classification (miRNA-gene arm interaction test)",
         subtitle = "Grey = not evaluable or tautological exclusion, NEVER a null (not-supported) result",
         x = "MCPcounter population", y = "Cancer type")
  save_plot(p_fdr, "03_fdr_supported_classification_heatmap", width = 16, height = max(8, ceiling(n_triplets / 6) * 3))
  log_fn("Figure written: figures/cell_context/03_fdr_supported_classification_heatmap.{pdf,png}")

  focused_pairs <- list(c("hsa-miR-146a-5p", "CD3D"), c("hsa-miR-155-5p", "CTLA4"))
  focused_present <- vapply(focused_pairs, function(pr) any(triplets$mirna_id == pr[1] & triplets$gene == pr[2]), logical(1))
  if (any(focused_present)) {
    focused_dt <- summary_dt[paste(mirna_id, gene) %in% vapply(focused_pairs[focused_present], function(pr) paste(pr[1], pr[2]), character(1))]
    p_focused <- ggplot(focused_dt, aes(x = population, y = cancer, fill = category)) +
      geom_tile(color = "white", linewidth = 0.3) +
      facet_wrap(~triplet_label, ncol = 2) +
      scale_fill_manual(values = cat_colors, drop = FALSE, name = "Cellular-context category") +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid = element_blank()) +
      labs(title = "Focused view: previously most composition-robust relationships",
           subtitle = "hsa-miR-146a-5p-CD3D and hsa-miR-155-5p-CTLA4 (if present among this run's primary triplets)",
           x = "MCPcounter population", y = "Cancer type")
    save_plot(p_focused, "04_focused_146a_CD3D_155_CTLA4", width = 10, height = max(6, n_cancers * 0.3))
    log_fn("Figure written: figures/cell_context/04_focused_146a_CD3D_155_CTLA4.{pdf,png} (%d of 2 focus pairs present in this run's primary triplets)", sum(focused_present))
  } else {
    log_fn("Focused figure skipped: neither hsa-miR-146a-5p/CD3D nor hsa-miR-155-5p/CTLA4 is among this run's own primary triplets -- this run's actual triplet set determines this, not assumed.")
  }

  # ==========================================================================
  # Report
  # ==========================================================================
  cat_counts <- table(summary_dt$category)
  n_triplet_cancer_pop <- nrow(summary_dt)
  strongest <- summary_dt[category %in% c("composition_enriched", "composition_depleted") & is.finite(min_interaction_fdr)]
  setorder(strongest, min_interaction_fdr)
  top_n <- min(15, nrow(strongest))

  report_lines <- c(
    "# Cancer-Specific Cellular-Composition Context Report",
    "",
    sprintf("Generated by `scripts/22_cancer_specific_cellular_context.R`, run `%s`, %s.", ctx$run_id, format(Sys.time(), "%Y-%m-%d %H:%M:%OS0")),
    "",
    "**What this analysis is and is not.** This stage tests whether each",
    "primary triplet's three-arm association is MODIFIED by a given",
    "MCPcounter cellular-composition population's estimated abundance,",
    "within each eligible cancer type, using a prespecified continuous",
    "interaction model as the primary significance test (population is",
    "never dichotomized for that test) and tertile strata as a purely",
    "descriptive, interpretive supplement. **MCPcounter estimates cellular",
    "composition in bulk tumor tissue; it does not, and cannot, localize a",
    "mature miRNA's expression to any individual cell. No result here",
    "proves the cellular origin of a mature miRNA, and no result here",
    "demonstrates direct regulation.** Population-composition association,",
    "shared upstream regulatory programs, and coordinated infiltration are",
    "equally consistent explanations for any composition-dependent pattern",
    "found below.",
    "",
    "## Scope",
    sprintf("- %d primary triplets (this run's own).", n_triplets),
    sprintf("- %d CPE-eligible cancer types.", n_cancers),
    sprintf("- %d MCPcounter populations, each analyzed separately (no composite immune score, no multivariable combination of overlapping populations).", n_populations),
    sprintf("- %d total (triplet, cancer, population) combinations evaluated.", n_triplet_cancer_pop),
    "",
    "## Classification counts (all triplets, cancers, populations -- no favorable selection)",
    paste(sprintf("- `%s`: %d", names(cat_counts), as.integer(cat_counts)), collapse = "\n"),
    "",
    "**Classification rule, stated explicitly:** `composition_enriched`/`composition_depleted`",
    sprintf("requires FDR<%.2f in AT LEAST ONE of the three arms (miRNA-CYT, gene-CYT,", fdr_alpha),
    "miRNA-gene), not all three, with the majority-agreeing-direction arms",
    "determining enriched-vs-depleted. A 1-vs-1 tie among significant arms,",
    "or a significant miRNA-CYT arm whose sibling arms are tautological for",
    "that population, instead falls through to `stable_across_composition`/",
    "`not_detected`/`tautological_exclusion` (the tautology check runs",
    "first, unconditionally, before significance is examined).",
    "",
    sprintf("**Of the %d `composition_enriched`/`composition_depleted` combinations above, only %d (%.1f%%) have ALL THREE arms simultaneously FDR-significant** (%d enriched, %d depleted; %.2f%% of all %d combinations) -- `tables/downstream/triplet_cell_context_all_three_arms_significant.tsv`. **A single significant arm should not be read as evidence that the complete miRNA-gene-CYT triplet is composition-dependent; the all-three-arms subset is the stricter, triplet-level result.**",
            sum(summary_dt$category %in% c("composition_enriched", "composition_depleted")),
            n_all_three_sig,
            100 * n_all_three_sig / max(1, sum(summary_dt$category %in% c("composition_enriched", "composition_depleted"))),
            sum(summary_dt$all_three_arms_fdr_significant & summary_dt$category == "composition_enriched"),
            sum(summary_dt$all_three_arms_fdr_significant & summary_dt$category == "composition_depleted"),
            100 * n_all_three_sig / nrow(summary_dt),
            nrow(summary_dt)),
    "",
    sprintf("## Strongest FDR-supported composition-dependent contexts (top %d by interaction FDR, ranked -- not cherry-picked for direction; includes combinations significant in as few as 1 of the 3 arms -- see the all-three-arms-significant table above for the stricter subset)", top_n),
    "",
    "| miRNA | Gene | Cancer | Population | Category | Interaction FDR | n |",
    "|---|---|---|---|---|---|---|"
  )
  if (top_n > 0) {
    for (i in seq_len(top_n)) {
      r <- strongest[i]
      report_lines <- c(report_lines, sprintf("| %s | %s | %s | %s | %s | %.4g | %d |",
                                                r$mirna_id, r$gene, r$cancer, r$population, r$category, r$min_interaction_fdr, r$n_samples))
    }
  } else {
    report_lines <- c(report_lines, "| (none) | | | | | | |", "", "No (triplet, cancer, population) combination reached FDR-supported composition-enriched or composition-depleted status in this run -- this negative result is retained and reported, not discarded.")
  }
  report_lines <- c(report_lines, "",
    "## Stability vs. composition-dependence",
    sprintf("%d (triplet, cancer, population) combinations show the association PRESENT but STABLE regardless of cellular composition (`stable_across_composition`) -- the baseline three-arm relationship is detected in that cancer, and no FDR-supported interaction with any tested population was found.",
            sum(summary_dt$category == "stable_across_composition")),
    sprintf("%d combinations show no association detected at all in that cancer (`not_detected`).", sum(summary_dt$category == "not_detected")),
    "",
    "## Tautology and biological-overlap exclusions",
    sprintf("%d (triplet, cancer, population) combinations were excluded as tautological -- the triplet's target gene is itself one of that population's own MCPcounter marker genes, making the population estimate partly a direct function of the gene's own expression. These are never reported as biological evidence, and never contribute to the interaction-test FDR family.",
            sum(summary_dt$category == "tautological_exclusion")),
    "",
    "## Not evaluable",
    sprintf("%d combinations could not be evaluated for at least one required arm (insufficient sample size after restricting to the population-estimate-available, CPE-purity-eligible cohort, or a non-estimable/singular model fit). These are retained and reported as `not_evaluable`, never silently dropped or conflated with a null (`not_detected`) result -- grey in every heatmap means not evaluable, never a null finding.",
            sum(summary_dt$category == "not_evaluable")),
    "",
    "## What can and cannot be concluded",
    "- A `composition_enriched`/`composition_depleted` classification means this run's data support a STATISTICAL interaction between the tested population's abundance and the triplet's association strength, at the stated FDR threshold, using bulk-tissue composition estimates -- it does NOT mean the miRNA or gene is expressed BY that cell population, and it does NOT establish a causal or mechanistic relationship.",
    "- A `stable_across_composition` classification means the association's strength does not detectably vary with the tested population's abundance in that cancer -- it does not itself validate the association as biologically meaningful beyond what the primary discovery pipeline and Stage 19's matched-null test already established (or did not establish).",
    "- Descriptive tertile-stratum correlations (full detail in `tables/downstream/triplet_cell_context_long.tsv`) are reported for interpretation only and were never used as this analysis's significance test.",
    "- Every classification, including every `not_evaluable`/`tautological_exclusion`, is reported for every triplet/cancer/population this run actually tested -- no favorable cancer or population was selected for presentation.",
    ""
  )
  writeLines(report_lines, file.path(ctx$run_root, "reports", "CANCER_SPECIFIC_CELLULAR_CONTEXT_REPORT.md"))
  log_fn("reports/CANCER_SPECIFIC_CELLULAR_CONTEXT_REPORT.md written.")

  record_stage_status(ctx, "23_cancer_specific_cellular_context", "COMPLETE")
  log_fn("Stage complete.")
  invisible(summary_dt)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=22 or source scripts/run.R in RStudio.")
}
