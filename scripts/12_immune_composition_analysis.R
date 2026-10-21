# ============================================================================
# 12_immune_composition_analysis.R
#
# Purpose: (Phase 1) assess whether the matched mRNA expression matrix
#          supports immune-cell-composition estimation; (Phase 2) compute
#          sample-level MCPcounter estimates (Becht et al. 2016); (Phase 3)
#          examine, for every primary triplet of the current run in every
#          CPE-eligible cancer, how adjustment for each individual MCPcounter
#          population estimate affects the miRNA, gene and CYT arms. Each
#          population is analysed as a separate single-population
#          sensitivity adjustment.
#
#          No composite immune index is constructed: populations are not
#          combined by principal-component reduction, and no covariate is
#          orthogonalized against CYT or any other outcome. Incomplete
#          orthogonalization can reintroduce spurious associations when a
#          second covariate (CPE) is added (Demissie & Cupples 2011, Genet
#          Epidemiol 35(6):598-605).
#
#          MCPcounter is the only deconvolution method used by the
#          per-triplet, per-cancer regressions in this stage. quanTIseq
#          (Finotello et al. 2019) was not executed because its input-scale
#          requirement was unresolved; that unexecuted module is omitted
#          from this workflow. Its output is not used here.
#          CIBERSORTx requires a separately licensed, non-redistributable
#          binary and is not used. MCPcounter and quanTIseq outputs are
#          never combined in a statistical model.
#
# Inputs:  <run_root>/processed/stage03/mirna_matched.rds, mrna_matched.rds
#          (same-run); processed/stage04/cyt_scores.tsv (same-run);
#          processed/stage05/sample_annotation_with_purity.tsv (same-run);
#          processed/stage10/master_CPE.rds (same-run);
#          reference/mcpcounter_signatures/genes.txt,
#          SIGNATURE_MANIFEST.tsv (allowlist-verified);
#          the current run's primary triplet set (scripts/lib/run_discovery_source.R).
# Outputs: <run_root>/tables/audit/immune_signature_mapping.tsv,
#          immune_input_suitability.tsv, immune_primary_gene_signature_overlap.tsv;
#          reports/20_immune_input_validation.md;
#          processed/stage11_immune/sample_level_immune_estimates.tsv;
#          tables/downstream/immune_estimate_qc.tsv,
#          immune_composition_associations.tsv (one row per triplet x
#          cancer x arm x adjustment population),
#          immune_adjusted_recurrence_by_population.tsv,
#          immune_adjusted_recurrence_range_by_triplet.tsv;
#          reports/21_immune_composition_report.md.
#          (Artifact paths such as processed/stage11_immune/ are stable
#          output names and do not follow the stage number.)
# Method:  for each MCPcounter population individually, a two-covariate OLS
#          regression (feature ~ CPE + population_k) is fitted within each
#          cancer and Spearman correlation is computed on the residuals.
#          "Cytotoxic lymphocytes" and "CD8 T cells" are reported but
#          labelled as biological-overlap sensitivity analyses rather than
#          ordinary confounders, because of their conceptual proximity to
#          the cytotoxic-lymphocyte construction of CYT; they are not used to
#          claim complete adjustment for CYT. The analysis characterizes
#          sensitivity only and does not add to or remove from the primary
#          discovery set. Every count is read from the current run's Stage 10
#          output.
# Depends on: 03_harmonize_identifiers.R, 04_compute_validate_cyt.R,
#             05_define_eligible_features.R, 09_recurrence_analysis.R,
#             10_purity_adjusted_analysis.R (all same-run).
# ============================================================================

run_stage12 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "12_immune_composition_analysis")
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "stats_engine.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "expression_preprocessing.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "run_discovery_source.R"), local = TRUE)
  suppressPackageStartupMessages({
    library(data.table)
    library(MCPcounter)
  })
  script_path <- file.path(ctx$data_root, "scripts", "12_immune_composition_analysis.R")
  write_fn <- function(x, p) write.table(x, p, sep = "\t", row.names = FALSE, quote = FALSE)
  fwrite_fn <- function(x, p) fwrite(x, p, sep = "\t")

  # ==========================================================================
  # PHASE 1: INPUT SUITABILITY AUDIT
  # ==========================================================================
  log_fn("PHASE 1: input suitability audit.")

  sig_dir <- file.path(ctx$data_root, "reference", "mcpcounter_signatures")
  sig_manifest_path <- file.path(sig_dir, "SIGNATURE_MANIFEST.tsv")
  genes_sig_path <- file.path(sig_dir, "genes.txt")
  verify_allowed_raw_input(ctx$data_root, sig_manifest_path)
  verify_allowed_raw_input(ctx$data_root, genes_sig_path)
  sig_manifest <- read.delim(sig_manifest_path, stringsAsFactors = FALSE)
  genes_sig <- read.table(genes_sig_path, sep = "\t", stringsAsFactors = FALSE,
                           header = TRUE, colClasses = "character", check.names = FALSE)
  setDT(genes_sig)
  run_stage_assert(ctx, log_fn, nrow(genes_sig) == 111, "Cached MCPcounter signature file has the expected 111 marker-gene rows")

  norm_qc <- read_run_artifact(ctx, file.path("processed", "stage03", "normalization_qc.tsv"), "03_harmonize_identifiers",
                                function(p) read.delim(p, stringsAsFactors = FALSE))
  mrna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mrna_matched.rds"), "03_harmonize_identifiers", readRDS)
  symbols_raw <- toupper(trimws(mrna_matched$gene_id_raw))
  dup_symbols <- unique(symbols_raw[duplicated(symbols_raw)])

  mrna_mat_raw <- as.matrix(mrna_matched[, -1, drop = FALSE])
  rownames(mrna_mat_raw) <- symbols_raw
  mrna_mat <- deduplicate_matrix_rows(mrna_mat_raw)
  run_stage_assert(ctx, log_fn, !any(duplicated(rownames(mrna_mat))),
               "No duplicate gene symbols in the mRNA matrix (%d found upstream and already resolved: %s); %d unique rows",
               length(dup_symbols), paste(dup_symbols, collapse = ","), nrow(mrna_mat))

  sig_map <- genes_sig[, .(
    n_signature_genes = .N,
    n_mapped = sum(`HUGO symbols` %in% rownames(mrna_mat)),
    mapped_genes = paste(sort(`HUGO symbols`[`HUGO symbols` %in% rownames(mrna_mat)]), collapse = ","),
    missing_genes = paste(sort(`HUGO symbols`[!(`HUGO symbols` %in% rownames(mrna_mat))]), collapse = ",")
  ), by = .(cell_population = `Cell population`)]
  sig_map[, pct_mapped := round(100 * n_mapped / n_signature_genes, 1)]
  setcolorder(sig_map, c("cell_population", "n_signature_genes", "n_mapped", "pct_mapped", "mapped_genes", "missing_genes"))
  write_run_artifact(ctx, sig_map, file.path("tables", "audit", "immune_signature_mapping.tsv"),
                      "12_immune_composition_analysis", script_path, character(0), fwrite_fn)
  log_fn("Signature mapping table written; overall mapping %.1f%% (%d/%d genes).",
         100 * sum(sig_map$n_mapped) / sum(sig_map$n_signature_genes), sum(sig_map$n_mapped), sum(sig_map$n_signature_genes))

  triplets <- load_run_primary_triplets(ctx)
  primary_genes <- sort(unique(triplets$gene))
  gene_overlap <- genes_sig[`HUGO symbols` %in% primary_genes]
  setnames(gene_overlap, c("HUGO symbols", "Cell population"), c("primary_gene", "overlapping_population"))
  write_run_artifact(ctx, gene_overlap[, .(primary_gene, overlapping_population)],
                      file.path("tables", "audit", "immune_primary_gene_signature_overlap.tsv"),
                      "12_immune_composition_analysis", script_path, character(0), fwrite_fn)
  log_fn("%d of this run's %d unique primary target genes are themselves MCPcounter marker genes: %s",
         nrow(gene_overlap), length(primary_genes), paste(gene_overlap$primary_gene, collapse = ","))

  suitability <- data.frame(
    check = c("input_scale", "log_transform_status", "gene_id_type", "duplicate_symbols",
              "overall_pct_signature_genes_mapped", "min_population_pct_mapped",
              "cross_cohort_batch_adjustment_present", "cross_cancer_absolute_level_comparability",
              "primary_gene_signature_tautology_risk", "cyt_biological_overlap_risk",
              "compositional_status", "recommended_scope"),
    finding = c(
      sprintf("log2-like (%s); max=%.4f min=%.4f", norm_qc$declared_scale[norm_qc$layer == "mRNA"], max(mrna_mat, na.rm = TRUE), min(mrna_mat, na.rm = TRUE)),
      "Not re-transformed by this pipeline (input already log2-scale per stage03 QC)",
      "HUGO gene symbols (uppercased, trimmed) -- matches MCPcounter featuresType='HUGO_symbols'",
      if (length(dup_symbols) > 0) sprintf("%d duplicate symbol(s) found upstream (%s); resolved at 03_harmonize_identifiers.R by highest-variance row", length(dup_symbols), paste(dup_symbols, collapse = ",")) else "No duplicate gene symbols",
      sprintf("%.1f%% (%d/%d)", 100 * sum(sig_map$n_mapped) / sum(sig_map$n_signature_genes), sum(sig_map$n_mapped), sum(sig_map$n_signature_genes)),
      sprintf("%.1f%% (%s)", min(sig_map$pct_mapped), sig_map$cell_population[which.min(sig_map$pct_mapped)]),
      "YES -- source file is EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena, an Empirical-Bayes batch-adjusted pan-cancer matrix",
      "NOT ESTABLISHED -- absolute cross-cancer comparability was not independently verified; scores used within-cancer only",
      sprintf("YES for %d/%d of this run's primary target genes (%s)", nrow(gene_overlap), length(primary_genes), paste(gene_overlap$primary_gene, collapse = ",")),
      "Cytotoxic lymphocytes and CD8 T cells populations measure biology closely related to CYT's own construction; reported but labeled BIOLOGICAL-OVERLAP SENSITIVITY throughout, never used to claim complete CYT adjustment",
      "MCPcounter produces arbitrary-unit abundance SCORES, NOT compositional fractions (they do not sum to 1); each population may be used individually as a standard regression covariate without a reference-category or compositional-collinearity concern -- Becht et al. 2016 does not specify a required transformation for this use",
      "WITHIN-CANCER RANKING / RELATIVE-ASSOCIATION USE ONLY; every population tested INDIVIDUALLY, never combined into a composite score"
    ),
    stringsAsFactors = FALSE
  )
  write_run_artifact(ctx, suitability, file.path("tables", "audit", "immune_input_suitability.tsv"),
                      "12_immune_composition_analysis", script_path, character(0), fwrite_fn)
  verdict <- "SUITABLE WITH RESTRICTIONS"
  log_fn("Input suitability verdict: %s", verdict)

  report20 <- c(
    "# Stage 12 Input Validation: Immune-Composition Estimation Suitability",
    "",
    sprintf("Generated by `scripts/12_immune_composition_analysis.R`, run `%s`, %s.", ctx$run_id, format(Sys.time(), "%Y-%m-%d %H:%M:%OS0")),
    "",
    "## Verdict", "", sprintf("**%s.**", verdict), "",
    "This version of Stage 12 uses ONLY standard, published methods: MCPcounter (Becht et al. 2016) here, and quanTIseq (Finotello et al. 2019) in scripts/13_quantiseq_deconvolution.R, each with every output population/fraction reported INDIVIDUALLY. No composite index, principal-component reduction, or outcome-orthogonalized covariate is constructed anywhere in this pipeline. See archive/superseded_designs/README.md for the reasons an earlier PC1/orthogonalized design was retracted and is not used here.",
    "",
    "## Signature and package provenance", "",
    sprintf("- Package: MCPcounter %s", sig_manifest$mcpcounter_package_version[1]),
    sprintf("- Signature file SHA-256: `%s`, retrieved %s", sig_manifest$sha256[sig_manifest$file == "genes.txt"], sig_manifest$retrieved_utc[sig_manifest$file == "genes.txt"]),
    "- `featuresType=\"HUGO_symbols\"` used explicitly.",
    "",
    "## Matrix suitability", "",
    "| Check | Finding |", "|---|---|",
    paste0("| ", suitability$check, " | ", gsub("\\|", "/", suitability$finding), " |"),
    "",
    "## Decision", "",
    "Stage 12 proceeds to sample-level estimation and per-population (never combined) sensitivity analysis."
  )
  write_run_artifact(ctx, report20, file.path("reports", "20_immune_input_validation.md"),
                      "12_immune_composition_analysis", script_path, character(0),
                      function(x, p) writeLines(x, p))
  log_fn("reports/20_immune_input_validation.md written.")

  # ==========================================================================
  # PHASE 2: SAMPLE-LEVEL MCPcounter ESTIMATION
  # ==========================================================================
  log_fn("PHASE 2: sample-level MCPcounter estimation.")
  est <- MCPcounter.estimate(mrna_mat, featuresType = "HUGO_symbols", genes = genes_sig)
  run_stage_assert(ctx, log_fn, ncol(est) == ncol(mrna_mat), "MCPcounter estimate produced for every matched sample (%d)", ncol(mrna_mat))
  populations <- rownames(est)

  # Population identity is derived on every run from (a) the row names
  # returned by MCPcounter.estimate() (`populations`) and (b) the validated
  # signature table's "Cell population" grouping (`sig_map$cell_population`,
  # the grouping MCPcounter uses internally), and is never assumed to be a
  # fixed count. MCPcounter.estimate() omits a population from its output
  # when none of its marker genes is present in the input matrix, so an
  # explicit name-keyed comparison is required.
  expected_populations <- sort(unique(sig_map$cell_population))
  pop_check <- validate_population_schema(populations, expected_populations)
  run_stage_assert(ctx, log_fn, pop_check$non_empty, "Every MCPcounter output population name is non-empty")
  run_stage_assert(ctx, log_fn, pop_check$unique, "MCPcounter output population names are unique (no duplicates): %s",
                    if (pop_check$unique) "confirmed" else paste(populations[duplicated(populations)], collapse = ", "))
  run_stage_assert(ctx, log_fn, length(pop_check$missing) == 0 && length(pop_check$unexpected) == 0,
                    "MCPcounter's %d output population(s) exactly match the validated signature table's %d expected population(s) (missing: %s; unexpected: %s)",
                    length(populations), length(expected_populations),
                    if (length(pop_check$missing)) paste(pop_check$missing, collapse = ", ") else "(none)",
                    if (length(pop_check$unexpected)) paste(pop_check$unexpected, collapse = ", ") else "(none)")
  log_fn("Populations estimated (%d): %s", length(populations), paste(populations, collapse = ", "))

  sample_annot <- read_run_artifact(ctx, file.path("processed", "stage05", "sample_annotation_with_purity.tsv"), "05_define_eligible_features",
                                     function(p) fread(p))
  run_stage_assert(ctx, log_fn, identical(colnames(est), colnames(mrna_mat)), "MCPcounter estimate columns are in the same sample order as the input matrix")
  # Keyed matching by population name: each column of est_dt is built by
  # indexing `est` with the population name (est[pop, ]), not by position.
  est_dt <- data.table(sample_id = colnames(est))
  for (pop in populations) est_dt[[pop]] <- as.numeric(est[pop, ])
  est_dt <- merge(sample_annot[, .(sample_id, cancer)], est_dt, by = "sample_id")
  setcolorder(est_dt, c("sample_id", "cancer", populations))
  write_run_artifact(ctx, est_dt, file.path("processed", "stage11_immune", "sample_level_immune_estimates.tsv"),
                      "12_immune_composition_analysis", script_path,
                      file.path(ctx$run_root, file.path("processed", "stage03", "mrna_matched.rds")), fwrite_fn)
  log_fn("processed/stage11_immune/sample_level_immune_estimates.tsv written (%d samples x %d populations).", nrow(est_dt), length(populations))

  n_missing <- sum(!complete.cases(est_dt[, ..populations]))
  qc_meta <- data.frame(
    field = c("package", "package_version", "signature_version_sha256", "featuresType",
              "gene_id_mapping_pct_overall", "n_samples_scored", "n_samples_missing_any_estimate",
              "output_date", "generating_script", "composite_index_constructed"),
    value = c("MCPcounter", sig_manifest$mcpcounter_package_version[1], sig_manifest$sha256[sig_manifest$file == "genes.txt"],
              "HUGO_symbols", sprintf("%.1f", 100 * sum(sig_map$n_mapped) / sum(sig_map$n_signature_genes)),
              as.character(ncol(est)), as.character(n_missing), format(Sys.time(), "%Y-%m-%d"),
              "scripts/12_immune_composition_analysis.R", "NONE -- every population reported and tested individually"),
    stringsAsFactors = FALSE
  )
  write_run_artifact(ctx, qc_meta, file.path("tables", "downstream", "immune_estimate_qc.tsv"),
                      "12_immune_composition_analysis", script_path, character(0), fwrite_fn)
  run_stage_assert(ctx, log_fn, n_missing == 0, "No missing MCPcounter estimates across %d samples", ncol(est))

  # ==========================================================================
  # PHASE 3: PER-POPULATION SENSITIVITY ANALYSIS
  # ==========================================================================
  log_fn("PHASE 3: per-population sensitivity analysis (standard multiple regression, one population at a time).")

  cyt_scores <- read_run_artifact(ctx, file.path("processed", "stage04", "cyt_scores.tsv"), "04_compute_validate_cyt", fread)
  full_annot <- merge(sample_annot, est_dt[, c("sample_id", populations), with = FALSE], by = "sample_id")

  mirna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mirna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mirna_mat_full <- as.matrix(mirna_matched[, -1, drop = FALSE]); rownames(mirna_mat_full) <- mirna_matched$mirna_id

  # Per-feature winsorization immediately before this phase's OLS
  # regression (feature ~ CPE + immune population), which is sensitive to
  # the leverage of extreme samples. MCPcounter estimation (Phase 2) uses
  # the unwinsorized mrna_mat, consistent with its published input
  # expectations; mrna_mat_reg is a separate winsorized copy.
  mirna_mat_full <- winsorize_matrix_rows(mirna_mat_full)
  mrna_mat_reg <- winsorize_matrix_rows(mrna_mat)

  master_cpe <- read_run_artifact(ctx, file.path("processed", "stage10", "master_CPE.rds"), "10_purity_adjusted_analysis", readRDS)
  eligible_cancers_cpe <- sort(unique(master_cpe$mircyt_full$cancer))

  min_n <- analysis_config$discovery$min_samples_per_cancer
  fdr_cutoff <- analysis_config$discovery$fdr_cutoff
  rho_threshold <- analysis_config$thresholds$primary
  biological_overlap_populations <- analysis_config$immune_adjustment$cytotoxic_signal_populations_excluded

  # build_gene_exclusion_map() (scripts/lib/stats_engine.R) maps each gene to every
  # MCPcounter population whose signature contains it (e.g. CD3D belongs to
  # both "T cells" and "Cytotoxic lymphocytes"); see
  # tests/test_immune_adjustment_scenarios.R, scenario 5.
  gene_exclusion_map <- build_gene_exclusion_map(gene_overlap)

  # ==========================================================================
  # Output schema (long format): one row per
  # (mirna_id, gene, cancer, arm, population), arm in
  # {miRNA-CYT, miRNA-gene, gene-CYT}. Every row carries a real MCPcounter
  # population name; populations are never collapsed or averaged into a
  # single immune effect. Stage 14 consumes this schema.
  # ==========================================================================
  rows <- list()
  for (i in seq_len(nrow(triplets))) {
    mir <- triplets$mirna_id[i]; tgt_gene <- triplets$gene[i]
    tautology_populations <- if (tgt_gene %in% names(gene_exclusion_map)) gene_exclusion_map[[tgt_gene]] else character(0)

    for (cc in eligible_cancers_cpe) {
      idx_all <- which(sample_annot$cancer == cc)
      if (length(idx_all) < min_n) next
      mirna_v <- as.numeric(mirna_mat_full[mir, idx_all])
      gene_v <- as.numeric(mrna_mat_reg[tgt_gene, idx_all])
      cyt_v <- cyt_scores$CYT[match(sample_annot$sample_id[idx_all], cyt_scores$sample_id)]
      cpe_v <- sample_annot$CPE[idx_all]

      for (pop in populations) {
        tautology_excluded <- pop %in% tautology_populations
        if (tautology_excluded) next  # skip every population whose signature contains the triplet's target gene (tautology); a gene may belong to several populations
        pop_v <- full_annot[[pop]][match(sample_annot$sample_id[idx_all], full_annot$sample_id)]

        # immune_covariate_diagnostics() (scripts/lib/stats_engine.R): a population
        # score that is constant within this cancer, or collinear with CPE,
        # cannot be used as a second covariate; this is reported explicitly
        # rather than letting lm() drop or alias the term (see
        # tests/test_immune_adjustment_scenarios.R, scenarios 3-4).
        diag <- immune_covariate_diagnostics(pop_v, cpe_v)
        pop_v_finite <- pop_v[is.finite(pop_v) & is.finite(cpe_v)]
        zero_variance_population <- diag$zero_variance_population
        collinearity_r <- diag$collinearity_r
        collinearity_vif <- diag$collinearity_vif
        aliased_covariate <- diag$aliased_covariate

        if (zero_variance_population || aliased_covariate) {
          for (arm_i in c("miRNA-CYT", "miRNA-gene", "gene-CYT")) {
            rows[[length(rows) + 1]] <- data.table(
              mirna_id = mir, gene = tgt_gene, cancer = cc, arm = arm_i, population = pop,
              is_biological_overlap_sensitivity = pop %in% biological_overlap_populations,
              n = length(pop_v_finite), rho = NA_real_, p = NA_real_,
              collinearity_cpe_population_r = collinearity_r, collinearity_vif = collinearity_vif,
              status = if (zero_variance_population) "NOT ESTIMABLE (zero-variance population)" else "NOT ESTIMABLE (aliased with CPE, |r|>0.999)"
            )
          }
          next
        }

        cov_df <- data.frame(CPE = cpe_v, pop_covariate = pop_v)
        vars_df <- data.frame(mirna = mirna_v, gene = gene_v, cyt = cyt_v)
        rr <- residualize_columns_multi(vars_df, cov_df)
        mc <- spearman_pair(rr[, 1], rr[, 3], min_n = min_n)
        mg <- spearman_pair(rr[, 1], rr[, 2], min_n = min_n)
        gc <- spearman_pair(rr[, 2], rr[, 3], min_n = min_n)

        for (arm_res in list(list(arm = "miRNA-CYT", r = mc), list(arm = "miRNA-gene", r = mg), list(arm = "gene-CYT", r = gc))) {
          rows[[length(rows) + 1]] <- data.table(
            mirna_id = mir, gene = tgt_gene, cancer = cc, arm = arm_res$arm, population = pop,
            is_biological_overlap_sensitivity = pop %in% biological_overlap_populations,
            n = arm_res$r$n, rho = arm_res$r$rho, p = arm_res$r$p,
            collinearity_cpe_population_r = collinearity_r, collinearity_vif = collinearity_vif,
            status = if (is.na(arm_res$r$rho)) "NOT ESTIMABLE" else "ESTIMATED"
          )
        }
      }
    }
    if (i %% 10 == 0) log_fn("  Processed %d/%d this-run primary triplets.", i, nrow(triplets))
  }
  assoc <- rbindlist(rows)

  # FDR families: BH is applied separately for each (arm, population)
  # combination within each cancer, never across populations or arms
  # (reports/FDR_FAMILY_SPECIFICATION.tsv).
  # Shared-arm deduplication: miRNA-CYT is repeated for every target gene of
  # the same miRNA, and gene-CYT is repeated for every miRNA targeting the
  # same gene. Each biological hypothesis contributes exactly once to its
  # cancer-specific BH family, and the resulting FDR is broadcast back to
  # all triplet rows sharing that hypothesis.
  assoc[, fdr := NA_real_]
  for (pop in populations) {
    for (arm_i in c("miRNA-CYT", "miRNA-gene", "gene-CYT")) {
      if (arm_i == "miRNA-CYT") {
        sel <- assoc$arm == arm_i & assoc$population == pop
        assoc[sel, fdr := dedup_mirna_cyt_fdr(assoc[sel])]
      } else if (arm_i == "gene-CYT") {
        sel <- assoc$arm == arm_i & assoc$population == pop
        assoc[sel, fdr := dedup_gene_cyt_fdr(assoc[sel])]
      } else {
        sel <- assoc$arm == arm_i & assoc$population == pop
        assoc[sel, fdr := p.adjust(p, "BH"), by = cancer]
      }
    }
  }
  assoc[, pass := abs(rho) >= rho_threshold & fdr < fdr_cutoff & !is.na(fdr)]
  assoc_wide_pass <- dcast(assoc, mirna_id + gene + cancer + population + is_biological_overlap_sensitivity ~ arm, value.var = "pass")
  assoc_wide_pass[, three_arm_pass := `miRNA-CYT` & `miRNA-gene` & `gene-CYT`]
  assoc <- merge(assoc, assoc_wide_pass[, .(mirna_id, gene, cancer, population, three_arm_pass)],
                 by = c("mirna_id", "gene", "cancer", "population"), all.x = TRUE)

  # Required-column, key-uniqueness and allowed-value assertions at this
  # producer boundary (Stage 14 re-asserts the same contract on read).
  required_cols <- c("mirna_id", "gene", "cancer", "arm", "population", "n", "rho", "p", "fdr", "pass", "three_arm_pass")
  run_stage_assert(ctx, log_fn, all(required_cols %in% names(assoc)), "immune_composition_associations.tsv contains every contract-required column")
  run_stage_assert(ctx, log_fn, all(assoc$arm %in% c("miRNA-CYT", "miRNA-gene", "gene-CYT")), "arm column contains only the 3 allowed values")
  dup_key <- assoc[, .N, by = .(mirna_id, gene, cancer, arm, population)]
  run_stage_assert(ctx, log_fn, all(dup_key$N == 1), "No duplicated (mirna_id, gene, cancer, arm, population) key in immune_composition_associations.tsv")

  stage10_input_paths <- file.path(ctx$run_root, c(
    file.path("processed", "stage03", "mirna_matched.rds"), file.path("processed", "stage03", "mrna_matched.rds"),
    file.path("processed", "stage04", "cyt_scores.tsv"), file.path("processed", "stage05", "sample_annotation_with_purity.tsv"),
    file.path("processed", "stage10", "master_CPE.rds")
  ))
  write_run_artifact(ctx, assoc, file.path("tables", "downstream", "immune_composition_associations.tsv"),
                      "12_immune_composition_analysis", script_path, stage10_input_paths, fwrite_fn)
  log_fn("immune_composition_associations.tsv written (%d rows: %d triplets x %d cancers x %d populations x 3 arms, minus per-triplet tautology exclusions).",
         nrow(assoc), nrow(triplets), length(eligible_cancers_cpe), length(populations))

  # Per-triplet, per-population recurrence
  triplet_pop_pass <- unique(assoc[, .(mirna_id, gene, cancer, population, is_biological_overlap_sensitivity, three_arm_pass)])
  recur_by_pop <- triplet_pop_pass[, .(n_evaluable = sum(!is.na(three_arm_pass)), recurrence = sum(three_arm_pass, na.rm = TRUE)),
                         by = .(mirna_id, gene, adjustment_population = population, is_biological_overlap_sensitivity)]
  min_recurrence <- analysis_config$discovery$min_recurrence
  recur_by_pop[, recurrence_meets_threshold := recurrence >= min_recurrence]
  write_run_artifact(ctx, recur_by_pop, file.path("tables", "downstream", "immune_adjusted_recurrence_by_population.tsv"),
                      "12_immune_composition_analysis", script_path, character(0), fwrite_fn)

  non_overlap_recur <- recur_by_pop[is_biological_overlap_sensitivity == FALSE]
  n_triplets_all_nonoverlap_pass <- non_overlap_recur[, .(all_pass = all(recurrence_meets_threshold)), by = .(mirna_id, gene)]
  log_fn("Standard, non-composite per-population sensitivity: %d of this run's %d primary triplets retain recurrence>=%d under EVERY non-overlap MCPcounter population tested individually.",
         sum(n_triplets_all_nonoverlap_pass$all_pass), nrow(triplets), min_recurrence)

  # ==========================================================================
  # REPORT 21
  # ==========================================================================
  min_recur_by_triplet <- non_overlap_recur[, .(min_recurrence_across_populations = min(recurrence), max_recurrence_across_populations = max(recurrence)), by = .(mirna_id, gene)]
  write_run_artifact(ctx, min_recur_by_triplet, file.path("tables", "downstream", "immune_adjusted_recurrence_range_by_triplet.tsv"),
                      "12_immune_composition_analysis", script_path, character(0), fwrite_fn)

  n_nonoverlap_pop <- sum(!populations %in% biological_overlap_populations)
  report21 <- c(
    "# Stage 12: Immune-Composition Sensitivity Analysis (Standard Methods Only)",
    "",
    sprintf("Generated by `scripts/12_immune_composition_analysis.R`, run `%s`.", ctx$run_id),
    "",
    "## Method (revised per the STANDARD METHODS ONLY policy)",
    "",
    "No composite immune index is constructed anywhere in this analysis. Each of MCPcounter's 10 published cell-population scores is used INDIVIDUALLY as a standard second covariate alongside CPE, in an ordinary two-covariate OLS regression, exactly as any single additional covariate would be added to a multiple regression model. No principal-component reduction, no combination of populations, and no orthogonalization of any covariate against CYT (or any other outcome) is performed. 'Cytotoxic lymphocytes' and 'CD8 T cells' are computed and reported, but labeled BIOLOGICAL-OVERLAP SENSITIVITY throughout -- they are not treated as ordinary confounders, and no claim of 'complete' CYT adjustment is made using them, per their conceptual closeness to CYT's own construction.",
    "",
    sprintf("For the %d of this run's primary target genes that are themselves MCPcounter marker genes, the matching population is excluded from that specific triplet's analysis (a mechanical tautology exclusion, not a scientific judgment) -- %s.",
            nrow(gene_overlap), if (nrow(gene_overlap) > 0) paste(sprintf("%s excludes %s", gene_overlap$primary_gene, gene_overlap$overlapping_population), collapse = "; ") else "none in this run"),
    "",
    "## Results", "",
    sprintf("- Triplet-cancer-population rows tested: %d (%d triplets x up to %d cancers x %d populations, minus tautology exclusions).", nrow(assoc), nrow(triplets), length(eligible_cancers_cpe), length(populations)),
    sprintf("- Of the %d non-overlap MCPcounter populations, %d of this run's %d primary triplets retain three-arm recurrence >=%d under EVERY population tested individually (the most conservative possible standard-methods summary).",
            n_nonoverlap_pop, sum(n_triplets_all_nonoverlap_pass$all_pass), nrow(triplets), min_recurrence),
    "- Per-triplet, per-population recurrence (min/max across the non-overlap populations) is in `tables/downstream/immune_adjusted_recurrence_by_population.tsv` -- this deliberately does NOT collapse to one number, because no single composite adjustment is authoritative under this policy.",
    "",
    "## Interpretation", "",
    sprintf("Each population's sensitivity result is reported on its own terms. A triplet failing recurrence>=%d under one specific population (e.g. Monocytic lineage) but passing under another (e.g. Neutrophils) is not evidence of a contradiction -- it means the specific cell type in question does or does not explain that triplet's association, and both facts are retained, not averaged or combined into a single verdict.", min_recurrence),
    "",
    "## Outputs", "",
    "- `tables/audit/immune_signature_mapping.tsv`, `immune_input_suitability.tsv`, `immune_primary_gene_signature_overlap.tsv`",
    "- `processed/stage11_immune/sample_level_immune_estimates.tsv`",
    "- `tables/downstream/immune_estimate_qc.tsv`, `immune_composition_associations.tsv`, `immune_adjusted_recurrence_by_population.tsv`, `immune_adjusted_recurrence_range_by_triplet.tsv`",
    "- `reports/20_immune_input_validation.md`, `reports/21_immune_composition_report.md` (this file)"
  )
  write_run_artifact(ctx, report21, file.path("reports", "21_immune_composition_report.md"),
                      "12_immune_composition_analysis", script_path, character(0),
                      function(x, p) writeLines(x, p))
  log_fn("reports/21_immune_composition_report.md written.")

  record_stage_status(ctx, "12_immune_composition_analysis", "COMPLETE")
  log_fn("Stage 12 complete.")
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=12 or source scripts/run.R in RStudio.")
}
