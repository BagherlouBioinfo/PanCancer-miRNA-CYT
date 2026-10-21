# ============================================================================
# 14_cross_cancer_meta_analysis.R
#
# Purpose: synthesize the cancer-specific evidence for each primary triplet
#          of the current run across every cancer type with adequate matched
#          data by random-effects meta-analysis. The synthesis universe is
#          not restricted to cancer-triplet records that pass the primary
#          discovery rule: eligibility depends only on prespecified data
#          adequacy (matched-sample minimum, purity-metric availability),
#          never on correlation magnitude, FDR or direction. Records that
#          pass the primary three-arm rule are flagged in a
#          `primary_support` column rather than used as a filter.
# Inputs:  the current run's primary triplet set (scripts/lib/run_discovery_source.R,
#          <run_root>/tables/primary/); <run_root>/tables/downstream/
#          immune_composition_associations.tsv (Stage 12: unadjusted,
#          CPE-adjusted and immune-adjusted per-cancer effects for the
#          primary triplets); processed/stage06/master_mirna_cyt.tsv;
#          processed/stage08/three_arm_threshold_*.tsv;
#          processed/stage10/master_CPE.rds and master_LUMP.rds;
#          processed/stage05/eligible_cancers_cpe.tsv and
#          eligible_cancers_lump.tsv; tables/audit/
#          cancer_level_evidence_CPE_rho030.tsv (source of the
#          primary_support flag). All inputs are same-run artifacts.
# Outputs: <run_root>/tables/downstream/meta_analysis_input_all_cancers.tsv,
#          meta_analysis_all_arms.tsv, meta_analysis_heterogeneity.tsv,
#          meta_analysis_influence.tsv, meta_analysis_leave_one_out.tsv,
#          meta_analysis_diagnostic_failures.tsv;
#          figures/main_candidates/meta_analysis/*.png;
#          reports/22_cross_cancer_meta_analysis_report.md.
# Method:  each cancer-specific Spearman rho is Fisher-z transformed
#          (z = atanh(rho)) with standard error SE_z = sqrt(1.06 / (n - 3))
#          (the Fieller, Hartley & Pearson 1957 variance correction for
#          Spearman correlations; scripts/lib/stats_engine.R::fisher_z_se_spearman())
#          and pooled by REML random-effects meta-analysis (metafor::rma)
#          separately for every (triplet, arm, adjustment) combination. A
#          model is fitted only where at least
#          analysis_config$meta_analysis$min_cancers_per_model cancers have a
#          valid estimate (a prespecified minimum independent of
#          significance). Pooled estimates are back-transformed to the rho
#          scale with metafor's transf.ztor. FDR is controlled separately
#          within each (arm, adjustment[, population]) family across the
#          run's triplets; these families are independent of the
#          discovery-stage FDR families of Stages 06, 08 and 10.
#          No expected triplet or cancer count is referenced; every count is
#          read from the current run's output.
# Depends on: 03-10 (same-run discovery pipeline),
#             12_immune_composition_analysis.R (same-run).
# ============================================================================

run_stage14 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "14_cross_cancer_meta_analysis")
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "run_discovery_source.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "relational_integrity.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "stats_engine.R"), local = TRUE)
  suppressPackageStartupMessages({
    library(data.table)
    library(metafor)
  })
  script_path <- file.path(ctx$data_root, "scripts", "13_cross_cancer_meta_analysis.R")
  fwrite_fn <- function(x, p) fwrite(x, p, sep = "\t")

  triplets <- load_run_primary_triplets(ctx)

  # ==========================================================================
  # CPE and LUMP eligible-cancer sets are read independently from their own
  # Stage 05 outputs and are never treated as one shared condition.
  # ==========================================================================
  eligible_cancers_cpe <- sort(read_run_artifact(ctx, file.path("processed", "stage05", "eligible_cancers_cpe.tsv"),
                                                  "05_define_eligible_features",
                                                  function(p) read.delim(p, stringsAsFactors = FALSE)$cancer))
  eligible_cancers_lump <- sort(read_run_artifact(ctx, file.path("processed", "stage05", "eligible_cancers_lump.tsv"),
                                                   "05_define_eligible_features",
                                                   function(p) read.delim(p, stringsAsFactors = FALSE)$cancer))
  run_stage_assert(ctx, log_fn, length(eligible_cancers_cpe) >= 1, "%d cancers have valid CPE purity data in this run", length(eligible_cancers_cpe))
  run_stage_assert(ctx, log_fn, length(eligible_cancers_lump) >= 1, "%d cancers have valid LUMP purity data in this run", length(eligible_cancers_lump))
  log_fn("CPE-eligible cancers: %d. LUMP-eligible cancers: %d. Sets identical: %s.",
         length(eligible_cancers_cpe), length(eligible_cancers_lump), setequal(eligible_cancers_cpe, eligible_cancers_lump))

  # ==========================================================================
  # Long-format meta-analysis input: one row per
  # (mirna_id, gene, cancer, arm, adjustment[, population]).
  #  - unadjusted: Stage 06 (miRNA-CYT arm, full eligible-miRNA family) and
  #    Stage 08 (miRNA-gene and gene-CYT arms, primary threshold).
  #  - CPE, LUMP: Stage 10's master_{CPE,LUMP}.rds (mircyt_full, mirgene,
  #    genecyt).
  #  - immune: Stage 12's immune_composition_associations.tsv, with
  #    population identity preserved; every MCPcounter population is a
  #    separate (adjustment = "immune", population = <name>) meta-analysis.
  # ==========================================================================
  master_mircyt_unadj <- read_run_artifact(ctx, file.path("processed", "stage06", "master_mirna_cyt.tsv"),
                                            "06_mirna_cyt_discovery", fread)
  primary_thr_tag <- sprintf("rho%03d", round(analysis_config$thresholds$primary * 100))
  three_arm_unadj <- read_run_artifact(ctx, file.path("processed", "stage08", sprintf("three_arm_threshold_%.2f.tsv", analysis_config$thresholds$primary)),
                                        "08_same_cancer_three_arm_testing", fread)
  master_cpe <- read_run_artifact(ctx, file.path("processed", "stage10", "master_CPE.rds"), "10_purity_adjusted_analysis", readRDS)
  master_lump <- read_run_artifact(ctx, file.path("processed", "stage10", "master_LUMP.rds"), "10_purity_adjusted_analysis", readRDS)
  immune_assoc <- read_run_artifact(ctx, file.path("tables", "downstream", "immune_composition_associations.tsv"),
                                     "12_immune_composition_analysis", fread)
  # Required-column contract check on the Stage 12 input (consumer-side
  # counterpart of the producer-side assertion in Stage 12).
  run_stage_assert(ctx, log_fn, all(c("mirna_id", "gene", "cancer", "arm", "population", "n", "rho", "p") %in% names(immune_assoc)),
                    "immune_composition_associations.tsv (Stage 12) contains every contract-required column this stage consumes")
  primary_evidence <- load_run_evidence(ctx, "CPE")

  fisher_z_correction <- analysis_config$meta_analysis$fisher_z_spearman_correction
  # Variance-scale Fisher-z correction for Spearman correlations,
  # SE = sqrt(correction / (n - 3)) (scripts/lib/stats_engine.R::fisher_z_se_spearman()).
  se_fisher <- function(n) fisher_z_se_spearman(n, fisher_z_correction)
  min_n_cfg <- analysis_config$discovery$min_samples_per_cancer

  finalize <- function(dt, adj_label, requires_purity, eligible_cancer_set) {
    dt[, adjustment := adj_label]
    dt[, eligible := !is.na(rho) & !is.na(n) & n >= min_n_cfg]
    dt[, exclusion_reason := NA_character_]
    if (requires_purity) {
      dt[!eligible & !(cancer %in% eligible_cancer_set), exclusion_reason := sprintf("no valid %s purity estimate available for this cancer type", adj_label)]
      dt[!eligible & (cancer %in% eligible_cancer_set) & adj_label == "immune", exclusion_reason := "immune-composition covariate could not be estimated or residualization failed (zero variance, aliasing with CPE, or insufficient complete observations)"]
      dt[!eligible & (cancer %in% eligible_cancer_set) & adj_label != "immune", exclusion_reason := "insufficient complete observations or zero variance after purity residualization"]
    } else {
      dt[eligible == FALSE, exclusion_reason := "insufficient complete observations or zero variance in this cancer"]
    }
    dt[, z := ifelse(eligible, atanh(rho), NA_real_)]
    dt[, se_z := ifelse(eligible, se_fisher(n), NA_real_)]
    dt[, direction := ifelse(eligible, sign(rho), NA_real_)]
    dt
  }

  # ==========================================================================
  # Relational grain: the cancer-level meta-analysis input has exactly one
  # row per valid (mirna_id, gene, cancer, arm, adjustment[, population]),
  # where (mirna_id, gene) ranges only over the current run's primary
  # triplets. Stage 08's three_arm_threshold_*.tsv and Stage 10's
  # master_*$mirgene cover the full pair universe and are therefore
  # restricted with restrict_to_keys() (an explicit keyed semi-join against
  # `triplet_pairs`, whose uniqueness is asserted first). No join uses
  # allow.cartesian = TRUE to silence a cardinality mismatch. The two
  # one-to-many broadcasts below (a miRNA-CYT or gene-CYT estimate shared
  # by every triplet containing that miRNA or gene) are verified with
  # assert_broadcast_cardinality().
  # ==========================================================================
  triplet_pairs <- unique(triplets[, .(mirna_id, gene)])
  run_stage_assert(ctx, log_fn, nrow(triplets[, .N, by = .(mirna_id, gene)][N > 1]) == 0,
                    "This run's %d primary triplets have a unique (mirna_id, gene) key (no triplet pair repeated)", nrow(triplets))
  log_fn("This run's primary triplet set: %d triplets, %d distinct (mirna_id, gene) pairs.", nrow(triplets), nrow(triplet_pairs))

  # --- unadjusted: one row per (mirna, gene, cancer, arm) from Stages 06 and
  # 08, restricted to the primary triplet pairs (three_arm_unadj covers the
  # full threshold-0.30 pair universe). ---
  three_arm_restricted <- restrict_to_keys(three_arm_unadj, triplet_pairs, by = c("mirna_id", "gene"))
  run_stage_assert(ctx, log_fn, three_arm_restricted$ok, "%s", if (three_arm_restricted$ok) "ok" else three_arm_restricted$problem)
  run_stage_assert(ctx, log_fn, three_arm_restricted$n_keys_matched == three_arm_restricted$n_keys_expected,
                    "Every one of this run's %d primary triplet pairs was found in Stage 08's three_arm_threshold_%.2f.tsv (%d matched)",
                    three_arm_restricted$n_keys_expected, analysis_config$thresholds$primary, three_arm_restricted$n_keys_matched)
  three_arm_unadj_r <- three_arm_restricted$data
  log_fn("three_arm_unadj restricted to this run's primary triplets: %d -> %d rows (%d -> %d unique pairs).",
         three_arm_restricted$n_before, three_arm_restricted$n_after, uniqueN(three_arm_unadj[, .(mirna_id, gene)]), three_arm_restricted$n_keys_matched)

  long_unadj_mc <- three_arm_unadj_r[, .(mirna_id, gene, cancer)]
  long_unadj_mc <- merge(long_unadj_mc, master_mircyt_unadj[, .(mirna_id, cancer, n = n_samples, rho = cor_mirna_cyt, p_raw = p_mirna_cyt)],
                          by = c("mirna_id", "cancer"), all.x = TRUE)
  long_unadj_mc[, arm := "miRNA-CYT"]
  long_unadj_mg <- three_arm_unadj_r[, .(mirna_id, gene, cancer, arm = "miRNA-gene",
                                          n = n_samples.x, rho = cor_mir_gene, p_raw = p_mir_gene)]
  long_unadj_gc <- three_arm_unadj_r[, .(mirna_id, gene, cancer, arm = "gene-CYT",
                                          n = n_samples.y, rho = cor_gene_cyt, p_raw = p_gene_cyt)]
  long_unadj <- rbindlist(list(long_unadj_mc[, .(mirna_id, gene, cancer, arm, n, rho, p_raw)], long_unadj_mg, long_unadj_gc), use.names = TRUE)
  long_unadj[, population := NA_character_]
  long_unadj <- finalize(long_unadj, "unadjusted", FALSE, character(0))

  # --- CPE / LUMP: built from Stage 10's master_{CPE,LUMP}.rds, every
  # component restricted to the primary triplets. ---
  build_purity_long <- function(master, adj_label, eligible_set) {
    mc_src <- master$mircyt_full[, .(mirna_id, cancer, arm = "miRNA-CYT", n = n_samples, rho = cor_mirna_cyt_adj, p_raw = p_mirna_cyt_adj)]
    # miRNA-CYT is estimated per (miRNA, cancer); it is broadcast
    # (one-to-many) to every primary (mirna, gene) pair containing that
    # miRNA, so the long table stays keyed at (mirna, gene, cancer, arm).
    # The expected fan-out is verified immediately below.
    mc <- merge(mc_src, triplet_pairs, by = "mirna_id", allow.cartesian = TRUE)
    bc_mc <- assert_broadcast_cardinality(mc_src, mc, by = "mirna_id", fanout_keys = triplet_pairs, fanout_by = "mirna_id",
                                           label = sprintf("%s miRNA-CYT broadcast to triplets", adj_label))
    run_stage_assert(ctx, log_fn, bc_mc$ok, "%s", if (bc_mc$ok) sprintf("[%s] miRNA-CYT broadcast to triplets matched its exact expected fan-out for all %d miRNA(s)", adj_label, bc_mc$n_groups_checked) else bc_mc$problem)

    # miRNA-gene: restricted to the primary triplet pairs (master$mirgene
    # covers the full union-0.20 pair universe).
    mg_restricted <- restrict_to_keys(master$mirgene[, .(mirna_id, gene, cancer, n = n_samples, rho = cor_mir_gene_adj, p_raw = p_mir_gene_adj)],
                                       triplet_pairs, by = c("mirna_id", "gene"))
    run_stage_assert(ctx, log_fn, mg_restricted$ok, "%s", if (mg_restricted$ok) "ok" else mg_restricted$problem)
    run_stage_assert(ctx, log_fn, mg_restricted$n_keys_matched == mg_restricted$n_keys_expected,
                      "[%s] Every one of this run's %d primary triplet pairs was found in Stage 10's master_%s.rds$mirgene (%d matched)",
                      adj_label, mg_restricted$n_keys_expected, adj_label, mg_restricted$n_keys_matched)
    log_fn("[%s] mirgene restricted to this run's primary triplets: %d -> %d rows.", adj_label, mg_restricted$n_before, mg_restricted$n_after)
    mg <- mg_restricted$data; mg[, arm := "miRNA-gene"]

    gc_src <- master$genecyt[, .(gene, cancer, arm = "gene-CYT", n = n_samples, rho = cor_gene_cyt_adj, p_raw = p_gene_cyt_adj)]
    # gene-CYT is estimated per gene; it is broadcast (one-to-many) to every
    # primary triplet sharing that gene, with the fan-out verified below.
    gc_expanded <- merge(gc_src, triplet_pairs, by = "gene", allow.cartesian = TRUE)
    bc_gc <- assert_broadcast_cardinality(gc_src, gc_expanded, by = "gene", fanout_keys = triplet_pairs, fanout_by = "gene",
                                           label = sprintf("%s gene-CYT broadcast to triplets", adj_label))
    run_stage_assert(ctx, log_fn, bc_gc$ok, "%s", if (bc_gc$ok) sprintf("[%s] gene-CYT broadcast to triplets matched its exact expected fan-out for all %d gene(s)", adj_label, bc_gc$n_groups_checked) else bc_gc$problem)

    long <- rbindlist(list(mc[, .(mirna_id, gene, cancer, arm, n, rho, p_raw)],
                            mg[, .(mirna_id, gene, cancer, arm, n, rho, p_raw)],
                            gc_expanded[, .(mirna_id, gene, cancer, arm, n, rho, p_raw)]), use.names = TRUE)
    long[, population := NA_character_]
    out <- finalize(long, adj_label, TRUE, eligible_set)
    cc <- assert_cardinality(out, c("mirna_id", "gene", "cancer", "arm"), max_expected = nrow(triplet_pairs) * 3 * length(eligible_set),
                              label = sprintf("[%s] purity-adjusted long table", adj_label))
    run_stage_assert(ctx, log_fn, cc$ok, "%s", if (cc$ok) sprintf("[%s] %d rows, %d unique (mirna,gene,cancer,arm) keys, 0 duplicates, within the dynamically derived maximum of %d", adj_label, cc$n_rows, cc$n_unique_keys, nrow(triplet_pairs) * 3 * length(eligible_set)) else cc$problem)
    out
  }
  long_cpe <- build_purity_long(master_cpe, "CPE", eligible_cancers_cpe)
  long_lump <- build_purity_long(master_lump, "LUMP", eligible_cancers_lump)

  # --- immune: one row per (mirna, gene, cancer, arm, population), with
  # population identity preserved; populations are separate meta-analyses.
  # Stage 12 already restricts itself to load_run_primary_triplets(ctx);
  # the triplet grain is nevertheless verified here.
  run_stage_assert(ctx, log_fn, uniqueN(immune_assoc[, .(mirna_id, gene)]) == nrow(triplet_pairs) &&
                      all(triplet_pairs[unique(immune_assoc[, .(mirna_id, gene)]), on = c("mirna_id", "gene"), nomatch = 0][, .N] == nrow(triplet_pairs)),
                    "Stage 12's immune_composition_associations.tsv already covers exactly this run's %d primary triplet pairs (no restriction needed, verified rather than assumed)", nrow(triplet_pairs))
  n_immune_populations <- uniqueN(immune_assoc$population)
  long_immune <- immune_assoc[, .(mirna_id, gene, cancer, arm, population, n, rho, p_raw = p)]
  long_immune <- finalize(long_immune, "immune", TRUE, eligible_cancers_cpe)

  meta_input <- rbindlist(list(long_unadj, long_cpe, long_lump, long_immune), use.names = TRUE, fill = TRUE)

  # ==========================================================================
  # Cardinality guard: the assembled row count cannot exceed the theoretical
  # maximum derived from the run's own triplets, arms, adjustments,
  # populations and eligible cancers. Each non-immune adjustment
  # (unadjusted, CPE, LUMP) contributes at most n_triplets x 3 arms x that
  # adjustment's cancer count; the immune adjustment contributes at most
  # n_triplets x 3 arms x n_populations x CPE-eligible cancers (tautology
  # exclusions can only reduce this).
  # ==========================================================================
  n_cancers_unadj <- uniqueN(three_arm_unadj_r$cancer)
  theoretical_max_rows <- nrow(triplet_pairs) * 3 * n_cancers_unadj +
    nrow(triplet_pairs) * 3 * length(eligible_cancers_cpe) +
    nrow(triplet_pairs) * 3 * length(eligible_cancers_lump) +
    nrow(triplet_pairs) * 3 * n_immune_populations * length(eligible_cancers_cpe)
  run_stage_assert(ctx, log_fn, nrow(meta_input) <= theoretical_max_rows,
                    "Cancer-level input row count (%d) does not exceed the theoretical maximum (%d) dynamically derived from this run's own %d triplet pairs x 3 arms x (unadjusted:%d + CPE:%d + LUMP:%d cancers + immune:%d populations x %d cancers)",
                    nrow(meta_input), theoretical_max_rows, nrow(triplet_pairs), n_cancers_unadj, length(eligible_cancers_cpe), length(eligible_cancers_lump), n_immune_populations, length(eligible_cancers_cpe))

  # --- primary_support flag: TRUE if (mirna_id, gene, cancer) is a
  # three-arm-passing supporting-cancer record of one of the current run's
  # primary triplets. primary_evidence (Stage 10's
  # cancer_level_evidence_CPE_rho030.tsv) covers the full threshold-0.30
  # pair universe, including pairs that pass in some cancers without
  # reaching analysis_config$discovery$min_recurrence; it is therefore
  # restricted to the primary triplet pairs before flagging. The flag has
  # (mirna_id, gene, cancer) grain, matching the CPE-specific three_arm_pass
  # column, and is broadcast through the (mirna_id, gene, cancer) join to
  # every arm/adjustment/population row of that triplet-cancer, so each
  # meta-analysis model can report how many of its contributing cancers
  # support the triplet in the primary analysis.
  primary_pass_restricted <- restrict_to_keys(primary_evidence[three_arm_pass == TRUE, .(mirna_id, gene, cancer)], triplet_pairs, by = c("mirna_id", "gene"))
  run_stage_assert(ctx, log_fn, primary_pass_restricted$ok, "%s", if (primary_pass_restricted$ok) "ok" else primary_pass_restricted$problem)
  primary_pass <- unique(primary_pass_restricted$data)
  primary_pass[, primary_support := TRUE]
  n_primary_pass_expected <- nrow(primary_pass)
  log_fn("This run's primary triplets' own supporting-cancer records: %d (restricted from %d total three-arm-passing (mirna,gene,cancer) records found across the full %d-pair threshold-0.30 universe -- the remainder belong to pairs that passed in some cancer(s) but never reached this run's own recurrence threshold).",
         n_primary_pass_expected, primary_pass_restricted$n_before, uniqueN(primary_evidence[, .(mirna_id, gene)]))
  meta_input <- merge(meta_input, primary_pass, by = c("mirna_id", "gene", "cancer"), all.x = TRUE)
  meta_input[is.na(primary_support), primary_support := FALSE]
  run_stage_assert(ctx, log_fn, nrow(unique(meta_input[primary_support == TRUE, .(mirna_id, gene, cancer)])) == n_primary_pass_expected,
               "primary_support flag correctly identifies all %d supporting-cancer records for this run's %d primary triplets",
               n_primary_pass_expected, nrow(triplets))

  setcolorder(meta_input, c("mirna_id", "gene", "arm", "adjustment", "population", "cancer", "n", "rho", "z", "se_z", "p_raw",
                             "direction", "primary_support", "eligible", "exclusion_reason"))

  # Key-integrity assertion: (mirna_id, gene, cancer, arm, adjustment,
  # population) must be unique.
  dup_key <- meta_input[, .N, by = .(mirna_id, gene, cancer, arm, adjustment, population)]
  run_stage_assert(ctx, log_fn, all(dup_key$N == 1), "No duplicated (mirna_id, gene, cancer, arm, adjustment, population) key in meta_analysis_input_all_cancers.tsv")
  run_stage_assert(ctx, log_fn, all(meta_input$arm %in% c("miRNA-CYT", "miRNA-gene", "gene-CYT")), "arm column contains only the 3 allowed values")
  run_stage_assert(ctx, log_fn, all(meta_input$adjustment %in% c("unadjusted", "CPE", "LUMP", "immune")), "adjustment column contains only the 4 allowed values")
  # Immune population labels must match Stage 12's output exactly (derived
  # from immune_assoc, never a hard-coded count or list).
  meta_populations <- sort(unique(meta_input[adjustment == "immune"]$population))
  stage11_populations <- sort(unique(immune_assoc$population))
  run_stage_assert(ctx, log_fn, identical(meta_populations, stage11_populations),
                    "The %d immune population label(s) in meta_analysis_input_all_cancers.tsv exactly match Stage 12's own %d output population(s)",
                    length(meta_populations), length(stage11_populations))
  # No non-triplet (mirna_id, gene) pair may appear in the assembled table.
  meta_pairs <- unique(meta_input[, .(mirna_id, gene)])
  illegitimate_pairs <- fsetdiff(meta_pairs, triplet_pairs)
  run_stage_assert(ctx, log_fn, nrow(illegitimate_pairs) == 0,
                    "Every (mirna_id, gene) pair in meta_analysis_input_all_cancers.tsv is one of this run's %d primary triplet pairs (0 illegitimate pairs found, out of %d unique pairs present)",
                    nrow(triplet_pairs), nrow(meta_pairs))

  base_input_paths <- file.path(ctx$run_root, c(
    file.path("tables", "primary", run_primary_triplet_filename(ctx)),
    file.path("processed", "stage06", "master_mirna_cyt.tsv"),
    file.path("processed", "stage08", sprintf("three_arm_threshold_%.2f.tsv", analysis_config$thresholds$primary)),
    file.path("processed", "stage10", "master_CPE.rds"), file.path("processed", "stage10", "master_LUMP.rds"),
    file.path("tables", "downstream", "immune_composition_associations.tsv"),
    file.path("processed", "stage05", "eligible_cancers_cpe.tsv"), file.path("processed", "stage05", "eligible_cancers_lump.tsv"),
    file.path("tables", "audit", run_primary_evidence_filename(ctx, "CPE"))
  ))
  write_run_artifact(ctx, meta_input, file.path("tables", "downstream", "meta_analysis_input_all_cancers.tsv"),
                      "14_cross_cancer_meta_analysis", script_path, base_input_paths, fwrite_fn)
  log_fn("meta_analysis_input_all_cancers.tsv written (%d rows: %d triplets x 3 arms x (3 non-immune adjustments + %d immune populations) x up to %d cancers; %d eligible).",
         nrow(meta_input), nrow(triplets), length(populations_in_immune <- unique(immune_assoc$population)), length(eligible_cancers_cpe), sum(meta_input$eligible))

  # ==========================================================================
  # Random-effects meta-analysis per (triplet, arm, adjustment[, population]).
  # Every MCPcounter population is a separate model.
  # ==========================================================================
  min_k <- analysis_config$meta_analysis$min_cancers_per_model
  meta_estimator <- analysis_config$meta_analysis$estimator
  log_fn("Fitting %s random-effects models (minimum %d cancers per model, prespecified).", meta_estimator, min_k)

  models <- unique(meta_input[, .(mirna_id, gene, arm, adjustment, population)])
  # Cardinality guard: the number of model groups cannot exceed its
  # theoretical maximum, derived as for the row-count guard above, and the
  # set of contributing (mirna_id, gene) pairs must equal the triplet set.
  theoretical_max_models <- nrow(triplet_pairs) * 3 * (1 + 1 + 1 + n_immune_populations)
  run_stage_assert(ctx, log_fn, nrow(models) <= theoretical_max_models,
                    "Model-group count (%d) does not exceed the theoretical maximum (%d) dynamically derived from %d triplet pairs x 3 arms x (unadjusted + CPE + LUMP + %d immune populations)",
                    nrow(models), theoretical_max_models, nrow(triplet_pairs), n_immune_populations)
  model_pairs_illegit <- fsetdiff(unique(models[, .(mirna_id, gene)]), triplet_pairs)
  run_stage_assert(ctx, log_fn, nrow(model_pairs_illegit) == 0,
                    "Every model group's (mirna_id, gene) pair is one of this run's %d primary triplet pairs (0 illegitimate pairs among %d model-contributing pairs)",
                    nrow(triplet_pairs), uniqueN(models[, .(mirna_id, gene)]))
  pooled_rows <- vector("list", nrow(models))
  loo_rows <- list()
  infl_rows <- list()
  fitted_objects <- list()
  diagnostic_failures <- list()
  pop_key <- function(p) ifelse(is.na(p), "NA", p)
  record_diagnostic_failure <- function(mir, tgt_gene, arm_i, adj_i, pop_i, diagnostic, err) {
    diagnostic_failures[[length(diagnostic_failures) + 1]] <<- data.table(
      mirna_id = mir, gene = tgt_gene, arm = arm_i, adjustment = adj_i, population = pop_i, diagnostic = diagnostic,
      error_message = conditionMessage(err), timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%OS"), run_id = ctx$run_id
    )
  }

  for (i in seq_len(nrow(models))) {
    mir <- models$mirna_id[i]; tgt_gene <- models$gene[i]; arm_i <- models$arm[i]; adj_i <- models$adjustment[i]; pop_i <- models$population[i]
    sub <- if (is.na(pop_i)) meta_input[mirna_id == mir & gene == tgt_gene & arm == arm_i & adjustment == adj_i & is.na(population) & eligible == TRUE]
           else meta_input[mirna_id == mir & gene == tgt_gene & arm == arm_i & adjustment == adj_i & population == pop_i & eligible == TRUE]
    k <- nrow(sub)
    # The number of contributing cancers k must equal the number of unique
    # cancers in `sub`, so that duplicated rows can never inflate it.
    k_unique_cancers <- uniqueN(sub$cancer)
    if (k != k_unique_cancers) {
      run_stage_fail(ctx, log_fn,
                      "Model (%s, %s, %s, %s, %s): row count (%d) does not equal unique contributing-cancer count (%d) -- a duplicated cancer row would silently inflate this model's evidence.",
                      mir, tgt_gene, arm_i, adj_i, pop_key(pop_i), k, k_unique_cancers)
    }
    model_key <- sprintf("%s|%s|%s|%s|%s", mir, tgt_gene, arm_i, adj_i, pop_key(pop_i))

    if (k < min_k) {
      pooled_rows[[i]] <- data.table(
        mirna_id = mir, gene = tgt_gene, arm = arm_i, adjustment = adj_i, population = pop_i,
        n_cancers = k, status = "NOT ESTIMABLE", reason = sprintf("fewer than %d cancers with a valid effect estimate (k=%d)", min_k, k),
        pooled_rho = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, raw_p = NA_real_,
        tau2 = NA_real_, QE = NA_real_, QEp = NA_real_, I2 = NA_real_,
        pi_lower = NA_real_, pi_upper = NA_real_, pct_same_direction = NA_real_,
        min_cancer_rho = NA_real_, max_cancer_rho = NA_real_, n_primary_qualifying_in_model = sum(sub$primary_support)
      )
      next
    }

    fit <- tryCatch(rma(yi = z, sei = se_z, data = sub, method = meta_estimator), error = function(e) NULL)
    if (is.null(fit)) {
      pooled_rows[[i]] <- data.table(
        mirna_id = mir, gene = tgt_gene, arm = arm_i, adjustment = adj_i, population = pop_i,
        n_cancers = k, status = "MODEL FIT FAILED", reason = "rma() did not converge",
        pooled_rho = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, raw_p = NA_real_,
        tau2 = NA_real_, QE = NA_real_, QEp = NA_real_, I2 = NA_real_,
        pi_lower = NA_real_, pi_upper = NA_real_, pct_same_direction = NA_real_,
        min_cancer_rho = NA_real_, max_cancer_rho = NA_real_, n_primary_qualifying_in_model = sum(sub$primary_support)
      )
      next
    }

    pr <- tryCatch(predict(fit, transf = transf.ztor), error = function(e) { record_diagnostic_failure(mir, tgt_gene, arm_i, adj_i, pop_i, "predict.rma", e); NULL })
    majority_sign <- sign(fit$beta[1])
    pooled_rows[[i]] <- data.table(
      mirna_id = mir, gene = tgt_gene, arm = arm_i, adjustment = adj_i, population = pop_i,
      n_cancers = k, status = "ESTIMATED", reason = NA_character_,
      pooled_rho = if (!is.null(pr)) pr$pred else tanh(fit$beta[1]),
      ci_lower = if (!is.null(pr)) pr$ci.lb else tanh(fit$ci.lb),
      ci_upper = if (!is.null(pr)) pr$ci.ub else tanh(fit$ci.ub),
      raw_p = fit$pval, tau2 = fit$tau2, QE = fit$QE, QEp = fit$QEp, I2 = fit$I2,
      pi_lower = if (!is.null(pr) && k >= 4) pr$pi.lb else NA_real_,
      pi_upper = if (!is.null(pr) && k >= 4) pr$pi.ub else NA_real_,
      pct_same_direction = round(100 * mean(sub$direction == majority_sign), 1),
      min_cancer_rho = min(sub$rho), max_cancer_rho = max(sub$rho),
      n_primary_qualifying_in_model = sum(sub$primary_support)
    )
    fitted_objects[[model_key]] <- list(fit = fit, sub = sub)

    # --- influence diagnostics (all models with k >= 4) ---
    if (k >= 4) {
      inf <- tryCatch(influence(fit), error = function(e) { record_diagnostic_failure(mir, tgt_gene, arm_i, adj_i, pop_i, "influence.rma", e); NULL })
      if (!is.null(inf)) {
        idf <- as.data.frame(inf$inf)
        idf$cancer <- sub$cancer
        idf$mirna_id <- mir; idf$gene <- tgt_gene; idf$arm <- arm_i; idf$adjustment <- adj_i; idf$population <- pop_i
        idf$is_influential <- inf$is.infl
        infl_rows[[length(infl_rows) + 1]] <- as.data.table(idf)
      }
    }

    # --- leave-one-out (primary CPE-adjusted meta-analysis, all three arms, all triplets with k >= 4) ---
    if (adj_i == "CPE" && k >= 4) {
      l1o <- tryCatch(leave1out(fit, transf = transf.ztor), error = function(e) { record_diagnostic_failure(mir, tgt_gene, arm_i, adj_i, pop_i, "leave1out.rma", e); NULL })
      if (!is.null(l1o)) {
        l1o_df <- as.data.frame(l1o)
        l1o_df$cancer_removed <- sub$cancer
        l1o_df$mirna_id <- mir; l1o_df$gene <- tgt_gene; l1o_df$arm <- arm_i
        loo_rows[[length(loo_rows) + 1]] <- as.data.table(l1o_df)
      }
    }
    if (i %% 40 == 0) log_fn("  Fitted %d/%d triplet-arm-adjustment(-population) models.", i, nrow(models))
  }

  pooled <- rbindlist(pooled_rows, use.names = TRUE, fill = TRUE)

  # --- BH-FDR, one prespecified family per (arm, adjustment, population)
  # combination. Each MCPcounter population is its own family and is never
  # pooled with another population or with the non-immune adjustments.
  pooled[, fdr := NA_real_]
  for (a in unique(pooled$arm)) for (adj in unique(pooled$adjustment)) for (pop in unique(pooled$population[pooled$adjustment == adj])) {
    sel <- if (is.na(pop)) pooled$arm == a & pooled$adjustment == adj & is.na(pooled$population) & pooled$status == "ESTIMATED"
           else pooled$arm == a & pooled$adjustment == adj & pooled$population == pop & pooled$status == "ESTIMATED"
    if (sum(sel) > 0) pooled$fdr[sel] <- p.adjust(pooled$raw_p[sel], method = "BH")
  }
  write_run_artifact(ctx, pooled, file.path("tables", "downstream", "meta_analysis_all_arms.tsv"),
                      "14_cross_cancer_meta_analysis", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "meta_analysis_input_all_cancers.tsv"), fwrite_fn)
  log_fn("meta_analysis_all_arms.tsv written (%d triplet-arm-adjustment-population models; %d estimated, %d not estimable).",
         nrow(pooled), sum(pooled$status == "ESTIMATED"), sum(pooled$status != "ESTIMATED"))

  heterogeneity <- pooled[status == "ESTIMATED", .(mirna_id, gene, arm, adjustment, population, n_cancers, tau2, QE, QEp, I2, pi_lower, pi_upper, pct_same_direction, min_cancer_rho, max_cancer_rho, n_primary_qualifying_in_model)]
  write_run_artifact(ctx, heterogeneity, file.path("tables", "downstream", "meta_analysis_heterogeneity.tsv"),
                      "14_cross_cancer_meta_analysis", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "meta_analysis_all_arms.tsv"), fwrite_fn)

  influence_dt <- if (length(infl_rows)) rbindlist(infl_rows, use.names = TRUE, fill = TRUE) else data.table()
  write_run_artifact(ctx, influence_dt, file.path("tables", "downstream", "meta_analysis_influence.tsv"),
                      "14_cross_cancer_meta_analysis", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "meta_analysis_all_arms.tsv"), fwrite_fn)
  log_fn("meta_analysis_influence.tsv written (%d cancer-level influence rows; %d flagged influential).",
         nrow(influence_dt), if (nrow(influence_dt)) sum(influence_dt$is_influential) else 0)

  loo_dt <- if (length(loo_rows)) rbindlist(loo_rows, use.names = TRUE, fill = TRUE) else data.table()
  write_run_artifact(ctx, loo_dt, file.path("tables", "downstream", "meta_analysis_leave_one_out.tsv"),
                      "14_cross_cancer_meta_analysis", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "meta_analysis_all_arms.tsv"), fwrite_fn)
  n_loo_triplets <- if (nrow(loo_dt)) uniqueN(loo_dt[, .(mirna_id, gene, arm)]) else 0
  log_fn("meta_analysis_leave_one_out.tsv written (%d rows, covering %d triplet-arm CPE-adjusted models).", nrow(loo_dt), n_loo_triplets)

  empty_diag_fail <- data.table(mirna_id = character(), gene = character(), arm = character(), adjustment = character(), population = character(),
                                 diagnostic = character(), error_message = character(), timestamp = character(), run_id = character())
  diag_fail_dt <- if (length(diagnostic_failures)) rbindlist(diagnostic_failures, use.names = TRUE, fill = TRUE) else empty_diag_fail
  write_run_artifact(ctx, diag_fail_dt, file.path("tables", "downstream", "meta_analysis_diagnostic_failures.tsv"),
                      "14_cross_cancer_meta_analysis", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "meta_analysis_all_arms.tsv"), fwrite_fn)
  log_fn("meta_analysis_diagnostic_failures.tsv written (%d recorded failures of predict()/influence()/leave1out() on an otherwise successfully estimated model).", nrow(diag_fail_dt))

  # ==========================================================================
  # Forest plots for highlighted candidates, selected by one data-driven rule:
  # the analysis_config$meta_analysis$n_highlighted_candidates triplets with
  # the highest strict_recurrence in the current run, ties broken
  # alphabetically by (mirna_id, gene). No miRNA or gene identity is
  # hard-coded.
  # ==========================================================================
  n_highlight <- analysis_config$meta_analysis$n_highlighted_candidates
  plot_targets <- if (nrow(triplets) > 0) {
    triplets[order(-strict_recurrence, mirna_id, gene)][seq_len(min(n_highlight, nrow(triplets))), .(mirna_id, gene)]
  } else {
    data.table(mirna_id = character(), gene = character())
  }
  fig_input_paths <- file.path(ctx$run_root, "tables", "downstream", "meta_analysis_all_arms.tsv")
  for (r in seq_len(nrow(plot_targets))) {
    mir <- plot_targets$mirna_id[r]; gene <- plot_targets$gene[r]
    for (arm_i in c("miRNA-CYT", "miRNA-gene", "gene-CYT")) {
      key <- sprintf("%s|%s|%s|CPE|NA", mir, gene, arm_i)
      obj <- fitted_objects[[key]]
      if (is.null(obj)) next
      fname_rel <- file.path("figures", "main_candidates", "meta_analysis",
                              sprintf("forest_%s_%s_%s_CPE.png", gsub("[^A-Za-z0-9]", "", mir), gene, gsub("[^A-Za-z0-9]", "", arm_i)))
      write_run_artifact(ctx, obj, fname_rel, "14_cross_cancer_meta_analysis", script_path, fig_input_paths,
        write_fn = function(o, p) {
          png(p, width = 900, height = 150 + 40 * nrow(o$sub), res = 120)
          tryCatch({
            forest(o$fit, slab = o$sub$cancer, transf = transf.ztor, xlab = "Spearman rho (CPE-adjusted)",
                   main = sprintf("%s - %s (%s, CPE-adjusted)", mir, gene, arm_i))
          }, error = function(e) log_fn("  Forest plot failed for %s: %s", key, conditionMessage(e)))
          dev.off()
        })
    }
  }
  log_fn("Forest plots written to figures/main_candidates/meta_analysis/ for %d triplets.", nrow(plot_targets))

  # ==========================================================================
  # REPORT 22
  # ==========================================================================
  # Highlighted-candidate detail block, applied identically to every entry
  # in plot_targets (the same selection as the forest plots above).
  fmt_highlighted_candidate <- function(mir, target_gene) {
    cand_pooled <- pooled[
      mirna_id == mir &
      gene == target_gene &
      adjustment == "CPE" &
      arm == "gene-CYT"
    ]
    cand_input <- meta_input[
      mirna_id == mir &
      gene == target_gene &
      adjustment == "CPE" &
      arm == "gene-CYT" &
      eligible == TRUE
    ]
    header <- sprintf("## %s - %s", mir, target_gene)
    body <- if (nrow(cand_pooled) && cand_pooled$status[1] == "ESTIMATED")
      sprintf("CPE-adjusted, gene-CYT arm: k=%d cancers (%s), pooled rho=%.4f (95%% CI %.4f-%.4f), raw p=%.4g, FDR=%.4g, tau^2=%.4f, I^2=%.1f%%, Q-p=%.4g, %.1f%% of cancers share the pooled effect's direction, range across cancers [%.4f, %.4f], %d of %d input cancers are among this run's own primary-qualifying records.",
              cand_pooled$n_cancers[1], paste(cand_input$cancer, collapse = ", "),
              cand_pooled$pooled_rho[1], cand_pooled$ci_lower[1], cand_pooled$ci_upper[1],
              cand_pooled$raw_p[1], cand_pooled$fdr[1], cand_pooled$tau2[1], cand_pooled$I2[1], cand_pooled$QEp[1],
              cand_pooled$pct_same_direction[1], cand_pooled$min_cancer_rho[1], cand_pooled$max_cancer_rho[1],
              cand_pooled$n_primary_qualifying_in_model[1], cand_pooled$n_cancers[1])
    else "Not estimable under one or more adjustment/arm combinations in this run, or this pair is not present in this run's matched-expression data; complete per-cancer inputs (where computable) remain in `tables/downstream/meta_analysis_input_all_cancers.tsv` and are not excluded."
    c(header, "", body, "")
  }
  highlighted_candidate_blocks <- unlist(lapply(seq_len(nrow(plot_targets)), function(i) {
    fmt_highlighted_candidate(plot_targets$mirna_id[i], plot_targets$gene[i])
  }))
  if (length(highlighted_candidate_blocks) == 0) {
    highlighted_candidate_blocks <- c("(No triplets were available to highlight in this run.)", "")
  }

  fmt_arm_table <- function(dt, adj) {
    sapply(c("miRNA-CYT", "miRNA-gene", "gene-CYT"), function(a) {
      d <- dt[arm == a & adjustment == adj]
      sprintf("| %s | %d | %d | %.1f%% | %.4f | %.4f-%.4f | %.4g | %.4g | %.1f%% |",
              a, nrow(d), sum(d$status == "ESTIMATED"),
              100 * sum(d$status == "ESTIMATED") / nrow(d),
              mean(d$pooled_rho, na.rm = TRUE), min(d$ci_lower, na.rm = TRUE), max(d$ci_upper, na.rm = TRUE),
              median(d$raw_p, na.rm = TRUE), median(d$fdr, na.rm = TRUE), mean(d$I2, na.rm = TRUE))
    })
  }

  report22 <- c(
    "# Stage 14: Cross-Cancer Random-Effects Meta-Analysis",
    "",
    sprintf("Generated by `scripts/13_cross_cancer_meta_analysis.R`, run `%s`.", ctx$run_id),
    "",
    "## Meta-analysis universe (essential design decision)",
    "",
    sprintf("Eligibility for meta-analysis depends only on prespecified data adequacy: at least %d matched samples for the unadjusted arm, and valid CPE/LUMP purity data (available in %d of this run's eligible cancers) for the CPE-, LUMP-, and immune-adjusted arms. **The synthesis is never restricted to this run's own primary-qualifying cancer-triplet records** -- every cancer meeting the adequacy rule is included regardless of correlation magnitude, FDR, direction, or primary-recurrence status. Each input row carries a `primary_support` flag (TRUE for this run's own primary-qualifying records) for transparency, not as a filter.",
            analysis_config$discovery$min_samples_per_cancer, length(eligible_cancers_cpe)),
    "",
    sprintf("Total candidate rows assembled: %d (%d triplets x 3 arms x (3 non-immune adjustments + %d immune populations) x up to %d cancers); %d are eligible and enter a pooled model.",
            nrow(meta_input), nrow(triplets), length(unique(meta_input$population[meta_input$adjustment == "immune"])), length(eligible_cancers_cpe), sum(meta_input$eligible)),
    "",
    "## Effect-size and pooling method",
    "",
    "Each cancer-specific Spearman rho is Fisher-z transformed (z = atanh(rho)); its standard error is approximated as SE_z = 1.06 / sqrt(n-3), the classical Spearman-specific correction to the Pearson Fisher-z SE (Fieller, Hartley & Pearson, 1957). Pooling uses `metafor::rma(method=\"REML\")`, a random-effects model, separately for every (triplet, arm, adjustment[, population]) combination -- EVERY MCPcounter population is its own model, never pooled with another population or averaged into a composite immune effect. A model is fit only when at least 3 cancers have a valid estimate for that combination; below that it is reported as `NOT ESTIMABLE` with a stated reason, never silently dropped from the output tables. Pooled estimates, confidence intervals, and prediction intervals are back-transformed to the interpretable rho scale via `transf.ztor`.",
    "",
    "## FDR family definition (prespecified)",
    "",
    sprintf("Benjamini-Hochberg correction is applied separately within each (arm, adjustment, population) combination -- %d families in total (%d arms x [3 non-immune adjustments + %d immune populations]), across that combination's triplets in this run. Arms, adjustment variants, and immune populations are never pooled into a single FDR family. This is a distinct set of FDR families from the discovery-stage miRNA-CYT arm's own family (Stages 06/08/10, documented in reports/FDR_FAMILY_METHOD_DECISION.md) -- this stage's synthesis operates downstream of, and does not recompute, that decision.",
            length(unique(pooled$arm)) * (3 + length(unique(pooled$population[pooled$adjustment == "immune"]))),
            length(unique(pooled$arm)), length(unique(pooled$population[pooled$adjustment == "immune"]))),
    "",
    "## Primary (CPE-adjusted) results by arm",
    "",
    "| Arm | Models | Estimated | % estimable | Mean pooled rho | Range of pooled CIs | Median raw p | Median FDR | Mean I^2 |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    fmt_arm_table(pooled, "CPE"),
    "",
    "## Sensitivity results by arm (unadjusted, LUMP)",
    "",
    "| Adjustment | Arm | Models | Estimated | % estimable | Mean pooled rho | Median raw p | Median FDR | Mean I^2 |",
    "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    unlist(lapply(c("unadjusted", "LUMP"), function(adj) {
      sapply(c("miRNA-CYT", "miRNA-gene", "gene-CYT"), function(a) {
        d <- pooled[arm == a & adjustment == adj]
        sprintf("| %s | %s | %d | %d | %.1f%% | %.4f | %.4g | %.4g | %.1f%% |",
                adj, a, nrow(d), sum(d$status == "ESTIMATED"), 100 * sum(d$status == "ESTIMATED") / nrow(d),
                mean(d$pooled_rho, na.rm = TRUE), median(d$raw_p, na.rm = TRUE), median(d$fdr, na.rm = TRUE), mean(d$I2, na.rm = TRUE))
      })
    })),
    "",
    "## Immune-composition sensitivity results, BY POPULATION (never averaged across populations)",
    "",
    "| Population | Arm | Models | Estimated | % estimable | Mean pooled rho | Median raw p | Median FDR | Mean I^2 |",
    "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    unlist(lapply(sort(unique(pooled$population[pooled$adjustment == "immune" & !is.na(pooled$population)])), function(pop) {
      sapply(c("miRNA-CYT", "miRNA-gene", "gene-CYT"), function(a) {
        d <- pooled[arm == a & adjustment == "immune" & population == pop]
        if (nrow(d) == 0) return(sprintf("| %s | %s | 0 | 0 | NA | NA | NA | NA | NA |", pop, a))
        sprintf("| %s | %s | %d | %d | %.1f%% | %.4f | %.4g | %.4g | %.1f%% |",
                pop, a, nrow(d), sum(d$status == "ESTIMATED"), 100 * sum(d$status == "ESTIMATED") / nrow(d),
                mean(d$pooled_rho, na.rm = TRUE), median(d$raw_p, na.rm = TRUE), median(d$fdr, na.rm = TRUE), mean(d$I2, na.rm = TRUE))
      })
    })),
    "",
    "## Heterogeneity and robustness",
    "",
    sprintf("- Prediction intervals reported for %d models with >=4 cancers.", sum(!is.na(pooled$pi_lower))),
    sprintf("- Influence diagnostics computed for %d models (k>=4); %d cancer-model combinations flagged influential by metafor's default rule (not removed, retained in all pooled estimates above).", uniqueN(influence_dt[, .(mirna_id, gene, arm, adjustment, population)]), if (nrow(influence_dt)) sum(influence_dt$is_influential) else 0),
    sprintf("- Leave-one-out diagnostics computed for all %d triplets, all 3 arms, primary CPE-adjusted models with k>=4 (%d triplet-arm combinations); models with k<4 are reported as not eligible for leave-one-out with the reason stated in `meta_analysis_all_arms.tsv`.", nrow(triplets), n_loo_triplets),
    "",
    sprintf("## Highlighted candidates (top %d by this run's own strict_recurrence, ties broken alphabetically -- the SAME explicit rule applied to every triplet; no specific miRNA or gene identity is hard-coded)", n_highlight),
    "",
    highlighted_candidate_blocks,
    "Full per-arm, per-adjustment results for every triplet (all cancers, all four adjustment variants) are retained in `meta_analysis_all_arms.tsv` and `meta_analysis_input_all_cancers.tsv`; no triplet is silently excluded from any table or figure in this stage regardless of whether it was selected for highlighting above.",
    "",
    "## Interpretation (restrained)",
    "",
    "High I^2 or a wide tau^2 in a given triplet-arm model reflects genuine cross-cancer heterogeneity in the underlying biology, technical differences in the EB++ pan-cancer batch adjustment across cohorts, or a combination; this report does not attribute heterogeneity to a single cause. A pooled estimate that changes materially between the CPE-adjusted and immune-adjusted sensitivity models is reported descriptively; it is not automatically interpreted as evidence that immune infiltration explains the association, for the same reasons given in `reports/21_immune_composition_report.md`.",
    "",
    "## Outputs",
    "",
    "- `tables/downstream/meta_analysis_input_all_cancers.tsv`, `meta_analysis_all_arms.tsv`, `meta_analysis_heterogeneity.tsv`, `meta_analysis_influence.tsv`, `meta_analysis_leave_one_out.tsv`",
    "- `figures/main_candidates/meta_analysis/*.png`",
    "- `reports/22_cross_cancer_meta_analysis_report.md` (this file)"
  )
  write_run_artifact(ctx, report22, file.path("reports", "22_cross_cancer_meta_analysis_report.md"),
                      "14_cross_cancer_meta_analysis", script_path, character(0), function(x, p) writeLines(x, p))
  log_fn("reports/22_cross_cancer_meta_analysis_report.md written.")

  record_stage_status(ctx, "14_cross_cancer_meta_analysis", "COMPLETE")
  log_fn("Stage 14 complete.")
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=13 or source scripts/run.R in RStudio.")
}
