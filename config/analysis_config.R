# ============================================================================
# analysis_config.R
#
# Single source of truth for every prespecified parameter, threshold, path
# and random seed used by the pipeline. No stage script defines or
# overrides a threshold, FDR cutoff, minimum sample size, purity metric or
# random seed. Paths under `paths` are relative to the project root and are
# resolved against ctx$data_root. Stage outputs are never located through
# this file; they are written under the current run directory through
# run_output_path() and write_run_artifact() (R/run_isolation.R).
# ============================================================================

analysis_config <- list(

  project_name = "miRNA-CYT pan-cancer triplet discovery",
  config_version = "1.0.0",

  thresholds = list(
    primary = 0.30,
    sensitivity = c(0.25, 0.20),
    all = c(0.30, 0.25, 0.20)
  ),

  discovery = list(
    cor_method = "spearman",
    fdr_cutoff = 0.05,
    min_samples_per_cancer = 30,
    min_recurrence = 3
  ),

  # Declared source scale of each raw expression matrix, based on the known
  # identity of the input files rather than inferred from an observed
  # maximum. The UCSC Xena Pan-Cancer EB++ batch-adjusted RNASeqV2
  # gene-expression matrix and the corresponding EB-adjusted miRNA matrix
  # are distributed as log2-transformed, batch-adjusted values. The
  # maximum-value check in
  # R/expression_preprocessing.R::validate_expression_scale() validates this
  # declaration and stops on a contradiction.
  expression = list(
    mrna_scale = "log2",
    mirna_scale = "log2",
    scale_sanity_max = 50,
    pseudo_count = 1,
    winsorize_probs = c(0.001, 0.999)
  ),

  purity = list(
    primary_metric = "CPE",
    sensitivity_metric = "LUMP"
  ),

  seeds = list(
    master_seed = 20260906,
    bootstrap = 20260906,
    kmeans = 20260906,
    permutation = 20260906,
    null_calibration = 20260906
  ),

  resampling = list(
    rf_n_boot = 200,
    rf_top_k_fraction = 0.5,
    mediation_n_boot = 1000,
    mediation_ci_method = "percentile",
    mediation_ci_level = 0.95
  ),

  meta_analysis = list(
    estimator = "REML",
    min_cancers_per_model = 3,
    heterogeneity_bins = c(0, 25, 50, 75, 100),
    q_multiplicity_method = "BH",
    fisher_z_spearman_correction = 1.06,
    # Number of triplets given a detailed forest plot and report subsection
    # in Stage 14's "Highlighted candidates" section, selected by a fixed
    # data-driven rule: highest strict_recurrence first, ties broken
    # alphabetically by (mirna_id, gene). This is a count; no triplet
    # identity is hard-coded.
    n_highlighted_candidates = 3
  ),

  survival = list(
    endpoint_time_primary = "OS.time",
    endpoint_event_primary = "OS",
    endpoint_time_sensitivity = c("DSS.time", "PFI.time"),
    endpoint_event_sensitivity = c("DSS", "PFI"),
    covariates = c("age_at_initial_pathologic_diagnosis", "gender", "stage_coalesced", "CYT"),
    stage_fields_coalesced = c("ajcc_pathologic_tumor_stage", "clinical_stage"),
    events_per_predictor_minimum = 10,
    stage_completeness_adequacy_pct = 86,
    grade_completeness_adequacy_pct = 95,
    ph_violation_alpha = 0.05,
    convergence_max_abs_coef = 10,
    tiers = list(
      Tier1 = c("feature_z", "age_z", "sex01"),
      Tier2 = c("feature_z", "age_z", "sex01", "stage_z"),
      Tier3 = c("feature_z", "cyt_z", "age_z", "sex01"),
      Tier4 = c("feature_z", "cyt_z", "age_z", "sex01", "stage_z")
    )
  ),

  immune_adjustment = list(
    leukocyte_populations = c("T cells", "B lineage", "NK cells", "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils"),
    cytotoxic_signal_populations_excluded = c("Cytotoxic lymphocytes", "CD8 T cells"),
    stromal_populations_excluded = c("Endothelial cells", "Fibroblasts")
  ),

  collinearity = list(
    vif_high_threshold = 10
  ),

  clustering = list(
    min_k = 2,
    max_k = 6,
    n_init = 200,
    adjust_cancer_effect_for_pca = TRUE,
    adjust_cancer_effect_for_clustering_matrix = TRUE,
    # Null-cluster-structure calibration (R/stats_engine.R::
    # run_null_cluster_calibration()): B independent per-feature-permuted
    # null matrices, each processed with the same within-cancer scaling and
    # k-selection/stability procedure as the observed data. With
    # n_null_iterations = 99 the smallest attainable empirical p-value is
    # 1/100 = 0.01. n_null_boot_stability is smaller than the observed
    # analysis's n_boot_stability = 100 for computational tractability; this
    # affects only the precision of each null draw's stability estimate,
    # not the validity of the comparison, which uses the same formula for
    # observed and null statistics.
    # n_init (above) is not reduced for null iterations, because k-selection
    # must use the same nstart as the observed data. The calibration
    # therefore performs about n_null_iterations x length(k_range) kmeans()
    # calls with nstart = 200 and is the longest-running step of Stage 18.
    n_null_iterations = 99,
    n_null_boot_stability = 20
  ),

  immune = list(
    method = "MCPcounter",
    fallback_method = "marker_gene_zscore",
    # quanTIseq (Stage 13) requires linear, TPM-like input. Whether the
    # EB++AdjustPANCAN batch-adjusted matrix, de-logged via 2^x - 1, satisfies
    # that requirement has not been established. Stage 13 refuses to run
    # while this flag is FALSE. It should be set to TRUE only after input-scale
    # compatibility has been established for this matrix's provenance (e.g.
    # from quanTIseq's published validation literature); the pipeline never
    # changes it.
    quantiseq_input_scale_validated = FALSE
  ),

  # Stage 21 (functional enrichment). `databases`: GO Biological Process
  # and Reactome, read from the versioned Bioconductor annotation packages
  # org.Hs.eg.db, reactome.db and GO.db (scripts/21_functional_enrichment.R).
  # KEGG is not used because its redistribution terms are restrictive.
  # `background`: the unique genes that could have been tested by the Stage
  # 08 three-arm rule, i.e. the gene column of Stage 07's loosest-threshold
  # (union) pair universe, rather than all expressed genes or the whole
  # genome. It is computed from each run's Stage 07 output; this string
  # documents the choice and is not used as a file path.
  enrichment = list(
    databases = c("GO_BP", "Reactome"),
    fdr_cutoff = 0.05,
    background = "stage07_union_pair_universe_genes"
  ),

  # Stage 22 (drug-gene annotation) through the DGIdb public GraphQL API
  # (no authentication key required). Genes are queried from the current
  # run's primary triplets.
  drug_annotation = list(
    database = "DGIdb",
    api_url = "https://dgidb.org/api/graphql"
  ),

  # Stage 23 (cancer-specific cellular-composition context). Primary test: a
  # continuous rank-based interaction model (predictor x standardized
  # MCPcounter population score) for every primary triplet x eligible
  # cancer x population x arm combination. The minimum sample size per
  # combination is analysis_config$discovery$min_samples_per_cancer.
  # min_samples_per_tertile is smaller because tertile strata are
  # descriptive only. interaction_fdr_cutoff is the significance threshold
  # for the interaction BH-FDR family defined in
  # scripts/23_cancer_specific_cellular_context.R.
  cell_context = list(
    min_samples_per_tertile = 10,
    interaction_fdr_cutoff = 0.05
  ),

  # Stage 24 (cellular-context conjunction analysis). All thresholds are
  # prespecified.
  cell_context_corrected = list(
    # Same minimum sample size as the primary pipeline; sufficient for the
    # five-parameter interaction model (four terms plus intercept).
    min_n = 30,
    # Significance threshold for the conjunction-test FDR family.
    conjunction_fdr_cutoff = 0.05,
    # Populations classified as CYT_biological_overlap. The MCPcounter
    # reference contains no literal marker-gene overlap with GZMA/PRF1, but
    # these lymphocyte populations are expected, by biological identity, to
    # correlate with cytotoxic effector gene expression.
    cyt_overlap_populations = c("T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells"),
    # |Spearman rho| above which two MCPcounter population scores are
    # flagged as correlated within a cancer.
    population_correlation_threshold = 0.70,
    # Permutation negative controls. Applied to every conjunction-supported
    # combination plus a fixed random sample of non-supported combinations,
    # rather than the full (triplet, cancer, population) grid, which would
    # require re-running the entire analysis n_permutations times. Monte
    # Carlo uncertainty is reported with every empirical p-value.
    n_permutations = 200,
    permutation_seed = 20260922,
    # Effect-size grid (standardized interaction coefficient) for the
    # simulation-based minimum-detectable-effect and power assessment.
    power_effect_grid = c(0.05, 0.10, 0.15, 0.20, 0.30, 0.40),
    power_n_simulations = 200
  ),

  # Raw and reference inputs are read from data/<file> and
  # reference/<subdirectory>/<file>; every input that is read is listed with
  # its SHA-256 in config/FROM_SCRATCH_ALLOWED_INPUTS.tsv. The miRTarBase,
  # drug_gene_interactions, raw_gene_id and *_root entries are not read by
  # any stage.
  paths = list(
    raw_mrna = "data/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena",
    raw_mirna = "data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena",
    raw_survival = "data/Survival_SupplementalTable_S1_20171025_xena_sp",
    raw_gene_id = "data/Gene.ID.csv",
    targetscan = "reference/targets/Predicted_Targets_Context_Scores.default_predictions.txt.zip",
    mirtarbase_mti = "data/miRTarBase_MTI.csv",
    mirtarbase_strong = "data/miRTarBase_SE_WR.csv",
    mir_family_info = "reference/targets/miR_Family_Info.txt.zip",
    drug_gene_interactions = "data/drug_gene_interactions.tsv",
    mcpcounter_genes = "reference/mcpcounter_signatures/genes.txt",
    mcpcounter_signature_manifest = "reference/mcpcounter_signatures/SIGNATURE_MANIFEST.tsv",
    processed_root = "processed",
    tables_root = "tables",
    figures_root = "figures",
    logs_root = "logs",
    reports_root = "reports"
  ),

  run_status_values = c("PENDING", "RUNNING", "COMPLETE", "FAILED", "SKIPPED")
)

validate_analysis_config <- function(cfg) {
  stopifnot(
    is.numeric(cfg$thresholds$primary),
    all(cfg$thresholds$primary %in% cfg$thresholds$all),
    cfg$discovery$fdr_cutoff > 0 && cfg$discovery$fdr_cutoff < 1,
    cfg$discovery$min_samples_per_cancer >= 1,
    cfg$discovery$min_recurrence >= 1,
    cfg$purity$primary_metric != cfg$purity$sensitivity_metric,
    cfg$discovery$cor_method == "spearman",
    cfg$expression$mrna_scale %in% c("log2", "linear"),
    cfg$expression$mirna_scale %in% c("log2", "linear")
  )
  invisible(TRUE)
}
