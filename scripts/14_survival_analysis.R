# ============================================================================
# 15_survival_analysis.R
#
# Purpose: test whether expression of each primary miRNA and gene feature of
#          the current run is associated with overall survival using
#          cancer-specific Cox models (rather than one pooled model with
#          cancer as an additive covariate, which would assume common
#          baseline hazards and covariate effects across cancer types), and
#          quantify whether each feature adds prognostic information beyond
#          CYT.
# Inputs:  <run_root>/processed/stage03/survival_matched.rds,
#          mirna_matched.rds, mrna_matched.rds (same-run);
#          processed/stage04/cyt_scores.tsv (same-run);
#          processed/stage05/sample_annotation_with_purity.tsv (same-run);
#          the current run's primary triplet set (scripts/lib/run_discovery_source.R).
# Outputs: <run_root>/tables/downstream/clinical_data_completeness_by_cancer.tsv,
#          survival_models.tsv, proportional_hazards_tests.tsv,
#          cyt_incremental_models.tsv;
#          reports/23_survival_input_audit.md,
#          24_survival_analysis_report.md.
# Method:  for each unique primary feature (miRNAs and genes, counted from
#          the run's triplet set) and each cancer with at least
#          analysis_config$discovery$min_samples_per_cancer matched samples
#          and at least one event, up to four prespecified Cox models are
#          fitted: Tier 1 (feature + age + sex), Tier 2 (Tier 1 + stage, only
#          in cancers with adequate stage completeness), Tier 3 (feature +
#          CYT + age + sex) and Tier 4 (Tier 3 + stage, same stage-adequate
#          cancers). Continuous covariates (feature expression, age, stage)
#          are standardized within each cancer; stage is a coalesced
#          AJCC/clinical ordinal (I-IV), not a data-driven cut point. A model
#          is fitted only if it satisfies a minimum events-per-predictor rule
#          (Peduzzi et al. 1996); otherwise it is reported as NOT ESTIMABLE
#          with a reason. The proportional-hazards assumption of every model
#          is tested (cox.zph) and global violations are flagged. Incremental
#          value beyond CYT is tested by a likelihood-ratio comparison of
#          Tier 3 against a CYT + age + sex baseline fitted on the identical
#          complete-case sample set, together with the change in
#          concordance index. No external validation is performed, and no
#          claim of validated prognostic performance is made.
# Depends on: 03_harmonize_identifiers.R, 04_compute_validate_cyt.R,
#             05_define_eligible_features.R, 09_recurrence_analysis.R,
#             10_purity_adjusted_analysis.R (all same-run).
# ============================================================================

run_stage15 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "15_survival_analysis")
  source(file.path(ctx$data_root, "scripts", "lib", "clinical_stage_mapping.R"), local = TRUE)
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "run_discovery_source.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "expression_preprocessing.R"), local = TRUE)
  suppressPackageStartupMessages({
    library(data.table)
    library(survival)
  })
  script_path <- file.path(ctx$data_root, "scripts", "14_survival_analysis.R")
  fwrite_fn <- function(x, p) fwrite(x, p, sep = "\t")

  triplets <- load_run_primary_triplets(ctx)
  features <- unique(rbind(
    data.table(feature = triplets$mirna_id, feature_type = "miRNA"),
    data.table(feature = triplets$gene, feature_type = "gene")
  ))
  n_mirna_expected <- uniqueN(triplets$mirna_id)
  n_gene_expected <- uniqueN(triplets$gene)
  run_stage_assert(ctx, log_fn, nrow(features) == n_mirna_expected + n_gene_expected &&
                        sum(features$feature_type == "miRNA") == n_mirna_expected &&
                        sum(features$feature_type == "gene") == n_gene_expected,
               "unique this-run primary feature count is internally consistent (%d miRNAs, %d genes)", n_mirna_expected, n_gene_expected)

  # ==========================================================================
  # PHASE 1: CLINICAL-DATA INPUT AUDIT
  # ==========================================================================
  log_fn("PHASE 1: clinical-data completeness audit.")

  surv <- as.data.table(read_run_artifact(ctx, file.path("processed", "stage03", "survival_matched.rds"), "03_harmonize_identifiers", readRDS))
  sample_annot <- read_run_artifact(ctx, file.path("processed", "stage05", "sample_annotation_with_purity.tsv"), "05_define_eligible_features", fread)
  cyt_scores <- read_run_artifact(ctx, file.path("processed", "stage04", "cyt_scores.tsv"), "04_compute_validate_cyt", fread)
  surv <- merge(surv, sample_annot[, .(sample_id, cancer, CPE, LUMP)], by = "sample_id")
  surv <- merge(surv, cyt_scores[, .(sample_id, CYT)], by = "sample_id")

  surv[, stage_coalesced_raw := ifelse(trimws(ajcc_pathologic_tumor_stage) != "", trimws(ajcc_pathologic_tumor_stage), trimws(clinical_stage))]
  surv[, stage_numeric := stage_to_numeric(stage_coalesced_raw)]
  surv[, grade_clean := trimws(histological_grade)]
  surv[, gender_clean := trimws(gender)]
  surv[, age := age_at_initial_pathologic_diagnosis]

  min_n <- analysis_config$discovery$min_samples_per_cancer
  stage_adequacy_cutoff <- analysis_config$survival$stage_completeness_adequacy_pct  # stage completeness is typically bimodal: most cancers lie well above the threshold, a few have near-zero completeness

  completeness <- surv[, .(
    N = .N,
    n_events = sum(OS == 1, na.rm = TRUE),
    pct_OS_complete = round(100 * mean(!is.na(OS) & !is.na(OS.time)), 1),
    pct_age_complete = round(100 * mean(!is.na(age)), 1),
    pct_sex_complete = round(100 * mean(gender_clean != "" & !is.na(gender_clean)), 1),
    pct_stage_complete_ajcc_only = round(100 * mean(trimws(ajcc_pathologic_tumor_stage) != ""), 1),
    pct_stage_complete_coalesced = round(100 * mean(!is.na(stage_numeric)), 1),
    pct_grade_complete = round(100 * mean(grade_clean != "" & !is.na(grade_clean)), 1),
    pct_CPE_complete = round(100 * mean(!is.na(CPE)), 1),
    pct_CYT_complete = round(100 * mean(!is.na(CYT)), 1)
  ), by = cancer][order(cancer)]
  completeness[, stage_adequate := pct_stage_complete_coalesced >= stage_adequacy_cutoff]
  completeness[, grade_adequate := pct_grade_complete >= analysis_config$survival$grade_completeness_adequacy_pct]
  completeness[, eligible_for_survival := N >= min_n & n_events >= 1]
  write_run_artifact(ctx, completeness, file.path("tables", "downstream", "clinical_data_completeness_by_cancer.tsv"),
                      "15_survival_analysis", script_path,
                      file.path(ctx$run_root, c(file.path("processed", "stage03", "survival_matched.rds"),
                                                 file.path("processed", "stage05", "sample_annotation_with_purity.tsv"),
                                                 file.path("processed", "stage04", "cyt_scores.tsv"))),
                      fwrite_fn)

  n_stage_adequate <- sum(completeness$stage_adequate)
  n_grade_zero <- sum(completeness$pct_grade_complete == 0)
  n_grade_adequate <- sum(completeness$grade_adequate)
  materially_corrected <- completeness[abs(pct_stage_complete_coalesced - pct_stage_complete_ajcc_only) > 10]
  log_fn("Stage-adequate cancers (coalesced field, >=%d%% complete): %d of %d.", stage_adequacy_cutoff, n_stage_adequate, nrow(completeness))
  log_fn("Cancers with 0%% grade completeness: %d of %d.", n_grade_zero, nrow(completeness))

  report23 <- c(
    "# Stage 15 Input Validation: Clinical-Covariate Completeness for Survival Analysis",
    "",
    sprintf("Generated by `scripts/14_survival_analysis.R`, run `%s`.", ctx$run_id),
    "",
    "## Tumor-stage completeness: coalesced field vs. raw AJCC field",
    "",
    sprintf("Tumor-stage completeness is computed here as the fraction of samples whose stage value actually parses to a usable ordinal Stage I-IV, using a coalesced field rather than the raw `ajcc_pathologic_tumor_stage` column alone, for two reasons: (1) TCGA records gynecological and lymphoid cancers (OV, CESC, UCEC, UCS, DLBC, THYM, and a modest fraction of HNSC) in `clinical_stage` instead of `ajcc_pathologic_tumor_stage`; coalescing the two fields recovers a usable stage value for these cancers. (2) A small number of cancers' `ajcc_pathologic_tumor_stage` field is mostly non-blank but contains some non-ordinal placeholder values (e.g. a TCGA-flagged conflicting-source record, or an in-situ code not on the invasive I-IV scale), which are correctly excluded from the parsed count, not coalescing artifacts. Combined, stage completeness (using the coalesced field, parsed to Stage I-IV) differs from the raw non-blank `ajcc_pathologic_tumor_stage` count by more than 10 percentage points for %d cancers:", nrow(materially_corrected)),
    "",
    "| Cancer | Old (non-blank ajcc field) % | Corrected (parses to Stage I-IV) % | Direction |",
    "|---|---:|---:|---|",
    sprintf("| %s | %.1f | %.1f | %s |", materially_corrected$cancer, materially_corrected$pct_stage_complete_ajcc_only, materially_corrected$pct_stage_complete_coalesced,
            ifelse(materially_corrected$pct_stage_complete_coalesced > materially_corrected$pct_stage_complete_ajcc_only, "recovered via clinical_stage", "excluded non-ordinal placeholder values")),
    "",
    sprintf("This stage uses the corrected, coalesced stage field throughout. %d of %d cancers now meet the prespecified %d%% stage-completeness threshold for Tiers 2 and 4 (cancers that remain structurally without a usable stage field in this data source are excluded from stage-adjusted tiers only, not from Tiers 1 and 3).", n_stage_adequate, nrow(completeness), stage_adequacy_cutoff),
    "",
    "## Reconciliation of the histological-grade completeness count",
    "",
    sprintf("Direct computation from this run's own `processed/stage03/survival_matched.rds` finds histological-grade completeness is exactly 0%% in %d of %d cancers. This pipeline uses its own directly-computed count as authoritative rather than any figure carried over from elsewhere. Grade is not used in any pan-cancer claim and is not included in any of the four prespecified survival tiers below, because it is missing entirely for a substantial fraction of cancer types; it would be usable only as a cancer-specific sensitivity covariate in the %d cancers with >=%d%% grade completeness, and no such sensitivity model is fit in this stage.",
            n_grade_zero, nrow(completeness), n_grade_adequate, analysis_config$survival$grade_completeness_adequacy_pct),
    "",
    "## Full per-cancer completeness table",
    "",
    "See `tables/downstream/clinical_data_completeness_by_cancer.tsv` for the complete table (N, events, OS/age/sex/stage/grade/CPE/CYT completeness, stage_adequate and eligible_for_survival flags) underlying every design decision in this report.",
    "",
    "## Prespecified tiered model design (fixed before any outcome was inspected)",
    "",
    "- Tier 1: feature + age + sex",
    "- Tier 2: feature + age + sex + stage (stage-adequate cancers only)",
    "- Tier 3: feature + CYT + age + sex",
    "- Tier 4: feature + CYT + age + sex + stage (stage-adequate cancers only)",
    "",
    "All continuous covariates (feature expression, age, and stage where included) are standardized (z-scored) within each cancer before fitting. Stage uses a coalesced-field ordinal encoding (I-IV), not a data-driven cutpoint. A nuisance covariate with no within-cancer variation is omitted because it cannot be estimated or confound a within-cancer association; the requested covariates, actual fitted covariates, and any invariant covariates omitted are recorded for every model. The feature itself is never omitted: zero feature variance makes that model `NOT ESTIMABLE`. The events-per-predictor requirement is calculated from the covariates actually fitted. Cancers/tiers failing a requirement are reported `NOT ESTIMABLE` with the specific reason, never silently omitted."
  )
  write_run_artifact(ctx, report23, file.path("reports", "23_survival_input_audit.md"),
                      "15_survival_analysis", script_path, character(0), function(x, p) writeLines(x, p))
  log_fn("reports/23_survival_input_audit.md written.")

  # ==========================================================================
  # PHASE 2: TIERED COX MODELS, PER FEATURE PER CANCER
  # ==========================================================================
  log_fn("PHASE 2: fitting tiered cancer-specific Cox models for %d features.", nrow(features))

  mirna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mirna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mirna_mat <- as.matrix(mirna_matched[, -1, drop = FALSE]); rownames(mirna_mat) <- mirna_matched$mirna_id
  mrna_matched <- read_run_artifact(ctx, file.path("processed", "stage03", "mrna_matched.rds"), "03_harmonize_identifiers", readRDS)
  mrna_mat <- as.matrix(mrna_matched[, -1, drop = FALSE]); rownames(mrna_mat) <- toupper(trimws(mrna_matched$gene_id_raw))
  if (any(duplicated(rownames(mrna_mat)))) {
    dup <- unique(rownames(mrna_mat)[duplicated(rownames(mrna_mat))])
    keep <- !logical(nrow(mrna_mat))
    for (s in dup) { idx <- which(rownames(mrna_mat) == s); v <- apply(mrna_mat[idx, , drop = FALSE], 1, var, na.rm = TRUE); keep[idx[-which.max(v)]] <- FALSE }
    mrna_mat <- mrna_mat[keep, , drop = FALSE]
  }
  # Per-feature winsorization before use as a Cox-model covariate: the Cox
  # partial likelihood is sensitive to the leverage of extreme expression
  # values, unlike the rank-based Spearman arms of Stages 06 and 08.
  mirna_mat <- winsorize_matrix_rows(mirna_mat)
  mrna_mat <- winsorize_matrix_rows(mrna_mat)

  eligible_cancers <- completeness[eligible_for_survival == TRUE]$cancer
  epp_min <- analysis_config$survival$events_per_predictor_minimum

  fit_one <- function(d, covars) {
    requested_covars <- covars
    variable_covariate <- vapply(covars, function(v) {
      values <- d[[v]]
      values <- values[is.finite(values)]
      length(unique(values)) >= 2
    }, logical(1))
    used_covars <- covars[variable_covariate]
    dropped_covars <- setdiff(covars, used_covars)

    if ("feature_z" %in% dropped_covars) {
      return(list(status = "NOT ESTIMABLE", reason = "zero variance in feature expression in the model-specific analysis set",
                  fit = NULL, n = nrow(d), n_events = sum(d$OS == 1, na.rm = TRUE),
                  requested_covars = requested_covars, used_covars = used_covars, dropped_covars = dropped_covars))
    }

    required_columns <- c("OS", "OS.time", used_covars)
    complete <- complete.cases(d[, ..required_columns])
    d_model <- d[complete]
    n_events <- sum(d_model$OS == 1, na.rm = TRUE)
    min_events_needed <- epp_min * length(used_covars)
    if (n_events < min_events_needed) {
      return(list(status = "NOT ESTIMABLE", reason = sprintf("insufficient events for EPP>=%d (needed >=%d, have %d)", epp_min, min_events_needed, n_events),
                  fit = NULL, n = nrow(d_model), n_events = n_events,
                  requested_covars = requested_covars, used_covars = used_covars, dropped_covars = dropped_covars))
    }
    form <- as.formula(paste("Surv(OS.time, OS) ~", paste(used_covars, collapse = " + ")))
    fit <- tryCatch(coxph(form, data = d_model, x = TRUE), error = function(e) NULL)
    if (is.null(fit) || !all(is.finite(coef(fit)))) {
      return(list(status = "MODEL FIT FAILED", reason = "coxph() did not converge or produced non-finite coefficients",
                  fit = NULL, n = nrow(d_model), n_events = n_events,
                  requested_covars = requested_covars, used_covars = used_covars, dropped_covars = dropped_covars))
    }
    converged <- all(abs(coef(fit)) < analysis_config$survival$convergence_max_abs_coef) && (is.null(fit$iter) || fit$iter[1] < 20)
    list(status = if (converged) "ESTIMATED" else "ESTIMATED (CONVERGENCE CAUTION: large |coef| or iteration limit reached)",
         reason = NA_character_, fit = fit, n = nrow(d_model), n_events = n_events,
         requested_covars = requested_covars, used_covars = used_covars, dropped_covars = dropped_covars)
  }

  model_rows <- list()
  ph_rows <- list()
  incr_rows <- list()

  for (fi in seq_len(nrow(features))) {
    feat <- features$feature[fi]; ftype <- features$feature_type[fi]
    expr_row <- if (ftype == "miRNA") mirna_mat[feat, ] else mrna_mat[feat, ]
    expr_dt <- data.table(sample_id = names(expr_row), feature_expr = as.numeric(expr_row))

    for (cc in eligible_cancers) {
      d <- surv[cancer == cc, .(sample_id, OS, OS.time, age, gender_clean, stage_numeric, CYT)]
      d <- merge(d, expr_dt, by = "sample_id")
      d <- d[!is.na(OS) & !is.na(OS.time) & OS.time > 0 & !is.na(age) & gender_clean %in% c("MALE", "FEMALE") & !is.na(feature_expr)]
      # Every attempted (feature, cancer) combination appears in the output,
      # even when no tier can be fitted: explicit NOT ESTIMABLE rows are
      # recorded for all four tiers.
      if (nrow(d) < min_n) {
        for (tname in names(analysis_config$survival$tiers)) {
          model_rows[[length(model_rows) + 1]] <- data.table(
            feature = feat, feature_type = ftype, cancer = cc, tier = tname,
            requested_covariates = paste(analysis_config$survival$tiers[[tname]], collapse = "+"),
            covariates = paste(analysis_config$survival$tiers[[tname]], collapse = "+"),
            fitted_covariates = NA_character_, dropped_invariant_covariates = NA_character_, n_predictors = NA_integer_,
            n = nrow(d), n_events = sum(d$OS == 1, na.rm = TRUE),
            status = "NOT ESTIMABLE", reason = sprintf("fewer than %d complete-case samples after requiring non-missing OS/age/sex/feature (n=%d)", min_n, nrow(d)),
            HR = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, raw_p = NA_real_, ph_global_p = NA_real_
          )
        }
        next
      }
      # Zero-variance guard: scale() on a constant vector divides by zero
      # and would produce NaN without an error or a recorded reason.
      if (length(unique(d$feature_expr)) < 2) {
        for (tname in names(analysis_config$survival$tiers)) {
          model_rows[[length(model_rows) + 1]] <- data.table(
            feature = feat, feature_type = ftype, cancer = cc, tier = tname,
            requested_covariates = paste(analysis_config$survival$tiers[[tname]], collapse = "+"),
            covariates = paste(analysis_config$survival$tiers[[tname]], collapse = "+"),
            fitted_covariates = NA_character_, dropped_invariant_covariates = NA_character_, n_predictors = NA_integer_,
            n = nrow(d), n_events = sum(d$OS == 1, na.rm = TRUE),
            status = "NOT ESTIMABLE", reason = "zero variance in feature expression after complete-case filtering in this cancer",
            HR = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, raw_p = NA_real_, ph_global_p = NA_real_
          )
        }
        next
      }
      d[, feature_z := as.numeric(scale(feature_expr))]
      d[, age_z := as.numeric(scale(age))]
      d[, sex01 := as.integer(gender_clean == "MALE")]
      d[, cyt_z := as.numeric(scale(CYT))]
      stage_adeq <- completeness[cancer == cc]$stage_adequate
      d_stage <- d[!is.na(stage_numeric)]
      if (nrow(d_stage) > 0) d_stage[, stage_z := as.numeric(scale(stage_numeric))]

      stage_usable <- isTRUE(stage_adeq) && nrow(d_stage) >= min_n
      stage_skip_reason <- if (isTRUE(stage_adeq)) sprintf("stage-adequate cancer but fewer than %d complete-case samples with a usable stage value", min_n) else sprintf("cancer type not stage-adequate (<%d%% coalesced stage completeness)", stage_adequacy_cutoff)

      cfg_tiers <- analysis_config$survival$tiers
      tier_data <- list(Tier1 = d, Tier2 = d_stage, Tier3 = d, Tier4 = d_stage)
      tier_applicable <- list(Tier1 = TRUE, Tier2 = stage_usable, Tier3 = TRUE, Tier4 = stage_usable)
      tier_defs <- setNames(lapply(names(cfg_tiers), function(tn) {
        list(covars = cfg_tiers[[tn]], data = tier_data[[tn]], applicable = tier_applicable[[tn]])
      }), names(cfg_tiers))

      for (tname in names(tier_defs)) {
        td <- tier_defs[[tname]]
        res <- if (td$applicable) fit_one(td$data, td$covars) else
          list(status = "NOT ESTIMABLE", reason = stage_skip_reason, fit = NULL, n = nrow(td$data), n_events = sum(td$data$OS == 1, na.rm = TRUE),
               requested_covars = td$covars, used_covars = character(), dropped_covars = character())
        row <- data.table(feature = feat, feature_type = ftype, cancer = cc, tier = tname,
                           requested_covariates = paste(res$requested_covars, collapse = "+"),
                           covariates = paste(res$requested_covars, collapse = "+"),
                           fitted_covariates = if (length(res$used_covars)) paste(res$used_covars, collapse = "+") else NA_character_,
                           dropped_invariant_covariates = if (length(res$dropped_covars)) paste(res$dropped_covars, collapse = "+") else "",
                           n_predictors = length(res$used_covars), n = res$n, n_events = res$n_events,
                           status = res$status, reason = res$reason,
                           HR = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, raw_p = NA_real_, ph_global_p = NA_real_)
        if (!is.null(res$fit)) {
          s <- summary(res$fit)
          coefrow <- s$coefficients["feature_z", ]
          confrow <- s$conf.int["feature_z", ]
          row$HR <- unname(coefrow["exp(coef)"]); row$raw_p <- unname(coefrow["Pr(>|z|)"])
          row$ci_lower <- unname(confrow["lower .95"]); row$ci_upper <- unname(confrow["upper .95"])
          zph <- tryCatch(cox.zph(res$fit), error = function(e) NULL)
          if (!is.null(zph)) {
            row$ph_global_p <- zph$table["GLOBAL", "p"]
            ph_rows[[length(ph_rows) + 1]] <- data.table(
              feature = feat, feature_type = ftype, cancer = cc, tier = tname,
              term = rownames(zph$table), chisq = zph$table[, "chisq"], df = zph$table[, "df"], p = zph$table[, "p"]
            )
          }
        }
        model_rows[[length(model_rows) + 1]] <- row
      }

      # --- incremental value beyond CYT: Tier 3 (feature + CYT + age + sex) vs. CYT + age + sex baseline on identical complete-case data ---
      if ("Tier3" %in% names(tier_defs)) {
        res3 <- fit_one(d, c("feature_z", "cyt_z", "age_z", "sex01"))
        res0 <- fit_one(d, c("cyt_z", "age_z", "sex01"))
        if (!is.null(res3$fit) && !is.null(res0$fit)) {
          lrt_stat <- 2 * (as.numeric(logLik(res3$fit)) - as.numeric(logLik(res0$fit)))
          lrt_p <- pchisq(lrt_stat, df = 1, lower.tail = FALSE)
          c3 <- summary(res3$fit)$concordance[1]; c0 <- summary(res0$fit)$concordance[1]
          incr_rows[[length(incr_rows) + 1]] <- data.table(
            feature = feat, feature_type = ftype, cancer = cc, n = res3$n, n_events = res3$n_events,
            baseline_model = paste(res0$used_covars, collapse = "+"), full_model = paste(res3$used_covars, collapse = "+"),
            dropped_invariant_covariates = paste(sort(unique(c(res3$dropped_covars, res0$dropped_covars))), collapse = "+"),
            concordance_baseline = c0, concordance_full = c3, concordance_delta = c3 - c0,
            lrt_chisq = lrt_stat, lrt_df = 1, lrt_p = lrt_p, status = "ESTIMATED"
          )
        } else {
          incr_rows[[length(incr_rows) + 1]] <- data.table(
            feature = feat, feature_type = ftype, cancer = cc, n = res3$n, n_events = res3$n_events,
            baseline_model = paste(res0$used_covars, collapse = "+"), full_model = paste(res3$used_covars, collapse = "+"),
            dropped_invariant_covariates = paste(sort(unique(c(res3$dropped_covars, res0$dropped_covars))), collapse = "+"),
            concordance_baseline = NA_real_, concordance_full = NA_real_, concordance_delta = NA_real_,
            lrt_chisq = NA_real_, lrt_df = 1, lrt_p = NA_real_,
            status = if (is.null(res3$fit)) sprintf("NOT ESTIMABLE (full model: %s)", res3$reason) else sprintf("NOT ESTIMABLE (baseline model: %s)", res0$reason)
          )
        }
      }
    }
    if (fi %% 10 == 0) log_fn("  Processed %d/%d features.", fi, nrow(features))
  }

  models_dt <- rbindlist(model_rows, use.names = TRUE, fill = TRUE)
  systematic_failures <- models_dt[status == "MODEL FIT FAILED", .N, by = .(cancer, tier)][N == nrow(features)]
  run_stage_assert(ctx, log_fn, nrow(systematic_failures) == 0,
                   "No cancer-tier combination has MODEL FIT FAILED for all %d tested features", nrow(features))
  models_dt[, fdr := NA_real_]
  for (tname in unique(models_dt$tier)) {
    sel <- models_dt$tier == tname & grepl("^ESTIMATED", models_dt$status)
    if (sum(sel) > 0) models_dt$fdr[sel] <- p.adjust(models_dt$raw_p[sel], method = "BH")
  }
  models_dt[, ph_violated := !is.na(ph_global_p) & ph_global_p < analysis_config$survival$ph_violation_alpha]
  survival_inputs <- file.path(ctx$run_root, c(
    file.path("processed", "stage03", "survival_matched.rds"), file.path("processed", "stage03", "mirna_matched.rds"),
    file.path("processed", "stage03", "mrna_matched.rds"), file.path("processed", "stage04", "cyt_scores.tsv"),
    file.path("processed", "stage05", "sample_annotation_with_purity.tsv")
  ))
  write_run_artifact(ctx, models_dt, file.path("tables", "downstream", "survival_models.tsv"),
                      "15_survival_analysis", script_path, survival_inputs, fwrite_fn)
  log_fn("survival_models.tsv written (%d rows; %d estimated, %d not estimable, %d PH-violated).",
         nrow(models_dt), sum(grepl("^ESTIMATED", models_dt$status)), sum(!grepl("^ESTIMATED", models_dt$status)), sum(models_dt$ph_violated, na.rm = TRUE))

  ph_dt <- rbindlist(ph_rows, use.names = TRUE, fill = TRUE)
  write_run_artifact(ctx, ph_dt, file.path("tables", "downstream", "proportional_hazards_tests.tsv"),
                      "15_survival_analysis", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "survival_models.tsv"), fwrite_fn)

  incr_dt <- rbindlist(incr_rows, use.names = TRUE, fill = TRUE)
  incr_dt[, lrt_fdr := NA_real_]
  sel <- incr_dt$status == "ESTIMATED"
  if (sum(sel) > 0) incr_dt$lrt_fdr[sel] <- p.adjust(incr_dt$lrt_p[sel], method = "BH")
  write_run_artifact(ctx, incr_dt, file.path("tables", "downstream", "cyt_incremental_models.tsv"),
                      "15_survival_analysis", script_path,
                      file.path(ctx$run_root, "tables", "downstream", "survival_models.tsv"), fwrite_fn)
  log_fn("cyt_incremental_models.tsv written (%d rows; %d estimated LRT comparisons, %d not estimable).",
         nrow(incr_dt), sum(incr_dt$status == "ESTIMATED"), sum(incr_dt$status != "ESTIMATED"))

  # ==========================================================================
  # REPORT 24
  # ==========================================================================
  n_total_models <- nrow(models_dt)
  n_estimated <- sum(grepl("^ESTIMATED", models_dt$status))
  n_not_estimable <- sum(models_dt$status == "NOT ESTIMABLE")
  n_fit_failed <- sum(models_dt$status == "MODEL FIT FAILED")
  n_invariant_adjustments <- sum(models_dt$dropped_invariant_covariates != "" & !is.na(models_dt$dropped_invariant_covariates))
  n_sig_raw <- sum(models_dt$raw_p < 0.05, na.rm = TRUE)
  n_sig_fdr <- sum(models_dt$fdr < 0.05, na.rm = TRUE)
  n_ph_violated <- sum(models_dt$ph_violated, na.rm = TRUE)
  n_incr_sig_fdr <- sum(incr_dt$lrt_fdr < 0.05, na.rm = TRUE)

  tier_summary <- models_dt[, .(
    n_models = .N, n_estimated = sum(grepl("^ESTIMATED", status)), n_not_estimable = sum(status == "NOT ESTIMABLE"),
    n_fit_failed = sum(status == "MODEL FIT FAILED"),
    n_sig_raw = sum(raw_p < 0.05, na.rm = TRUE), n_sig_fdr = sum(fdr < 0.05, na.rm = TRUE), n_ph_violated = sum(ph_violated, na.rm = TRUE)
  ), by = tier][order(tier)]

  report24 <- c(
    "# Stage 15: Cancer-Specific Survival Analysis",
    "",
    sprintf("Generated by `scripts/14_survival_analysis.R`, run `%s`.", ctx$run_id),
    "",
    "## Design",
    "",
    sprintf("Every model is fit within a single cancer type; no model pools samples across cancer types. %d unique this-run primary features (%d miRNAs, %d genes) were tested against overall survival in each of %d cancers with >=%d matched samples and at least one event, using the four prespecified tiers in `reports/23_survival_input_audit.md`. A minimum of %d events per predictor (Peduzzi et al., 1996) was required before any model was fit; models failing this rule are reported `NOT ESTIMABLE`, not omitted.",
            nrow(features), n_mirna_expected, n_gene_expected, length(eligible_cancers), min_n, epp_min),
    "",
    "## Overall model counts",
    "",
    sprintf("- Total (feature x cancer x tier) models attempted: %d", n_total_models),
    sprintf("- Estimated: %d; Not estimable: %d; Model-fit failures: %d", n_estimated, n_not_estimable, n_fit_failed),
    sprintf("- Models fitted after transparently omitting at least one invariant nuisance covariate: %d", n_invariant_adjustments),
    sprintf("- Raw p < 0.05: %d; BH-FDR < 0.05 (within each tier's own family): %d", n_sig_raw, n_sig_fdr),
    sprintf("- Global proportional-hazards test violated (p < 0.05): %d of %d estimated models", n_ph_violated, n_estimated),
    "",
    "## Results by tier",
    "",
    "| Tier | Models | Estimated | Not estimable | Fit failed | Raw p<0.05 | FDR<0.05 | PH violated |",
    "|---|---:|---:|---:|---:|---:|---:|---:|",
    sprintf("| %s | %d | %d | %d | %d | %d | %d | %d |", tier_summary$tier, tier_summary$n_models, tier_summary$n_estimated,
            tier_summary$n_not_estimable, tier_summary$n_fit_failed, tier_summary$n_sig_raw, tier_summary$n_sig_fdr, tier_summary$n_ph_violated),
    "",
    "PH violations are reported per model in `tables/downstream/survival_models.tsv` (`ph_violated` column) and per-term in `proportional_hazards_tests.tsv`; a violated model's HR is not withdrawn, but is not presented as a time-constant effect either -- interpret its `raw_p`/`HR` as a summary of the average effect over follow-up only, per the caution recorded here.",
    "",
    "## Incremental value beyond CYT",
    "",
    sprintf("Likelihood-ratio comparison of `feature+CYT+age+sex` against a `CYT+age+sex`-only baseline, fit on identical complete-case samples: %d of %d (feature, cancer) combinations were estimable; %d reach BH-FDR<0.05 for the incremental LRT (family: all estimable feature-cancer LRTs). This tests statistical evidence for incremental fit only -- it is not a discrimination or validated-prediction claim; the concordance-index delta (`concordance_delta` in `cyt_incremental_models.tsv`) is reported alongside as a discrimination-scale quantity, and no external validation of any model was performed in this or any prior stage.", sum(incr_dt$status == "ESTIMATED"), nrow(incr_dt), n_incr_sig_fdr),
    "",
    "## Scope and limitations",
    "",
    "- This is a TCGA-only, exploratory, cancer-specific association analysis. No claim of validated prognostic performance or clinical utility is made.",
    "- Sex is omitted automatically in a cancer when it has no within-cancer variation (for example, a sex-specific cancer); such omission is recorded in `dropped_invariant_covariates`, and the EPP threshold uses the number of predictors actually fitted.",
    "- Grade is excluded from every tier because it is entirely missing (0% complete) in a substantial number of cancers; see `reports/23_survival_input_audit.md` for the exact per-cancer count.",
    sprintf("- Stage is included only in the subset of cancers meeting the prespecified %d%% completeness threshold; cancers structurally lacking a usable stage field in this data source are represented only in Tiers 1 and 3.", stage_adequacy_cutoff),
    "- No composite or multi-feature survival score is constructed in this stage; each feature is tested individually, so outcome leakage from a derived composite score does not arise here.",
    "",
    "## Outputs",
    "",
    "- `tables/downstream/clinical_data_completeness_by_cancer.tsv`, `survival_models.tsv`, `proportional_hazards_tests.tsv`, `cyt_incremental_models.tsv`",
    "- `reports/23_survival_input_audit.md`, `reports/24_survival_analysis_report.md` (this file)"
  )
  write_run_artifact(ctx, report24, file.path("reports", "24_survival_analysis_report.md"),
                      "15_survival_analysis", script_path, character(0), function(x, p) writeLines(x, p))
  log_fn("reports/24_survival_analysis_report.md written.")

  record_stage_status(ctx, "15_survival_analysis", "COMPLETE")
  log_fn("Stage 15 complete.")
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=14 or source scripts/run.R in RStudio.")
}
