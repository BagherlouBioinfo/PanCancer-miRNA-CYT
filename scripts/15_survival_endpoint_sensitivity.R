# ============================================================================
# 16_survival_endpoint_sensitivity.R
#
# Purpose: survival-endpoint sensitivity analysis. Stage 15 fits its tiered
#          Cox design using overall survival (OS). Because the event process
#          for death from any cause can differ from disease-specific
#          mortality or progression, and late-separating Kaplan-Meier curves
#          can indicate non-proportional hazards, this stage re-fits the
#          identical prespecified design with disease-specific survival
#          (DSS) and the progression-free interval (PFI) as endpoints. Both
#          endpoints are present in the harmonized clinical data
#          (processed/stage03/survival_matched.rds columns DSS/DSS.time and
#          PFI/PFI.time; Liu et al. 2018, Cell 173:400-416, TCGA Pan-Cancer
#          Clinical Data Resource). No new data, covariate, threshold or
#          candidate feature is introduced, and the primary feature and
#          triplet sets are unchanged.
#
# Design (identical to Stage 15 except for the endpoint):
#   - The same primary features (miRNAs and genes, counted from the run's
#     triplet set), the same four tiers (Tier 1: feature + age + sex;
#     Tier 2: + stage; Tier 3: feature + CYT + age + sex; Tier 4: + stage),
#     the same per-cancer standardization, coalesced stage field and
#     stage-adequacy gate, the same events-per-predictor rule (Peduzzi et
#     al. 1996), the same cox.zph test per model, and the same Tier 3 vs.
#     (CYT + age + sex) likelihood-ratio test and concordance change for
#     incremental value beyond CYT.
#   - A nuisance covariate with no within-cancer variation is omitted from
#     that cancer's model and recorded explicitly; the primary feature is
#     never omitted in this way. The events-per-predictor rule uses the
#     number of predictors actually fitted, as in Stage 15.
#   - Per-cancer eligibility is recomputed for each endpoint
#     (N >= min_samples_per_cancer and at least one event on that
#     endpoint), because DSS and PFI completeness and event counts differ
#     from OS.
#   - Features, thresholds and tiers are not re-selected: all four tiers x
#     both endpoints x all primary features x all eligible cancers are
#     fitted and reported.
# Inputs:  <run_root>/processed/stage03/survival_matched.rds,
#          mirna_matched.rds, mrna_matched.rds (same-run);
#          processed/stage04/cyt_scores.tsv (same-run);
#          processed/stage05/sample_annotation_with_purity.tsv (same-run);
#          tables/downstream/survival_models.tsv (Stage 15, same-run; used
#          only for the OS vs. DSS/PFI direction-concordance comparison);
#          the current run's primary triplet set (scripts/lib/run_discovery_source.R).
# Outputs: <run_root>/tables/downstream/survival_models_DSS.tsv,
#          survival_models_PFI.tsv, proportional_hazards_tests_DSS.tsv,
#          proportional_hazards_tests_PFI.tsv, cyt_incremental_models_DSS.tsv,
#          cyt_incremental_models_PFI.tsv, survival_endpoint_completeness.tsv,
#          survival_OS_vs_DSS_concordance.tsv, survival_OS_vs_PFI_concordance.tsv;
#          reports/SURVIVAL_ENDPOINT_SENSITIVITY_REPORT.md.
# Depends on: 03_harmonize_identifiers.R, 04_compute_validate_cyt.R,
#             05_define_eligible_features.R, 09_recurrence_analysis.R,
#             10_purity_adjusted_analysis.R, 15_survival_analysis.R
#             (all same-run).
# ============================================================================

run_stage16 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "16_survival_endpoint_sensitivity")
  source(file.path(ctx$data_root, "scripts", "lib", "clinical_stage_mapping.R"), local = TRUE)
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "run_discovery_source.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "expression_preprocessing.R"), local = TRUE)
  suppressPackageStartupMessages({
    library(data.table)
    library(survival)
  })
  script_path <- file.path(ctx$data_root, "scripts", "15_survival_endpoint_sensitivity.R")
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
               "unique this-run primary feature count is internally consistent (%d miRNAs, %d genes) -- matches Stage 15", n_mirna_expected, n_gene_expected)

  surv <- as.data.table(read_run_artifact(ctx, file.path("processed", "stage03", "survival_matched.rds"), "03_harmonize_identifiers", readRDS))
  sample_annot <- read_run_artifact(ctx, file.path("processed", "stage05", "sample_annotation_with_purity.tsv"), "05_define_eligible_features", fread)
  cyt_scores <- read_run_artifact(ctx, file.path("processed", "stage04", "cyt_scores.tsv"), "04_compute_validate_cyt", fread)
  surv <- merge(surv, sample_annot[, .(sample_id, cancer, CPE, LUMP)], by = "sample_id")
  surv <- merge(surv, cyt_scores[, .(sample_id, CYT)], by = "sample_id")
  surv[, stage_coalesced_raw := ifelse(trimws(ajcc_pathologic_tumor_stage) != "", trimws(ajcc_pathologic_tumor_stage), trimws(clinical_stage))]
  surv[, stage_numeric := stage_to_numeric(stage_coalesced_raw)]
  surv[, gender_clean := trimws(gender)]
  surv[, age := age_at_initial_pathologic_diagnosis]

  min_n <- analysis_config$discovery$min_samples_per_cancer
  stage_adequacy_cutoff <- analysis_config$survival$stage_completeness_adequacy_pct
  epp_min <- analysis_config$survival$events_per_predictor_minimum

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
  # Per-feature winsorization before use as a Cox-model covariate, as in
  # scripts/14_survival_analysis.R.
  mirna_mat <- winsorize_matrix_rows(mirna_mat)
  mrna_mat <- winsorize_matrix_rows(mrna_mat)

  fit_one <- function(d, covars, time_col, event_col) {
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
                  fit = NULL, n = nrow(d), n_events = sum(d[[event_col]] == 1, na.rm = TRUE),
                  requested_covars = requested_covars, used_covars = used_covars, dropped_covars = dropped_covars))
    }

    required_columns <- c(event_col, time_col, used_covars)
    complete <- complete.cases(d[, ..required_columns])
    d_model <- d[complete]
    n_events <- sum(d_model[[event_col]] == 1, na.rm = TRUE)
    min_events_needed <- epp_min * length(used_covars)
    if (n_events < min_events_needed) {
      return(list(status = "NOT ESTIMABLE", reason = sprintf("insufficient events for EPP>=%d (needed >=%d, have %d)", epp_min, min_events_needed, n_events),
                  fit = NULL, n = nrow(d_model), n_events = n_events,
                  requested_covars = requested_covars, used_covars = used_covars, dropped_covars = dropped_covars))
    }
    form <- as.formula(paste(sprintf("Surv(%s, %s) ~", time_col, event_col), paste(used_covars, collapse = " + ")))
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

  endpoints <- list(
    DSS = list(time = "DSS.time", event = "DSS", label = "Disease-Specific Survival"),
    PFI = list(time = "PFI.time", event = "PFI", label = "Progression-Free Interval")
  )

  completeness_rows <- list()
  models_dt_by_endpoint <- list()  # in-memory models_dt per endpoint, used by the concordance summary below
  base_inputs <- file.path(ctx$run_root, c(
    file.path("processed", "stage03", "survival_matched.rds"), file.path("processed", "stage03", "mirna_matched.rds"),
    file.path("processed", "stage03", "mrna_matched.rds"), file.path("processed", "stage04", "cyt_scores.tsv"),
    file.path("processed", "stage05", "sample_annotation_with_purity.tsv")
  ))

  for (ep_name in names(endpoints)) {
    ep <- endpoints[[ep_name]]
    time_col <- ep$time; event_col <- ep$event
    log_fn("=== Endpoint %s (%s): starting. ===", ep_name, ep$label)

    ep_completeness <- surv[, .(
      N = .N,
      n_events = sum(get(event_col) == 1, na.rm = TRUE),
      pct_complete = round(100 * mean(!is.na(get(event_col)) & !is.na(get(time_col))), 1),
      pct_stage_complete_coalesced = round(100 * mean(!is.na(stage_numeric)), 1)
    ), by = cancer][order(cancer)]
    ep_completeness[, stage_adequate := pct_stage_complete_coalesced >= stage_adequacy_cutoff]
    ep_completeness[, eligible_for_survival := N >= min_n & n_events >= 1]
    ep_completeness[, endpoint := ep_name]
    completeness_rows[[length(completeness_rows) + 1]] <- ep_completeness
    eligible_cancers <- ep_completeness[eligible_for_survival == TRUE]$cancer
    log_fn("  %d/%d cancers eligible for %s (N>=%d and >=1 event).", length(eligible_cancers), nrow(ep_completeness), ep_name, min_n)

    model_rows <- list(); ph_rows <- list(); incr_rows <- list()

    for (fi in seq_len(nrow(features))) {
      feat <- features$feature[fi]; ftype <- features$feature_type[fi]
      expr_row <- if (ftype == "miRNA") mirna_mat[feat, ] else mrna_mat[feat, ]
      expr_dt <- data.table(sample_id = names(expr_row), feature_expr = as.numeric(expr_row))

      for (cc in eligible_cancers) {
        d <- surv[cancer == cc, c("sample_id", event_col, time_col, "age", "gender_clean", "stage_numeric", "CYT"), with = FALSE]
        d <- merge(d, expr_dt, by = "sample_id")
        d <- d[!is.na(get(event_col)) & !is.na(get(time_col)) & get(time_col) > 0 & !is.na(age) & gender_clean %in% c("MALE", "FEMALE") & !is.na(feature_expr)]
        # Every attempted (feature, cancer) combination appears in the output
        # for every tier, even when it cannot be fitted (as in
        # scripts/14_survival_analysis.R).
        if (nrow(d) < min_n) {
          for (tname in names(analysis_config$survival$tiers)) {
            model_rows[[length(model_rows) + 1]] <- data.table(
              feature = feat, feature_type = ftype, cancer = cc, tier = tname,
              requested_covariates = paste(analysis_config$survival$tiers[[tname]], collapse = "+"),
              covariates = paste(analysis_config$survival$tiers[[tname]], collapse = "+"),
              fitted_covariates = NA_character_, dropped_invariant_covariates = NA_character_, n_predictors = NA_integer_,
              n = nrow(d), n_events = sum(d[[event_col]] == 1, na.rm = TRUE),
              status = "NOT ESTIMABLE", reason = sprintf("fewer than %d complete-case samples after requiring non-missing %s/age/sex/feature (n=%d)", min_n, ep_name, nrow(d)),
              HR = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, raw_p = NA_real_, ph_global_p = NA_real_
            )
          }
          next
        }
        if (length(unique(d$feature_expr)) < 2) {
          for (tname in names(analysis_config$survival$tiers)) {
            model_rows[[length(model_rows) + 1]] <- data.table(
              feature = feat, feature_type = ftype, cancer = cc, tier = tname,
              requested_covariates = paste(analysis_config$survival$tiers[[tname]], collapse = "+"),
              covariates = paste(analysis_config$survival$tiers[[tname]], collapse = "+"),
              fitted_covariates = NA_character_, dropped_invariant_covariates = NA_character_, n_predictors = NA_integer_,
              n = nrow(d), n_events = sum(d[[event_col]] == 1, na.rm = TRUE),
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
        stage_adeq <- ep_completeness[cancer == cc]$stage_adequate
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
          res <- if (td$applicable) fit_one(td$data, td$covars, time_col, event_col) else
            list(status = "NOT ESTIMABLE", reason = stage_skip_reason, fit = NULL, n = nrow(td$data), n_events = sum(td$data[[event_col]] == 1, na.rm = TRUE),
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

        if ("Tier3" %in% names(tier_defs)) {
          res3 <- fit_one(d, c("feature_z", "cyt_z", "age_z", "sex01"), time_col, event_col)
          res0 <- fit_one(d, c("cyt_z", "age_z", "sex01"), time_col, event_col)
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
      if (fi %% 10 == 0) log_fn("  [%s] Processed %d/%d features.", ep_name, fi, nrow(features))
    }

    models_dt <- rbindlist(model_rows, use.names = TRUE, fill = TRUE)
    systematic_failures <- models_dt[status == "MODEL FIT FAILED", .N, by = .(cancer, tier)][N == nrow(features)]
    run_stage_assert(ctx, log_fn, nrow(systematic_failures) == 0,
                     "[%s] No cancer-tier combination has MODEL FIT FAILED for all %d tested features", ep_name, nrow(features))
    models_dt[, fdr := NA_real_]
    for (tname in unique(models_dt$tier)) {
      sel <- models_dt$tier == tname & grepl("^ESTIMATED", models_dt$status)
      if (sum(sel) > 0) models_dt$fdr[sel] <- p.adjust(models_dt$raw_p[sel], method = "BH")
    }
    models_dt[, ph_violated := !is.na(ph_global_p) & ph_global_p < analysis_config$survival$ph_violation_alpha]
    write_run_artifact(ctx, models_dt, file.path("tables", "downstream", sprintf("survival_models_%s.tsv", ep_name)),
                        "16_survival_endpoint_sensitivity", script_path, base_inputs, fwrite_fn)
    models_dt_by_endpoint[[ep_name]] <- models_dt  # retained in memory for the concordance summary below (not re-read from disk)

    ph_dt <- rbindlist(ph_rows, use.names = TRUE, fill = TRUE)
    write_run_artifact(ctx, ph_dt, file.path("tables", "downstream", sprintf("proportional_hazards_tests_%s.tsv", ep_name)),
                        "16_survival_endpoint_sensitivity", script_path,
                        file.path(ctx$run_root, "tables", "downstream", sprintf("survival_models_%s.tsv", ep_name)), fwrite_fn)

    incr_dt <- rbindlist(incr_rows, use.names = TRUE, fill = TRUE)
    incr_dt[, lrt_fdr := NA_real_]
    sel <- incr_dt$status == "ESTIMATED"
    if (sum(sel) > 0) incr_dt$lrt_fdr[sel] <- p.adjust(incr_dt$lrt_p[sel], method = "BH")
    write_run_artifact(ctx, incr_dt, file.path("tables", "downstream", sprintf("cyt_incremental_models_%s.tsv", ep_name)),
                        "16_survival_endpoint_sensitivity", script_path,
                        file.path(ctx$run_root, "tables", "downstream", sprintf("survival_models_%s.tsv", ep_name)), fwrite_fn)

    log_fn("=== Endpoint %s done: %d rows, %d estimated, %d not estimable, %d model-fit failures, %d PH-violated. ===",
           ep_name, nrow(models_dt), sum(grepl("^ESTIMATED", models_dt$status)),
           sum(models_dt$status == "NOT ESTIMABLE"), sum(models_dt$status == "MODEL FIT FAILED"),
           sum(models_dt$ph_violated, na.rm = TRUE))
  }

  completeness_all <- rbindlist(completeness_rows, use.names = TRUE, fill = TRUE)
  write_run_artifact(ctx, completeness_all, file.path("tables", "downstream", "survival_endpoint_completeness.tsv"),
                      "16_survival_endpoint_sensitivity", script_path, base_inputs, fwrite_fn)

  # --- Concordance-direction summary vs. the primary OS result (Stage 15, same-run) ---
  # os_models is read with read_run_artifact() from Stage 15, which is
  # already COMPLETE in this run. The DSS and PFI models are taken from
  # memory (models_dt_by_endpoint, populated in the per-endpoint loop)
  # rather than re-read from disk, because read_run_artifact() requires the
  # generating stage to be recorded COMPLETE, and this stage records
  # completion only after all of its outputs, including this summary, have
  # been written.
  os_models <- read_run_artifact(ctx, file.path("tables", "downstream", "survival_models.tsv"), "15_survival_analysis", fread)
  dss_models <- models_dt_by_endpoint[["DSS"]]
  pfi_models <- models_dt_by_endpoint[["PFI"]]

  compare_tier3 <- function(alt_dt, alt_name) {
    os3 <- os_models[tier == "Tier3" & grepl("^ESTIMATED", status), .(feature, cancer, OS_HR = HR, OS_p = raw_p, OS_fdr = fdr)]
    alt3 <- alt_dt[tier == "Tier3" & grepl("^ESTIMATED", status), .(feature, cancer, HR, raw_p, fdr)]
    setnames(alt3, c("HR", "raw_p", "fdr"), paste0(alt_name, c("_HR", "_p", "_fdr")))
    m <- merge(os3, alt3, by = c("feature", "cancer"))
    m[, direction_concordant := sign(log(OS_HR)) == sign(log(get(paste0(alt_name, "_HR"))))]
    m
  }
  dss_compare <- compare_tier3(dss_models, "DSS")
  pfi_compare <- compare_tier3(pfi_models, "PFI")
  endpoint_summary <- function(x) {
    list(
      estimated = sum(grepl("^ESTIMATED", x$status)),
      not_estimable = sum(x$status == "NOT ESTIMABLE"),
      fit_failed = sum(x$status == "MODEL FIT FAILED"),
      invariant_adjustments = sum(x$dropped_invariant_covariates != "" & !is.na(x$dropped_invariant_covariates))
    )
  }
  dss_summary <- endpoint_summary(dss_models)
  pfi_summary <- endpoint_summary(pfi_models)
  concordance_inputs <- file.path(ctx$run_root, c(file.path("tables", "downstream", "survival_models.tsv"),
                                                    file.path("tables", "downstream", "survival_models_DSS.tsv"),
                                                    file.path("tables", "downstream", "survival_models_PFI.tsv")))
  write_run_artifact(ctx, dss_compare, file.path("tables", "downstream", "survival_OS_vs_DSS_concordance.tsv"),
                      "16_survival_endpoint_sensitivity", script_path, concordance_inputs, fwrite_fn)
  write_run_artifact(ctx, pfi_compare, file.path("tables", "downstream", "survival_OS_vs_PFI_concordance.tsv"),
                      "16_survival_endpoint_sensitivity", script_path, concordance_inputs, fwrite_fn)

  report <- c(
    "# Survival Endpoint Sensitivity: DSS and PFI",
    "",
    sprintf("Generated by `scripts/15_survival_endpoint_sensitivity.R`, run `%s`. Reuses Stage 15's exact tiered Cox design (Tier1: feature+age+sex; Tier2: +stage; Tier3: feature+CYT+age+sex; Tier4: +stage), same covariates, same standardization, same events-per-predictor rule, same `cox.zph` PH test, same Tier3-vs-baseline LRT for incremental value beyond CYT -- substituting Disease-Specific Survival (DSS) and Progression-Free Interval (PFI) for Overall Survival (OS). No feature, triplet, or threshold was re-selected; this run's own %d primary features (%d miRNAs, %d genes) are unchanged from Stage 15.", ctx$run_id, nrow(features), n_mirna_expected, n_gene_expected),
    "",
    "## 1. Per-cancer eligibility, by endpoint",
    "",
    sprintf("OS (Stage 15, for reference): eligibility computed there. DSS: %d/%d cancers eligible (N>=%d and >=1 DSS event). PFI: %d/%d cancers eligible (N>=%d and >=1 PFI event). Full table: `tables/downstream/survival_endpoint_completeness.tsv`.",
            sum(completeness_all[endpoint == "DSS"]$eligible_for_survival), nrow(completeness_all[endpoint == "DSS"]), min_n,
            sum(completeness_all[endpoint == "PFI"]$eligible_for_survival), nrow(completeness_all[endpoint == "PFI"]), min_n),
    "",
    "## 2. Tier3 (feature+CYT+age+sex) results, by endpoint",
    "",
    sprintf("DSS: %d/%d Tier3 models ESTIMATED; %d/%d PH-violated (global cox.zph p<%.2f); %d/%d reach BH-FDR<0.05.",
            sum(dss_models$tier == "Tier3" & grepl("^ESTIMATED", dss_models$status)),
            sum(dss_models$tier == "Tier3"),
            sum(dss_models$tier == "Tier3" & dss_models$ph_violated, na.rm = TRUE),
            sum(dss_models$tier == "Tier3" & grepl("^ESTIMATED", dss_models$status)),
            analysis_config$survival$ph_violation_alpha,
            sum(dss_models$tier == "Tier3" & dss_models$fdr < 0.05, na.rm = TRUE),
            sum(dss_models$tier == "Tier3" & grepl("^ESTIMATED", dss_models$status))),
    sprintf("PFI: %d/%d Tier3 models ESTIMATED; %d/%d PH-violated; %d/%d reach BH-FDR<0.05.",
            sum(pfi_models$tier == "Tier3" & grepl("^ESTIMATED", pfi_models$status)),
            sum(pfi_models$tier == "Tier3"),
            sum(pfi_models$tier == "Tier3" & pfi_models$ph_violated, na.rm = TRUE),
            sum(pfi_models$tier == "Tier3" & grepl("^ESTIMATED", pfi_models$status)),
            sum(pfi_models$tier == "Tier3" & pfi_models$fdr < 0.05, na.rm = TRUE),
            sum(pfi_models$tier == "Tier3" & grepl("^ESTIMATED", pfi_models$status))),
    "",
    "## 3. Model accounting and invariant nuisance covariates",
    "",
    sprintf("DSS: %d estimated, %d not estimable, %d model-fit failures; %d model rows were fitted after explicitly dropping at least one invariant nuisance covariate. PFI: %d estimated, %d not estimable, %d model-fit failures; %d such invariant-covariate adjustments. Requested, fitted, and dropped covariates are recorded row-by-row in the endpoint model tables; EPP uses the number of predictors actually fitted.",
            dss_summary$estimated, dss_summary$not_estimable, dss_summary$fit_failed, dss_summary$invariant_adjustments,
            pfi_summary$estimated, pfi_summary$not_estimable, pfi_summary$fit_failed, pfi_summary$invariant_adjustments),
    "",
    "## 4. Direction concordance against the primary OS result (Tier3, per feature-cancer pair estimable under both endpoints)",
    "",
    sprintf("DSS vs. OS: %d/%d estimable feature-cancer pairs have a concordant hazard-ratio direction (same side of HR=1). PFI vs. OS: %d/%d concordant. Full pairwise tables: `tables/downstream/survival_OS_vs_DSS_concordance.tsv`, `survival_OS_vs_PFI_concordance.tsv`.",
            sum(dss_compare$direction_concordant), nrow(dss_compare),
            sum(pfi_compare$direction_concordant), nrow(pfi_compare)),
    "",
    "## 5. Interpretation constraints",
    "",
    "This is a sensitivity analysis for alternative survival endpoints; it is not a claim of validated prognostic performance (no external validation exists anywhere in this project), and a PH violation under any endpoint means that endpoint's unadjusted HR should not be reported as a single time-constant effect without qualification. Direction-concordance with OS is reported descriptively; it is not itself a formal test of endpoint equivalence. Sex is omitted automatically when it has no within-cancer variation (for example, in a sex-specific cancer); that omission is transparent and does not remove the feature under study."
  )
  write_run_artifact(ctx, report, file.path("reports", "SURVIVAL_ENDPOINT_SENSITIVITY_REPORT.md"),
                      "16_survival_endpoint_sensitivity", script_path, character(0), function(x, p) writeLines(x, p))
  log_fn("reports/SURVIVAL_ENDPOINT_SENSITIVITY_REPORT.md written.")

  record_stage_status(ctx, "16_survival_endpoint_sensitivity", "COMPLETE")
  log_fn("Stage 16 complete.")
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=15 or source scripts/run.R in RStudio.")
}
