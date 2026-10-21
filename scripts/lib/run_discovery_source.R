# ============================================================================
# run_discovery_source.R
#
# Single access point through which downstream stages (11 onward) obtain
# the discovery result of the current run. Every function routes through
# read_run_artifact() (scripts/lib/runtime.R), which verifies run identity, the
# generating stage's COMPLETE status in this run, the artifact hash, parent-input hashes and current configuration. No frozen checkpoint,
# project-level table or other run directory is ever read.
#
# This file contains no hard-coded expected triplet count, miRNA count,
# gene count or triplet identity. The primary discovery set is whatever
# Stage 10 (or Stage 09, for unadjusted evidence) produced in the current
# run, with no comparison against any external number.
# ============================================================================

# Constructs the primary-threshold triplet filename from the configuration,
# using the naming convention of scripts/10_purity_adjusted_analysis.R
# (sprintf("rho%03d", round(T * 100))).
run_primary_triplet_filename <- function(ctx) {
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = (cfg_env <- new.env()))
  cfg <- cfg_env$analysis_config
  thr_tag <- sprintf("rho%03d", round(cfg$thresholds$primary * 100))
  sprintf("triplets_%s_%s.tsv", thr_tag, cfg$purity$primary_metric)
}

run_primary_evidence_filename <- function(ctx, adjustment = "CPE") {
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = (cfg_env <- new.env()))
  cfg <- cfg_env$analysis_config
  thr_tag <- sprintf("rho%03d", round(cfg$thresholds$primary * 100))
  sprintf("cancer_level_evidence_%s_%s.tsv", adjustment, thr_tag)
}

# The primary triplet table (triplets_<threshold>_<metric>.tsv) is written
# by Stage 10 when purity$primary_metric is CPE or LUMP, the only supported
# configuration. Stage 09 also writes an unadjusted table to tables/primary/
# (triplets_<threshold>_unadjusted.tsv), which is not the primary set. The
# function stops if the configured metric is neither CPE nor LUMP.
.run_primary_generating_stage <- function(ctx) {
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = (cfg_env <- new.env()))
  metric <- cfg_env$analysis_config$purity$primary_metric
  if (!metric %in% c("CPE", "LUMP")) {
    stop(sprintf(
      "run_discovery_source.R: analysis_config$purity$primary_metric is '%s', neither CPE nor LUMP. This module assumes the primary triplet table is always a Stage 10 (purity-adjusted) output; an unadjusted primary metric is not currently supported by any downstream (11+) stage and must be resolved before proceeding.",
      metric
    ), call. = FALSE)
  }
  "10_purity_adjusted_analysis"
}

# Loads the primary triplet set (qualifies == TRUE) from this run's
# tables/primary/ directory via read_run_artifact().
load_run_primary_triplets <- function(ctx) {
  suppressPackageStartupMessages(library(data.table))
  rel <- file.path("tables", "primary", run_primary_triplet_filename(ctx))
  dt <- read_run_artifact(ctx, rel, .run_primary_generating_stage(ctx), fread)
  dt[qualifies == TRUE]
}

# Loads the corresponding cancer-level evidence table from this run:
# adjustment = "CPE" or "LUMP" reads Stage 10 output; "unadjusted" reads
# Stage 09 output. Both are read from tables/audit/ under ctx$run_root.
load_run_evidence <- function(ctx, adjustment = c("CPE", "LUMP", "unadjusted")) {
  adjustment <- match.arg(adjustment)
  suppressPackageStartupMessages(library(data.table))
  rel <- file.path("tables", "audit", run_primary_evidence_filename(ctx, adjustment))
  gen_stage <- if (adjustment == "unadjusted") "09_recurrence_analysis" else "10_purity_adjusted_analysis"
  read_run_artifact(ctx, rel, gen_stage, fread)
}

# Rejects a candidate table larger than this run's primary set, so that a
# sensitivity-threshold candidate pool cannot be mistaken for the primary
# set.
assert_using_run_primary_set <- function(ctx, candidate_dt) {
  if (!all(c("mirna_id", "gene") %in% names(candidate_dt))) {
    stop("Candidate table passed to assert_using_run_primary_set() is missing mirna_id/gene columns.", call. = FALSE)
  }
  primary <- load_run_primary_triplets(ctx)
  if (nrow(candidate_dt) > nrow(primary)) {
    stop(sprintf(
      paste(
        "Candidate table has %d rows, more than this run's own primary set's %d --",
        "a downstream stage must never enlarge the primary discovery set",
        "(sensitivity-threshold candidates belong only in clearly-labeled sensitivity outputs)."
      ),
      nrow(candidate_dt), nrow(primary)
    ), call. = FALSE)
  }
  invisible(TRUE)
}
