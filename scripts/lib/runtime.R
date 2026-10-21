# Stable results-directory runtime. Scientific estimators are defined in the stage scripts.
suppressPackageStartupMessages(library(digest))
.current_config_hash <- function(data_root) {
  cfg_env <- new.env()
  source(file.path(data_root, "config", "analysis_config.R"), local = cfg_env)
  digest(cfg_env$analysis_config)
}

.script_hash <- function(script_path) {
  if (!file.exists(script_path)) return(NA_character_)
  digest(script_path, algo = "sha256", file = TRUE)
}

# Writes `obj` to `rel_path` (relative to ctx$run_root) using `write_fn`
# (e.g. function(x, path) fwrite(x, path, sep = "\t") or saveRDS), then
# appends one row to the run's artifact_provenance.tsv recording what
# produced it: run ID, stage, generation time, the SHA-256 of the
# generating script (computed at write time), the configuration hash, and
# the SHA-256 of every declared input file. `input_paths` must list the
# files actually read to produce `obj` (same-run artifacts or allowlisted
# raw/reference inputs); an empty vector is valid only for metadata that
# reads no scientific input.
write_run_artifact <- function(ctx, obj, rel_path, stage, script_path, input_paths = character(0), write_fn) {
  out_path <- run_output_path(ctx, rel_path)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  write_fn(obj, out_path)

  input_hashes <- vapply(input_paths, function(p) {
    if (!file.exists(p)) stop(sprintf("write_run_artifact(): declared input '%s' does not exist at write time.", p), call. = FALSE)
    digest(p, algo = "sha256", file = TRUE)
  }, character(1))

  row <- data.frame(
    run_id = ctx$run_id,
    stage = stage,
    relative_path = rel_path,
    generated_at_utc = format(Sys.time(), tz = "UTC", "%Y-%m-%d %H:%M:%OS3"),
    script_sha256 = .script_hash(script_path),
    config_hash = .current_config_hash(ctx$data_root),
    input_paths = paste(input_paths, collapse = ";"),
    input_hashes = paste(input_hashes, collapse = ";"),
    artifact_sha256 = digest(out_path, algo = "sha256", file = TRUE),
    stringsAsFactors = FALSE
  )
  write.table(row, ctx$manifest_path, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
  invisible(out_path)
}

# --- Run-scoped logging (written under ctx$run_root/logs/) ------------------

new_run_logger <- function(ctx, stage_name) {
  log_file <- file.path(ctx$run_root, "logs", paste0(stage_name, ".log"))
  if (file.exists(log_file)) file.remove(log_file)
  log_fn <- function(...) {
    line <- sprintf("[%s] [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"), stage_name, sprintf(...))
    cat(line, "\n")
    cat(line, "\n", file = log_file, append = TRUE)
  }
  attr(log_fn, "log_file") <- log_file
  log_fn
}

run_stage_fail <- function(ctx, log_fn, ...) {
  msg <- sprintf(...)
  log_fn("FATAL: %s", msg)
  stop(msg, call. = FALSE)
}

run_stage_assert <- function(ctx, log_fn, cond, ...) {
  if (!isTRUE(cond)) run_stage_fail(ctx, log_fn, ...)
  log_fn("ASSERT OK: %s", sprintf(...))
}

record_stage_status <- function(ctx, stage, status, detail = "") {
  row <- data.frame(run_id = ctx$run_id, stage = stage, status = status,
                     timestamp_utc = format(Sys.time(), tz = "UTC", "%Y-%m-%d %H:%M:%OS3"), detail = detail,
                     stringsAsFactors = FALSE)
  write.table(row, ctx$status_path, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
}


verify_allowed_raw_input <- function(data_root, path) {
  allowlist_path <- file.path(data_root, "config", "FROM_SCRATCH_ALLOWED_INPUTS.tsv")
  if (!file.exists(allowlist_path)) {
    stop(sprintf("config/FROM_SCRATCH_ALLOWED_INPUTS.tsv not found at %s. A from-scratch run cannot verify any raw input without it.", allowlist_path), call. = FALSE)
  }
  allowlist <- read.delim(allowlist_path, stringsAsFactors = FALSE)
  norm_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  # Resolved from data_root + relative_path so the check is portable across
  # clones; the allowlist's absolute_path column is informational only.
  norm_abs <- file.path(data_root, allowlist$relative_path)
  match_idx <- which(normalizePath(norm_abs, winslash = "/", mustWork = FALSE) == norm_path)
  if (length(match_idx) == 0) {
    stop(sprintf("NOT AN ALLOWED INPUT: '%s' does not appear in config/FROM_SCRATCH_ALLOWED_INPUTS.tsv. A from-scratch run may only read the explicitly declared raw/reference inputs on that allowlist.", path), call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("ALLOWED INPUT MISSING: '%s' is on the allowlist but does not exist on disk.", path), call. = FALSE)
  }
  actual_hash <- digest(path, algo = "sha256", file = TRUE)
  expected_hash <- allowlist$sha256[match_idx[1]]
  if (!identical(actual_hash, expected_hash)) {
    stop(sprintf(
      "ALLOWED INPUT HASH MISMATCH: '%s' has SHA-256 %s now, but the allowlist recorded %s. The file has changed since the allowlist was built; a from-scratch run must not proceed against a raw input whose integrity cannot be confirmed.",
      path, actual_hash, expected_hash
    ), call. = FALSE)
  }
  invisible(TRUE)
}

# Load the resumable execution adapter relative to this sourced helper.
.runtime_files <- unlist(lapply(sys.frames(),function(f) {
  x <- if(!is.null(f$ofile)) f$ofile else f$file
  if(is.character(x)&&length(x)==1L&&grepl('runtime\\.R$',x)) x else NULL
}))
if(!length(.runtime_files)) stop('Cannot locate scripts/lib/runtime.R')
source(file.path(dirname(tail(.runtime_files,1)),'execution.R'),local=TRUE)
rm(.runtime_files)
