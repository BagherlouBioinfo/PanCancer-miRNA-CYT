# Stage 24 resolves project/result paths; local_external_dir may be set for offline input.
STAGE25_CONFIG <- list(
  project_root = ".",
  discovery_run_id = "analysis",
  output_parent = "results/external_validation",
  expected_primary_pairs = 47L,
  min_n = 30L,
  rho_threshold = 0.30,
  alpha = 0.05,
  winsor_probs = c(0.001, 0.999),
  bootstrap_B = 1000L, # paired-patient bootstrap CIs; descriptive, not test p-values
  seed = 20261017L,
  allow_download = TRUE,
  # Optional folder containing the four exact CPTAC filenames listed in README.
  # Set this for offline use; source files are copied and hashed, never edited.
  local_external_dir = NULL,
  write_figures = TRUE
)
