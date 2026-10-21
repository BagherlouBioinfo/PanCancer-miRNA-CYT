All generated analysis outputs are written here, using stable filenames.
This directory is ignored by git except for this README.

`processed/`, `tables/`, `figures/`, `reports/`, `logs/` and `manifests/`
hold the core analyses. `external_validation/`, `cluster_survival/` and
`mir146/` hold the integrated additional analyses. Final multipanel figures,
tables, legends and count summaries are written to `manuscript/`.

A failed stage can be rerun in this same directory. Completed outputs are
verified before reuse. Partial outputs from a failed stage are not accepted
by downstream readers.
