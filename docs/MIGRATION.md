# Repository reorganisation

The source of the 23 retained original analysis stages is the supplied
Ubuntu code archive. The 13 scientific helper files were supplied
separately. The external-validation, cluster-survival, conditional-miR-146
and figure/table modules were recovered from their earlier code archives.
No analysis is reconstructed from manuscript prose or plotted values.

## Preserved

- The complete scientific `analysis_config.R`, sample exclusions and raw
  input allowlist are retained without changing their contents.
- Statistical estimators, thresholds, hypothesis families, seeds,
  resampling counts, cohort/feature rules and output table schemas remain
  in their original scientific stage functions.
- The numerical helpers of all three added analyses are unchanged.
- Main and supplementary R drawing functions retain the original panel
  construction, styles, dimensions and labels. The original workflow
  drawing was recovered separately; its drawing instructions are retained.

## Changed

- Executable files follow a single numbered sequence under `scripts/`.
  The unexecuted original quanTIseq stage is not included. The new script
  13 is meta-analysis, not quanTIseq.
- Shared helpers are under `scripts/lib/`, drawing/table functions under
  `scripts/figures/`. Stage source paths and script-hash metadata paths
  point to those new locations.
- Stage 03's expression-layer exclusion summary now has its own intermediate
  filename. Stage 05 reads that file and writes the unchanged combined-summary
  filename. This prevents stage 05 from overwriting an input that it must read
  again on restart; the summary calculations and final columns are unchanged.
- The orchestrator uses a stable `results/` directory. It can rerun a
  failed stage without discarding completed earlier stages. Latest status,
  content hashes, source/configuration fingerprints and dependency
  invalidation replace the old unique-run and file-age restrictions.
- Importing the original expensive clustering result is explicit and
  guarded. The importer preserves its source provenance and does not
  describe imported values as newly generated.
- Additional analyses write to their own stable result subdirectories.
  Restarting one can replace that analysis's generated subdirectory, but
  path checks prevent its output argument from targeting the project or
  arbitrary directories.
- Figure/table readers use newly generated additional-analysis manifests
  when an additional analysis has been rerun. They no longer demand that
  a new table be byte-identical to an archived run solely because it was
  bundled with the earlier figure code.
- Failed figure/table tasks now stop stage 27 instead of leaving an
  apparently complete build with only a warning.

Historical producer IDs and result basenames are intentionally retained;
they are compatibility metadata, not a requirement for unique run folders.

## Audit records

`SOURCE_MAP.tsv` records supplied source digests and destination paths.
`SCIENTIFIC_SOURCE_COMPARISON.tsv` records executable-expression equality
for all retained core scripts and scientific helpers, after the stated
mechanical path substitutions and direct-launch and historical-stage diagnostic changes.
`DEPENDENCY_CHECK.tsv` checks literal upstream producer references against
the new execution dependency graph. `VALIDATION.md` records the tests that
were actually run and their limits.

This is a source-preserving reorganisation, not a full-data numerical
revalidation. Byte-identical figures across operating systems are not
guaranteed: font engines, R/ggplot2 versions and graphics devices can affect
rendering even when the plotting instructions are unchanged.
