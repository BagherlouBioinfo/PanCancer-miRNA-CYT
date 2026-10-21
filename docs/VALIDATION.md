# Validation record

Date: 19 October 2026.

Environment: Ubuntu 24.04, R 4.3.3. This differs from the original reported
R 4.3.1 environment; the original reported package versions are recorded
separately in `REFERENCE_PACKAGE_VERSIONS.tsv`.

## Completed checks

1. All R files parsed successfully.
2. The complete `tests/run_tests.R` suite passed: 13 numerical, synthetic,
   runtime, integration and structural test files. The execution transcript
   is retained in `TEST_OUTPUT.txt`. Optional full-data spot checks were
   reported as skipped, not as successful biological validation.
3. All 23 retained core stage scripts and 13 supplied scientific helpers
   passed executable-expression comparisons after the documented path
   substitutions, separate exclusion-summary intermediate path and
   direct-launch and historical-stage diagnostic changes: 36 comparisons.
4. Eighteen named numerical helpers from the three added analyses retained
   identical R function formals and bodies: four cluster-survival helpers,
   three conditional-miR-146 helpers and eleven external-validation helpers.
5. All 13 main/supplementary R figure-building functions retained identical
   function formals and bodies. Original figure configuration values for
   style, dimensions, layout, thresholds and displayed content were unchanged.
6. The recovered original Figure 1 drawing was rerendered at 3600 × 3600
   pixels. Its regenerated PNG was pixel-identical to the supplied original.
   Only the drawing script's output-directory handling was changed.
7. The actual added-analysis entry points ran against synthetic,
   provenance-registered input files: cluster-survival and conditional-miR-146
   outputs, figures, completion markers and output-table checksums passed.
8. The actual runner passed a deliberate failure/restart exercise. Earlier
   output was retained; the repaired failed stage completed; verified
   completed stages were skipped; repeating a producer invalidated consumers.
9. Structural checks covered all 27 registered stages, 111 distinct static
   output paths, known literal producer dependencies, unchanged scientific
   configuration bytes and rejection of a report-only clustering cache.
10. The included MCP-counter signatures and signature manifest matched the
    original input allowlist SHA-256 digests.

## Not established by these checks

- A complete rerun of the full TCGA matrices and all 27 steps was not
  performed during this repository reorganisation. The large source
  matrices were not supplied in the code archive.
- No full original Stage 18 clustering cache, including sample membership
  and original provenance, was available in the recovered archives. It has
  not been fabricated or declared imported. The guarded import tool still
  requires the original files, or script 17 must be executed.
- The complete CPTAC reanalysis was not rerun. Its numerical primitives
  passed their original synthetic checks; the original reported CPTAC tables
  are retained separately with their original checksums.
- The final complete set of numerical manuscript figures was not regenerated
  from a new full-data run. Drawing-code preservation and the Figure 1 pixel
  comparison do not establish equality of results from different data or
  annotation-package versions.
- Existing null-clustering convergence warnings were not repaired by changing
  the algorithm. Their scientific interpretation still requires the
  corresponding original numerical diagnostics.
- Live DGIdb and updated annotation packages can return different records
  from those obtained on the original retrieval date.

The checks support code preservation, numerical helper behaviour and restart
integration. They do not constitute a new validation of every manuscript
claim or a guarantee of identical full-data results on another environment.
