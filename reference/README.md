The original MCP-counter signatures and signature manifest are included and
were verified against the supplied analysis allowlist.

Place the two original TargetScan 7.2 ZIP files in `reference/targets/`.
Their exact names and digests are in `config/FROM_SCRATCH_ALLOWED_INPUTS.tsv`.

`reported_results/` contains the previously generated CPTAC and additional
analysis tables recovered from the supplied figure-code archive. These are
archival reported results, not newly computed outputs. The main runner does
not silently substitute them for a missing current analysis. Their original
checksum files are retained; new outputs are checked against their own
generated manifests.

Third-party data and signatures retain the terms of their original sources.
