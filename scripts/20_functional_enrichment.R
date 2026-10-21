# ============================================================================
# 21_functional_enrichment.R
#
# Purpose: over-representation analysis (ORA) of the unique target genes of
#          the current run's primary triplets against GO Biological Process
#          and Reactome. The analysis is descriptive; it makes no claim of
#          pathway activation, causal regulation or independent validation
#          of miRNA-mediated regulation.
# Inputs:  <run_root>/tables/primary/triplets_rho030_CPE.tsv (primary
#          triplets, via load_run_primary_triplets(ctx));
#          <run_root>/processed/stage07/pair_universe_threshold_0.20.tsv
#          (background universe, see below); <run_root>/tables/audit/
#          immune_primary_gene_signature_overlap.tsv (Stage 12, same-run;
#          gene list for the marker-overlap sensitivity check).
# Outputs: <run_root>/tables/downstream/enrichment_gene_universe_provenance.tsv,
#          enrichment_gene_id_mapping.tsv, enrichment_go_bp_results.tsv,
#          enrichment_reactome_results.tsv,
#          enrichment_sensitivity_marker_removed_go_bp.tsv,
#          enrichment_sensitivity_marker_removed_reactome.tsv,
#          enrichment_sensitivity_leave_one_out.tsv,
#          enrichment_terms_robust_to_leave_one_out.tsv,
#          enrichment_terms_single_gene_driven.tsv;
#          reports/25_functional_enrichment_report.md.
# Method:  gene symbols are mapped to Entrez Gene IDs via org.Hs.eg.db with
#          an explicit mapping table classifying each symbol as mapped,
#          unmapped, ambiguous or duplicated
#          (R/enrichment_analysis.R::map_gene_symbols_to_entrez()). GO
#          Biological Process terms use org.Hs.eg.db's propagated
#          (GOALL/ONTOLOGYALL) annotation; Reactome pathways come from
#          reactome.db, restricted to Homo sapiens (R-HSA-) pathways. Both
#          gene-set collections are restricted to the background before
#          testing, and terms with no background genes are not tested.
#          Significance is assessed by the hypergeometric test
#          (R/enrichment_analysis.R::hypergeometric_ora(), equivalent to a
#          one-sided Fisher's exact test) for every term with at least one
#          background gene, and the complete result table is written.
#          BH-FDR is applied separately within the two prespecified
#          families (GO_BP, Reactome; analysis_config$enrichment). KEGG is
#          not queried because of its restrictive redistribution terms.
#
#          Background universe: the unique genes that could have been
#          tested by the Stage 08 three-arm rule, i.e. the `gene` column of
#          Stage 07's loosest-threshold (0.20, union) pair universe. Genes
#          outside this set could never have entered a primary triplet,
#          because Stage 08 evaluates only pairs in this universe; neither
#          all expressed genes nor the whole genome is used.
# Depends on: 07_target_universe_construction.R, 09_recurrence_analysis.R,
#             10_purity_adjusted_analysis.R, 12_immune_composition_analysis.R
#             (all same-run).
# ============================================================================

run_stage21 <- function(ctx) {
  log_fn <- new_run_logger(ctx, "21_functional_enrichment")
  source(file.path(ctx$data_root, "config", "analysis_config.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "run_discovery_source.R"), local = TRUE)
  source(file.path(ctx$data_root, "scripts", "lib", "enrichment_analysis.R"), local = TRUE)
  suppressPackageStartupMessages({
    library(data.table)
    library(org.Hs.eg.db)
    library(reactome.db)
    library(GO.db)
    library(AnnotationDbi)
  })
  script_path <- file.path(ctx$data_root, "scripts", "20_functional_enrichment.R")
  fwrite_fn <- function(x, p) fwrite(x, p, sep = "\t")

  # ==========================================================================
  # 1. Selected genes: the current run's primary triplets.
  # ==========================================================================
  triplets <- load_run_primary_triplets(ctx)
  selected_genes_raw <- triplets$gene
  selected_genes <- sort(unique(selected_genes_raw))
  run_stage_assert(ctx, log_fn, length(selected_genes) > 0, "This run's primary triplets contain at least one target gene (%d found)", length(selected_genes))
  log_fn("Selected gene set: %d unique genes from %d primary triplets (%d unique miRNAs).",
         length(selected_genes), nrow(triplets), uniqueN(triplets$mirna_id))

  # ==========================================================================
  # 2. Background universe: unique genes in Stage 07's loosest-threshold
  #    (union) pair universe (see the file header).
  # ==========================================================================
  loosest_threshold <- min(analysis_config$thresholds$all)
  union_rel <- file.path("processed", "stage07", sprintf("pair_universe_threshold_%.2f.tsv", loosest_threshold))
  union_pairs <- read_run_artifact(ctx, union_rel, "07_target_universe_construction", fread)
  background_genes <- sort(unique(union_pairs$gene))
  log_fn("Background universe (Stage 07 threshold-%.2f union pair universe): %d unique genes.", loosest_threshold, length(background_genes))

  missing_from_bg <- setdiff(selected_genes, background_genes)
  run_stage_assert(ctx, log_fn, length(missing_from_bg) == 0,
                    "Every one of this run's %d selected genes is a member of the %d-gene background universe (0 missing)",
                    length(selected_genes), length(background_genes))

  provenance_dt <- data.table(
    gene_set = c("selected", "background"),
    n_genes = c(length(selected_genes), length(background_genes)),
    source = c(
      sprintf("this run's %d primary triplets (tables/primary/%s), unique gene column", nrow(triplets), run_primary_triplet_filename(ctx)),
      sprintf("Stage 07 threshold-%.2f (loosest/union) pair universe, unique gene column -- the exact set of genes that could ever have been tested by Stage 08's three-arm rule", loosest_threshold)
    )
  )
  provenance_inputs <- file.path(ctx$run_root, c(file.path("tables", "primary", run_primary_triplet_filename(ctx)), union_rel))
  write_run_artifact(ctx, provenance_dt, file.path("tables", "downstream", "enrichment_gene_universe_provenance.tsv"),
                      "21_functional_enrichment", script_path, provenance_inputs, fwrite_fn)

  # ==========================================================================
  # 3. Gene-ID mapping (explicit, tabulated -- mapped/unmapped/ambiguous/duplicated)
  # ==========================================================================
  org_version <- tryCatch(as.character(packageVersion("org.Hs.eg.db")), error = function(e) NA_character_)
  reactome_version <- tryCatch(as.character(packageVersion("reactome.db")), error = function(e) NA_character_)
  godb_version <- tryCatch(as.character(packageVersion("GO.db")), error = function(e) NA_character_)
  log_fn("Annotation package versions: org.Hs.eg.db %s, reactome.db %s, GO.db %s.", org_version, reactome_version, godb_version)

  bg_map <- as.data.table(map_gene_symbols_to_entrez(background_genes, org.Hs.eg.db))
  bg_map[, in_selected_set := symbol %in% selected_genes]
  write_run_artifact(ctx, bg_map, file.path("tables", "downstream", "enrichment_gene_id_mapping.tsv"),
                      "21_functional_enrichment", script_path, file.path(ctx$run_root, union_rel), fwrite_fn)

  n_unmapped_sel <- sum(bg_map$status == "unmapped" & bg_map$in_selected_set)
  n_ambiguous_sel <- sum(bg_map$status == "ambiguous" & bg_map$in_selected_set)
  log_fn("Gene-ID mapping: %d/%d background genes mapped, %d unmapped, %d ambiguous (selected-set genes: %d unmapped, %d ambiguous).",
         sum(bg_map$status == "mapped"), nrow(bg_map), sum(bg_map$status == "unmapped"), sum(bg_map$status == "ambiguous"),
         n_unmapped_sel, n_ambiguous_sel)

  background_entrez <- unique(bg_map$entrez_id[bg_map$status == "mapped"])
  selected_entrez <- unique(bg_map$entrez_id[bg_map$status == "mapped" & bg_map$in_selected_set])
  run_stage_assert(ctx, log_fn, all(selected_entrez %in% background_entrez),
                    "Every mapped selected-gene Entrez ID is a member of the mapped background Entrez ID universe")

  # ==========================================================================
  # 4-9. GO BP and Reactome ORA (every tested term written; BH-FDR applied
  # separately within each prespecified family).
  # ==========================================================================
  go_sets <- build_go_bp_gene_sets(background_entrez, org.Hs.eg.db)
  log_fn("GO Biological Process gene sets built: %d terms with >=1 background gene.", length(go_sets))
  go_res <- as.data.table(hypergeometric_ora(selected_entrez, background_entrez, go_sets))
  go_res[, fdr := p.adjust(p_value, method = "BH")]
  go_term_names <- suppressMessages(AnnotationDbi::select(GO.db, keys = go_res$term_id, keytype = "GOID", columns = "TERM"))
  go_res[, term_name := go_term_names$TERM[match(term_id, go_term_names$GOID)]]
  setorder(go_res, p_value)
  write_run_artifact(ctx, go_res, file.path("tables", "downstream", "enrichment_go_bp_results.tsv"),
                      "21_functional_enrichment", script_path, file.path(ctx$run_root, union_rel), fwrite_fn)
  log_fn("GO BP ORA: %d terms tested; %d significant at FDR<%.2f.", nrow(go_res), sum(go_res$fdr < analysis_config$enrichment$fdr_cutoff), analysis_config$enrichment$fdr_cutoff)

  reactome_extid2pathid <- as.list(reactomeEXTID2PATHID)
  reactome_pathid2name <- as.list(reactomePATHID2NAME)
  reactome_sets <- build_reactome_gene_sets(background_entrez, reactome_extid2pathid)
  log_fn("Reactome gene sets built: %d Homo sapiens pathways with >=1 background gene.", length(reactome_sets))
  reactome_res <- as.data.table(hypergeometric_ora(selected_entrez, background_entrez, reactome_sets))
  reactome_res[, fdr := p.adjust(p_value, method = "BH")]
  reactome_res[, term_name := vapply(term_id, function(x) { nm <- reactome_pathid2name[[x]]; if (is.null(nm)) NA_character_ else nm }, character(1))]
  setorder(reactome_res, p_value)
  write_run_artifact(ctx, reactome_res, file.path("tables", "downstream", "enrichment_reactome_results.tsv"),
                      "21_functional_enrichment", script_path, file.path(ctx$run_root, union_rel), fwrite_fn)
  log_fn("Reactome ORA: %d pathways tested; %d significant at FDR<%.2f.", nrow(reactome_res), sum(reactome_res$fdr < analysis_config$enrichment$fdr_cutoff), analysis_config$enrichment$fdr_cutoff)

  if (nrow(go_res) == 0 && nrow(reactome_res) == 0) {
    log_fn("No GO BP or Reactome term had any background gene -- 0 terms tested. This negative result is retained, not treated as an error.")
  }
  if (sum(go_res$fdr < analysis_config$enrichment$fdr_cutoff) == 0 && sum(reactome_res$fdr < analysis_config$enrichment$fdr_cutoff) == 0) {
    log_fn("No pathway survives FDR<%.2f in either family -- this negative result is retained and reported, not discarded.", analysis_config$enrichment$fdr_cutoff)
  }

  # ==========================================================================
  # 10. Sensitivity checks
  # ==========================================================================
  # (a) Remove the MCPcounter marker-overlap genes (read from the current
  # run's Stage 12 output).
  marker_overlap <- read_run_artifact(ctx, file.path("tables", "audit", "immune_primary_gene_signature_overlap.tsv"), "12_immune_composition_analysis", fread)
  marker_genes <- sort(unique(marker_overlap$primary_gene))
  log_fn("Sensitivity check (a): removing %d MCPcounter marker-overlap gene(s) found in THIS run's own Stage 12 output: %s.", length(marker_genes), paste(marker_genes, collapse = ", "))
  selected_no_markers <- setdiff(selected_genes, marker_genes)
  selected_entrez_no_markers <- unique(bg_map$entrez_id[bg_map$status == "mapped" & bg_map$symbol %in% selected_no_markers])

  go_res_nomarkers <- as.data.table(hypergeometric_ora(selected_entrez_no_markers, background_entrez, go_sets))
  go_res_nomarkers[, fdr := p.adjust(p_value, method = "BH")]
  setorder(go_res_nomarkers, p_value)
  write_run_artifact(ctx, go_res_nomarkers, file.path("tables", "downstream", "enrichment_sensitivity_marker_removed_go_bp.tsv"),
                      "21_functional_enrichment", script_path, file.path(ctx$run_root, union_rel), fwrite_fn)

  reactome_res_nomarkers <- as.data.table(hypergeometric_ora(selected_entrez_no_markers, background_entrez, reactome_sets))
  reactome_res_nomarkers[, fdr := p.adjust(p_value, method = "BH")]
  reactome_res_nomarkers[, term_name := vapply(term_id, function(x) { nm <- reactome_pathid2name[[x]]; if (is.null(nm)) NA_character_ else nm }, character(1))]
  setorder(reactome_res_nomarkers, p_value)
  write_run_artifact(ctx, reactome_res_nomarkers, file.path("tables", "downstream", "enrichment_sensitivity_marker_removed_reactome.tsv"),
                      "21_functional_enrichment", script_path, file.path(ctx$run_root, union_rel), fwrite_fn)

  # (b) Genes present in several triplets: ORA operates on set membership
  # and is not weighted by recurrence. The deduplicated Entrez set is
  # reconstructed independently from the raw (duplicate-containing) triplet
  # gene column, ORA is rerun on it, and the result is required to be
  # numerically identical to the production result (`go_res`,
  # `reactome_res`, built from `selected_entrez`).
  raw_symbol_status <- bg_map$status[match(selected_genes_raw, bg_map$symbol)]
  raw_mapped_entrez <- bg_map$entrez_id[match(selected_genes_raw, bg_map$symbol)]
  raw_mapped_entrez <- raw_mapped_entrez[!is.na(raw_symbol_status) & raw_symbol_status == "mapped"]
  independent_dedup_entrez <- sort(unique(raw_mapped_entrez))
  run_stage_assert(ctx, log_fn, setequal(independent_dedup_entrez, selected_entrez),
                    "Independently reconstructing the deduplicated selected-gene Entrez set directly from this run's raw (pre-dedup) triplet gene column (%d raw entries -> %d unique mapped IDs) exactly matches the Entrez set actually used for the ORA above (%d IDs)",
                    length(selected_genes_raw), length(independent_dedup_entrez), length(selected_entrez))

  verify_ora_matches_independent_recompute <- function(production_res, term_sets, family_label) {
    independent_res <- as.data.table(hypergeometric_ora(independent_dedup_entrez, background_entrez, term_sets))
    independent_res[, fdr := p.adjust(p_value, method = "BH")]
    cmp <- merge(production_res[, .(term_id, n_selected, n_selected_in_term, p_value, fdr)],
                 independent_res[, .(term_id, n_selected, n_selected_in_term, p_value, fdr)],
                 by = "term_id", suffixes = c("_production", "_independent"))
    ok <- nrow(cmp) == nrow(production_res) &&
      all(cmp$n_selected_production == cmp$n_selected_independent) &&
      all(cmp$n_selected_in_term_production == cmp$n_selected_in_term_independent) &&
      isTRUE(all.equal(cmp$p_value_production, cmp$p_value_independent, tolerance = 1e-10)) &&
      isTRUE(all.equal(cmp$fdr_production, cmp$fdr_independent, tolerance = 1e-10))
    run_stage_assert(ctx, log_fn, ok,
                      "%s: the production ORA result (selected_entrez, n=%d) is numerically identical (term count, n, k, p, FDR, tolerance 1e-10) to an independently-reconstructed ORA run computed fresh from the raw triplet gene column (n=%d) for all %d terms",
                      family_label, length(selected_entrez), length(independent_dedup_entrez), nrow(cmp))
  }
  verify_ora_matches_independent_recompute(go_res, go_sets, "GO BP")
  verify_ora_matches_independent_recompute(reactome_res, reactome_sets, "Reactome")

  # (c) Leave-one-gene-out: tests whether any term significant in the
  # primary (all-genes) result depends on a single gene, reporting an
  # explicit significance-flip column
  # (R/enrichment_analysis.R::leave_one_gene_out_check()).
  sig_go_terms <- go_res$term_id[go_res$fdr < analysis_config$enrichment$fdr_cutoff]
  sig_reactome_terms <- reactome_res$term_id[reactome_res$fdr < analysis_config$enrichment$fdr_cutoff]
  go_p_by_term <- setNames(go_res$p_value, go_res$term_id); go_fdr_by_term <- setNames(go_res$fdr, go_res$term_id)
  reactome_p_by_term <- setNames(reactome_res$p_value, reactome_res$term_id); reactome_fdr_by_term <- setNames(reactome_res$fdr, reactome_res$term_id)
  loo_go <- as.data.table(leave_one_gene_out_check(selected_entrez, background_entrez, go_sets, sig_go_terms,
                                                     go_p_by_term, go_fdr_by_term, analysis_config$enrichment$fdr_cutoff))
  loo_reactome <- as.data.table(leave_one_gene_out_check(selected_entrez, background_entrez, reactome_sets, sig_reactome_terms,
                                                           reactome_p_by_term, reactome_fdr_by_term, analysis_config$enrichment$fdr_cutoff))
  loo_go[, family := "GO_BP"]; loo_reactome[, family := "Reactome"]
  loo_all <- rbindlist(list(loo_go, loo_reactome), use.names = TRUE, fill = TRUE)
  write_run_artifact(ctx, loo_all, file.path("tables", "downstream", "enrichment_sensitivity_leave_one_out.tsv"),
                      "21_functional_enrichment", script_path, file.path(ctx$run_root, union_rel), fwrite_fn)

  n_sig_terms_total <- length(sig_go_terms) + length(sig_reactome_terms)
  terms_with_any_flip <- if (nrow(loo_all) > 0) unique(loo_all$term_id[loo_all$significance_flip %in% TRUE]) else character(0)
  terms_robust_to_every_removal <- setdiff(c(sig_go_terms, sig_reactome_terms), terms_with_any_flip)
  log_fn("Sensitivity check (c): leave-one-gene-out recomputed for %d significant term(s) (%d GO_BP, %d Reactome; %d single-gene-removal comparisons total). %d of %d significant terms show a significance flip on at least one single-gene removal (single-gene-driven, not robust pathway-level evidence); %d remain significant under EVERY leave-one-gene-out removal (robust).",
         n_sig_terms_total, length(sig_go_terms), length(sig_reactome_terms), nrow(loo_all),
         length(terms_with_any_flip), n_sig_terms_total, length(terms_robust_to_every_removal))

  robust_terms_dt <- rbindlist(list(
    if (length(sig_go_terms) > 0) go_res[term_id %in% intersect(sig_go_terms, terms_robust_to_every_removal), .(term_id, term_name, family = "GO_BP", p_value, fdr)] else NULL,
    if (length(sig_reactome_terms) > 0) reactome_res[term_id %in% intersect(sig_reactome_terms, terms_robust_to_every_removal), .(term_id, term_name, family = "Reactome", p_value, fdr)] else NULL
  ), use.names = TRUE, fill = TRUE)
  single_gene_driven_dt <- rbindlist(list(
    if (length(sig_go_terms) > 0) go_res[term_id %in% intersect(sig_go_terms, terms_with_any_flip), .(term_id, term_name, family = "GO_BP", p_value, fdr)] else NULL,
    if (length(sig_reactome_terms) > 0) reactome_res[term_id %in% intersect(sig_reactome_terms, terms_with_any_flip), .(term_id, term_name, family = "Reactome", p_value, fdr)] else NULL
  ), use.names = TRUE, fill = TRUE)
  if (nrow(single_gene_driven_dt) > 0) {
    driving_genes <- loo_all[significance_flip %in% TRUE, .(driving_genes_removed = paste(sort(unique(gene_removed)), collapse = ";")), by = term_id]
    single_gene_driven_dt <- merge(single_gene_driven_dt, driving_genes, by = "term_id", all.x = TRUE)
  } else {
    single_gene_driven_dt$driving_genes_removed <- character(0)
  }
  write_run_artifact(ctx, robust_terms_dt, file.path("tables", "downstream", "enrichment_terms_robust_to_leave_one_out.tsv"),
                      "21_functional_enrichment", script_path, file.path(ctx$run_root, union_rel), fwrite_fn)
  write_run_artifact(ctx, single_gene_driven_dt, file.path("tables", "downstream", "enrichment_terms_single_gene_driven.tsv"),
                      "21_functional_enrichment", script_path, file.path(ctx$run_root, union_rel), fwrite_fn)

  # ==========================================================================
  # Concise report
  # ==========================================================================
  n_sig_go <- sum(go_res$fdr < analysis_config$enrichment$fdr_cutoff)
  n_sig_reactome <- sum(reactome_res$fdr < analysis_config$enrichment$fdr_cutoff)
  report <- c(
    "# Stage 21: Functional Enrichment (Over-Representation Analysis)",
    "",
    sprintf("Generated by `scripts/20_functional_enrichment.R`, run `%s`, %s.", ctx$run_id, format(Sys.time(), "%Y-%m-%d %H:%M:%OS0")),
    "",
    "**This is descriptive over-representation analysis only.** No pathway",
    "activation, causal regulation, or therapeutic action is claimed. No",
    "result here is presented as independent validation of the discovery",
    "pipeline's miRNA-regulatory claim -- the tested genes were themselves",
    "selected via the same CYT-linked structure the enrichment result would",
    "otherwise appear to \"confirm\".",
    "",
    "## Gene sets",
    sprintf("- Selected: %d unique genes from this run's %d primary triplets.", length(selected_genes), nrow(triplets)),
    sprintf("- Background: %d unique genes (Stage 07 threshold-%.2f union pair universe).", length(background_genes), loosest_threshold),
    sprintf("- Gene-ID mapping: %d/%d background genes mapped to a unique Entrez ID (%d unmapped, %d ambiguous).",
            sum(bg_map$status == "mapped"), nrow(bg_map), sum(bg_map$status == "unmapped"), sum(bg_map$status == "ambiguous")),
    "",
    "## GO Biological Process (full result table: tables/downstream/enrichment_go_bp_results.tsv)",
    sprintf("%d terms tested; %d significant at FDR<%.2f.", nrow(go_res), n_sig_go, analysis_config$enrichment$fdr_cutoff),
    "",
    "## Reactome (full result table: tables/downstream/enrichment_reactome_results.tsv)",
    sprintf("%d pathways tested; %d significant at FDR<%.2f. KEGG was not queried (see limitations).", nrow(reactome_res), n_sig_reactome, analysis_config$enrichment$fdr_cutoff),
    "",
    "## Sensitivity checks",
    sprintf("- Marker-overlap removal (%s): see enrichment_sensitivity_marker_removed_*.tsv.", paste(marker_genes, collapse = ", ")),
    "- Duplicate-collapse invariance: the production ORA result was verified numerically identical (term count, n, k, p, FDR, tolerance 1e-10) to an independently-reconstructed ORA run computed fresh from this run's raw, pre-dedup triplet gene column -- not merely an identity comparison of the same object to itself.",
    sprintf("- Leave-one-gene-out: %d recomputation(s) for the %d term(s) significant in the primary result.", nrow(loo_all), n_sig_terms_total),
    sprintf("  **%d of %d significant terms (%s) show a significance flip when at least one contributing gene is removed -- these are SINGLE-GENE-DRIVEN and are NOT described anywhere in this report as robust pathway-level evidence** (full list, with the specific driving gene(s): tables/downstream/enrichment_terms_single_gene_driven.tsv).",
            length(terms_with_any_flip), n_sig_terms_total,
            if (n_sig_terms_total > 0) sprintf("%.0f%%", 100 * length(terms_with_any_flip) / n_sig_terms_total) else "N/A"),
    sprintf("  %d of %d significant terms remain significant under EVERY single-gene removal (tables/downstream/enrichment_terms_robust_to_leave_one_out.tsv) -- these are the only terms this report treats as having withstood the leave-one-out check; even these are still descriptive over-representation, not independent validation (see disclaimer above).",
            length(terms_robust_to_every_removal), n_sig_terms_total),
    "",
    "## Limitations",
    "- KEGG was not queried because of its restrictive redistribution terms; only GO Biological Process and Reactome are reported.",
    "- Gene-set membership was read directly from the versioned Bioconductor annotation packages org.Hs.eg.db and reactome.db, and significance was assessed with the standard hypergeometric test.",
    "- This analysis does not, and cannot, distinguish whether an enriched term reflects direct miRNA-target regulation, shared cellular origin, or coordinated immune-cell infiltration; no such distinction is claimed anywhere in this output."
  )
  write_run_artifact(ctx, report, file.path("reports", "25_functional_enrichment_report.md"),
                      "21_functional_enrichment", script_path, character(0), function(x, p) writeLines(x, p))
  log_fn("reports/25_functional_enrichment_report.md written.")

  record_stage_status(ctx, "21_functional_enrichment", "COMPLETE")
  log_fn("Stage complete.")
  invisible(ctx)
}

if (sys.nframe() == 0) {
  stop("Use Rscript scripts/run.R --stages=20 or source scripts/run.R in RStudio.")
}
