# ============================================================================
# enrichment_analysis.R
#
# Helpers for Stage 21 (functional enrichment). Every function is a pure
# function of its arguments, with no file path, run context or
# project-specific literal, so each can be unit-tested with synthetic gene
# and annotation data (tests/test_enrichment_and_drug_annotation.R).
#
# Method: hypergeometric over-representation analysis (equivalent to a
# one-sided Fisher's exact test on the 2x2 table {selected, not selected} x
# {in term, not in term}), computed with base R's phyper(). Gene-set
# membership is read from the versioned Bioconductor annotation packages
# org.Hs.eg.db (GO Biological Process) and reactome.db (Reactome pathways).
# ============================================================================

# Maps a (possibly duplicated) vector of gene symbols to Entrez Gene IDs via
# org.Hs.eg.db and classifies every unique input symbol as:
#   - "mapped": exactly one non-NA ENTREZID
#   - "unmapped": no ENTREZID
#   - "ambiguous": more than one distinct ENTREZID
# `duplicated_in_input` and `n_occurrences_in_input` are reported separately
# from mapping status; a symbol can be mapped and duplicated (e.g. a gene
# present in more than one triplet).
#
# org.Hs.eg.db is passed explicitly as `orgdb`, so the function has no
# hidden library() dependency and can be tested with a synthetic stand-in.
map_gene_symbols_to_entrez <- function(symbols_raw, orgdb) {
  symbols_unique <- unique(symbols_raw)
  raw <- suppressMessages(AnnotationDbi::select(orgdb, keys = symbols_unique, keytype = "SYMBOL", columns = "ENTREZID"))
  by_symbol <- split(raw$ENTREZID, raw$SYMBOL)

  occurrence_counts <- table(symbols_raw)

  rows <- lapply(symbols_unique, function(s) {
    ids <- by_symbol[[s]]
    ids <- unique(ids[!is.na(ids)])
    status <- if (length(ids) == 0) "unmapped" else if (length(ids) > 1) "ambiguous" else "mapped"
    entrez_id <- if (length(ids) == 1) ids else NA_character_
    n_occ <- as.integer(occurrence_counts[[s]])
    data.frame(symbol = s, entrez_id = entrez_id, status = status,
               n_occurrences_in_input = n_occ, duplicated_in_input = n_occ > 1,
               all_candidate_entrez_ids = paste(ids, collapse = ";"),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Builds a named list of term_id -> Entrez Gene IDs for GO Biological
# Process terms, restricted to `background_entrez` (terms with no background
# genes are dropped). Uses the propagated (true-path) annotation, in which a
# gene annotated to a child term also counts for every ancestor term.
build_go_bp_gene_sets <- function(background_entrez, orgdb) {
  ann <- suppressMessages(AnnotationDbi::select(orgdb, keys = background_entrez, keytype = "ENTREZID", columns = c("GOALL", "ONTOLOGYALL")))
  ann <- ann[!is.na(ann$GOALL) & ann$ONTOLOGYALL == "BP", c("ENTREZID", "GOALL")]
  ann <- unique(ann)
  split(ann$ENTREZID, ann$GOALL)
}

# Builds a named list of pathway_id -> Entrez Gene IDs for Homo sapiens
# Reactome pathways (prefix "R-HSA-"), restricted to `background_entrez`.
# `extid2pathid` is as.list(reactome.db::reactomeEXTID2PATHID), passed
# explicitly so the function can be tested with a synthetic list.
build_reactome_gene_sets <- function(background_entrez, extid2pathid) {
  hits <- extid2pathid[background_entrez]
  hits <- hits[!vapply(hits, is.null, logical(1))]
  out <- list()
  for (eg in names(hits)) {
    paths <- hits[[eg]]
    paths <- paths[grepl("^R-HSA-", paths)]
    for (p in paths) out[[p]] <- c(out[[p]], eg)
  }
  lapply(out, unique)
}

# Hypergeometric over-representation test, one row per term with at least
# one background gene. Every tested term is returned; FDR is not applied
# here. The caller applies p.adjust() once per prespecified family (all GO
# BP terms; all Reactome pathways).
#
# Urn model: N = |background|, K = |term genes in background|,
# n = |selected|, k = |selected in term|;
# p = P(X >= k) = phyper(k - 1, K, N - K, n, lower.tail = FALSE).
hypergeometric_ora <- function(selected, background, term_gene_sets) {
  not_in_bg <- setdiff(selected, background)
  if (length(not_in_bg) > 0) {
    stop(sprintf("hypergeometric_ora(): %d selected gene(s) are not in the background universe (%s) -- every selected gene must be a member of its own background by construction.",
                 length(not_in_bg), paste(not_in_bg, collapse = ",")), call. = FALSE)
  }
  N <- length(background)
  n <- length(selected)
  term_ids <- names(term_gene_sets)
  rows <- vector("list", length(term_ids))
  for (i in seq_along(term_ids)) {
    term_genes <- intersect(term_gene_sets[[i]], background)
    K <- length(term_genes)
    if (K == 0) next
    hits <- intersect(selected, term_genes)
    k <- length(hits)
    p <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
    rows[[i]] <- data.frame(
      term_id = term_ids[i], n_background_in_term = K, n_selected_in_term = k,
      n_selected = n, n_background = N,
      expected_by_chance = n * K / N,
      fold_enrichment = if (n > 0 && K > 0) (k / n) / (K / N) else NA_real_,
      p_value = p, gene_hits = paste(sort(hits), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

# Leave-one-gene-out sensitivity: for every gene in `selected`, recomputes
# hypergeometric_ora() with that gene removed, for the supplied `term_ids`
# (typically the terms significant in the primary result). Returns one row
# per (term_id, gene_removed) with the recomputed p-value and a
# significance-flip indicator.
#
# `original_p_by_term` and `original_fdr_by_term` are named vectors (or
# lists), term_id -> value, giving each term's raw p-value and BH-adjusted
# FDR in the primary result. `fdr_cutoff` is the threshold used for the
# primary result.
#
# The leave-one-out p-value is compared against a per-term raw-p threshold
# obtained by inverting that term's original BH scale factor,
# `loo_p_threshold = original_p_value * (fdr_cutoff / original_fdr)`,
# rather than by re-running BH on the significant-terms subset, which would
# understate the correction because that subset was selected for
# significance. This is a linear approximation of a full family refit that
# assumes single-gene removal does not reorder the family's p-value ranks;
# the assumption should be stated alongside any reported result.
#
# `significance_flip` is TRUE for a term that was significant originally
# and is no longer below its threshold with one gene removed, i.e. its
# significance depends on that single gene.
#
# An empty `term_ids` returns a zero-row data.frame.
leave_one_gene_out_check <- function(selected, background, term_gene_sets, term_ids,
                                      original_p_by_term, original_fdr_by_term, fdr_cutoff = 0.05) {
  empty_cols <- function() data.frame(
    term_id = character(0), gene_removed = character(0), p_value_without_gene = numeric(0),
    original_p_value = numeric(0), original_fdr = numeric(0), loo_p_threshold = numeric(0),
    loo_significant = logical(0), significance_flip = logical(0), stringsAsFactors = FALSE
  )
  if (length(term_ids) == 0) return(empty_cols())
  sets_subset <- term_gene_sets[term_ids]
  rows <- list()
  for (g in selected) {
    reduced_selected <- setdiff(selected, g)
    res <- hypergeometric_ora(reduced_selected, background, sets_subset)
    for (tid in term_ids) {
      p_new <- if (tid %in% res$term_id) res$p_value[res$term_id == tid] else NA_real_
      orig_p <- unname(original_p_by_term[[tid]])
      orig_fdr <- unname(original_fdr_by_term[[tid]])
      loo_thr <- orig_p * (fdr_cutoff / orig_fdr)
      loo_sig <- !is.na(p_new) && p_new < loo_thr
      flip <- isTRUE(!is.na(orig_fdr) && orig_fdr < fdr_cutoff) && !loo_sig
      rows[[length(rows) + 1]] <- data.frame(
        term_id = tid, gene_removed = g, p_value_without_gene = p_new,
        original_p_value = orig_p, original_fdr = orig_fdr, loo_p_threshold = loo_thr,
        loo_significant = loo_sig, significance_flip = flip, stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
