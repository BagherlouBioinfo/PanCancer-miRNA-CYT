# Independent post-run extension. It never modifies the discovery run or its status.
# Primary question: external BRCA association support for a FIXED discovery pair set.
# It is not a replication of CPE-adjusted discovery, recurrence across >=3 cancers,
# causal targeting, prognostic prediction, or the cellular-context interaction model.

s25_bool <- function(x) {
  y <- toupper(trimws(as.character(x)))
  if (any(!is.na(y) & !y %in% c("TRUE", "FALSE", "1", "0"))) stop("Invalid logical flag.")
  y %in% c("TRUE", "1") & !is.na(y)
}
s25_need <- function(x, cols, label) {
  missing <- setdiff(cols, names(x))
  if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ", "))
}
s25_key <- function(x, cols) do.call(paste, c(x[cols], sep = "\034"))
s25_unique <- function(x, cols, label) {
  if (anyDuplicated(s25_key(x, cols))) stop(label, " has duplicate keys: ", paste(cols, collapse = ", "))
}
s25_read <- function(path) as.data.frame(data.table::fread(path, check.names = FALSE, na.strings = c("NA", "NaN", "")))
s25_write <- function(x, path) data.table::fwrite(x, path, sep = "\t", na = "NA", quote = "auto")
s25_sha <- function(path) digest::digest(file = path, algo = "sha256")
s25_winsor <- function(x, probs) {
  ok <- is.finite(x)
  if (sum(ok) > 1L) {
    limits <- stats::quantile(x[ok], probs, names = FALSE, type = 7)
    x[ok] <- pmin(pmax(x[ok], limits[1L]), limits[2L])
  }
  x
}
s25_family_adjust <- function(p, planned_n, method = "BH") {
  out <- rep(NA_real_, length(p)); ok <- is.finite(p)
  if (any(ok)) out[ok] <- stats::p.adjust(p[ok], method = method, n = planned_n)
  out
}
s25_cor <- function(x, y, min_n, z = NULL, probs = c(.001, .999)) {
  # Complete observations selected BEFORE ranking; no missing value receives a rank.
  if (is.null(z)) ok <- is.finite(x) & is.finite(y)
  else ok <- is.finite(x) & is.finite(y) & is.finite(z)
  x <- x[ok]; y <- y[ok]; n <- length(x)
  empty <- function(reason) list(n = n, rho = NA_real_, p_raw = NA_real_, status = reason)
  if (n < min_n) return(empty("insufficient_complete_samples"))
  if (stats::sd(x) == 0 || stats::sd(y) == 0) return(empty("zero_variance"))
  if (is.null(z)) {
    rho <- stats::cor(x, y, method = "spearman")
    df <- n - 2L
  } else {
    # Sensitivity only: partial Spearman, ranks residualised against ranked WES purity.
    z <- z[ok]
    if (stats::sd(z) == 0) return(empty("zero_variance_covariate"))
    xr <- rank(s25_winsor(x, probs), ties.method = "average")
    yr <- rank(s25_winsor(y, probs), ties.method = "average")
    zr <- rank(s25_winsor(z, probs), ties.method = "average")
    fit <- qr(cbind(1, zr))
    rx <- qr.resid(fit, xr); ry <- qr.resid(fit, yr)
    if (stats::sd(rx) < 1e-12 || stats::sd(ry) < 1e-12) return(empty("zero_residual_variance"))
    rho <- stats::cor(rx, ry); df <- n - 3L
  }
  rho <- max(-1, min(1, rho))
  p <- if (abs(rho) >= 1 - 1e-15) 0 else 2 * stats::pt(-abs(rho * sqrt(df / (1 - rho^2))), df)
  list(n = n, rho = rho, p_raw = p, status = "estimable")
}
s25_ci <- function(x, y, cfg) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < cfg$min_n || stats::sd(x) == 0 || stats::sd(y) == 0) return(c(NA_real_, NA_real_, 0L))
  values <- replicate(cfg$bootstrap_B, {
    ii <- sample.int(length(x), replace = TRUE)
    if (stats::sd(x[ii]) == 0 || stats::sd(y[ii]) == 0) NA_real_
    else stats::cor(x[ii], y[ii], method = "spearman")
  })
  valid <- is.finite(values)
  if (sum(valid) < ceiling(.95 * cfg$bootstrap_B)) return(c(NA_real_, NA_real_, sum(valid)))
  c(stats::quantile(values[valid], c(.025, .975), names = FALSE), sum(valid))
}
s25_expression <- function(path) {
  dt <- s25_read(path)
  if (ncol(dt) < 3L || nrow(dt) < 2L) stop("Invalid expression matrix: ", path)
  ids <- as.character(dt[[1L]]); sample_ids <- names(dt)[-1L]
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("Missing or duplicate feature IDs: ", path)
  if (anyNA(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) stop("Invalid sample IDs: ", path)
  if (!all(vapply(dt[-1L], is.numeric, logical(1)))) stop("Non-numeric expression column: ", path)
  mat <- as.matrix(dt[-1L]); rownames(mat) <- ids
  if (any(is.infinite(mat))) stop("Infinite expression value: ", path)
  if (any(grepl("TCGA[-.]", sample_ids, ignore.case = TRUE))) stop("External expression includes a TCGA-like sample ID.")
  mat
}
s25_map_genes <- function(mat, genes) {
  raw <- rownames(mat); stable <- sub("\\.[0-9]+$", "", raw)
  wanted <- sort(unique(genes))
  ann <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
    keys = wanted, keytype = "SYMBOL", columns = "ENSEMBL")
  ann <- unique(ann[!is.na(ann$ENSEMBL), c("SYMBOL", "ENSEMBL")])
  out <- matrix(NA_real_, nrow = length(wanted), ncol = ncol(mat), dimnames = list(wanted, colnames(mat)))
  audit <- vector("list", length(wanted))
  for (i in seq_along(wanted)) {
    gene <- wanted[i]
    direct <- which(raw == gene)
    ix <- if (length(direct)) direct else which(stable %in% ann$ENSEMBL[ann$SYMBOL == gene])
    status <- if (!length(ix)) "absent" else if (length(ix) != 1L) "ambiguous_multiple_expression_rows" else "mapped"
    if (status == "mapped") {
      # A selected Ensembl ID mapping to multiple symbols is not silently accepted.
      owners <- unique(ann$SYMBOL[ann$ENSEMBL == stable[ix]])
      if (!length(direct) && length(owners) > 1L) status <- "ambiguous_shared_ensembl"
      else out[i, ] <- mat[ix, ]
    }
    audit[[i]] <- data.frame(gene = gene, status = status, n_expression_rows = length(ix),
      source_ids = paste(raw[ix], collapse = ";"), stringsAsFactors = FALSE)
  }
  list(matrix = out, audit = do.call(rbind, audit), annotation = ann)
}
s25_primary <- function(path, cfg) {
  x <- s25_read(path)
  s25_need(x, c("mirna_id", "gene", "qualifies", "strict_recurrence", "purity_setting", "threshold", "run_id"), "Primary table")
  x <- x[s25_bool(x$qualifies), , drop = FALSE]
  if (nrow(x) != cfg$expected_primary_pairs) stop("Qualified primary count differs from the configured frozen count.")
  s25_unique(x, c("mirna_id", "gene"), "Qualified primary set")
  if (anyNA(x$mirna_id) || anyNA(x$gene) || any(x$strict_recurrence < 3L) || anyNA(x$strict_recurrence)) stop("Invalid primary keys or recurrence.")
  if (!all(x$purity_setting == "CPE") || !all(x$threshold == cfg$rho_threshold) || !all(x$run_id == cfg$discovery_run_id)) stop("Primary table provenance/threshold mismatch.")
  x[order(x$mirna_id, x$gene), , drop = FALSE]
}
s25_reference <- function(path, primary) {
  ref <- s25_read(path)
  s25_need(ref, c("mirna_id", "gene", "cancer", "arm", "adjustment", "rho", "eligible", "primary_support"), "Reference effects")
  ref <- ref[ref$cancer == "BRCA" & ref$adjustment %in% c("unadjusted", "CPE"), , drop = FALSE]
  ref <- ref[s25_key(ref, c("mirna_id", "gene")) %in% s25_key(primary, c("mirna_id", "gene")), , drop = FALSE]
  ref$arm <- c("miRNA-CYT" = "mirna_cyt", "miRNA-gene" = "mirna_gene", "gene-CYT" = "gene_cyt")[ref$arm]
  if (anyNA(ref$arm)) stop("Unknown reference arm label.")
  s25_unique(ref, c("mirna_id", "gene", "arm", "adjustment"), "BRCA reference effects")
  ref$rho[!s25_bool(ref$eligible)] <- NA_real_
  ref
}
s25_fetch <- function(name, url, cfg, dest) {
  old_timeout <- getOption("timeout")
  options(timeout = max(600, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  path <- file.path(dest, name)
  if (!is.null(cfg$local_external_dir)) {
    src <- file.path(cfg$local_external_dir, name)
    if (!file.exists(src)) stop("Missing offline file: ", src)
    if (!file.copy(src, path, overwrite = FALSE)) stop("Cannot copy offline file: ", name)
    method <- "copied_from_local_external_dir"
  } else {
    if (!isTRUE(cfg$allow_download)) stop("Downloads disabled and no local_external_dir configured.")
    message("Downloading ", name)
    tmp <- paste0(path, ".part")
    on.exit(unlink(tmp), add = TRUE)
    status <- utils::download.file(url, tmp, mode = "wb", quiet = FALSE, method = "libcurl")
    if (status != 0L || !file.exists(tmp) || file.info(tmp)$size < 1000) stop("Download failed: ", name)
    first <- readLines(tmp, n = 1L, warn = FALSE)
    if (!length(first) || grepl("<html|<!DOCTYPE", first, ignore.case = TRUE)) stop("Downloaded an HTML response: ", name)
    if (!file.rename(tmp, path)) stop("Cannot commit download: ", name)
    method <- "downloaded_https"
  }
  data.frame(file = name, url = url, acquisition = method, timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    bytes = file.info(path)$size, sha256 = s25_sha(path), stringsAsFactors = FALSE)
}

run_stage25 <- function(cfg, package_dir) {
  required <- c("data.table", "digest", "AnnotationDbi", "org.Hs.eg.db")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
  if (cfg$bootstrap_B < 1000L || cfg$min_n < 30L || cfg$alpha != .05 || cfg$rho_threshold != .30) stop("Unexpected validation settings; edit the analysis plan deliberately before execution.")
  project <- normalizePath(cfg$project_root, winslash = "/", mustWork = TRUE)
  discovery <- normalizePath(file.path(project, "results"), winslash = "/", mustWork = TRUE)
  parent <- normalizePath(cfg$output_parent, winslash = "/", mustWork = FALSE)
  expected_out <- normalizePath(file.path(project,'results/external_validation'),winslash='/',mustWork=FALSE)
  if(!identical(parent,expected_out)) stop('External outputs must use results/external_validation')
  status_path <- file.path(discovery, "logs", "run_status.tsv")
  status <- s25_read(status_path); s25_need(status, c("stage", "status"), "Run status")
  rows <- grepl("^24_", status$stage)
  if (!any(rows) || tail(status$status[rows], 1L) != "COMPLETE") stop("Discovery Stage 24 is not recorded COMPLETE.")
  primary_path <- file.path(discovery, "tables", "primary", "triplets_rho030_CPE.tsv")
  reference_path <- file.path(discovery, "tables", "downstream", "meta_analysis_input_all_cancers.tsv")
  annotation_path <- file.path(discovery, "processed", "stage05", "sample_annotation_with_purity.tsv")
  source_paths <- c(primary_path, reference_path, annotation_path, status_path)
  if (any(!file.exists(source_paths))) stop("Required same-run source missing: ", paste(source_paths[!file.exists(source_paths)], collapse = ", "))
  original_hashes <- vapply(source_paths, s25_sha, character(1))
  primary <- s25_primary(primary_path, cfg); reference <- s25_reference(reference_path, primary)
  ann <- s25_read(annotation_path); s25_need(ann, "sample_id", "Discovery sample annotation")
  out <- parent
  if (dir.exists(out)) unlink(out, recursive=TRUE)
  dir.create(out, recursive = TRUE)
  for (sub in c("inputs", "tables", "figures", "reports", "logs", "code")) dir.create(file.path(out, sub))
  done <- FALSE
  on.exit({
    unchanged <- identical(vapply(source_paths, s25_sha, character(1)), original_hashes)
    writeLines(c(paste0("status=", if (done && unchanged) "COMPLETE" else "FAILED"), paste0("discovery_inputs_unchanged=", unchanged)), file.path(out, "logs", "status.txt"))
    if (!unchanged) warning("Discovery inputs changed during validation; do not use output.")
  }, add = TRUE)
  files_code <- list.files(package_dir, pattern = "\\.(R|md)$", full.names = TRUE)
  if (!all(file.copy(files_code, file.path(out, "code"), overwrite = FALSE))) stop("Cannot archive validation code.")
  writeLines(capture.output(dput(cfg)), file.path(out, "code", "resolved_config.R.txt"))
  s25_write(primary, file.path(out, "inputs", "frozen_primary_pairs.tsv"))
  s25_write(reference, file.path(out, "inputs", "frozen_BRCA_reference_effects.tsv"))
  s25_write(data.frame(path = source_paths, sha256 = unname(original_hashes)), file.path(out, "inputs", "discovery_source_hashes.tsv"))
  plan <- c("# Locked external validation analysis plan", paste("Created UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Discovery run:", cfg$discovery_run_id), paste("Fixed pairs:", nrow(primary)),
    "External source: CPTAC prospective BRCA, LinkedOmics pan-cancer freeze v1.2.",
    "No new candidate selection, no refitting of discovery, no threshold tuning after results.",
    "Primary: two-sided Spearman tests on log-scale data without additional log transform or winsorization.",
    "CYT=(GZMA+PRF1)/2 on the provided log2 expression scale; cross-platform absolute scores are not comparable.",
    "One unique test per miRNA-CYT, gene-CYT, and miRNA-gene hypothesis; BH separately per arm.",
    "The planned family sizes include missing/ambiguous/non-estimable hypotheses.",
    "Triplet association support: max of the three raw two-sided p-values; BH and BY across all fixed pairs.",
    "Effect threshold: absolute rho>=0.30 on all three arms. Direction agreement with frozen TCGA-BRCA effects is reported separately.",
    "WES-purity sensitivity: 0.1/99.9 percentile winsorization then partial Spearman via residualized ranks, df=n-3.",
    "WES purity is not CPE. No CPE replication claim or new pan-cancer recurrence claim is made.",
    "Paired-patient percentile bootstrap CIs for primary unique arms, B=1000 by default, seed fixed; no bootstrapped p-value.",
    "No survival or interaction-model validation is attempted in this module.")
  writeLines(plan, file.path(out, "reports", "ANALYSIS_PLAN.md"))
  base_url <- "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/BRCA/"
  names_external <- c(mrna = "BRCA_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt", mirna = "BRCA_miRNAseq_mature_miRNA_RPM_log2_Tumor.txt", meta = "BRCA_meta.txt", phenotype = "BRCA_phenotype.txt")
  manifests <- lapply(names_external, function(name) s25_fetch(name, paste0(base_url, name), cfg, file.path(out, "inputs")))
  s25_write(do.call(rbind, manifests), file.path(out, "inputs", "external_download_manifest.tsv"))
  mrna <- s25_expression(file.path(out, "inputs", names_external[["mrna"]]))
  mirna <- s25_expression(file.path(out, "inputs", names_external[["mirna"]]))
  meta <- s25_read(file.path(out, "inputs", names_external[["meta"]])); s25_need(meta, "case_id", "External metadata")
  meta <- meta[meta$case_id != "data_type" & !is.na(meta$case_id), , drop = FALSE]
  s25_unique(meta, "case_id", "External metadata")
  samples <- sort(Reduce(intersect, list(colnames(mrna), colnames(mirna), as.character(meta$case_id))))
  if (length(samples) < cfg$min_n) stop("Too few exactly matched external samples.")
  if (length(intersect(samples, as.character(ann$sample_id)))) stop("External/discovery sample-ID overlap detected.")
  # Exact IDs only. No substring/fuzzy matching or pairing by row order.
  mrna <- mrna[, samples, drop = FALSE]; mirna <- mirna[, samples, drop = FALSE]
  mapping <- s25_map_genes(mrna, c(primary$gene, "GZMA", "PRF1"))
  gm <- mapping$matrix
  s25_write(mapping$audit, file.path(out, "tables", "gene_mapping_audit.tsv"))
  s25_write(mapping$annotation, file.path(out, "tables", "ensembl_symbol_mapping.tsv"))
  if (!all(mapping$audit$status[match(c("GZMA", "PRF1"), mapping$audit$gene)] == "mapped")) stop("CYT genes absent or ambiguously mapped; see gene_mapping_audit.tsv.")
  cyt <- (gm["GZMA", ] + gm["PRF1", ]) / 2
  if (sum(is.finite(cyt)) < cfg$min_n) stop("Too few finite CYT scores.")
  mi_map <- data.frame(mirna_id = sort(unique(primary$mirna_id)), stringsAsFactors = FALSE)
  mi_map$status <- ifelse(mi_map$mirna_id %in% rownames(mirna), "exact_mature_id_match", "absent")
  s25_write(mi_map, file.path(out, "tables", "mirna_mapping_audit.tsv"))
  phenotype <- s25_read(file.path(out, "inputs", names_external[["phenotype"]]))
  s25_need(phenotype, c("idx", "WES_purity"), "External phenotype")
  s25_unique(phenotype, "idx", "External phenotype")
  purity <- suppressWarnings(as.numeric(phenotype$WES_purity[match(samples, as.character(phenotype$idx))]))
  invalid_purity <- is.finite(purity) & (purity < 0 | purity > 1)
  if (any(invalid_purity)) stop("WES purity outside [0,1]; no automatic percent conversion.")
  sample_table <- data.frame(sample_id = samples, CYT = unname(cyt), GZMA = gm["GZMA", ], PRF1 = gm["PRF1", ], WES_purity = purity, stringsAsFactors = FALSE)
  s25_write(sample_table, file.path(out, "tables", "external_sample_annotation.tsv"))
  tests <- rbind(
    data.frame(arm = "mirna_cyt", entity = sort(unique(primary$mirna_id)), stringsAsFactors = FALSE),
    data.frame(arm = "gene_cyt", entity = sort(unique(primary$gene)), stringsAsFactors = FALSE),
    data.frame(arm = "mirna_gene", entity = s25_key(primary, c("mirna_id", "gene")), stringsAsFactors = FALSE))
  # Internal compound keys use a control delimiter, replaced in exported tables below.
  get_vector <- function(id, mat) if (id %in% rownames(mat)) as.numeric(mat[id, ]) else rep(NA_real_, length(samples))
  arms_all <- list(); set.seed(cfg$seed)
  for (adjustment in c("unadjusted", "WES_purity")) {
    rows_out <- vector("list", nrow(tests))
    for (i in seq_len(nrow(tests))) {
      arm <- tests$arm[i]; entity <- tests$entity[i]
      mi <- gene <- NA_character_
      if (arm == "mirna_cyt") { mi <- entity; x <- get_vector(mi, mirna); y <- cyt }
      else if (arm == "gene_cyt") { gene <- entity; x <- get_vector(gene, gm); y <- cyt }
      else {
        pair <- strsplit(entity, "\034", fixed = TRUE)[[1L]]; mi <- pair[1L]; gene <- pair[2L]
        x <- get_vector(mi, mirna); y <- get_vector(gene, gm)
      }
      feature_available <- (is.na(mi) || mi %in% rownames(mirna)) && (is.na(gene) || mapping$audit$status[match(gene, mapping$audit$gene)] == "mapped")
      fit <- s25_cor(x, y, cfg$min_n, if (adjustment == "WES_purity") purity else NULL, cfg$winsor_probs)
      if (!isTRUE(feature_available)) fit$status <- "missing_or_ambiguous_feature"
      ci <- if (adjustment == "unadjusted" && fit$status == "estimable") s25_ci(x, y, cfg) else c(NA_real_, NA_real_, 0L)
      rows_out[[i]] <- data.frame(arm = arm, mirna_id = mi, gene = gene, adjustment = adjustment,
        n = fit$n, rho = fit$rho, p_raw = fit$p_raw, status = fit$status,
        ci_low = ci[1L], ci_high = ci[2L], bootstrap_valid = as.integer(ci[3L]), stringsAsFactors = FALSE)
    }
    arms <- do.call(rbind, rows_out); arms$fdr_BH <- arms$fdr_BY <- NA_real_
    for (arm in unique(arms$arm)) {
      ix <- which(arms$arm == arm)
      arms$fdr_BH[ix] <- s25_family_adjust(arms$p_raw[ix], length(ix))
      arms$fdr_BY[ix] <- s25_family_adjust(arms$p_raw[ix], length(ix), "BY")
    }
    arms_all[[adjustment]] <- arms
  }
  unique_arms <- do.call(rbind, arms_all); rownames(unique_arms) <- NULL
  s25_write(unique_arms, file.path(out, "tables", "external_unique_arm_tests.tsv"))
  long <- list(); result <- list()
  for (adjustment in names(arms_all)) {
    long_adj <- list()
    for (i in seq_len(nrow(primary))) for (arm in c("mirna_cyt", "mirna_gene", "gene_cyt")) {
      mi <- primary$mirna_id[i]; gene <- primary$gene[i]; arms <- arms_all[[adjustment]]
      take <- arms$arm == arm
      if (arm != "gene_cyt") take <- take & !is.na(arms$mirna_id) & arms$mirna_id == mi
      if (arm != "mirna_cyt") take <- take & !is.na(arms$gene) & arms$gene == gene
      row <- arms[which(take), , drop = FALSE]
      if (nrow(row) != 1L) stop("Non-unique external arm join.")
      row$mirna_id <- mi; row$gene <- gene
      for (ref_adj in c("CPE", "unadjusted")) {
        rr <- reference[reference$mirna_id == mi & reference$gene == gene & reference$arm == arm & reference$adjustment == ref_adj, , drop = FALSE]
        ref_rho <- if (nrow(rr)) rr$rho else NA_real_
        row[[paste0("TCGA_BRCA_", ref_adj, "_rho")]] <- ref_rho
        row[[paste0("direction_agrees_", ref_adj)]] <- if (!is.finite(ref_rho) || ref_rho == 0 || !is.finite(row$rho) || row$rho == 0) NA else sign(row$rho) == sign(ref_rho)
      }
      long_adj[[length(long_adj) + 1L]] <- row
    }
    long_adj <- do.call(rbind, long_adj)
    pair_rows <- lapply(seq_len(nrow(primary)), function(i) {
      q <- long_adj[long_adj$mirna_id == primary$mirna_id[i] & long_adj$gene == primary$gene[i], , drop = FALSE]
      good <- all(q$status == "estimable") && all(is.finite(q$p_raw))
      all_known_true <- function(v) if (anyNA(v)) NA else all(v)
      data.frame(mirna_id = primary$mirna_id[i], gene = primary$gene[i], adjustment = adjustment,
        discovery_recurrence = primary$strict_recurrence[i], n_min = min(q$n), evaluable = good,
        p_conjunction = if (good) max(q$p_raw) else NA_real_,
        minimum_abs_rho = if (good) min(abs(q$rho)) else NA_real_,
        all_arms_BH_supported = if (good) all(q$fdr_BH < cfg$alpha) else NA,
        all_directions_agree_CPE = all_known_true(q$direction_agrees_CPE),
        all_directions_agree_unadjusted = all_known_true(q$direction_agrees_unadjusted), stringsAsFactors = FALSE)
    })
    pair_rows <- do.call(rbind, pair_rows)
    pair_rows$conjunction_fdr_BH <- s25_family_adjust(pair_rows$p_conjunction, nrow(primary))
    pair_rows$conjunction_fdr_BY <- s25_family_adjust(pair_rows$p_conjunction, nrow(primary), "BY")
    pair_rows$association_support_BH <- ifelse(pair_rows$evaluable,
      pair_rows$conjunction_fdr_BH < cfg$alpha & pair_rows$minimum_abs_rho >= cfg$rho_threshold, NA)
    pair_rows$directional_replication_vs_CPE_BRCA <- pair_rows$association_support_BH & pair_rows$all_directions_agree_CPE
    pair_rows$directional_replication_vs_unadjusted_BRCA <- pair_rows$association_support_BH & pair_rows$all_directions_agree_unadjusted
    pair_rows$strict_three_arm_BH_effect_direction <- pair_rows$all_arms_BH_supported &
      pair_rows$minimum_abs_rho >= cfg$rho_threshold & pair_rows$all_directions_agree_CPE
    result[[adjustment]] <- pair_rows; long[[adjustment]] <- long_adj
  }
  pair_results <- do.call(rbind, result); long_results <- do.call(rbind, long)
  rownames(pair_results) <- rownames(long_results) <- NULL
  s25_write(pair_results, file.path(out, "tables", "external_pair_validation.tsv"))
  s25_write(long_results, file.path(out, "tables", "external_pair_arm_comparison.tsv"))
  coverage <- data.frame(metric = c("fixed_pairs", "fixed_miRNAs", "fixed_genes", "matched_external_samples", "finite_CYT", "valid_WES_purity", "mapped_selected_genes", "matched_mature_miRNAs", "exact_sample_id_overlap_with_discovery"),
    value = c(nrow(primary), length(unique(primary$mirna_id)), length(unique(primary$gene)), length(samples), sum(is.finite(cyt)), sum(is.finite(purity)),
      sum(mapping$audit$gene %in% primary$gene & mapping$audit$status == "mapped"), sum(mi_map$status == "exact_mature_id_match"), 0L))
  s25_write(coverage, file.path(out, "tables", "coverage_summary.tsv"))
  families <- aggregate(p_raw ~ adjustment + arm, unique_arms, function(x) sum(is.finite(x)), na.action = na.pass)
  names(families)[3L] <- "n_evaluable"
  families$planned_n <- vapply(families$arm, function(a) sum(tests$arm == a), integer(1))
  s25_write(families, file.path(out, "tables", "testing_families.tsv"))
  checks <- c(fixed_pair_count = nrow(primary) == cfg$expected_primary_pairs,
    all_pairs_retained_each_adjustment = all(vapply(result, nrow, integer(1)) == nrow(primary)),
    unique_pair_keys = !anyDuplicated(s25_key(pair_results, c("mirna_id", "gene", "adjustment"))),
    three_arms_per_pair = nrow(long_results) == 3L * nrow(pair_results),
    unique_tests_not_repeated_in_BH = nrow(unique_arms) == 2L * nrow(tests),
    valid_p_values = all(unique_arms$p_raw[is.finite(unique_arms$p_raw)] >= 0 & unique_arms$p_raw[is.finite(unique_arms$p_raw)] <= 1),
    original_inputs_unchanged = identical(vapply(source_paths, s25_sha, character(1)), original_hashes))
  s25_write(data.frame(check = names(checks), pass = unname(checks)), file.path(out, "logs", "validation_checks.tsv"))
  if (!all(checks)) stop("Output consistency check failed.")
  figure_status <- "not_requested"
  if (isTRUE(cfg$write_figures)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) figure_status <- "SKIPPED_missing_ggplot2"
    else {
      primary_plot <- long_results[long_results$adjustment == "unadjusted", , drop = FALSE]
      primary_plot$pair <- paste(primary_plot$mirna_id, primary_plot$gene, sep = " / ")
      primary_plot$pair <- factor(primary_plot$pair, levels = rev(paste(primary$mirna_id, primary$gene, sep = " / ")))
      primary_plot$display_rho <- ifelse(is.finite(primary_plot$rho), sprintf("%.2f", primary_plot$rho), "NA")
      fig <- ggplot2::ggplot(primary_plot, ggplot2::aes(x = arm, y = pair, fill = rho)) +
        ggplot2::geom_tile(colour = "white") + ggplot2::geom_text(ggplot2::aes(label = display_rho), size = 2.4) +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1), na.value = "grey85", name = "Spearman rho") +
        ggplot2::labs(x = NULL, y = NULL) + ggplot2::theme_minimal(base_size = 9) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom")
      ggplot2::ggsave(file.path(out, "figures", "external_BRCA_associations.pdf"), fig, width = 7.5, height = 12)
      ggplot2::ggsave(file.path(out, "figures", "external_BRCA_associations.png"), fig, width = 7.5, height = 12, dpi = 400)
      figure_status <- "written"
    }
  }
  summary <- do.call(rbind, lapply(names(result), function(a) {
    r <- result[[a]]
    data.frame(adjustment = a, fixed_pairs = nrow(r), evaluable = sum(r$evaluable),
      missing_or_nonestimable = sum(!r$evaluable), conjunction_BH_and_effect = sum(r$association_support_BH %in% TRUE),
      directional_vs_CPE_BRCA = sum(r$directional_replication_vs_CPE_BRCA %in% TRUE),
      directional_vs_unadjusted_BRCA = sum(r$directional_replication_vs_unadjusted_BRCA %in% TRUE),
      strict_arm_BH_effect_CPE_direction = sum(r$strict_three_arm_BH_effect_direction %in% TRUE))
  }))
  s25_write(summary, file.path(out, "tables", "validation_summary.tsv"))
  report <- c("# External BRCA validation", paste("Discovery run:", cfg$discovery_run_id),
    paste("Exactly matched external patients:", length(samples)), paste("Figure status:", figure_status),
    "", capture.output(print(summary, row.names = FALSE)), "",
    "Interpretation: these are independent-cohort BRCA association tests of fixed discovery pairs.",
    "The primary external analysis is unadjusted. WES-purity sensitivity is distinct from CPE adjustment.",
    "Direction agreement is assessed against each pair's frozen TCGA-BRCA estimates, not assumed positive.",
    "No significant result is guaranteed. Missing/ambiguous features remain in the denominator and are labelled non-evaluable.",
    "Pair-level BH conjunction and unique-arm BH criteria are different summaries; neither is substituted for the other.",
    "All bootstrap CIs are pointwise descriptive intervals, not multiplicity-adjusted intervals.",
    "Independence evidence: prospectively collected CPTAC BRCA cohort (Krug et al., Cell 2020), non-TCGA sample IDs and no exact overlap with discovery identifiers.",
    "Identifier comparison alone cannot prove donor independence when aliases differ; retain the cohort publication and accession lineage in the manuscript.",
    "This does not validate recurrence in >=3 external cancers, survival prediction, causal targeting or population-specific interaction effects.",
    "Reference: Krug K et al. Cell. 2020;183:1436-1456.e31. doi:10.1016/j.cell.2020.10.036.",
    "Data portal: https://www.linkedomics.org/data_download/CPTAC-pancan-BRCA/", "",
    "Figure legend: External CPTAC-BRCA Spearman correlations for all fixed discovery pairs. Cell labels show rho; grey cells labelled NA are non-estimable. Effect sizes alone do not indicate FDR support; consult external_pair_validation.tsv.")
  writeLines(report, file.path(out, "reports", "EXTERNAL_VALIDATION_REPORT.md"))
  writeLines(capture.output(sessionInfo()), file.path(out, "logs", "sessionInfo.txt"))
  packages <- unique(c(required, "ggplot2"))
  s25_write(data.frame(package = packages, version = vapply(packages, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_, character(1))), file.path(out, "logs", "package_versions.tsv"))
  all_files <- list.files(out, recursive = TRUE, full.names = TRUE)
  s25_write(data.frame(path = substring(all_files, nchar(out) + 2L), sha256 = vapply(all_files, s25_sha, character(1))), file.path(out, "logs", "output_sha256.tsv"))
  done <- TRUE
  message("Stage 25 COMPLETE: ", out)
  out
}
