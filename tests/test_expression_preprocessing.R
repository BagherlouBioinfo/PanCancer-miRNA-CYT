# ============================================================================
# test_expression_preprocessing.R
#
# Synthetic-data regression tests for R/expression_preprocessing.R, covering
# every scenario required by the Section H (expression preprocessing) audit:
# already-log data, linear data, extreme outliers, NA/Inf/non-numeric
# values, duplicate standardized barcodes (technical duplicate and
# conflicting record), aliquot duplicates, duplicated gene symbols,
# zero-variance features, shuffled sample order, and distinct miRNA/mRNA
# scale policies. Also covers two specific R data.frame subsetting defects
# found and fixed during this audit (setdiff()+match() silently skipping
# later occurrences of a duplicated raw column name during barcode
# standardization; `[.data.frame` silently applying make.unique()-style
# renaming on any multi-column selection when the result contains a
# repeated name) -- see reports/EXPRESSION_PREPROCESSING_AUDIT_AND_DECISION.md
# section 5 for the full account.
#
# SYNTHETIC DATA ONLY. No project data is loaded.
#
# Usage: Rscript tests/test_expression_preprocessing.R <project_root>
# ============================================================================

run_expression_preprocessing_tests <- function(project_root) {
  source(file.path(project_root, "scripts", "lib", "tcga_identifiers.R"))
  source(file.path(project_root, "scripts", "lib", "expression_preprocessing.R"))

  results <- list()
  add <- function(id, desc, pass, detail) {
    results[[length(results) + 1]] <<- data.frame(test_id = id, description = desc,
      result = if (isTRUE(pass)) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE)
    cat(sprintf("[%s] %s -- %s\n", if (isTRUE(pass)) "PASS" else "FAIL", id, detail))
  }

  # EP01: already-log2 data is NOT re-transformed
  df1 <- data.frame(gene = c("A", "B", "C"),
    "TCGA.02.0001.01A" = c(5.1, 3.2, 8.8), "TCGA.02.0002.01A" = c(4.9, 3.5, 8.1), "TCGA.02.0003.01A" = c(5.5, 3.1, 8.9),
    check.names = FALSE)
  res1 <- preprocess_expression_matrix(df1, "gene_id", declared_scale = "log2", pseudo_count = 1)
  add("EP01_already_log2", "Already-log2 data is not re-transformed and scale check is consistent",
      res1$ok && !res1$transformed_by_pipeline && isTRUE(res1$scale_check$consistent),
      sprintf("ok=%s transformed=%s consistent=%s", res1$ok, res1$transformed_by_pipeline, res1$scale_check$consistent))

  # EP02: linear data IS log2-transformed
  df2 <- data.frame(gene = c("A", "B", "C"),
    "TCGA.02.0001.01A" = c(120, 5000, 30), "TCGA.02.0002.01A" = c(200, 100, 60), "TCGA.02.0003.01A" = c(80, 300, 45),
    check.names = FALSE)
  res2 <- preprocess_expression_matrix(df2, "gene_id", declared_scale = "linear", pseudo_count = 1)
  add("EP02_linear_transformed", "Linear data is log2(x+pseudo_count)-transformed and matches by hand",
      res2$ok && res2$transformed_by_pipeline && isTRUE(res2$scale_check$consistent) && abs(res2$data[1, 2] - log2(121)) < 1e-9,
      sprintf("ok=%s transformed=%s value_match=%s", res2$ok, res2$transformed_by_pipeline, abs(res2$data[1, 2] - log2(121)) < 1e-9))

  # EP03: declared scale contradicted by observed range (linear declared, but bounded [0,1])
  df3 <- data.frame(gene = c("A", "B"), "TCGA.02.0001.01A" = c(0.1, 0.5), "TCGA.02.0002.01A" = c(0.2, 0.6), check.names = FALSE)
  res3 <- preprocess_expression_matrix(df3, "gene_id", declared_scale = "linear", pseudo_count = 1)
  add("EP03_scale_misdeclared_linear", "A [0,1]-bounded matrix declared 'linear' is flagged inconsistent, not silently trusted",
      res3$ok && !isTRUE(res3$scale_check$consistent), sprintf("consistent=%s", res3$scale_check$consistent))

  # EP04: extreme outlier beyond the log2 sanity bound is flagged, not silently accepted
  df4 <- data.frame(gene = c("A", "B"), "TCGA.02.0001.01A" = c(5, 999999), "TCGA.02.0002.01A" = c(4, 3), check.names = FALSE)
  res4 <- preprocess_expression_matrix(df4, "gene_id", declared_scale = "log2", pseudo_count = 1)
  add("EP04_extreme_outlier_flagged", "A value far exceeding the log2 sanity bound is flagged inconsistent",
      res4$ok && !isTRUE(res4$scale_check$consistent), sprintf("consistent=%s max=%s", res4$scale_check$consistent, res4$scale_check$empirical_max))

  # EP05: duplicate standardized barcode, IDENTICAL values -> technical_duplicate, auto-resolved
  df5 <- data.frame(gene = c("A", "B"), "TCGA.02.0001.01A" = c(5, 3), "TCGA.02.0001.01A" = c(5, 3), "TCGA.02.0002.01A" = c(6, 4), check.names = FALSE)
  res5 <- preprocess_expression_matrix(df5, "gene_id", declared_scale = "log2", pseudo_count = 1)
  add("EP05_technical_duplicate_resolved", "Identical-valued duplicate barcode columns are classified technical_duplicate and resolved to one column",
      res5$ok && nrow(res5$duplicate_column_report) == 1 && res5$duplicate_column_report$classification[1] == "technical_duplicate" && ncol(res5$data) == 3,
      sprintf("ok=%s n_dup_groups=%d classification=%s n_out_cols=%d", res5$ok, nrow(res5$duplicate_column_report),
              if (nrow(res5$duplicate_column_report)) res5$duplicate_column_report$classification[1] else NA, ncol(res5$data)))

  # EP06: duplicate standardized barcode, CONFLICTING values -> stops, never auto-resolved
  df6 <- data.frame(gene = c("A", "B"), "TCGA.02.0001.01A" = c(5, 3), "TCGA.02.0001.01A" = c(9, 1), "TCGA.02.0002.01A" = c(6, 4), check.names = FALSE)
  res6 <- preprocess_expression_matrix(df6, "gene_id", declared_scale = "log2", pseudo_count = 1)
  add("EP06_conflicting_duplicate_stops", "Conflicting-valued duplicate barcode columns stop the stage (ok=FALSE), never silently resolved",
      !res6$ok && identical(res6$stage, "duplicate_column_conflict"), sprintf("ok=%s stage=%s", res6$ok, res6$stage))

  # EP07: aliquot duplicates (distinct full barcodes, same 15-char prefix) -- resolved to the more-complete aliquot
  df7 <- data.frame(gene = c("A", "B"),
    "TCGA.02.0001.01A.11.123.07" = c(5, NA), "TCGA.02.0001.01B.22.456.07" = c(4.9, 3.0), "TCGA.02.0002.01A.11.789.07" = c(6, 4),
    check.names = FALSE)
  res7 <- preprocess_expression_matrix(df7, "gene_id", declared_scale = "log2", pseudo_count = 1)
  add("EP07_aliquot_duplicate_resolved", "Two aliquots of the same 15-char sample are detected and resolved to the more-complete one (fewer missing values)",
      res7$ok && nrow(res7$aliquot_report) == 1 && res7$aliquot_report$n_aliquots[1] == 2 && ncol(res7$data) == 3 && res7$aliquot_dropped == "TCGA-02-0001-01A-11-123-07",
      sprintf("ok=%s n_aliquot_groups=%d dropped=%s", res7$ok, nrow(res7$aliquot_report), paste(res7$aliquot_dropped, collapse = ";")))

  # EP08: non-numeric value coerced to NA and explicitly counted, not silently absorbed
  df8 <- data.frame(gene = c("A", "B"), "TCGA.02.0001.01A" = c("5", "not_a_number"), "TCGA.02.0002.01A" = c("4", "3"), check.names = FALSE)
  res8 <- preprocess_expression_matrix(df8, "gene_id", declared_scale = "log2", pseudo_count = 1)
  add("EP08_non_numeric_reported", "A non-numeric, non-missing-token cell is coerced to NA and explicitly counted",
      res8$ok && res8$n_newly_na_from_coercion == 1, sprintf("n_newly_na_from_coercion=%d", res8$n_newly_na_from_coercion))

  # EP08b: recognized missing tokens (blank, "NA") do NOT count as newly-NA (expected missingness)
  df8b <- data.frame(gene = c("A", "B"), "TCGA.02.0001.01A" = c("5", "NA"), "TCGA.02.0002.01A" = c("4", ""), check.names = FALSE)
  res8b <- preprocess_expression_matrix(df8b, "gene_id", declared_scale = "log2", pseudo_count = 1)
  add("EP08b_expected_missingness_not_overcounted", "Recognized missing tokens (blank, 'NA') are not counted as newly-NA unparseable values",
      res8b$ok && res8b$n_newly_na_from_coercion == 0, sprintf("n_newly_na_from_coercion=%d", res8b$n_newly_na_from_coercion))

  # EP09: duplicated gene symbols resolved by highest-variance row
  m9 <- matrix(c(1, 1, 1, 5, 1, 9), nrow = 2, byrow = TRUE); rownames(m9) <- c("G1", "G1")
  d9 <- deduplicate_matrix_rows(m9)
  add("EP09_duplicate_gene_symbol_highest_variance", "Duplicated gene-symbol rows resolve to the highest-variance row",
      nrow(d9) == 1 && all(d9[1, ] == c(5, 1, 9)), sprintf("n_rows_out=%d", nrow(d9)))

  # EP10: zero-variance feature removed and counted
  df10 <- data.frame(gene = c("A", "B"), "TCGA.02.0001.01A" = c(5, 5), "TCGA.02.0002.01A" = c(5, 3), "TCGA.02.0003.01A" = c(5, 7), check.names = FALSE)
  res10 <- preprocess_expression_matrix(df10, "gene_id", declared_scale = "log2", pseudo_count = 1, remove_zero_variance = TRUE)
  add("EP10_zero_variance_removed", "A zero-variance feature is removed and the removal count is reported",
      res10$ok && nrow(res10$data) == 1 && res10$n_zero_variance_removed == 1,
      sprintf("n_features_out=%d n_zero_variance_removed=%d", nrow(res10$data), res10$n_zero_variance_removed))

  # EP11: winsorize_matrix_rows clips a per-row extreme value, leaves a non-extreme row materially unaffected
  m11 <- matrix(c(1, 2, 3, 4, 100, 1, 2, 3, 4, 5), nrow = 2, byrow = TRUE)
  w11 <- winsorize_matrix_rows(m11, probs = c(0.1, 0.9))
  add("EP11_winsorize_matrix_rows", "winsorize_matrix_rows() clips an extreme per-feature value and preserves dimensions",
      w11[1, 5] < 100 && all(dim(w11) == dim(m11)), sprintf("clipped_value=%.2f", w11[1, 5]))

  # EP12: shuffled sample column order produces the same result set (order-independence)
  df12a <- data.frame(gene = c("A", "B"), "TCGA.02.0001.01A" = c(5, 3), "TCGA.02.0002.01A" = c(6, 4), "TCGA.02.0003.01A" = c(7, 5), check.names = FALSE)
  df12b <- df12a[, c(1, 3, 2, 4)]
  res12a <- preprocess_expression_matrix(df12a, "gene_id", declared_scale = "log2", pseudo_count = 1)
  res12b <- preprocess_expression_matrix(df12b, "gene_id", declared_scale = "log2", pseudo_count = 1)
  order_independent <- res12a$ok && res12b$ok && setequal(colnames(res12a$data), colnames(res12b$data))
  if (order_independent) {
    o <- match(colnames(res12a$data), colnames(res12b$data))
    order_independent <- all(res12a$data == res12b$data[, o])
  }
  add("EP12_shuffled_sample_order", "Shuffling the input sample-column order produces an identical result set",
      order_independent, sprintf("order_independent=%s", order_independent))

  # EP13: miRNA and mRNA scale policies are declared and applied independently
  res_mirna <- preprocess_expression_matrix(df1, "mirna_id", declared_scale = "log2", pseudo_count = 1)
  res_mrna <- preprocess_expression_matrix(df2, "gene_id_raw", declared_scale = "linear", pseudo_count = 1)
  add("EP13_distinct_mirna_mrna_policy", "miRNA (log2, not transformed) and mRNA (linear, transformed) scale policies are independently declared and applied",
      res_mirna$declared_scale == "log2" && res_mrna$declared_scale == "linear" &&
        !res_mirna$transformed_by_pipeline && res_mrna$transformed_by_pipeline,
      sprintf("mirna: declared=%s transformed=%s | mrna: declared=%s transformed=%s",
              res_mirna$declared_scale, res_mirna$transformed_by_pipeline, res_mrna$declared_scale, res_mrna$transformed_by_pipeline))

  # EP14: standardization applies to EVERY sample-column position, not only the first
  # occurrence of each unique raw name (regression test for the setdiff()+match() defect
  # found and fixed this pass -- see file header).
  df14 <- data.frame(gene = c("A", "B"),
    "TCGA.02.0001.01A" = c(5, 3), "TCGA.02.0001.01A" = c(5, 3), "TCGA.02.0002.01A" = c(6, 4),
    check.names = FALSE)
  lst14 <- unclass(df14)
  standardized14 <- standardize_tcga_id(names(lst14)[-1])
  add("EP14_standardization_applies_to_every_position", "standardize_tcga_id() is applied to every sample-column position, not only the first occurrence of a duplicated raw name",
      length(unique(standardized14)) == 2 && standardized14[1] == standardized14[2],
      sprintf("standardized_names=%s", paste(standardized14, collapse = ",")))

  # EP15: clinical/survival duplicate-row classification (technical vs conflicting)
  clin_a <- data.frame(sample_id = c("S1", "S1", "S2"), OS = c(1, 1, 0), OS.time = c(100, 100, 200), stringsAsFactors = FALSE)
  dup_a <- classify_duplicate_rows(clin_a, "sample_id")
  clin_b <- data.frame(sample_id = c("S1", "S1", "S2"), OS = c(1, 0, 0), OS.time = c(100, 999, 200), stringsAsFactors = FALSE)
  dup_b <- classify_duplicate_rows(clin_b, "sample_id")
  add("EP15_clinical_duplicate_rows_classified", "Clinical duplicate rows are classified technical_duplicate (identical) vs conflicting_record (differing) correctly",
      nrow(dup_a) == 1 && dup_a$classification[1] == "technical_duplicate" &&
        nrow(dup_b) == 1 && dup_b$classification[1] == "conflicting_record",
      sprintf("identical_case=%s differing_case=%s", dup_a$classification[1], dup_b$classification[1]))

  # ==========================================================================
  # EP16-EP19: deduplicate_matrix_rows()/highest_variance_row()'s
  # all-NA-variance silent group-loss (Phase D finding, corrected). When
  # EVERY candidate row in a duplicate-rowname group has an uncomputable
  # (NA) variance, the entire group must NOT silently vanish.
  # ==========================================================================
  mat_all_na <- matrix(c(NA, NA, NA, NA, NA, NA), nrow = 2, byrow = TRUE)
  rownames(mat_all_na) <- c("GENE", "GENE"); colnames(mat_all_na) <- c("S1", "S2", "S3")
  out_all_na <- deduplicate_matrix_rows(mat_all_na)
  add("EP16_all_na_variance_duplicate_group_not_silently_dropped",
      "A duplicate-rowname group where EVERY candidate row is fully NA (variance uncomputable for all) does not silently vanish -- highest_variance_row() falls back to the first candidate rather than which.max() returning integer(0) and dropping the whole group",
      nrow(out_all_na) == 1 && "GENE" %in% rownames(out_all_na),
      sprintf("nrow=%d GENE present=%s", nrow(out_all_na), "GENE" %in% rownames(out_all_na)))

  mat_sparse_na <- matrix(c(5, NA, NA, 7, NA, NA), nrow = 2, byrow = TRUE)
  rownames(mat_sparse_na) <- c("GENE2", "GENE2"); colnames(mat_sparse_na) <- c("S1", "S2", "S3")
  out_sparse_na <- deduplicate_matrix_rows(mat_sparse_na)
  add("EP17_single_non_na_value_per_row_duplicate_group_not_silently_dropped",
      "A duplicate-rowname group where every candidate row has only 1 non-NA value (var() returns NA for n<2) does not silently vanish",
      nrow(out_sparse_na) == 1 && "GENE2" %in% rownames(out_sparse_na),
      sprintf("nrow=%d GENE2 present=%s", nrow(out_sparse_na), "GENE2" %in% rownames(out_sparse_na)))

  mat_normal_dup <- rbind(c(1, 1, 1), c(1, 2, 3), c(1, 1, 100))
  rownames(mat_normal_dup) <- c("G", "G", "G")
  out_normal_dup <- deduplicate_matrix_rows(mat_normal_dup)
  add("EP18_normal_highest_variance_tiebreak_unaffected",
      "The fix does not change behavior for the normal case: a 3-way duplicate with genuinely distinct, computable variances still keeps the highest-variance row",
      nrow(out_normal_dup) == 1 && isTRUE(all.equal(as.numeric(out_normal_dup["G", ]), c(1, 1, 100))),
      sprintf("kept row values: %s (expected 1,1,100 -- the highest-variance row)", paste(out_normal_dup["G", ], collapse = ",")))

  mat_tied_dup <- rbind(c(1, 2, 3), c(10, 11, 12))
  rownames(mat_tied_dup) <- c("G", "G")
  out_tied_dup <- deduplicate_matrix_rows(mat_tied_dup)
  add("EP19_tied_variance_deterministic_first_row_unaffected",
      "The fix does not change behavior for tied-variance duplicates: still deterministically keeps the first-occurring row",
      nrow(out_tied_dup) == 1 && isTRUE(all.equal(as.numeric(out_tied_dup["G", ]), c(1, 2, 3))),
      sprintf("kept row values: %s (expected 1,2,3 -- the first-occurring row)", paste(out_tied_dup["G", ], collapse = ",")))

  df <- do.call(rbind, results)
  out_path <- file.path(project_root, "tests", "output", "test_expression_preprocessing_results.tsv")
  write.table(df, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
  n_fail <- sum(df$result == "FAIL")
  cat(sprintf("\n%d/%d PASS, %d FAIL\n", sum(df$result == "PASS"), nrow(df), n_fail))
  if (n_fail > 0) stop(sprintf("%d expression-preprocessing test(s) FAILED. See %s.", n_fail, out_path), call. = FALSE)
  invisible(df)
}

if (sys.nframe() == 0) {
  cli_args <- commandArgs(trailingOnly = TRUE)
  if (length(cli_args) < 1) stop("Usage: Rscript test_expression_preprocessing.R <project_root>")
  invisible(run_expression_preprocessing_tests(normalizePath(cli_args[[1]], mustWork = TRUE)))
}
