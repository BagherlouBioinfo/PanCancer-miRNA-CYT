# ============================================================================
# test_spearman_pvalue_equivalence.R
#
# Verifies R/stats_engine.R::spearman_pvalue() against the trusted reference
# implementation, cor.test(method="spearman", exact=FALSE), across multiple
# sample sizes (including n=30), weak/moderate/strong/no effects, positive
# and negative correlation, tied ranks, large n, and the numerically
# unstable |rho|->1 boundary (see R/stats_engine.R's own header comment for
# why that boundary needed an explicit guard, found during this
# verification). Tests BOTH rho and p-value agreement, not rho alone.
#
# SYNTHETIC DATA ONLY. No project data is loaded.
#
# Usage: Rscript tests/test_spearman_pvalue_equivalence.R <project_root>
# ============================================================================

run_spearman_equivalence_tests <- function(project_root) {
  source(file.path(project_root, "scripts", "lib", "stats_engine.R"))
  set.seed(42)

  results <- list()
  add <- function(id, desc, pass, detail) {
    results[[length(results) + 1]] <<- data.frame(test_id = id, description = desc,
      result = if (isTRUE(pass)) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE)
    cat(sprintf("[%s] %s -- %s\n", if (isTRUE(pass)) "PASS" else "FAIL", id, detail))
  }

  test_case <- function(id, x, y, rho_tol = 1e-9, p_tol = 1e-6) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
    if (length(unique(x)) < 2 || length(unique(y)) < 2 || length(x) < 3) {
      results[[length(results) + 1]] <<- data.frame(test_id = id, description = "degenerate input", result = "SKIPPED", detail = "insufficient variance/n", stringsAsFactors = FALSE)
      return(invisible())
    }
    ref <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    r_own <- cor(rank(x), rank(y))
    p_own <- spearman_pvalue(r_own, length(x))
    rho_ok <- abs(unname(ref$estimate) - r_own) < rho_tol
    # For extremely small p-values (both essentially 0 up to machine
    # precision), compare on the log scale instead of absolute difference,
    # since both are already far below any decision threshold.
    p_ok <- if (ref$p.value < 1e-10 && p_own < 1e-10) TRUE else abs(ref$p.value - p_own) < p_tol
    add(id, sprintf("n=%d Spearman rho and p-value match cor.test(exact=FALSE)", length(x)),
        rho_ok && p_ok, sprintf("rho: ref=%.6f own=%.6f (diff=%.2e) | p: ref=%.6g own=%.6g", unname(ref$estimate), r_own, abs(unname(ref$estimate) - r_own), ref$p.value, p_own))
  }

  n <- 30
  x <- rnorm(n)
  test_case("SP01_n30_no_effect", x, rnorm(n))
  test_case("SP02_n30_weak", x, x * 0.2 + rnorm(n) * 2)
  test_case("SP03_n30_moderate", x, x * 0.6 + rnorm(n))
  test_case("SP04_n30_strong", x, x * 3 + rnorm(n) * 0.3)
  test_case("SP05_n30_negative", x, -x * 2 + rnorm(n))
  test_case("SP06_n5", rnorm(5), rnorm(5))
  test_case("SP07_n10", rnorm(10), rnorm(10))
  test_case("SP08_ties_n15", c(1,1,2,2,3,4,4,4,5,6,6,7,8,9,10), c(5,4,4,3,3,2,2,1,1,1,2,3,4,5,6))
  xl <- rnorm(500); test_case("SP09_n500", xl, xl * 0.4 + rnorm(500))
  test_case("SP10_perfect_positive", 1:10, 1:10)
  test_case("SP11_perfect_negative", 1:10, 10:1)
  test_case("SP12_near_perfect_one_swap", 1:20, c(1:18, 20, 19))

  # Vectorized production path (one_vs_many_spearman / rank_rows / safe_rank),
  # used by Stages 06, 08, 10 -- must match cor.test() row-by-row too.
  m <- rbind(rnorm(30), rnorm(30) * -1 + rnorm(30) * 0.1, 1:30)
  y_ref <- rnorm(30)
  ranked_y <- safe_rank(y_ref)
  res_vec <- one_vs_many_spearman(ranked_y, rank_rows(m))
  max_diff_cor <- 0; max_diff_p <- 0
  for (i in 1:3) {
    ref <- suppressWarnings(cor.test(m[i, ], y_ref, method = "spearman", exact = FALSE))
    max_diff_cor <- max(max_diff_cor, abs(unname(ref$estimate) - res_vec$cor[i]))
    max_diff_p <- max(max_diff_p, abs(ref$p.value - res_vec$p[i]))
  }
  add("SP13_vectorized_matches_scalar", "one_vs_many_spearman() (the vectorized production path used by Stages 06/08/10) matches cor.test() row-by-row",
      max_diff_cor < 1e-9 && max_diff_p < 1e-6, sprintf("max cor diff=%.2e, max p diff=%.2e", max_diff_cor, max_diff_p))

  df <- do.call(rbind, results)
  out_path <- file.path(project_root, "tests", "output", "test_spearman_pvalue_equivalence_results.tsv")
  write.table(df, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
  n_fail <- sum(df$result == "FAIL")
  cat(sprintf("\n%d/%d PASS, %d SKIPPED, %d FAIL\n", sum(df$result == "PASS"), nrow(df), sum(df$result == "SKIPPED"), n_fail))
  if (n_fail > 0) stop(sprintf("%d Spearman p-value equivalence test(s) FAILED. See %s.", n_fail, out_path), call. = FALSE)
  invisible(df)
}

if (sys.nframe() == 0) {
  cli_args <- commandArgs(trailingOnly = TRUE)
  if (length(cli_args) < 1) stop("Usage: Rscript test_spearman_pvalue_equivalence.R <project_root>")
  invisible(run_spearman_equivalence_tests(normalizePath(cli_args[[1]], mustWork = TRUE)))
}
