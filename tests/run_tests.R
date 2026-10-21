# Synthetic/numerical checks only; does not require TCGA data.
run_tests <- function(project_root='.') {
  project_root <- normalizePath(project_root,winslash='/',mustWork=TRUE)
  dir.create(file.path(project_root,'tests/output'),recursive=TRUE,showWarnings=FALSE)
  all_r <- list.files(project_root,pattern='\\.R$',recursive=TRUE,full.names=TRUE)
  all_r <- all_r[!grepl('/results/',all_r,fixed=TRUE)]
  for(f in all_r) parse(f)
  cat('[PASS] All R files parse\n')
  tests <- c(test_expression_preprocessing.R='run_expression_preprocessing_tests',
    test_pairwise_complete_spearman.R='run_pairwise_complete_tests',
    test_spearman_pvalue_equivalence.R='run_spearman_equivalence_tests',
    test_cell_context_helpers.R='run_cell_context_helper_tests',
    test_cell_context_corrected_helpers.R='run_cell_context_corrected_helper_tests',
    test_mediation_triplet_rollup.R='run_mediation_triplet_rollup_tests',
    test_stage16_kmeans_convergence.R='run_stage16_kmeans_convergence_tests',
    test_runtime.R='test_runtime',
    test_integrated_addons.R='test_integrated_addons',
    test_structure.R='test_structure')
  failed <- character()
  for(f in names(tests)) {
    cat('\n=== ',f,' ===\n',sep='')
    e <- new.env(parent=globalenv());e$project_root <- project_root
    ok <- tryCatch({sys.source(file.path(project_root,'tests',f),envir=e);get(tests[[f]],envir=e)(project_root);TRUE},
      error=function(err){cat('[FAIL] ',conditionMessage(err),'\n',sep='');FALSE})
    if(!ok) failed <- c(failed,f)
  }
  for(f in c('test_addon.R','test_mir146.R','test_stage25.R')) {
    cat('\n=== ',f,' ===\n',sep='')
    e <- new.env(parent=globalenv());e$project_root <- project_root;e$code_root <- file.path(project_root,'scripts/lib')
    ok <- tryCatch({sys.source(file.path(project_root,'tests',f),envir=e);TRUE},
      error=function(err){cat('[FAIL] ',conditionMessage(err),'\n',sep='');FALSE})
    if(!ok) failed <- c(failed,f)
  }
  if(length(failed)) stop('Failed test files: ',paste(failed,collapse=', '),call.=FALSE)
  cat('\nAll 13 synthetic/numerical/integration test files passed.\n')
  invisible(TRUE)
}
if(sys.nframe()==0L) {
  a <- commandArgs(trailingOnly=TRUE)
  f <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
  root <- if(length(a)) a[1] else dirname(dirname(normalizePath(f,winslash='/',mustWork=TRUE)))
  run_tests(root)
}
