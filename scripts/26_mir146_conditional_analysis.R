run_mir146_conditional <- function(ctx) {
  e <- new.env(parent=globalenv());lib <- file.path(ctx$data_root,'scripts/lib')
  sys.source(file.path(lib,'context_reader.R'),envir=e)
  sys.source(file.path(lib,'mir146_analysis.R'),envir=e)
  out <- e$run_mir146(ctx$data_root,ctx$run_id,file.path(ctx$run_root,'mir146'),lib)
  register_extra_outputs(ctx,'mir146_conditional_analysis',out,
    file.path(ctx$data_root,'scripts/26_mir146_conditional_analysis.R'),
    file.path(ctx$run_root,'tables/primary/triplets_rho030_CPE.tsv'))
}
