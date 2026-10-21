cs_require <- function(x,cols,label) {m<-setdiff(cols,names(x));if(length(m)) stop(label,': missing ',paste(m,collapse=', '))}
cs_existing_context <- function(root, run_id) {
  e <- new.env(parent=globalenv())
  sys.source(file.path(root,'scripts/lib/runtime.R'), envir=e)
  ctx <- e$analysis_context(root)
  if (!identical(ctx$run_id, run_id)) stop('Analysis identity mismatch')
  attr(ctx,'start_time_evidence') <- file.path(ctx$run_root,'manifests/analysis_identity.txt')
  ctx
}

