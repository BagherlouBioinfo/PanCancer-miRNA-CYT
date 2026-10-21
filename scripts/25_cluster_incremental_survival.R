run_cluster_incremental <- function(ctx) {
  e <- new.env(parent=globalenv())
  lib <- file.path(ctx$data_root,'scripts/lib')
  sys.source(file.path(lib,'cluster_survival.R'),envir=e)
  out <- e$run_cluster_survival(ctx$data_root,ctx$run_id,
    file.path(ctx$run_root,'cluster_survival'),code_root=lib)
  register_extra_outputs(ctx,'cluster_incremental_survival',out,
    file.path(ctx$data_root,'scripts/25_cluster_incremental_survival.R'),
    file.path(ctx$run_root,'tables/downstream/cluster_membership.tsv'))
}
