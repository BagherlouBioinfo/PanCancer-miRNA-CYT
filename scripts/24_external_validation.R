run_external_validation <- function(ctx) {
  e <- new.env(parent=globalenv())
  sys.source(file.path(ctx$data_root,'config/external_validation.R'),envir=e)
  sys.source(file.path(ctx$data_root,'scripts/lib/stage25_external_validation.R'),envir=e)
  cfg <- e$STAGE25_CONFIG
  cfg$project_root <- ctx$data_root;cfg$discovery_run_id <- ctx$run_id
  cfg$output_parent <- file.path(ctx$run_root,'external_validation')
  e$run_stage25(cfg,file.path(ctx$data_root,'scripts/lib'))
  register_extra_outputs(ctx,'external_validation',cfg$output_parent,
    file.path(ctx$data_root,'scripts/24_external_validation.R'),
    file.path(ctx$run_root,c('tables/primary/triplets_rho030_CPE.tsv',
      'tables/downstream/meta_analysis_input_all_cancers.tsv')))
}
