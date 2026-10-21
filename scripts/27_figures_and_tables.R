# Final main/supplementary figures, tables and editable legends.
# Original plotting functions, dimensions, layouts, palettes and labels retained.
run_figures_tables <- function(ctx) {
  project_root <- ctx$data_root
  code_root <- file.path(project_root,"scripts/figures")
  env <- new.env(parent=environment())
  sys.source(file.path(project_root,"config/figures.R"),envir=env)
  env$CONFIG <- utils::modifyList(env$CONFIG,list(
    run_root=ctx$run_root,run_id=ctx$run_id,output_root=file.path(ctx$run_root,"manuscript"),
    external_root=file.path(ctx$run_root,"external_validation"),
    cluster_survival_root=file.path(ctx$run_root,"cluster_survival"),
    mir146_root=file.path(ctx$run_root,"mir146")))
  overrides <- getOption("mirna.figure.config")
  if(!is.null(overrides)) env$CONFIG <- utils::modifyList(env$CONFIG,overrides)
  if(any(vapply(env$CONFIG[c("external_root","cluster_survival_root","mir146_root")],is.null,logical(1))))
    stop("Current analysis outputs are required. Archived reported results are not a fallback for a new analysis.")
  for(pkg in c("data.table","ggplot2","digest")) if(!requireNamespace(pkg,quietly=TRUE)) stop("Missing package: ",pkg)
  env$OUT <- env$CONFIG$output_root
  expected_out <- file.path(ctx$run_root,'manuscript')
  if(!identical(normalizePath(env$OUT,winslash='/',mustWork=FALSE),normalizePath(expected_out,winslash='/',mustWork=FALSE)))
    stop('Figures and tables must use results/manuscript')
  if(dir.exists(env$OUT)) unlink(env$OUT,recursive=TRUE)
  for(d in c("figures","tables","source_data","text","logs")) dir.create(file.path(env$OUT,d),recursive=TRUE,showWarnings=FALSE)
  env$sources <- list();env$legends <- list();env$build_status <- list();env$table_issues <- list();env$build_warnings <- list()
  env$code_root <- code_root;env$project_root <- project_root
  for(f in c("common.R","validate_outputs.R","figures_main.R","figures_supplementary.R","tables_and_text.R","external_validation.R","review_addons.R"))
    sys.source(file.path(code_root,f),envir=env)
  with(env,{
    st<-read_run("logs/run_status.tsv",c("stage","status"))
    if(!any(grepl("^24_",st$stage)&st$status=="COMPLETE")) stop("No completed final stage in this run")
    validate_outputs()
    validate_external_outputs()
    validate_review_addons()
    if(isTRUE(CONFIG$write_workflow)) {
      py <- Sys.which(CONFIG$python)
      if(!nzchar(py)) py <- Sys.which('python')
      if(!nzchar(py)) stop('Figure 1 needs Python with matplotlib. Install it or set write_workflow=FALSE in config/figures.R.')
      status <- system2(py,c(shQuote(file.path(code_root,'draw_workflow.py')),
        '--output-dir',shQuote(file.path(OUT,'figures'))))
      if(status!=0L) stop('Workflow drawing failed; check Python/matplotlib.')
      legends[['Figure_1']] <- paste(readLines(file.path(code_root,'Figure_1_legend.txt'),warn=FALSE),collapse=' ')
    }
    jobs<-list(Figure_2=figure_discovery,Figure_3=figure_composition,Figure_4=figure_survival,
      Figure_5=figure_cell_context,Figure_S1=figure_meta_sensitivity,Figure_S2=figure_mediation,
      Figure_S3=figure_enrichment,Figure_S4=figure_clustering,Figure_S5=figure_specificity,
      Figure_6=figure_external,Figure_S6=figure_external_sensitivity,
      Figure_S7=figure_cluster_incremental,Figure_S8=figure_mir146_conditional,External_tables=write_external_tables,
      Supplementary_tables=write_tables,Results_draft=write_results_draft)
    if(length(setdiff(CONFIG$modules,names(jobs)))) stop("Unknown configured publication module")
    for(nm in CONFIG$modules) task(nm,jobs[[nm]])
    write_external_handoff()
    write_review_handoff()
    writeLines(unlist(lapply(names(legends),function(nm)c(paste0("## ",nm),"",legends[[nm]],""))),
      file.path(OUT,"text","Figure_legends.md"),useBytes=TRUE)
    status<-do.call(rbind,build_status)
    data.table::fwrite(status,file.path(OUT,"logs","BUILD_STATUS.tsv"),sep="\t")
    data.table::fwrite(do.call(rbind,sources),file.path(OUT,"logs","INPUT_TABLES.tsv"),sep="\t")
    if(length(table_issues)) data.table::fwrite(do.call(rbind,table_issues),
      file.path(OUT,"logs","MISSING_TABLES.tsv"),sep="\t")
    if(length(build_warnings)) data.table::fwrite(do.call(rbind,build_warnings),
      file.path(OUT,"logs","BUILD_WARNINGS.tsv"),sep="\t")
    writeLines(capture.output(sessionInfo()),file.path(OUT,"logs","sessionInfo.txt"))
    file.copy(file.path(project_root,"config/figures.R"),file.path(OUT,"logs","config_used.R"))
    saveRDS(CONFIG,file.path(OUT,"logs","effective_config.rds"))
    for(f in list.files(code_root,pattern="\\.(R|md)$",full.names=TRUE))
      file.copy(f,file.path(OUT,"logs",basename(f)),overwrite=TRUE)
    print(status,row.names=FALSE)
    cat("\nOutput directory:",OUT,"\n")
    if(any(status$status!="OK")||length(table_issues))
      stop("Some outputs are incomplete. Review BUILD_STATUS.tsv and MISSING_TABLES.tsv.")
  })
  register_extra_outputs(ctx,"figures_and_tables",env$OUT,
    file.path(ctx$data_root,"scripts/27_figures_and_tables.R"))
  invisible(env$OUT)
}
