# Command line: Rscript scripts/run.R --stages=1-27
# RStudio: source('scripts/run.R'); run_pipeline(stages=1:27)
.runner_files <- unlist(lapply(sys.frames(),function(f) {
  x <- if(!is.null(f$ofile)) f$ofile else f$file
  if(is.character(x)&&length(x)==1L&&grepl('(^|[/\\\\])run\\.R$',x)&&file.exists(x)) x else NULL
}))
.runner_cli <- grep('^--file=',commandArgs(),value=TRUE)
.runner_cli <- .runner_cli[file.exists(sub('^--file=','',.runner_cli))]
.runner_path <- if(length(.runner_files)) tail(.runner_files,1) else if(length(.runner_cli)) sub('^--file=','',.runner_cli[1]) else 'scripts/run.R'
.default_root <- dirname(dirname(normalizePath(.runner_path,winslash='/',mustWork=TRUE)))

pipeline_registry <- function(root) {
  r <- read.delim(file.path(root,'config/stages.tsv'),stringsAsFactors=FALSE,check.names=FALSE)
  extra <- data.frame(number=24:27,file=c('24_external_validation.R','25_cluster_incremental_survival.R',
    '26_mir146_conditional_analysis.R','27_figures_and_tables.R'),
    function_name=c('run_external_validation','run_cluster_incremental','run_mir146_conditional','run_figures_tables'),
    stage_id=c('external_validation','cluster_incremental_survival','mir146_conditional_analysis','figures_and_tables'),
    source_number=c('25','addon','addon','figures'),stringsAsFactors=FALSE)
  names(extra)[3] <- 'function'
  rbind(r,extra)
}
pipeline_dependencies <- function() list(
  integer(),1L,2L,3L,c(3L,4L),c(3L,4L,5L),c(3L,5L,6L),c(3L,4L,5L,6L,7L),8L,
  c(3L,4L,5L,7L,9L),c(3L,4L,5L,7L,10L),c(3L,4L,5L,10L),
  c(3L,4L,5L,10L,12L),c(3L,4L,5L,10L),c(3L,4L,5L,10L,14L),
  c(3L,4L,5L,10L,13L),c(3L,4L,5L,10L,12L),c(3L,5L,7L,10L,12L),
  c(3L,7L,10L),c(7L,8L,9L,10L,12L),10L,c(3L,4L,5L,10L,12L),
  c(3L,4L,5L,10L,12L,22L),c(5L,10L,13L,23L),c(3L,4L,5L,14L,17L),
  c(3L,4L,5L,10L),1:26)

run_pipeline <- function(project_root=.default_root,stages=1:27,force=FALSE,dry_run=FALSE) {
  root <- normalizePath(project_root,winslash='/',mustWork=TRUE)
  stages <- as.integer(stages)
  if(!length(stages)||anyNA(stages)||any(!stages%in%1:27)||anyDuplicated(stages)||is.unsorted(stages))
    stop('Stages must be distinct increasing integers between 1 and 27')
  r <- pipeline_registry(root);deps <- pipeline_dependencies()
  if(dry_run) {print(r[r$number%in%stages,c('number','file','stage_id')],row.names=FALSE);return(invisible(r))}
  env <- new.env(parent=globalenv())
  sys.source(file.path(root,'scripts/lib/runtime.R'),envir=env)
  ctx <- env$analysis_context(root)
  lock <- file.path(ctx$run_root,'manifests/execution.lock')
  if(!dir.create(lock,showWarnings=FALSE)) stop('Another execution may be active. Inspect results/manifests/execution.lock before removing it.')
  writeLines(c(paste('PID',Sys.getpid()),paste('Host',Sys.info()[['nodename']]),paste('UTC',Sys.time())),file.path(lock,'owner.txt'))
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  writeLines(capture.output(sessionInfo()),file.path(ctx$run_root,'logs/sessionInfo.txt'))
  for(n in stages) {
    rm(list=ls(env$.sha_cache,all.names=TRUE),envir=env$.sha_cache)
    row <- r[r$number==n,,drop=FALSE];sid <- row$stage_id
    meta <- file.path(ctx$run_root,'manifests',sprintf('stage_%02d.rds',n))
    fingerprint <- env$stage_code_fingerprint(root,row$file)
    valid <- env$stage_completed_this_run(ctx,sid) && file.exists(meta)
    if(valid) {
      saved <- readRDS(meta)
      valid <- identical(saved$fingerprint,fingerprint)
      if(valid) valid <- tryCatch({
        m <- read.delim(ctx$manifest_path,stringsAsFactors=FALSE)
        latest <- m[!duplicated(m$relative_path,fromLast=TRUE),,drop=FALSE]
        own <- latest[latest$stage==sid,,drop=FALSE]
        if(!nrow(own)) stop('No registered outputs')
        for(p in own$relative_path) env$read_run_artifact(ctx,p,sid,function(p) TRUE)
        TRUE
      },error=function(e){message('Stale stage ',n,': ',conditionMessage(e));FALSE})
    }
    if(valid&&!force) {cat(sprintf('[SKIP] %02d %s: verified complete\n',n,row$file));next}
    # A failed/repeated producer invalidates transitive consumers. Unrelated
    # branches (including clustering after a survival-only change) stay valid.
    affected <- n
    repeat {
      nextset <- union(affected,which(vapply(deps,function(d) any(d%in%affected),logical(1))))
      if(setequal(nextset,affected)) break
      affected <- nextset
    }
    for(j in setdiff(affected,n)) if(env$stage_completed_this_run(ctx,r$stage_id[j]))
      env$record_stage_status(ctx,r$stage_id[j],'STALE',paste('Upstream stage restarted:',n))
    missing <- deps[[n]][!vapply(deps[[n]],function(j) env$stage_completed_this_run(ctx,r$stage_id[j]),logical(1))]
    if(length(missing)) stop('Stage ',n,' needs completed stages: ',paste(missing,collapse=', '),'. Run those stages first.')
    cat(sprintf('\n[START] %02d %s\n',n,row$file))
    env$record_stage_status(ctx,sid,'RUNNING')
    warning_rows <- list()
    result <- tryCatch(withCallingHandlers({
      stage_env <- new.env(parent=env)
      sys.source(file.path(root,'scripts',row$file),envir=stage_env)
      get(row[['function']],envir=stage_env)(ctx)
      TRUE
    },warning=function(w) {
      warning_rows[[length(warning_rows)+1L]] <<- data.frame(message=conditionMessage(w),call=paste(deparse(conditionCall(w)),collapse=' '))
    }),error=function(e) e)
    if(length(warning_rows)) write.table(do.call(rbind,warning_rows),
      file.path(ctx$run_root,'logs',sprintf('stage_%02d_warnings.tsv',n)),sep='\t',row.names=FALSE)
    if(inherits(result,'error')) {
      env$record_stage_status(ctx,sid,'FAILED',gsub('[\t\n\r]',' ',conditionMessage(result)))
      stop('Stage ',n,' failed: ',conditionMessage(result),'\nFix the cause and rerun this stage; earlier completed outputs are retained.',call.=FALSE)
    }
    env$record_stage_status(ctx,sid,'COMPLETE')
    saveRDS(list(fingerprint=fingerprint,completed_utc=Sys.time(),script=row$file),meta)
    cat(sprintf('[COMPLETE] %02d\n',n))
  }
  cat('\nResults:',ctx$run_root,'\n')
  invisible(ctx)
}

if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  root_arg <- grep('^--project=',args,value=TRUE)
  stage_arg <- grep('^--stages=',args,value=TRUE)
  allowed <- grepl('^--(project|stages)=',args)|args%in%c('--force','--dry-run','--help')
  if(any(!allowed)) stop('Unknown arguments: ',paste(args[!allowed],collapse=' '))
  if('--help'%in%args) {
    cat('Rscript scripts/run.R [--project=/path/to/project] [--stages=1-27|17|24,25,26,27] [--force] [--dry-run]\n')
  } else {
    spec <- if(length(stage_arg)) sub('^--stages=','',stage_arg[1]) else '1-27'
    if(grepl('^[0-9]+-[0-9]+$',spec)) {bounds <- as.integer(strsplit(spec,'-',fixed=TRUE)[[1]]);selection <- seq.int(bounds[1],bounds[2])}
    else selection <- as.integer(strsplit(spec,',',fixed=TRUE)[[1]])
    run_pipeline(if(length(root_arg)) sub('^--project=','',root_arg[1]) else .default_root,
      selection,force='--force'%in%args,dry_run='--dry-run'%in%args)
  }
}
