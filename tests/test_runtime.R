test_runtime <- function(project_root) {
  e <- new.env(parent=globalenv());sys.source(file.path(project_root,'scripts/lib/runtime.R'),envir=e)
  temp <- tempfile('runtime_test_');dir.create(temp)
  on.exit(unlink(temp,recursive=TRUE),add=TRUE)
  dir.create(file.path(temp,'config'));dir.create(file.path(temp,'scripts/lib'),recursive=TRUE)
  for(p in list.files(file.path(project_root,'scripts/lib'),full.names=TRUE)) file.copy(p,file.path(temp,'scripts/lib'))
  for(p in list.files(file.path(project_root,'config'),full.names=TRUE)) file.copy(p,file.path(temp,'config'))
  file.copy(file.path(project_root,'scripts/run.R'),file.path(temp,'scripts/run.R'))
  check <- function(ok,label) {if(!isTRUE(ok)) stop(label);cat('[PASS] ',label,'\n',sep='')}
  raises <- function(expr) inherits(tryCatch({force(expr);NULL},error=identity),'error')
  ctx <- e$analysis_context(temp)
  check(identical(ctx$run_root,file.path(normalizePath(temp,winslash='/'),'results')),'stable results directory')
  check(identical(ctx,e$analysis_context(temp)),'context can reconnect without creating a new run')
  input <- file.path(temp,'data.txt');writeLines('unchanged',input)
  script <- file.path(temp,'scripts/test.R');writeLines('# test',script)
  rel <- 'tables/test.tsv'
  e$write_run_artifact(ctx,data.frame(x=1),rel,'test',script,input,
    write_fn=function(x,p) write.table(x,p,sep='\t',row.names=FALSE))
  check(raises(e$read_run_artifact(ctx,rel,'test',read.delim)),'incomplete producer rejected')
  e$record_stage_status(ctx,'test','COMPLETE')
  check(e$read_run_artifact(ctx,rel,'test',read.delim)$x==1,'verified artifact accepted')
  check(raises(e$read_run_artifact(ctx,rel,'wrong',read.delim)),'wrong generating stage rejected')
  writeLines('modified',input)
  check(raises(e$read_run_artifact(ctx,rel,'test',read.delim)),'changed upstream input rejected')
  writeLines('unchanged',input)
  cat('\n',file=file.path(ctx$run_root,rel),append=TRUE)
  check(raises(e$read_run_artifact(ctx,rel,'test',read.delim)),'modified artifact rejected')
  e$record_stage_status(ctx,'test','FAILED')
  check(!e$stage_completed_this_run(ctx,'test'),'latest failed status overrides earlier completion')
  check(raises(e$run_output_path(ctx,'../outside.tsv')),'path traversal rejected')
  if(.Platform$OS.type!='windows') {
    file.symlink(temp,file.path(ctx$run_root,'escape'))
    check(raises(e$run_output_path(ctx,'escape/uncreated.tsv')),'nonexistent leaf through escaping symlink rejected')
  }
  # Test the actual runner with tiny synthetic stages rather than invoking
  # any biological analysis or stubbing the restart implementation.
  registry <- data.frame(number=1:2,file=c('01_fake.R','02_fake.R'),fun=c('fake1','fake2'),
    stage_id=c('first','second'),source_number=c('test','test'))
  names(registry)[3] <- 'function'
  write.table(registry,file.path(temp,'config/stages.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
  one <- "fake1 <- function(ctx) {write_run_artifact(ctx,data.frame(x=1),'tables/first.tsv','first',file.path(ctx$data_root,'scripts/01_fake.R'),write_fn=function(x,p) write.table(x,p,sep='\\t',row.names=FALSE))}"
  two <- "fake2 <- function(ctx) {x<-read_run_artifact(ctx,'tables/first.tsv','first',read.delim);write_run_artifact(ctx,x,'tables/second.tsv','second',file.path(ctx$data_root,'scripts/02_fake.R'),file.path(ctx$run_root,'tables/first.tsv'),write_fn=function(x,p) write.table(x,p,sep='\\t',row.names=FALSE))}"
  writeLines(one,file.path(temp,'scripts/01_fake.R'))
  writeLines("fake2 <- function(ctx) stop('deliberate synthetic failure')",file.path(temp,'scripts/02_fake.R'))
  runner <- new.env(parent=globalenv());sys.source(file.path(temp,'scripts/run.R'),envir=runner)
  runner$pipeline_registry <- function(root) registry
  runner$pipeline_dependencies <- function() list(integer(),1L)
  check(raises(runner$run_pipeline(temp,1:2)),'runner stops on failure')
  old_hash <- digest::digest(file=file.path(ctx$run_root,'tables/first.tsv'),algo='sha256')
  writeLines(two,file.path(temp,'scripts/02_fake.R'))
  runner$run_pipeline(temp,2L)
  check(digest::digest(file=file.path(ctx$run_root,'tables/first.tsv'),algo='sha256')==old_hash,'failed stage resumes without rewriting earlier output')
  check(e$stage_completed_this_run(ctx,'second'),'resumed stage completes')
  before <- readRDS(file.path(ctx$run_root,'manifests/stage_02.rds'))$completed_utc
  runner$run_pipeline(temp,1:2)
  check(identical(before,readRDS(file.path(ctx$run_root,'manifests/stage_02.rds'))$completed_utc),'verified completed stages are skipped')
  runner$run_pipeline(temp,1L,force=TRUE)
  check(!e$stage_completed_this_run(ctx,'second'),'rerun producer invalidates downstream completion')
  check(!dir.exists(file.path(ctx$run_root,'manifests/execution.lock')),'execution lock removed after success or failure')
  invisible(TRUE)
}
