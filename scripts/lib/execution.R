# Execution, restart and provenance checks for a stable results/ directory.
# No thresholds or statistical estimators are defined here.
.sha_cache <- new.env(parent=emptyenv())
file_sha256_cached <- function(path) {
  p <- normalizePath(path,winslash='/',mustWork=TRUE)
  info <- file.info(p)
  stamp <- c(as.numeric(info$size),as.numeric(info$mtime),as.numeric(info$ctime))
  old <- .sha_cache[[p]]
  if(!is.null(old)&&identical(old$stamp,stamp)) return(old$hash)
  hash <- digest::digest(file=p,algo='sha256')
  .sha_cache[[p]] <- list(stamp=stamp,hash=hash)
  hash
}
analysis_context <- function(project_root='.') {
  root <- normalizePath(project_root,winslash='/',mustWork=TRUE)
  out <- file.path(root,'results')
  for(d in c('processed','tables','figures','logs','reports','checkpoints','manifests'))
    dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
  id <- file.path(out,'manifests/analysis_identity.txt')
  if(!file.exists(id)) writeLines(c('analysis',format(Sys.time(),tz='UTC','%Y-%m-%d %H:%M:%OS6')),id)
  identity <- readLines(id,warn=FALSE)
  if(length(identity)!=2L || !nzchar(identity[1])) stop('Invalid analysis identity file')
  start <- as.POSIXct(identity[2],tz='UTC')
  if(is.na(start)) stop('Invalid recorded analysis start time')
  manifest <- file.path(out,'manifests/artifact_provenance.tsv')
  status <- file.path(out,'logs/run_status.tsv')
  if(!file.exists(manifest)) writeLines(paste(c('run_id','stage','relative_path','generated_at_utc',
    'script_sha256','config_hash','input_paths','input_hashes','artifact_sha256'),collapse='\t'),manifest)
  if(!file.exists(status)) writeLines('run_id\tstage\tstatus\ttimestamp_utc\tdetail',status)
  structure(list(run_id=identity[1],run_root=out,data_root=root,start_time_utc=start,
    manifest_path=manifest,status_path=status),class='analysis_context')
}

forbid_legacy_path <- function(path,data_root,ctx=NULL) {
  if(is.null(ctx)) return(invisible(TRUE))
  root <- normalizePath(ctx$run_root,winslash='/',mustWork=TRUE)
  # Resolve existing parents too: normalizePath(mustWork=FALSE) alone does
  # not canonicalise a nonexistent leaf or an ancestor symlink consistently.
  tail <- character();p <- path
  while(!file.exists(p) && !dir.exists(p)) {
    par <- dirname(p)
    if(identical(par,p)) stop('Cannot resolve output path: ',path)
    tail <- c(basename(p),tail);p <- par
  }
  p <- normalizePath(p,winslash='/',mustWork=TRUE)
  if(length(tail)) for(x in tail) {
    if(x%in%c('.','..') || grepl('[/\\\\]',x)) stop('Unsafe output path: ',path)
    p <- paste0(p,'/',x)
  }
  if(!(identical(p,root)||startsWith(p,paste0(root,'/')))) stop('PATH ESCAPES RESULTS ROOT: ',path)
  invisible(TRUE)
}
run_output_path <- function(ctx,...) {
  p <- file.path(ctx$run_root,...);forbid_legacy_path(p,ctx$data_root,ctx);p
}
stage_completed_this_run <- function(ctx,stage) {
  s <- read.delim(ctx$status_path,stringsAsFactors=FALSE)
  z <- s[s$run_id==ctx$run_id & s$stage==stage,,drop=FALSE]
  nrow(z)>0L && identical(tail(z$status,1),'COMPLETE')
}

read_run_artifact <- function(ctx,rel_path,expected_stage,reader,.seen=character()) {
  if(rel_path%in%.seen) stop('Cyclic artifact provenance: ',rel_path)
  p <- run_output_path(ctx,rel_path)
  if(!file.exists(p)) stop('Missing upstream artifact: ',rel_path)
  m <- read.delim(ctx$manifest_path,stringsAsFactors=FALSE,check.names=FALSE)
  z <- m[m$relative_path==rel_path,,drop=FALSE]
  if(!nrow(z)) stop('No provenance record: ',rel_path)
  z <- tail(z,1)
  if(!identical(z$run_id,ctx$run_id)) stop('Analysis identity mismatch: ',rel_path)
  if(!identical(z$stage,expected_stage)) stop('Generating-stage mismatch: ',rel_path)
  if(!stage_completed_this_run(ctx,expected_stage)) stop('Upstream stage is incomplete or stale: ',expected_stage)
  registry <- read.delim(file.path(ctx$data_root,'config/stages.tsv'),stringsAsFactors=FALSE)
  idx <- match(expected_stage,registry$stage_id)
  if(!is.na(idx)) {
    meta <- file.path(ctx$run_root,'manifests',sprintf('stage_%02d.rds',registry$number[idx]))
    if(file.exists(meta) && !identical(readRDS(meta)$fingerprint,stage_code_fingerprint(ctx$data_root,registry$file[idx])))
      stop('Producing code changed; rerun stage ',registry$number[idx],': ',rel_path)
  }
  if(!identical(file_sha256_cached(p),z$artifact_sha256)) stop('Artifact checksum mismatch: ',rel_path)
  if(!identical(z$config_hash,.current_config_hash(ctx$data_root))) stop('Configuration changed; rerun producing stage: ',expected_stage)
  # Parent input digests catch upstream replacement without requiring unique
  # run folders or rejecting an intact output solely because it is old.
  if(!is.na(z$input_paths) && nzchar(z$input_paths)) {
    inputs <- strsplit(z$input_paths,';',fixed=TRUE)[[1]]
    hashes <- strsplit(z$input_hashes,';',fixed=TRUE)[[1]]
    if(length(inputs)!=length(hashes)) stop('Malformed input provenance: ',rel_path)
    for(i in seq_along(inputs)) {
      if(!file.exists(inputs[i]) || file_sha256_cached(inputs[i])!=hashes[i])
        stop('Upstream input changed or missing; rerun ',expected_stage,': ',inputs[i])
      canonical <- normalizePath(inputs[i],winslash='/',mustWork=TRUE)
      if(startsWith(canonical,paste0(ctx$run_root,'/'))) {
        parent_rel <- substring(canonical,nchar(ctx$run_root)+2)
        parent_row <- tail(m[m$relative_path==parent_rel,,drop=FALSE],1)
        if(nrow(parent_row) && parent_rel!=rel_path)
          read_run_artifact(ctx,parent_rel,parent_row$stage,function(p) TRUE,.seen=c(.seen,rel_path))
      }
    }
  }
  reader(p)
}

stage_code_fingerprint <- function(root,file) {
  n <- suppressWarnings(as.integer(substr(file,1,2)))
  addons <- c('cluster_survival.R','context_reader.R','mir146_analysis.R','plot_mir146.R','stage25_external_validation.R')
  libs <- list.files(file.path(root,'scripts/lib'),pattern='\\.R$',full.names=TRUE)
  libs <- libs[!basename(libs)%in%addons]
  selected <- if(n==24) 'stage25_external_validation.R' else if(n==25) 'cluster_survival.R' else if(n==26)
    c('context_reader.R','mir146_analysis.R','plot_mir146.R') else character()
  paths <- c(file.path(root,'scripts',file),
    libs,file.path(root,'scripts/lib',selected),file.path(root,'config',c('analysis_config.R','sample_exclusions.tsv','FROM_SCRATCH_ALLOWED_INPUTS.tsv','stages.tsv')))
  if(n==24) paths <- c(paths,file.path(root,'config/external_validation.R'))
  if(n==27) paths <- c(paths,file.path(root,'config/figures.R'),
    list.files(file.path(root,'scripts/figures'),pattern='\\.(R|py|txt)$',full.names=TRUE))
  paths <- sort(unique(paths[file.exists(paths)]))
  digest::digest(vapply(paths,file_sha256_cached,character(1)),algo='sha256')
}

register_extra_outputs <- function(ctx,stage,dir,script_path,input_paths=character()) {
  paths <- list.files(dir,recursive=TRUE,full.names=TRUE)
  paths <- paths[!dir.exists(paths)]
  for(p in paths) {
    rel <- substring(normalizePath(p,winslash='/'),nchar(ctx$run_root)+2)
    # Register an existing file without rewriting its scientific payload.
    write_run_artifact(ctx,NULL,rel,stage,script_path,input_paths,
      write_fn=function(x,p) invisible(NULL))
  }
}

verify_table_manifest <- function(path,root) {
  f <- file.path(root,'output_sha256.tsv')
  if(!file.exists(f)) stop('Missing output manifest: ',f)
  m <- data.table::fread(f,data.table=FALSE)
  rel <- substring(normalizePath(path,winslash='/'),nchar(normalizePath(root,winslash='/'))+2)
  i <- match(rel,m$file)
  if(is.na(i) || digest::digest(file=path,algo='sha256')!=m$sha256[i]) stop('Result checksum mismatch: ',path)
  invisible(TRUE)
}

# Allow scripts/lib/runtime.R to expose this implementation to all readers.
