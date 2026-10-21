# Explicit reuse of completed ORIGINAL Stage 18 (new numbered Stage 17).
# Requires the full original run directory and its original project code.
# Values are copied byte-for-byte; imported provenance is retained separately.
import_clustering <- function(project_root,source_run) {
  root <- normalizePath(project_root,winslash='/',mustWork=TRUE)
  old <- normalizePath(source_run,winslash='/',mustWork=TRUE)
  e <- new.env(parent=globalenv());sys.source(file.path(root,'scripts/lib/runtime.R'),envir=e)
  ctx <- e$analysis_context(root)
  if(identical(old,ctx$run_root)) stop('Source and destination are identical')
  if(dir.exists(file.path(ctx$run_root,'manifests/execution.lock'))) stop('Stop the pipeline before importing clustering')
  sm <- file.path(old,'manifests/artifact_provenance.tsv');ss <- file.path(old,'logs/run_status.tsv')
  if(!all(file.exists(c(sm,ss)))) stop('Original full provenance and status files are required; summary tables/images are insufficient')
  m <- read.delim(sm,stringsAsFactors=FALSE);s <- read.delim(ss,stringsAsFactors=FALSE)
  sid <- '18_clustering_and_state_analysis'
  z <- s[s$stage==sid,,drop=FALSE]
  if(!nrow(z)||tail(z$status,1)!='COMPLETE') stop('Original clustering did not complete')
  rows <- m[m$stage==sid,,drop=FALSE]
  rows <- rows[!duplicated(rows$relative_path,fromLast=TRUE),,drop=FALSE]
  required <- paste0('tables/downstream/',c('cluster_validity','cluster_stability','cluster_null_calibration',
    'cluster_membership','cluster_characteristics','cluster_pca_coordinates'),'.tsv')
  if(!all(required%in%rows$relative_path)) stop('Clustering cache is incomplete')
  sources <- read.delim(file.path(root,'docs/SOURCE_MAP.tsv'),stringsAsFactors=FALSE)
  expected <- sources$source_sha256[match('scripts/17_clustering_and_state_analysis.R',sources$destination)]
  if(length(expected)!=1L||is.na(expected)||any(rows$script_sha256!=expected)) stop('Original clustering script differs from the supplied scientific source')
  if(any(rows$config_hash!=e$.current_config_hash(root))) stop('Original analysis configuration differs')
  # Check the original scientific helpers too, rather than trusting only
  # the stage filename recorded in the older provenance format.
  old_project <- dirname(dirname(old))
  for(h in c('stats_engine.R','expression_preprocessing.R','run_discovery_source.R')) {
    source <- file.path(old_project,'R',h)
    target <- paste0('scripts/lib/',h)
    sha <- sources$source_sha256[match(target,sources$destination)]
    if(!file.exists(source)||is.na(sha)||digest::digest(file=source,algo='sha256')!=sha)
      stop('Original scientific helper missing or changed: ',source)
  }
  # Verify every original artifact against its original recorded hash.
  for(i in seq_len(nrow(rows))) {
    p <- file.path(old,rows$relative_path[i])
    if(!file.exists(p)||digest::digest(file=p,algo='sha256')!=rows$artifact_sha256[i])
      stop('Original cache checksum mismatch: ',p)
  }
  inputs <- c('tables/primary/triplets_rho030_CPE.tsv','processed/stage03/mirna_matched.rds',
    'processed/stage03/mrna_matched.rds','processed/stage04/cyt_scores.tsv',
    'processed/stage05/sample_annotation_with_purity.tsv','processed/stage11_immune/sample_level_immune_estimates.tsv')
  stages <- c('10_purity_adjusted_analysis',rep('03_harmonize_identifiers',2),
    '04_compute_validate_cyt','05_define_eligible_features','12_immune_composition_analysis')
  clean <- function(x) {
    if(inherits(x,'data.table')) x <- as.data.frame(x)
    if(is.data.frame(x)&&'run_id'%in%names(x)) x$run_id <- NULL
    x
  }
  for(i in seq_along(inputs)) {
    rel <- inputs[i];p <- file.path(old,rel)
    record <- tail(m[m$relative_path==rel,,drop=FALSE],1)
    if(!file.exists(p)||!nrow(record)||digest::digest(file=p,algo='sha256')!=record$artifact_sha256)
      stop('Original input unverifiable: ',rel)
    reader <- if(grepl('\\.rds$',rel)) readRDS else function(p) as.data.frame(data.table::fread(p))
    current <- e$read_run_artifact(ctx,rel,stages[i],reader)
    previous <- reader(p)
    if(!isTRUE(all.equal(clean(current),clean(previous),tolerance=0,check.attributes=TRUE)))
      stop('Clustering input differs from the original: ',rel)
  }
  # Do not import before validation is complete. Imported figures and report
  # keep their original embedded dates/identity and their source SHA digests.
  importdir <- file.path(ctx$run_root,'manifests/clustering_import')
  dir.create(importdir,recursive=TRUE,showWarnings=FALSE)
  write.table(rows,file.path(importdir,'original_artifact_provenance.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
  writeLines(c(paste('Source:',old),paste('Imported UTC:',format(Sys.time(),tz='UTC',usetz=TRUE)),
    'Mode: imported original clustering; not recomputed',
    'Six inputs matched exactly after removal of analysis-identity metadata only.',
    'Original stage script, scientific helpers, configuration and cache checksums verified.'),
    file.path(importdir,'IMPORT_RECORD.txt'))
  e$record_stage_status(ctx,sid,'RUNNING','Explicit original clustering import')
  for(i in seq_len(nrow(rows))) {
    rel <- rows$relative_path[i];dest <- e$run_output_path(ctx,rel)
    dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
    if(!file.copy(file.path(old,rel),dest,overwrite=TRUE)) stop('Copy failed: ',rel)
    e$write_run_artifact(ctx,NULL,rel,sid,file.path(root,'scripts/17_clustering_and_state_analysis.R'),
      input_paths=file.path(ctx$run_root,inputs),write_fn=function(x,p) invisible(NULL))
  }
  e$record_stage_status(ctx,sid,'COMPLETE','Imported original Stage 18; see manifests/clustering_import')
  saveRDS(list(fingerprint=e$stage_code_fingerprint(root,'17_clustering_and_state_analysis.R'),
    completed_utc=Sys.time(),mode='IMPORTED_ORIGINAL'),file.path(ctx$run_root,'manifests/stage_17.rds'))
  cat('Original clustering imported without recomputation. Continue with stages 18-27.\n')
  invisible(ctx)
}
if(sys.nframe()==0L) {
  a <- commandArgs(trailingOnly=TRUE)
  if(length(a)!=2L) stop('Usage: Rscript scripts/import_clustering.R /new/project /original/project/runs/original_run')
  import_clustering(a[1],a[2])
}
