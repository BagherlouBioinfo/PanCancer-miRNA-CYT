test_structure <- function(project_root) {
  e <- new.env(parent=globalenv());sys.source(file.path(project_root,'scripts/run.R'),e)
  r <- e$pipeline_registry(project_root);deps <- e$pipeline_dependencies()
  stopifnot(identical(as.integer(r$number),1:27),all(file.exists(file.path(project_root,'scripts',r$file))))
  constant_path <- function(x) {
    if(is.character(x)&&length(x)==1L) return(x)
    if(is.call(x)&&identical(x[[1]],as.name('file.path'))) {
      args <- as.list(x)[-1]
      if(all(vapply(args,function(a) is.character(a)&&length(a)==1L,logical(1)))) return(do.call(file.path,args))
    }
    NULL
  }
  writes <- list();reads <- vector('list',23)
  for(n in 1:23) {
    found <- character()
    walk <- function(x) {
      if(missing(x)) return(NULL)
      if(is.call(x)&&identical(x[[1]],as.name('write_run_artifact'))&&length(x)>=4) {
        p <- constant_path(x[[4]])
        if(!is.null(p)) writes[[length(writes)+1]] <<- data.frame(stage=n,path=p)
      }
      if(is.call(x)&&identical(x[[1]],as.name('read_run_artifact'))&&length(x)>=4&&is.character(x[[4]]))
        found <<- c(found,x[[4]])
      if(is.call(x)||is.expression(x)||is.pairlist(x)) for(el in as.list(x)) walk(el)
    }
    walk(parse(file.path(project_root,'scripts',r$file[n])))
    reads[[n]] <- setdiff(na.omit(match(unique(found),r$stage_id)),n)
  }
  w <- unique(do.call(rbind,writes))
  if(anyDuplicated(w$path)) stop('Different stages overwrite the same static output path: ',paste(unique(w$path[duplicated(w$path)]),collapse=', '))
  for(n in 1:23) {
    closure <- deps[[n]]
    repeat {
      more <- union(closure,unlist(deps[closure]))
      if(setequal(more,closure)) break
      closure <- more
    }
    if(length(setdiff(reads[[n]],closure))) stop('Dependency graph omits a known producer for stage ',n)
  }
  original_hashes <- read.delim(file.path(project_root,'docs/REFERENCE_CONFIG_HASHES.tsv'),stringsAsFactors=FALSE)
  for(i in seq_len(nrow(original_hashes))) stopifnot(digest::digest(file=file.path(project_root,original_hashes$file[i]),algo='sha256')==original_hashes$sha256[i])
  # A report-only directory must never be accepted as a clustering cache.
  temp <- tempfile('incomplete_clustering_');dir.create(temp)
  on.exit(unlink(temp,recursive=TRUE),add=TRUE)
  writeLines('synthetic report',file.path(temp,'report.md'))
  sys.source(file.path(project_root,'scripts/import_clustering.R'),envir=e)
  copy <- tempfile('import_destination_');dir.create(copy)
  on.exit(unlink(copy,recursive=TRUE),add=TRUE)
  for(d in c('scripts','config','docs')) file.copy(file.path(project_root,d),copy,recursive=TRUE)
  err <- tryCatch({e$import_clustering(copy,temp);NULL},error=identity)
  stopifnot(inherits(err,'error'),grepl('full provenance',conditionMessage(err),fixed=TRUE))
  cat('[PASS] 27-stage registry, ',nrow(w),' distinct static output paths, known producer dependencies, unchanged scientific configuration and incomplete-cache rejection\n',sep='')
  invisible(TRUE)
}
