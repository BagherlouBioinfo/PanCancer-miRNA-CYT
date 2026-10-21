# Read-only package/input inventory. Does not install or run any analysis.
check_environment <- function(project_root='.') {
  cran <- c('data.table','digest','survival','cluster','metafor','httr2','jsonlite','ggplot2','sandwich')
  bioc <- c('TCGAbiolinks','org.Hs.eg.db','GO.db','reactome.db','AnnotationDbi')
  packages <- c(cran,bioc,'MCPcounter','openxlsx')
  required <- !packages%in%c('MCPcounter','openxlsx')
  installed <- vapply(packages,requireNamespace,logical(1),quietly=TRUE)
  tab <- data.frame(package=packages,required=required,installed=installed,
    version=vapply(packages,function(p) if(requireNamespace(p,quietly=TRUE)) as.character(packageVersion(p)) else NA_character_,character(1)))
  print(tab,row.names=FALSE)
  cat('\nMCPcounter: recommended; the original documented marker-score fallback is retained.\n')
  cat('openxlsx: needed for XLSX export; TSV exports remain available without it.\n')
  e <- new.env();sys.source(file.path(project_root,'config/analysis_config.R'),envir=e)
  keys <- c('raw_mrna','raw_mirna','raw_survival','targetscan','mir_family_info','mcpcounter_genes','mcpcounter_signature_manifest')
  paths <- unlist(e$analysis_config$paths[keys],use.names=FALSE)
  full <- file.path(project_root,paths)
  print(data.frame(path=paths,exists=file.exists(full)),row.names=FALSE)
  cat('\nThis inventory does not replace the SHA-256 input validation in stage 01.\n')
  invisible(all(installed[required])&&all(file.exists(full)))
}
if(sys.nframe()==0L) {
  f <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
  root <- dirname(dirname(normalizePath(f,winslash='/',mustWork=TRUE)))
  if(!check_environment(root)) quit(status=1L)
}
