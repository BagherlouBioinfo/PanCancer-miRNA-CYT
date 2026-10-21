# Explicit installation only: source() defines the function without installing.
# Rscript scripts/install_dependencies.R --install
install_dependencies <- function(tests_only=FALSE) {
  cran <- c('data.table','digest','survival','cluster','metafor','httr2','jsonlite','ggplot2','sandwich','openxlsx','BiocManager','remotes')
  if(isTRUE(tests_only)) cran <- c('data.table','digest','survival','cluster','metafor','ggplot2','sandwich')
  missing <- cran[!vapply(cran,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) install.packages(missing,repos='https://cloud.r-project.org',dependencies=NA)
  if(isTRUE(tests_only)) {cat('Test dependencies installed. Run Rscript tests/run_tests.R.\n');return(invisible(TRUE))}
  # R 4.3 uses Bioconductor 3.18; newer R uses its compatible release.
  if(getRversion()>='4.3'&&getRversion()<'4.4'&&as.character(BiocManager::version())!='3.18')
    BiocManager::install(version='3.18',ask=FALSE,update=FALSE)
  bioc <- c('TCGAbiolinks','org.Hs.eg.db','GO.db','reactome.db','AnnotationDbi')
  missing <- bioc[!vapply(bioc,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) BiocManager::install(missing,ask=FALSE,update=FALSE)
  if(!requireNamespace('MCPcounter',quietly=TRUE))
    remotes::install_github('ebecht/MCPcounter',ref='master',subdir='Source',upgrade='never',dependencies=NA)
  cat('Installation finished. Run scripts/check_environment.R before analysis.\n')
}
if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if(identical(args,'--tests-only')) install_dependencies(tests_only=TRUE)
  else if(identical(args,'--install')) install_dependencies()
  else stop('Use --install for all dependencies or --tests-only for synthetic checks')
}
