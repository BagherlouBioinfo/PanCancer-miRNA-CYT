# Exercise the real add-on entry points against synthetic, provenance-registered
# files. No TCGA/CPTAC data or published numerical results are used here.
test_integrated_addons <- function(project_root) {
  suppressPackageStartupMessages(library(data.table))
  temp <- tempfile('integrated_addons_');dir.create(temp)
  on.exit(unlink(temp,recursive=TRUE),add=TRUE)
  for(d in c('scripts','config','docs')) file.copy(file.path(project_root,d),temp,recursive=TRUE)
  e <- new.env(parent=globalenv());sys.source(file.path(temp,'scripts/lib/runtime.R'),envir=e)
  ctx <- e$analysis_context(temp)
  set.seed(1198);n <- 480L;ids <- sprintf('sample_%04d',seq_len(n))
  cancer <- rep(c('BRCA','LUAD'),each=n/2)
  a <- rnorm(n);b <- .4*a+rnorm(n);cy <- .4*a+.3*b+rnorm(n)
  ann <- data.table(sample_id=ids,cancer=cancer,CPE=runif(n,.3,.95),LUMP=runif(n,.3,.95))
  cyt <- data.table(sample_id=ids,CYT=cy)
  cluster <- rep(1:2,length.out=n)
  surv <- data.table(sample_id=ids,OS=as.integer(runif(n)<.9),OS.time=rexp(n,exp(.4*cluster)/1000)+1,
    age_at_initial_pathologic_diagnosis=runif(n,30,80),gender=rep(c('MALE','FEMALE'),length.out=n),
    ajcc_pathologic_tumor_stage=rep(c('Stage I','Stage II','Stage III','Stage IV'),length.out=n),clinical_stage='')
  for(ep in c('DSS','PFI')) {surv[,(ep):=OS];surv[,(paste0(ep,'.time')):=OS.time]}
  member <- data.table(sample_id=ids,cancer=cancer,cluster=cluster)
  primary <- CJ(mirna_id=c('hsa-miR-146a-5p','hsa-miR-146b-5p'),gene=c('CD3D','SPI1'))
  primary[,`:=`(qualifies=TRUE,strict_recurrence=3L,run_id=ctx$run_id)]
  mi <- as.data.table(rbind(a,b));setnames(mi,ids);mi[,mirna_id:=c('hsa-miR-146a-5p','hsa-miR-146b-5p')]
  ge <- as.data.table(rbind(.6*a+.4*b+rnorm(n),.3*a+.5*b+rnorm(n)));setnames(ge,ids);ge[,gene_id_raw:=c('CD3D','SPI1')]
  save <- function(x,rel,stage,rds=FALSE) e$write_run_artifact(ctx,x,rel,stage,
    file.path(temp,'scripts/01_configure_and_validate_inputs.R'),
    write_fn=if(rds) saveRDS else function(x,p) fwrite(x,p,sep='\t'))
  save(surv,'processed/stage03/survival_matched.rds','03_harmonize_identifiers',TRUE)
  save(mi,'processed/stage03/mirna_matched.rds','03_harmonize_identifiers',TRUE)
  save(ge,'processed/stage03/mrna_matched.rds','03_harmonize_identifiers',TRUE)
  save(ann,'processed/stage05/sample_annotation_with_purity.tsv','05_define_eligible_features')
  save(cyt,'processed/stage04/cyt_scores.tsv','04_compute_validate_cyt')
  save(member,'tables/downstream/cluster_membership.tsv','18_clustering_and_state_analysis')
  save(data.table(cancer=unique(cancer),stage_adequate=TRUE),'tables/downstream/clinical_data_completeness_by_cancer.tsv','15_survival_analysis')
  save(primary,'tables/primary/triplets_rho030_CPE.tsv','10_purity_adjusted_analysis')
  for(setting in c('cpe','lump')) save(data.table(cancer=unique(cancer)),paste0('processed/stage05/eligible_cancers_',setting,'.tsv'),'05_define_eligible_features')
  for(stage in c('03_harmonize_identifiers','04_compute_validate_cyt','05_define_eligible_features',
    '10_purity_adjusted_analysis','15_survival_analysis','18_clustering_and_state_analysis')) e$record_stage_status(ctx,stage,'COMPLETE')
  for(script in c('25_cluster_incremental_survival.R','26_mir146_conditional_analysis.R')) sys.source(file.path(temp,'scripts',script),envir=e)
  e$run_cluster_incremental(ctx)
  e$run_mir146_conditional(ctx)
  stopifnot(identical(readLines(file.path(ctx$run_root,'cluster_survival/STATUS.txt')),'COMPLETE'),
    identical(readLines(file.path(ctx$run_root,'mir146/STATUS.txt')),'COMPLETE'))
  for(kind in c('cluster_survival','mir146')) {
    dir <- file.path(ctx$run_root,kind)
    for(p in list.files(file.path(dir,'tables'),pattern='\\.tsv$',full.names=TRUE)) e$verify_table_manifest(p,dir)
  }
  stopifnot(nrow(fread(file.path(ctx$run_root,'cluster_survival/tables/cluster_incremental_survival.tsv')))==12L)
  cat('[PASS] Integrated cluster-survival and miR-146 entry points, context readers, paths, figures and output manifests\n')
  invisible(TRUE)
}
