# Frozen-label, cancer-specific incremental survival analysis. No reclustering.
suppressPackageStartupMessages({library(data.table); library(survival); library(digest)})
cs_require <- function(x, cols, label) {
  miss <- setdiff(cols, names(x)); if(length(miss)) stop(label, ': missing ', paste(miss, collapse=', '))
}
cs_unique <- function(x, label) {
  cs_require(x, 'sample_id', label)
  if(anyNA(x$sample_id) || anyDuplicated(x$sample_id)) stop(label, ': missing/duplicate sample_id')
}
cs_fit <- function(d, covars, tv=FALSE) {
  messages <- character()
  f <- as.formula(paste('Surv(time, event) ~', if(length(covars)) paste(covars, collapse=' + ') else '1',
                 if(tv && length(covars)) paste('+', paste(paste0('tt(',covars,')'),collapse=' + ')) else ''))
  ans <- tryCatch(withCallingHandlers(
    coxph(f, data=d, ties='efron', x=!tv, y=TRUE, singular.ok=FALSE,
          control=coxph.control(iter.max=50),
          tt=if(tv) function(x,t,...) x*log((t+1)/366) else NULL),
    warning=function(w){messages <<- c(messages, conditionMessage(w)); invokeRestart('muffleWarning')}),
    error=function(e) e)
  if(inherits(ans,'error')) return(list(fit=NULL, reason=conditionMessage(ans)))
  if(length(messages) || any(!is.finite(coef(ans))) || any(!is.finite(vcov(ans))) ||
     any(diag(vcov(ans))<=0) || ans$iter>=50)
    return(list(fit=NULL, reason=paste(c('fit rejected',messages),collapse='; ')))
  list(fit=ans,reason='')
}
cs_ph <- function(fit) {
  z <- tryCatch(cox.zph(fit, transform='km'),error=function(e)e)
  if(inherits(z,'error')) return(data.table(term='GLOBAL',chisq=NA_real_,df=NA_real_,p=NA_real_,diagnostic=conditionMessage(z)))
  data.table(term=rownames(z$table),chisq=z$table[,'chisq'],df=z$table[,'df'],p=z$table[,'p'],diagnostic='')
}
cs_compare <- function(d, min_n=30, epp=10, tv=FALSE, max_coef=10) {
  vars <- c('cyt_z','age_z','sex01',if('stage_z'%in%names(d)) 'stage_z')
  vars <- vars[vapply(vars,function(v) length(unique(d[[v]]))>1,logical(1))]
  if(length(unique(d$cluster01))!=2) return(list(status='NOT_ESTIMABLE',reason='both frozen clusters required'))
  if(nrow(d)<min_n) return(list(status='NOT_ESTIMABLE',reason='too few complete cases'))
  np <- (length(vars)+1L)*if(tv) 2L else 1L
  if(sum(d$event)<epp*np) return(list(status='NOT_ESTIMABLE',reason=sprintf('events %d < required %d',sum(d$event),epp*np)))
  # No-event cluster groups can cause an infinite Cox coefficient.
  if(any(tapply(d$event,d$cluster01,sum)==0)) return(list(status='NOT_ESTIMABLE',reason='zero events in one cluster'))
  b <- cs_fit(d,vars,tv); f <- cs_fit(d,c(vars,'cluster01'),tv)
  if(is.null(b$fit)||is.null(f$fit)) return(list(status='FIT_FAILED',reason=paste(b$reason,f$reason)))
  if(any(abs(coef(b$fit))>max_coef)||any(abs(coef(f$fit))>max_coef))
    return(list(status='FIT_CAUTION',reason='coefficient exceeds configured convergence bound'))
  stopifnot(b$fit$n==f$fit$n, b$fit$n==nrow(d), b$fit$nevent==f$fit$nevent)
  lr <- 2*(as.numeric(logLik(f$fit))-as.numeric(logLik(b$fit)))
  df <- length(coef(f$fit))-length(coef(b$fit))
  if(!is.finite(lr)||lr < -1e-6) stop('Invalid nested likelihood comparison')
  list(status='ESTIMATED',reason='', baseline=b$fit,full=f$fit,
       covariates=paste(vars,collapse='+'),lrt=max(0,lr),df=df,
       p=pchisq(max(0,lr),df,lower.tail=FALSE))
}
cs_prepare <- function(surv, ann, cyt, membership, stage_fun) {
  for(nm in c('surv','ann','cyt','membership')) cs_unique(get(nm),nm)
  cs_require(surv,c('OS','OS.time','age_at_initial_pathologic_diagnosis','gender',
                   'ajcc_pathologic_tumor_stage','clinical_stage'),'survival')
  cs_require(ann,'cancer','annotation');cs_require(cyt,'CYT','CYT');cs_require(membership,c('cluster','cancer'),'membership')
  if(!setequal(surv$sample_id,ann$sample_id)||!setequal(ann$sample_id,cyt$sample_id)||
     !setequal(ann$sample_id,membership$sample_id)) stop('Input sample sets do not match; refusing silent inner-join loss')
  m <- match(ann$sample_id,membership$sample_id)
  if(anyNA(membership$cancer[m]) || any(ann$cancer!=membership$cancer[m])) stop('Cancer assignments disagree')
  d <- copy(surv); if('cancer'%in%names(d)) d[,cancer:=NULL]; d <- merge(d,ann[,.(sample_id,cancer)],by='sample_id',all=TRUE)
  if('CYT'%in%names(d)) d[,CYT:=NULL]
  d <- merge(d,cyt[,.(sample_id,CYT)],by='sample_id',all=TRUE)
  d <- merge(d,membership[,.(sample_id,cluster)],by='sample_id',all=TRUE)
  levels <- sort(unique(as.character(d$cluster)))
  if(length(levels)!=2 || anyNA(d$cluster)) stop('This add-on requires exactly two frozen clusters')
  a <- trimws(as.character(d$ajcc_pathologic_tumor_stage)); b <- trimws(as.character(d$clinical_stage))
  # Mirrors Stage 15 coalescing: nonblank AJCC first; no new clinical recoding.
  d[,stage_numeric:=stage_fun(ifelse(a!='',a,b))]
  d[,age:=as.numeric(age_at_initial_pathologic_diagnosis)]
  d[,sex01:=fifelse(trimws(as.character(gender))=='MALE',1,
                     fifelse(trimws(as.character(gender))=='FEMALE',0,NA_real_))]
  d[,cluster01:=as.integer(as.character(cluster)==levels[2])]
  attr(d,'cluster_levels') <- levels;d
}
# Reconnect read-only using persisted evidence; never use file ctime or Sys.time().
cs_existing_context <- function(root, run_id) {
  e <- new.env(parent=globalenv())
  sys.source(file.path(root,'scripts/lib/runtime.R'), envir=e)
  ctx <- e$analysis_context(root)
  if (!identical(ctx$run_id, run_id)) stop('Analysis identity mismatch')
  attr(ctx,'start_time_evidence') <- file.path(ctx$run_root,'manifests/analysis_identity.txt')
  ctx
}

run_cluster_survival <- function(project_root, run_id, output_root, code_root=NULL) {
  for(p in c('data.table','survival','digest','ggplot2')) if(!requireNamespace(p,quietly=TRUE)) stop('Missing package: ',p)
  root <- normalizePath(project_root,winslash='/',mustWork=TRUE)
  run <- normalizePath(file.path(root,'results'),winslash='/',mustWork=TRUE)
  outroot <- normalizePath(output_root,winslash='/',mustWork=FALSE)
  expected_out <- normalizePath(file.path(root,'results/cluster_survival'),winslash='/',mustWork=FALSE)
  if(!identical(outroot,expected_out)) stop('Cluster outputs must use results/cluster_survival')
  # Add-on files are separate from discovery artifacts.
  env <- new.env(parent=globalenv())
  helpers <- file.path(root,c('scripts/lib/runtime.R','scripts/lib/clinical_stage_mapping.R','config/analysis_config.R'))
  if(any(!file.exists(helpers))) stop('Missing project helper: ',paste(helpers[!file.exists(helpers)],collapse=', '))
  for(h in helpers) sys.source(h,envir=env)
  cfg <- env$analysis_config
  ctx <- cs_existing_context(root,run_id)
  # Reuse the project's same-run, provenance and checksum guard.
  inputs <- c('processed/stage03/survival_matched.rds','processed/stage05/sample_annotation_with_purity.tsv',
              'processed/stage04/cyt_scores.tsv','tables/downstream/cluster_membership.tsv',
              'tables/downstream/clinical_data_completeness_by_cancer.tsv')
  stages <- c('03_harmonize_identifiers','05_define_eligible_features','04_compute_validate_cyt',
              '18_clustering_and_state_analysis','15_survival_analysis')
  objects <- lapply(seq_along(inputs),function(i) as.data.table(env$read_run_artifact(ctx,inputs[i],stages[i],
    if(grepl('rds$',inputs[i])) readRDS else fread)))
  before <- vapply(unique(c(file.path(run,inputs),helpers,ctx$manifest_path,ctx$status_path,attr(ctx,'start_time_evidence'))),function(p) digest::digest(file=p,algo='sha256'),character(1))
  d <- cs_prepare(objects[[1]],objects[[2]],objects[[3]],objects[[4]],env$stage_to_numeric)
  completeness <- objects[[5]];cs_require(completeness,c('cancer','stage_adequate'),'saved completeness')
  if(anyDuplicated(completeness$cancer)) stop('Duplicate cancer in saved completeness')
  levels <- attr(d,'cluster_levels')
  out <- outroot
  if(dir.exists(out)) unlink(out,recursive=TRUE)
  dir.create(out,recursive=TRUE);for(x in c('tables','figures','logs','inputs')) dir.create(file.path(out,x))
  if(!is.null(code_root)) {
    code_files <- file.path(code_root,c('cluster_survival.R'))
    if(any(!file.exists(code_files))) stop('Add-on source files missing')
    dir.create(file.path(out,'code_snapshot'))
    file.copy(code_files,file.path(out,'code_snapshot'),overwrite=FALSE)
  }
  writeLines('RUNNING',file.path(out,'STATUS.txt'))
  writeLines(c(paste('Start UTC:',format(ctx$start_time_utc,tz='UTC',usetz=TRUE)),
    paste('Evidence:',attr(ctx,'start_time_evidence'))),file.path(out,'inputs/run_context.txt'))
  completed <- FALSE
  on.exit(if(!completed) writeLines('FAILED - inspect error; do not use partial outputs',file.path(out,'STATUS.txt')),add=TRUE)
  writeLines(capture.output(sessionInfo()),file.path(out,'logs/sessionInfo.txt'))
  fwrite(data.table(path=names(before),sha256=unname(before)),file.path(out,'inputs/input_sha256.tsv'),sep='\t')
  settings <- data.table(setting=c('min_n','events_per_parameter','stage_completeness_pct','PH_alpha','time_transform','status'),
    value=c(cfg$discovery$min_samples_per_cancer,cfg$survival$events_per_predictor_minimum,
      cfg$survival$stage_completeness_adequacy_pct,cfg$survival$ph_violation_alpha,
      'x * log((time_days + 1)/366); reference day 365','post-review exploratory; not original prespecification'))
  fwrite(settings,file.path(out,'inputs/settings.tsv'),sep='\t')
  fwrite(data.table(cluster=levels,cluster01=0:1),file.path(out,'inputs/cluster_contrast.tsv'),sep='\t')
  rows<-list();phrows<-list();effects<-list();samples<-list();audit<-list()
  min_n<-cfg$discovery$min_samples_per_cancer;epp<-cfg$survival$events_per_predictor_minimum
  for(ep in c('OS','DSS','PFI')) for(cc in sort(unique(d$cancer))) for(tier in c('CYT_age_sex','CYT_age_sex_stage')) {
    message(ep,' / ',cc,' / ',tier)
    row <- data.table(endpoint=ep,cancer=cc,tier=tier,n=0L,events=0L,n_reference=0L,n_comparison=0L,
      status='NOT_ESTIMABLE',reason='',covariates=NA_character_,HR=NA_real_,lower=NA_real_,upper=NA_real_,
      lrt=NA_real_,lrt_df=NA_integer_,lrt_p=NA_real_,c_baseline=NA_real_,c_full=NA_real_,c_delta=NA_real_,
      ph_global_p=NA_real_,ph_cluster_p=NA_real_,ph_baseline_p=NA_real_,
      tv_status='NOT_RUN',tv_reason='OS sensitivity only',tv_lrt=NA_real_,tv_df=NA_integer_,tv_p=NA_real_)
    z <- copy(d[cancer==cc]); row[,n:=nrow(z)]
    if(!all(c(ep,paste0(ep,'.time'))%in%names(z))) {row[,reason:='endpoint unavailable'];rows[[length(rows)+1]]<-row;next}
    if(tier=='CYT_age_sex_stage' && !isTRUE(completeness[cancer==cc]$stage_adequate)) {
      row[,reason:='not stage-adequate according to saved Stage 15 table'];rows[[length(rows)+1]]<-row;next}
    z[,time:=as.numeric(get(paste0(ep,'.time')))];z[,event:=as.numeric(get(ep))]
    needed <- c('time','event','age','sex01','CYT','cluster01',if(tier=='CYT_age_sex_stage') 'stage_numeric')
    keep <- complete.cases(z[,..needed]) & Reduce(`&`,lapply(z[,..needed],is.finite)) & z$time>0 & z$event%in%c(0,1)
    audit[[length(audit)+1]]<-data.table(endpoint=ep,cancer=cc,tier=tier,n_source=nrow(z),n_complete=sum(keep),n_excluded=sum(!keep))
    z<-z[keep]
    # Scale within the identical model-specific complete-case set. Constant nuisance variables omitted.
    scale0 <- function(x) if(length(unique(x))>1) as.numeric(scale(x)) else rep(0,length(x))
    z[,cyt_z:=scale0(CYT)];z[,age_z:=scale0(age)]
    if(tier=='CYT_age_sex_stage') z[,stage_z:=scale0(stage_numeric)]
    row[,`:=`(n=nrow(z),events=sum(z$event),n_reference=sum(z$cluster01==0),n_comparison=sum(z$cluster01==1))]
    samples[[length(samples)+1]]<-z[,.(sample_id,endpoint=ep,cancer,tier=tier)]
    a<-cs_compare(z,min_n,epp,max_coef=cfg$survival$convergence_max_abs_coef)
    row[,`:=`(status=a$status,reason=a$reason)]
    if(a$status=='ESTIMATED') {
      beta<-coef(a$full)['cluster01'];se<-sqrt(vcov(a$full)['cluster01','cluster01'])
      row[,`:=`(covariates=a$covariates,HR=exp(beta),lower=exp(beta-1.96*se),upper=exp(beta+1.96*se),
        lrt=a$lrt,lrt_df=a$df,lrt_p=a$p,c_baseline=unname(summary(a$baseline)$concordance[1]),
        c_full=unname(summary(a$full)$concordance[1]))]
      row[,c_delta:=c_full-c_baseline]
      for(kind in c('baseline','full')) {
        ph<-cs_ph(a[[kind]]);ph[,`:=`(endpoint=ep,cancer=cc,tier=tier,model=kind)]
        phrows[[length(phrows)+1]]<-ph
        if(kind=='full') row[,`:=`(ph_global_p=ph[term=='GLOBAL']$p,ph_cluster_p=if(any(ph$term=='cluster01')) ph[term=='cluster01']$p else NA_real_)]
        else row[,ph_baseline_p:=ph[term=='GLOBAL']$p]
      }
      # Fixed sensitivity for ALL estimable OS models, not selected by PH or significance.
      if(ep=='OS') {
        v<-cs_compare(z,min_n,epp,tv=TRUE,max_coef=cfg$survival$convergence_max_abs_coef)
        row[,`:=`(tv_status=v$status,tv_reason=v$reason)]
        if(v$status=='ESTIMATED') {
          row[,`:=`(tv_lrt=v$lrt,tv_df=v$df,tv_p=v$p)]
          for(day in c(365,1095,1825)) {
            w<-c(1,log((day+1)/366));nm<-c('cluster01','tt(cluster01)')
            bb<-sum(coef(v$full)[nm]*w);ss<-sqrt(drop(t(w)%*%vcov(v$full)[nm,nm]%*%w))
            nr<-sum(z$time>=day);supported<-nr>=10 && day<=max(z$time)
            effects[[length(effects)+1]]<-data.table(endpoint=ep,cancer=cc,tier=tier,day=day,n_at_risk=nr,
              supported=supported,HR=if(supported) exp(bb) else NA_real_,
              lower=if(supported) exp(bb-1.96*ss) else NA_real_,upper=if(supported) exp(bb+1.96*ss) else NA_real_)
          }
        }
      }
    }
    rows[[length(rows)+1]]<-row
  }
  results<-rbindlist(rows,fill=TRUE)
  # Six separate primary/sensitivity families: endpoint x clinical tier, across cancers.
  results[,lrt_fdr_BH:=p.adjust(lrt_p,'BH',n=.N),by=.(endpoint,tier)]
  results[,tv_fdr_BH:=p.adjust(tv_p,'BH',n=.N),by=.(endpoint,tier)]
  alpha<-cfg$survival$ph_violation_alpha
  results[,ph_global_flag:=!is.na(ph_global_p)&ph_global_p<alpha]
  results[,ph_cluster_flag:=!is.na(ph_cluster_p)&ph_cluster_p<alpha]
  results[,ph_status:=fcase(is.na(ph_global_p),'Unavailable',ph_global_p<alpha,'Flagged',default='Not flagged')]
  write_table<-function(x,name) fwrite(x,file.path(out,'tables',name),sep='\t',na='NA')
  write_table(results,'cluster_incremental_survival.tsv')
  write_table(if(length(phrows)) rbindlist(phrows,fill=TRUE) else data.table(term=character(),p=numeric()),'PH_diagnostics.tsv')
  write_table(if(length(effects)) rbindlist(effects) else data.table(day=numeric(),HR=numeric()),'OS_time_varying_cluster_HR.tsv')
  write_table(rbindlist(audit,fill=TRUE),'sample_flow.tsv')
  write_table(rbindlist(samples,fill=TRUE),'model_sample_ids.tsv')
  summ<-results[,.(attempted=.N,estimated=sum(status=='ESTIMATED'),BH_significant=sum(lrt_fdr_BH<0.05,na.rm=TRUE),
    PH_flagged=sum(ph_global_flag),TV_estimated=sum(tv_status=='ESTIMATED'),TV_BH_significant=sum(tv_fdr_BH<0.05,na.rm=TRUE)),by=.(endpoint,tier)]
  write_table(summ,'summary.tsv')
  # Editable ggplot objects plus vector PDF and 600 dpi PNG, no embedded figure titles.
  plotdata<-results[status=='ESTIMATED'&is.finite(HR)&is.finite(lower)&is.finite(upper)&lower>0&upper>0]
  if(nrow(plotdata)) {
    library(ggplot2)
    plotdata[,cancer:=factor(cancer,levels=rev(sort(unique(cancer))))]
    plotdata[,endpoint:=factor(endpoint,levels=c('OS','DSS','PFI'))]
    plotdata[,tier:=factor(tier,levels=c('CYT_age_sex','CYT_age_sex_stage'),
      labels=c('CYT + age + sex','CYT + age + sex + tumour stage'))]
    p1<-ggplot(plotdata,aes(x=HR,y=cancer))+geom_vline(xintercept=1,colour='grey65',linetype=2)+
      geom_segment(aes(x=lower,xend=upper,yend=cancer),linewidth=.35)+geom_point(aes(shape=ph_status),size=1.8)+
      scale_x_log10()+facet_grid(tier~endpoint,scales='free_y',space='free_y')+
      labs(x=paste0('Hazard ratio: cluster ',levels[2],' vs ',levels[1]),y=NULL,shape='Global PH flag')+theme_bw(base_size=10)+theme(legend.position='bottom')
    p2<-ggplot(plotdata,aes(x=c_delta,y=cancer))+geom_vline(xintercept=0,colour='grey65',linetype=2)+
      geom_point(size=1.8)+facet_grid(tier~endpoint,scales='free_y',space='free_y')+
      labs(x='Change in apparent concordance',y=NULL)+theme_bw(base_size=10)
    for(i in 1:2) {
      p<-if(i==1) p1 else p2
      saveRDS(p,file.path(out,'figures',paste0('panel_',LETTERS[i],'.rds')))
      ggsave(file.path(out,'figures',paste0('panel_',LETTERS[i],'.pdf')),p,width=11,height=10)
      ggsave(file.path(out,'figures',paste0('panel_',LETTERS[i],'.png')),p,width=11,height=10,dpi=600)
    }
    draw_combined <- function() {
    grid::grid.newpage();grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,2)))
    print(p1,vp=grid::viewport(layout.pos.row=1,layout.pos.col=1));print(p2,vp=grid::viewport(layout.pos.row=1,layout.pos.col=2))
    grid::grid.text('A',x=.01,y=.99,just=c('left','top'),gp=grid::gpar(fontface='bold'))
    grid::grid.text('B',x=.51,y=.99,just=c('left','top'),gp=grid::gpar(fontface='bold'))
    }
    grDevices::pdf(file.path(out,'figures','Figure_cluster_survival_AB.pdf'),width=15,height=10)
    draw_combined();grDevices::dev.off()
    grDevices::png(file.path(out,'figures','Figure_cluster_survival_AB.png'),width=15,height=10,units='in',res=600)
    draw_combined();grDevices::dev.off()
  }
  writeLines(c('# Frozen-cluster incremental survival analysis','',paste('Discovery run:',run_id),
    'Exploratory analysis added after peer review. Frozen labels; no reclustering, feature reselection, outcome cutpoints, or external prognostic validation.',
    'Within each endpoint/cancer/tier, baseline and full models use identical complete cases. Age, CYT and ordinal stage are standardised within that set; constant nuisance variables are omitted.',
    'Primary endpoint: OS. DSS and PFI are exploratory endpoint sensitivities, not universally recommended clinical endpoints across all cancers.',
    'Primary comparisons: CYT + age + sex versus the same model plus fixed cluster; stage-adjusted comparisons use saved Stage 15 eligibility.',
    'BH correction is separate by endpoint and tier across ALL attempted cancer comparisons, retaining non-estimable tests in the denominator.',
    'Global and term-specific Schoenfeld PH tests are diagnostic, nominal and unadjusted. A non-significant test does not prove PH.',
    'For every estimable OS comparison, a fixed time-varying sensitivity is attempted: all retained nuisance covariates have main and x*log((t+1)/366) terms in BOTH models; the full model adds cluster and tt(cluster), with a joint 2-df LRT.',
    'The time-varying family has a separate BH correction per tier. EPP is applied to the doubled parameter count. This model is a sensitivity analysis, not a guaranteed correction of all non-PH forms. cox.zph is not applied to tt models.',
    'HR(t) is reported at 1, 3 and 5 years only where at least 10 participants remain at risk. These are pointwise, not simultaneous confidence intervals.',
    'Changes in concordance are in-sample apparent values, not validated prediction gains. Fixed cluster contrasts preserve original numeric label order.',
    'Warning-generating, singular, non-finite or convergence-caution fits are excluded from inference and reported with reasons.',
    '',capture.output(print(summ))),file.path(out,'ANALYSIS_REPORT.md'))
  writeLines(c('Panel A: Cancer-specific cluster hazard ratios (comparison and reference labels in inputs/cluster_contrast.tsv), adjusted for continuous CYT, age, sex and, in the second row, ordinal tumour stage. Bars show pointwise 95% confidence intervals. Shapes identify nominal global PH violations; PH-flagged estimates require the accompanying sensitivity analysis.',
    'Panel B: Change in apparent Harrell concordance after adding the frozen cluster to the identical complete-case clinical/CYT model. Values describe internal fit and do not establish external prognostic performance.',
    'OS is primary; DSS and PFI are exploratory sensitivities. Missing cancer/panel combinations were not estimable; see full tables. Survival associations do not validate cluster biology or causal effects.'),file.path(out,'FIGURE_LEGEND.txt'))
  after<-vapply(names(before),function(p) digest::digest(file=p,algo='sha256'),character(1))
  if(!identical(before,after)) stop('Input changed during analysis')
  writeLines(capture.output(sessionInfo()),file.path(out,'logs/sessionInfo_end.txt'))
  files<-list.files(out,recursive=TRUE,full.names=TRUE);files<-files[basename(files)!='STATUS.txt']
  fwrite(data.table(file=substring(files,nchar(out)+2),sha256=vapply(files,function(p) digest::digest(file=p,algo='sha256'),character(1))),file.path(out,'output_sha256.tsv'),sep='\t')
  writeLines('COMPLETE',file.path(out,'STATUS.txt'));completed<-TRUE
  message('COMPLETE: ',normalizePath(out,winslash='/'));invisible(out)
}
