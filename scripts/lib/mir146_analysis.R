suppressPackageStartupMessages({library(data.table);library(digest)})
m146_rank <- function(x) {
  x <- pmin(pmax(x,quantile(x,.001,names=FALSE,type=7)),quantile(x,.999,names=FALSE,type=7))
  r <- rank(x,ties.method='average');if(sd(r)==0) return(rep(0,length(x)))
  as.numeric(scale(r))
}
m146_fit <- function(y,a,b,p,min_n=30) {
  if(length(unique(c(length(y),length(a),length(b),length(p))))!=1) stop('Unequal vector lengths')
  ok <- is.finite(y)&is.finite(a)&is.finite(b)&is.finite(p)
  d<-data.frame(y=y[ok],a=a[ok],b=b[ok],p=p[ok]);n<-nrow(d)
  bad <- function(why) list(status='NOT_ESTIMABLE',reason=why,n=n,keep=which(ok))
  if(n<min_n) return(bad('fewer than minimum complete observations'))
  if(any(vapply(d,function(x) length(unique(x))<2,logical(1)))) return(bad('constant outcome, miRNA or purity'))
  d[]<-lapply(d,m146_rank)
  X<-model.matrix(~a+b+p,d);q<-qr(X)
  if(q$rank<ncol(X)) return(bad('rank-deficient joint design; effects cannot be separated'))
  f<-lm.fit(X,d$y);h<-rowSums(qr.Q(q)^2)
  if(any(h>=.99)) return(bad('maximum leverage >= 0.99'))
  if(sum(f$residuals^2)<1e-12) return(bad('near-perfect outcome fit; inference degenerate'))
  inv<-tryCatch(solve(crossprod(X)),error=function(e)NULL)
  if(is.null(inv)) return(bad('singular crossproduct'))
  # HC3 covariance; normal-reference Wald inference is an asymptotic sensitivity.
  w<-(f$residuals/(1-h))^2
  V<-inv %*% crossprod(X,X*w) %*% inv
  fullR2<-1-sum(f$residuals^2)/sum((d$y-mean(d$y))^2)
  stats<-lapply(c('a','b'),function(v){
    other<-if(v=='a') 'b' else 'a'
    Z<-model.matrix(reformulate(c(other,'p')),d)
    ry<-lm.fit(Z,d$y)$residuals;rx<-lm.fit(Z,d[[v]])$residuals
    if(sum(rx^2)<1e-12) return(NULL)
    rho<-cor(rx,ry);se<-sqrt(V[v,v]);beta<-unname(f$coefficients[v])
    if(!is.finite(se)||se<=0) return(NULL)
    r2reduced<-1-sum(ry^2)/sum((d$y-mean(d$y))^2)
    data.table(focal=v,partial_rho=rho,beta_rank=beta,HC3_se=se,
      beta_lower=beta-qnorm(.975)*se,beta_upper=beta+qnorm(.975)*se,
      p_HC3=2*pnorm(-abs(beta/se)),R2_joint=fullR2,R2_reduced=r2reduced,
      delta_R2=max(0,fullR2-r2reduced),partial_R2=rho^2,
      VIF=sum((d[[v]]-mean(d[[v]]))^2)/sum(rx^2),max_leverage=max(h))
  })
  if(any(vapply(stats,is.null,logical(1)))) return(bad('degenerate conditional variance'))
  list(status='ESTIMATED',reason='',n=n,keep=which(ok),effects=rbindlist(stats),
       mirna_rho=cor(d$a,d$b),mirna_partial_purity=cor(lm.fit(model.matrix(~p,d),d$a)$residuals,
                                                      lm.fit(model.matrix(~p,d),d$b)$residuals))
}
m146_extract <- function(tab,idcol,requested,ids,uppercase=FALSE) {
  tab<-as.data.table(tab);if(!idcol%in%names(tab)) stop('Missing expression identifier: ',idcol)
  keys<-trimws(as.character(tab[[idcol]]));if(uppercase) keys<-toupper(keys)
  if(anyNA(match(requested,keys))) stop('Missing expression features: ',paste(setdiff(requested,keys),collapse=', '))
  if(any(vapply(requested,function(k) sum(keys==k)!=1L,logical(1)))) stop('Duplicate requested expression feature')
  if(!all(ids%in%names(tab))) stop('Expression sample columns do not cover annotation')
  mat<-as.matrix(tab[match(requested,keys),..ids])
  if(!is.numeric(mat)) stop('Expression values are not numeric')
  rownames(mat)<-requested;mat
}
run_mir146 <- function(project_root,run_id,output_root,code_root) {
  root<-normalizePath(project_root,winslash='/',mustWork=TRUE)
  outroot<-normalizePath(output_root,winslash='/',mustWork=FALSE)
  expected_out<-normalizePath(file.path(root,'results/mir146'),winslash='/',mustWork=FALSE)
  if(!identical(outroot,expected_out)) stop('Conditional outputs must use results/mir146')
  e<-new.env(parent=globalenv())
  helper_paths<-file.path(root,c('scripts/lib/runtime.R','scripts/lib/run_discovery_source.R','config/analysis_config.R'))
  if(any(!file.exists(helper_paths))) stop('Missing project helper')
  for(h in helper_paths) sys.source(h,envir=e)
  ctx<-cs_existing_context(root,run_id)
  cfg<-e$analysis_config;min_n<-cfg$discovery$min_samples_per_cancer
  trip<-as.data.table(e$load_run_primary_triplets(ctx))
  if('qualifies'%in%names(trip)) trip<-trip[qualifies%in%TRUE]
  cs_require(trip,c('mirna_id','gene'),'primary pairs')
  if('run_id'%in%names(trip) && any(is.na(trip$run_id)|trip$run_id!=run_id)) stop('Primary table run mismatch')
  if(anyDuplicated(trip[,.(mirna_id,gene)])) stop('Duplicate primary pair')
  aa<-'hsa-miR-146a-5p';bb<-'hsa-miR-146b-5p'
  if(!all(c(aa,bb)%in%trip$mirna_id)) stop('Both mature miRNAs must occur in frozen primary set')
  shared<-sort(intersect(trip[mirna_id==aa]$gene,trip[mirna_id==bb]$gene))
  if(!length(shared)) stop('No shared primary genes')
  if(any(shared%in%c('GZMA','PRF1'))) stop('CYT-defining gene in shared set; explicit tautology handling needed')
  paths<-c('processed/stage03/mirna_matched.rds','processed/stage03/mrna_matched.rds',
    'processed/stage05/sample_annotation_with_purity.tsv','processed/stage04/cyt_scores.tsv',
    'processed/stage05/eligible_cancers_cpe.tsv','processed/stage05/eligible_cancers_lump.tsv')
  stages<-c(rep('03_harmonize_identifiers',2),'05_define_eligible_features','04_compute_validate_cyt',rep('05_define_eligible_features',2))
  readi<-function(i) e$read_run_artifact(ctx,paths[i],stages[i],if(i<=2) readRDS else fread)
  ann<-as.data.table(readi(3));cy<-as.data.table(readi(4))
  for(x in list(ann,cy)) if(anyNA(x$sample_id)||anyDuplicated(x$sample_id)) stop('Missing/duplicate annotation sample_id')
  cs_require(ann,c('sample_id','cancer','CPE','LUMP'),'annotation');cs_require(cy,c('sample_id','CYT'),'CYT')
  if(!setequal(ann$sample_id,cy$sample_id)) stop('Annotation/CYT sample mismatch')
  ids<-ann$sample_id;ann[,CYT:=cy$CYT[match(ids,cy$sample_id)]]
  for(v in c('CPE','LUMP')) if(any(is.finite(ann[[v]]) & (ann[[v]]<0|ann[[v]]>1))) stop('Purity outside [0,1]: ',v)
  mi<-m146_extract(readi(1),'mirna_id',c(aa,bb),ids)
  ge<-m146_extract(readi(2),'gene_id_raw',shared,ids,uppercase=TRUE)
  elig<-list(CPE=as.data.table(readi(5)),LUMP=as.data.table(readi(6)))
  for(v in names(elig)) cs_require(elig[[v]],'cancer',paste(v,'eligibility'))
  # Hash every input actually read by the same-run helper, including primary source records.
  # Track any primary artifact loaded by the source helper using the manifest's stage 10 primary outputs.
  manifest<-fread(ctx$manifest_path)
  primary_rel<-manifest[stage=='10_purity_adjusted_analysis' & grepl('^tables/primary/',gsub('\\\\','/',relative_path))]$relative_path
  input_files<-unique(c(file.path(ctx$run_root,c(paths,primary_rel)),helper_paths,ctx$manifest_path,ctx$status_path,
                         attr(ctx,'start_time_evidence')))
  before<-vapply(input_files,function(p) digest(file=p,algo='sha256'),character(1))
  out<-outroot
  if(dir.exists(out)) unlink(out,recursive=TRUE)
  for(v in c('tables','figures','inputs','code_snapshot','logs')) dir.create(file.path(out,v),recursive=TRUE)
  complete<-FALSE;on.exit(if(!complete) writeLines('FAILED; partial outputs must not be used',file.path(out,'STATUS.txt')),add=TRUE)
  writeLines('RUNNING',file.path(out,'STATUS.txt'))
  file.copy(list.files(code_root,pattern='\\.R$',full.names=TRUE),file.path(out,'code_snapshot'))
  fwrite(data.table(path=names(before),sha256=unname(before)),file.path(out,'inputs/input_sha256.tsv'),sep='\t')
  fwrite(trip[mirna_id%in%c(aa,bb)],file.path(out,'inputs/frozen_mir146_primary_pairs.tsv'),sep='\t')
  fwrite(data.table(outcome=c('CYT',shared),type=c('CYT',rep('shared_primary_gene',length(shared)))),
    file.path(out,'inputs/fixed_outcomes.tsv'),sep='\t')
  results<-list();redundancy<-list();sample_rows<-list()
  for(setting in c('CPE','LUMP')) for(cc in sort(unique(elig[[setting]]$cancer))) for(outcome in c('CYT',shared)) {
    message(setting,' / ',cc,' / ',outcome)
    idx<-which(ann$cancer==cc)
    y<-if(outcome=='CYT') ann$CYT[idx] else as.numeric(ge[outcome,idx])
    a<-m146_fit(y,as.numeric(mi[aa,idx]),as.numeric(mi[bb,idx]),ann[[setting]][idx],min_n)
    base<-data.table(setting=setting,cancer=cc,outcome=outcome,outcome_type=if(outcome=='CYT') 'CYT' else 'shared_gene',
      n_source=length(idx),n=a$n,status=a$status,reason=a$reason)
    if(a$status=='ESTIMATED') {
      r<-cbind(base[rep(1,2)],a$effects)
      r[,mirna_id:=ifelse(focal=='a',aa,bb)]
      r[,other_mirna:=ifelse(focal=='a',bb,aa)]
      redundancy[[length(redundancy)+1]]<-data.table(setting=setting,cancer=cc,outcome=outcome,n=a$n,
        mirna_rho=a$mirna_rho,mirna_partial_purity=a$mirna_partial_purity)
    } else r<-cbind(base[rep(1,2)],data.table(focal=c('a','b'),mirna_id=c(aa,bb),other_mirna=c(bb,aa)))
    results[[length(results)+1]]<-r
    if(length(a$keep)) sample_rows[[length(sample_rows)+1]]<-data.table(setting=setting,cancer=cc,outcome=outcome,sample_id=ids[idx[a$keep]])
  }
  r<-rbindlist(results,fill=TRUE)
  for(v in c('partial_rho','beta_rank','HC3_se','beta_lower','beta_upper','p_HC3','R2_joint','R2_reduced','delta_R2','partial_R2','VIF','max_leverage'))
    if(!v%in%names(r)) r[,(v):=NA_real_]
  # Planned families include both miRNAs and all eligible cancers; no significant-only selection.
  r[,FDR_BH:=p.adjust(p_HC3,'BH',n=.N),by=.(setting,outcome_type)]
  r[,FDR_BY:=p.adjust(p_HC3,'BY',n=.N),by=.(setting,outcome_type)]
  r[,VIF_flag:=!is.na(VIF)&VIF>10]
  summary<-r[,.(attempted=.N,estimated=sum(status=='ESTIMATED'),BH_supported=sum(FDR_BH<.05,na.rm=TRUE),
    BY_supported=sum(FDR_BY<.05,na.rm=TRUE),VIF_over_10=sum(VIF_flag)),by=.(setting,outcome_type,mirna_id)]
  joint<-r[,.(both_estimated=all(status=='ESTIMATED'),both_BH_supported=all(!is.na(FDR_BH)&FDR_BH<.05)),by=.(setting,cancer,outcome)]
  save<-function(x,f) fwrite(x,file.path(out,'tables',f),sep='\t',na='NA')
  save(r,'conditional_mir146_associations.tsv');save(summary,'summary.tsv');save(joint,'both_mirnas_support.tsv')
  save(if(length(redundancy)) rbindlist(redundancy,fill=TRUE) else data.table(setting=character(),cancer=character(),outcome=character(),n=integer(),mirna_rho=numeric(),mirna_partial_purity=numeric()),'mirna_expression_redundancy.tsv')
  save(if(length(sample_rows)) rbindlist(sample_rows,fill=TRUE) else data.table(setting=character(),cancer=character(),outcome=character(),sample_id=character()),'model_sample_ids.tsv')
  # Shared axes, no data-driven thresholds; editable plotting script uses these TSVs.
  source(file.path(code_root,'plot_mir146.R'),local=TRUE)
  plot_mir146(out)
  report<-c('# miR-146a/b complementary association analysis','',paste('Discovery run:',run_id),
    paste('Fixed shared primary genes:',paste(shared,collapse=', ')),
    'Exploratory analysis added after review; not part of original prespecification. No change to discovery set.',
    'CPE is primary; LUMP is sensitivity. Eligible cancers are read from the final run, not selected by these outcomes.',
    'Within each outcome/cancer/setting, both miRNAs, the outcome and purity must be finite. Identical complete cases are used for both conditional effects and reduced-model comparisons.',
    'Variables are winsorised within that complete-case set at 0.1/99.9 percentiles (R quantile type 7), average-ranked, then standardised. This is a new rank-based conditional sensitivity, not an exact reproduction of the discovery residualisation.',
    'Full model: standardised rank(outcome) ~ standardised rank(miR-146a-5p) + standardised rank(miR-146b-5p) + standardised rank(purity).',
    'partial_rho is the correlation of outcome and focal-miRNA rank residuals after regressing both on the other miRNA and purity. It is distinct from beta_rank.',
    'Inference: two-sided HC3 Wald z tests and pointwise normal-reference 95% intervals for beta_rank. These are approximate large-sample tests; no confidence interval for partial_rho is claimed.',
    'Delta R2 compares the full model with the identical-sample reduced model omitting the focal miRNA. It is apparent in-sample fit improvement, not external predictive value.',
    'BH and BY are separate by purity setting and outcome class (CYT vs shared genes), pooling BOTH miRNAs and ALL eligible cancer/outcome tests, with non-estimable planned tests retained in the denominator.',
    'VIF >10 is flagged, not silently removed. Rank deficiency, constant purity or variables, leverage >=0.99, fewer than the configured minimum samples, or degenerate fits are non-estimable.',
    'Conditional associations do not establish biological independence, causal regulation, direct targeting, or independence from immune-cell composition. Selected shared primary genes are not an unbiased transcriptome-wide validation sample.',
    'Absence of significance does not prove redundancy; multicollinearity and power must be considered.',
    '',capture.output(print(summary)))
  writeLines(report,file.path(out,'ANALYSIS_REPORT.md'))
  writeLines(capture.output(sessionInfo()),file.path(out,'logs/sessionInfo.txt'))
  after<-vapply(names(before),function(p) digest(file=p,algo='sha256'),character(1))
  if(!identical(before,after)) stop('Source inputs changed during analysis')
  files<-list.files(out,full.names=TRUE,recursive=TRUE);files<-files[basename(files)!='STATUS.txt']
  fwrite(data.table(file=substring(files,nchar(out)+2),sha256=vapply(files,function(p) digest(file=p,algo='sha256'),character(1))),
    file.path(out,'output_sha256.tsv'),sep='\t')
  writeLines('COMPLETE',file.path(out,'STATUS.txt'));complete<-TRUE
  message('COMPLETE: ',out);invisible(out)
}
