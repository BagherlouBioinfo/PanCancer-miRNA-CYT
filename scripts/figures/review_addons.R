review_section_labels <- c(
  Cluster_Summary="Summary of frozen-cluster incremental survival analyses",
  Cluster_Models="Cancer-specific incremental survival models for frozen clusters",
  Cluster_PH="Proportional-hazards diagnostics for frozen-cluster models",
  Cluster_Time_HR="Time-varying overall-survival cluster effects",
  Cluster_SampleFlow="Complete-case sample flow for cluster survival analyses",
  Mir146_Summary="Summary of conditional miR-146a/b associations",
  Mir146_Conditional="Conditional associations of miR-146a-5p and miR-146b-5p",
  Mir146_Both="Contexts supporting both miRNAs after mutual adjustment",
  Mir146_Expression="Expression similarity of miR-146a-5p and miR-146b-5p")
# Publication-only integration of stored post-review analyses; no inferential refits.
review_bundle <- function() {
  if(exists('REVIEW_DATA',envir=environment(review_bundle),inherits=FALSE))
    return(get('REVIEW_DATA',envir=environment(review_bundle),inherits=FALSE))
  base<-file.path(project_root,'reference/reported_results/addon_data')
  manifest<-data.table::fread(file.path(base,'INPUT_SHA256.tsv'),data.table=FALSE)
  for(i in seq_len(nrow(manifest))) {
    path<-file.path(base,manifest$file[i])
    if(!file.exists(path)||digest::digest(file=path,algo='sha256')!=manifest$sha256[i]) stop('Bundled add-on input hash mismatch: ',manifest$file[i])
  }
  cluster_names<-c('cluster_incremental_survival','PH_diagnostics','OS_time_varying_cluster_HR','sample_flow','summary')
  mir_names<-c('conditional_mir146_associations','both_mirnas_support','mirna_expression_redundancy','summary')
  result<-list();records<-list()
  for(kind in c('cluster','mir146')) {
    custom<-if(kind=='cluster') CONFIG$cluster_survival_root else CONFIG$mir146_root
    root<-if(is.null(custom)) file.path(base,kind) else normalizePath(custom,winslash='/',mustWork=TRUE)
    if(!is.null(custom)) {
      if(!file.exists(file.path(root,'STATUS.txt')) || !identical(trimws(readLines(file.path(root,'STATUS.txt'),warn=FALSE)),'COMPLETE')) stop('Add-on is not COMPLETE: ',root)
      report<-readLines(file.path(root,'ANALYSIS_REPORT.md'),warn=FALSE)
      if(!any(grepl(CONFIG$run_id,report,fixed=TRUE))) stop('Add-on discovery run mismatch')
      root<-if(dir.exists(file.path(root,'tables'))) file.path(root,'tables') else root
    }
    nms<-if(kind=='cluster') cluster_names else mir_names
    result[[kind]]<-list()
    for(nm in nms) {
      path<-file.path(root,paste0(nm,'.tsv'))
      if(!is.null(custom)) verify_table_manifest(path, custom)
      d<-data.table::fread(path,na.strings=c('NA','NaN',''),data.table=FALSE,check.names=FALSE)
      result[[kind]][[nm]]<-d
      records[[paste(kind,nm)]]<-data.frame(path=paste0('review_addons/',kind,'/',nm,'.tsv'),rows=nrow(d),columns=ncol(d))
      source_table(d,paste0('Original_',kind,'_',nm))
    }
  }
  sources<<-c(sources,records)
  data.table::fwrite(manifest,file.path(OUT,'logs','ADDON_BUNDLE_SHA256.tsv'),sep='\t')
  assign('REVIEW_DATA',result,envir=environment(review_bundle));result
}
review_cluster_models <- function() {
  d<-review_bundle()$cluster$cluster_incremental_survival
  # Original add-on initialised excluded stage rows with zero events and source cohort N.
  # Correct presentation only; original raw table stays in source_data.
  skipped<-grepl('not stage-adequate',ifelse(is.na(d$reason),'',d$reason),fixed=TRUE)
  d$n_source_if_stage_ineligible<-ifelse(skipped,d$n,NA_integer_)
  d$n[skipped]<-NA_integer_;d$events[skipped]<-NA_integer_
  corrections<-data.frame(correction='Stage-ineligible placeholder events=0 and source N presented as NA model counts',rows=sum(skipped),inferential_values_changed=FALSE)
  data.table::fwrite(corrections,file.path(OUT,'logs','REPORTING_CORRECTIONS.tsv'),sep='\t')
  d
}
validate_review_addons <- function() {
  z<-review_bundle();c<-z$cluster$cluster_incremental_survival;m<-z$mir146$conditional_mir146_associations
  need(c,c('endpoint','cancer','tier','status','lrt_p','lrt_fdr_BH','ph_global_p','tv_fdr_BH'),'Cluster add-on')
  need(m,c('setting','cancer','outcome','mirna_id','status','p_HC3','FDR_BH','FDR_BY','VIF','partial_rho'),'miR-146 add-on')
  unique_key(c,c('endpoint','cancer','tier'),'cluster comparisons')
  unique_key(m,c('setting','cancer','outcome','mirna_id'),'conditional miR-146 hypotheses')
  stopifnot(nrow(c)==186L,nrow(m)==390L,all(m$status=='ESTIMATED'))
  stopifnot(setequal(unique(m$outcome),c('CYT','CD3D','CD80','LCP2','SPI1')))
  stopifnot(setequal(unique(m$mirna_id),c('hsa-miR-146a-5p','hsa-miR-146b-5p')))
  # Verify stored corrections, never overwrite them.
  for(indices in split(seq_len(nrow(c)),interaction(c$endpoint,c$tier,drop=TRUE))) {
    stopifnot(isTRUE(all.equal(p.adjust(c$lrt_p[indices],'BH',n=length(indices)),c$lrt_fdr_BH[indices],tolerance=1e-10)))
  }
  for(indices in split(seq_len(nrow(m)),interaction(m$setting,m$outcome_type,drop=TRUE))) for(method in c('BH','BY'))
    stopifnot(isTRUE(all.equal(p.adjust(m$p_HC3[indices],method,n=length(indices)),m[[paste0('FDR_',method)]][indices],tolerance=1e-10)))
  stopifnot(all(abs(m$partial_R2-m$partial_rho^2)<1e-10),all(m$VIF<10))
  cm<-c[c$status=='ESTIMATED',,drop=FALSE]
  stopifnot(all(cm$n==cm$n_reference+cm$n_comparison),all(cm$events<=cm$n))
  for(ep in c('OS','DSS','PFI')) for(ti in c('CYT_age_sex','CYT_age_sex_stage')) {
    idx<-c$endpoint==ep&c$tier==ti
    s<-z$cluster$summary[z$cluster$summary$endpoint==ep&z$cluster$summary$tier==ti,,drop=FALSE]
    stopifnot(nrow(s)==1,sum(c$status[idx]=='ESTIMATED')==s$estimated,
      sum(c$lrt_fdr_BH[idx]<.05,na.rm=TRUE)==s$BH_significant)
  }
  message('Stored cluster-survival and conditional miR-146 consistency checks passed.')
  invisible(TRUE)
}
review_table_sections <- function() {
  z<-review_bundle()
  list(Cluster_Summary=z$cluster$summary,Cluster_Models=review_cluster_models(),Cluster_PH=z$cluster$PH_diagnostics,
    Cluster_Time_HR=z$cluster$OS_time_varying_cluster_HR,Cluster_SampleFlow=z$cluster$sample_flow,
    Mir146_Summary=z$mir146$summary,Mir146_Conditional=z$mir146$conditional_mir146_associations,
    Mir146_Both=z$mir146$both_mirnas_support,Mir146_Expression=z$mir146$mirna_expression_redundancy)
}
review_tier <- function(x) factor(x,levels=c('CYT_age_sex','CYT_age_sex_stage'),
  labels=c('CYT + age + sex','CYT + age + sex + tumour stage'))
figure_cluster_incremental <- function() {
  library(ggplot2)
  z<-review_bundle()$cluster;d<-review_cluster_models()
  d<-d[d$endpoint=='OS'&d$status=='ESTIMATED',,drop=FALSE]
  d$cancer<-factor(d$cancer,levels=rev(sort(unique(d$cancer))))
  d$adjustment<-review_tier(d$tier)
  d$PH<-ifelse(is.na(d$ph_global_p),'Unavailable',ifelse(d$ph_global_p<.05,'Flagged','Not flagged'))
  source_table(d,'Figure_S7_OS_models')
  A<-ggplot(d,aes(x=HR,y=cancer))+geom_vline(xintercept=1,linetype=2,colour='grey60')+
    geom_segment(aes(x=lower,xend=upper,yend=cancer),linewidth=.35)+geom_point(aes(shape=PH),size=1.8)+
    scale_x_log10()+facet_grid(adjustment~.,scales='free_y',space='free_y')+
    labs(x='Cluster 2 versus cluster 1: hazard ratio',y=NULL,shape='Global PH diagnostic')
  B<-ggplot(d,aes(x=c_delta,y=cancer))+geom_vline(xintercept=0,linetype=2,colour='grey60')+
    geom_point(aes(shape=PH),size=1.8)+facet_grid(adjustment~.,scales='free_y',space='free_y')+
    labs(x='Change in apparent concordance',y=NULL,shape='Global PH diagnostic')
  t<-z$OS_time_varying_cluster_HR
  t<-t[t$cancer%in%CONFIG$cluster_time_example_cancers & truth(t$supported),,drop=FALSE]
  t$adjustment<-review_tier(t$tier);t$years<-t$day/365
  source_table(t,'Figure_S7_time_varying_examples')
  C<-ggplot(t,aes(x=years,y=HR,colour=adjustment,group=adjustment))+
    geom_hline(yintercept=1,linetype=2,colour='grey60')+geom_line(linewidth=.45)+
    geom_errorbar(aes(ymin=lower,ymax=upper),width=.08,position=position_dodge(.12))+
    geom_point(size=1.7,position=position_dodge(.12))+scale_y_log10()+
    facet_wrap(~cancer,nrow=1)+scale_colour_manual(values=c('#2166AC','#B35806'))+scale_x_continuous(breaks=c(1,3,5))+
    labs(x='Years after baseline',y='Time-varying cluster hazard ratio',colour='Adjustment')
  save_panels('Figure_S7',list(A,B,C),paste(
    'Supplementary Figure S7. Internal survival associations of frozen molecular clusters beyond CYT.',
    '(A) Cancer-specific OS cluster hazard ratios and pointwise 95% confidence intervals from models adjusted for CYT, age and sex, with a separate stage-adjusted analysis where eligible. Global PH diagnostic flags are nominal and do not by themselves identify which covariate violates PH.',
    '(B) Change in apparent Harrell concordance after adding the frozen cluster on identical complete cases. These are internal fit differences, not externally validated prediction gains.',
    '(C) Illustrative OS time-varying sensitivities for LGG and LUAD, selected for display after reviewing the associations; all attempted models remain in Supplementary Table S5. All retained covariates receive main and log-time terms in both nested models; cluster contribution is tested jointly with two degrees of freedom. Displayed pointwise intervals are available only where at least ten patients remain at risk. Connecting segments guide the eye between displayed time points. LGG is not stage-adequate.',
    'Static-model HRs in PH-flagged models require caution. Stored BH corrections are separate by endpoint and adjustment tier, across all attempted cancers; time-varying OS tests form separate sensitivity families. DSS/PFI results and non-estimable comparisons are retained in Supplementary Table S5. Clusters were not reselected and this is not independent biological or external prognostic validation.'))
}
figure_mir146_conditional <- function() {
  library(ggplot2)
  d<-review_bundle()$mir146$conditional_mir146_associations
  source_table(d,'Figure_S8_conditional_miR146')
  panels<-list()
  for(st in c('CPE','LUMP')) for(mi in c('hsa-miR-146a-5p','hsa-miR-146b-5p')) {
    z<-d[d$setting==st&d$mirna_id==mi,,drop=FALSE]
    z$outcome<-factor(z$outcome,levels=rev(c('CYT','CD3D','CD80','LCP2','SPI1')))
    z$cancer<-factor(z$cancer,levels=sort(unique(z$cancer)))
    z$star<-ifelse(!is.na(z$FDR_BH)&z$FDR_BH<.05,'*','')
    panels[[length(panels)+1]]<-ggplot(z,aes(x=cancer,y=outcome,fill=partial_rho))+
      geom_tile(colour='white',linewidth=.2)+geom_text(aes(label=star),size=3.4)+
      scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',midpoint=0,limits=c(-1,1),
        breaks=c(-1,-.5,0,.5,1),na.value='grey85',name='Conditional rank correlation')+
      labs(x=NULL,y=NULL)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=7.5),axis.ticks=element_blank())
  }
  save_panels('Figure_S8',panels,paste(
    'Supplementary Figure S8. Complementary conditional associations of miR-146a-5p and miR-146b-5p.',
    '(A/B) CPE-adjusted results for miR-146a-5p/miR-146b-5p, respectively; (C/D) corresponding LUMP sensitivities.',
    'Each cell shows the correlation between residualised outcome and focal-miRNA ranks after controlling both for the other miRNA and for purity. Variables were winsorised at the within-analysis 0.1st and 99.9th percentiles, average-ranked and standardised. A common -1 to 1 colour scale is used.',
    'Asterisks identify stored BH FDR <0.05 for the HC3 Wald coefficient test, not a separate correlation test. Families pool both miRNAs and all eligible cancers separately for CYT and the four fixed shared primary genes, and separately for CPE and LUMP. Shared genes are CD3D, CD80, LCP2 and SPI1.',
    'Supplementary Table S8 retains coefficients, pointwise confidence intervals, partial correlations, VIF, apparent incremental R2, BH/BY values and all planned tests. Conditional association does not establish direct targeting, causal effects or biological independence, and the shared genes were selected in the same discovery cohort.'))
}
write_review_handoff <- function() {
  z<-review_bundle()
  map<-data.frame(asset=c('Figure_S7','Figure_S8','Supplementary_Table_S5','Supplementary_Table_S8'),
    citation=c('Supplementary Fig. S7','Supplementary Fig. S8','Supplementary Table S5','Supplementary Table S8'),
    placement=c('After cancer-specific survival and incremental-information results','After seed-family/miRNA-redundancy results',
      'Survival section, including frozen-cluster add-on','Target-specificity and miRNA-family section, including conditional miR-146 add-on'),
    contents=c('OS static HR, apparent concordance, time-varying examples','Conditional rank effects, CPE/LUMP, both miRNAs',
      'Full cluster models, PH diagnostics, time-varying HR, sample flow and summary',
      'Full conditional results, expression redundancy, simultaneous support and summary'))
  data.table::fwrite(map,file.path(OUT,'text','NEW_ANALYSIS_CITATION_MAP.tsv'),sep='\t')
  # Keep data ledger separate from later full manuscript editing.
  writeLines(c('# Post-review analysis integration map','',
    'Existing main figures 2-6, supplementary figures S1-S6, main tables 1-2 and supplementary tables S1-S11 retain their numbering.',
    'Insert frozen-cluster survival results after the existing feature-level survival/incremental analysis. Cite Supplementary Fig. S7 and Supplementary Table S5.',
    'Insert conditional miR-146 results after the seed-family redundancy analysis. Cite Supplementary Fig. S8 and Supplementary Table S8.',
    'Methods: add identical-complete-case nested cluster models, PH diagnostics and fixed OS time-varying sensitivity; separately add standardised-rank conditional miRNA regression with HC3 and its own BH/BY families.',
    'Discussion: retain context-dependent/internal nature of cluster survival results and distinguish complementary statistical miRNA information from biological independence.',
    'Reviewer response: cluster incremental survival addresses the subtype-beyond-CYT request; conditional miR-146 addresses complementary information beyond seed-family similarity.',
    '', '## Stored cluster-survival summary',capture.output(print(z$cluster$summary,row.names=FALSE)),
    '', '## Stored conditional-miRNA summary',capture.output(print(z$mir146$summary,row.names=FALSE)),
    '', 'Do not replace the original feature-level survival results with cluster-level results. Both are different analyses.',
    'For stage-ineligible cluster rows, apparent event count 0 and cohort N were initialisation placeholders; publication model counts are NA. Original rows are retained as Original_cluster_cluster_incremental_survival.tsv; inferential values are unchanged.'),
    file.path(OUT,'text','Post_review_integration_handoff.md'))
  captions<-c(
    'Table 1. Distribution of the fixed primary miRNA-gene pairs across miRNAs. Recurrence is the number of cancer types meeting all three within-cancer discovery criteria.',
    'Table 2. External validation of the fixed 47-pair discovery set in CPTAC breast cancer. Stored BH/BY support combines the conjunction criterion with absolute correlation of at least 0.30 in each arm. Direction concordance is reported against both TCGA-BRCA references. WES purity is distinct from CPE.',
    'Supplementary Table S1. Cohort composition, purity availability, analysis eligibility and clinical endpoint completeness by cancer type. Availability counts and analysis-specific eligible sets are distinct.',
    'Supplementary Table S2. Primary discovery pairs, threshold and purity sensitivities, and purity-conflict exclusion sensitivity. Qualifying lists are shown in the grouped table; full source listings, including nonqualifying records where available, remain in the technical archive. These analyses do not establish direct targeting.',
    'Supplementary Table S3. Immune-composition sensitivity of recurrent associations. Population-specific adjustment results must be interpreted with marker overlap and non-estimable comparisons; MCP-counter scores are not cell fractions.',
    'Supplementary Table S4. Cancer-specific association meta-analysis estimates and heterogeneity diagnostics. Shared association arms are not independent replications; pooled estimates do not imply uniform effects across cancers.',
    'Supplementary Table S5. Feature-level survival analyses and internal incremental survival analyses of the frozen clusters. Added sections retain every attempted cancer/endpoint/tier combination, PH diagnostics, OS time-varying sensitivities and sample-flow counts. Model-ineligible counts are NA. Apparent concordance is not external predictive validation.',
    'Supplementary Table S6. Exploratory within-cancer mediation models. Cross-sectional sampling and selection in the same cohort preclude independent causal validation. Indirect-to-total effect ratios exceeding one are retained as ratios, not fractions or percentages mediated.',
    'Supplementary Table S7. Clustering validity, stability, null calibration and internal molecular/immune characteristics. These diagnostics do not establish externally validated biological subtypes; the feature-wise permutation null removes coexpression structure.',
    'Supplementary Table S8. Matched-null target-specificity, seed-family redundancy and conditional miR-146a/b analyses. Added conditional sections retain both miRNAs and all tested cancers/outcomes, HC3 coefficient inference, BH/BY corrections, VIF and incremental in-sample R2. These conditional associations do not establish direct regulation or biological independence.',
    'Supplementary Table S9. Functional enrichment and drug-gene annotations, including identifier resolution and leave-one-gene-out sensitivities. Database-listed interactions, resolved genes with no interactions and unresolved genes are distinct record classes. Drug approval annotations do not establish efficacy for the cancers analysed.',
    'Supplementary Table S10. Corrected cancer-specific cellular-context conjunction analysis, sequential primary evidence funnel and diagnostic sensitivities. Results identify triplet-cancer-population contexts; they do not establish cell of origin or a causal mechanism.',
    'Supplementary Table S11. Complete external-validation results, unique-arm tests, coverage, mapping and testing-family definitions. Stored FDR estimates are retained. WES sensitivity uses winsorised rank residualisation and is not an identical replication of CPE adjustment. Shared arms are not independent repetitions.')
  writeLines(c('# Table legends','',as.vector(rbind(captions,rep('',length(captions))))),file.path(OUT,'text','Table_legends.md'))
  ledger<-file.path(OUT,'text','Results_counts_summary.md')
  if(file.exists(ledger)) {
    tx<-readLines(ledger,warn=FALSE)
    cluster_text<-paste('Frozen-cluster membership added OS model information beyond CYT, age and sex in LGG and LUAD after BH correction (2/19 estimable cancers); both models had global PH diagnostic flags.',
      'No static stage-adjusted OS comparison met BH FDR<0.05 (0/15). In the separate time-varying sensitivity, LGG and LUAD retained support without stage adjustment and LUAD retained support with stage adjustment; LGG was not stage-adequate',
      '(Supplementary Fig. S7; Supplementary Table S5). These are internal exploratory associations, not validated prognostic subtypes.')
    mirna_text<-paste('After controlling for the other miRNA and CPE, miR-146a-5p and miR-146b-5p retained CYT associations in 15/20 and 8/20 cancers, respectively, under the stored HC3/BH criterion; corresponding LUMP counts were 12/19 and 8/19.',
      'No conditional model had VIF>10. The fixed shared-gene conditional results are reported separately; these findings support complementary statistical information without establishing biological independence',
      '(Supplementary Fig. S8; Supplementary Table S8).')
    pos<-grep('^Across all four OS model tiers',tx)
    if(length(pos)==1L) tx<-append(tx,c('',cluster_text),after=pos)
    pos<-grep('^In the matched-null comparison',tx)
    if(length(pos)==1L) tx<-append(tx,c('',mirna_text),after=pos)
    tx<-c(tx,'',paste('In independent CPTAC breast tumours, the fixed 47-pair set was assessed without reselection in 101 samples, with WES-purity sensitivity in 96.',
      'All 47 pairs met the stored three-arm criteria in both settings; all-three-arm direction agreement with the CPE-adjusted TCGA-BRCA reference was 46/47.',
      'This is association replication in a single breast-cancer cohort, not validation of pan-cancer recurrence, clusters or survival prediction (Fig. 6; Table 2; Supplementary Fig. S6; Supplementary Table S11).'))
    writeLines(tx,ledger)
  }
  allfig<-paste0('Figure_',c(2:6,paste0('S',1:8)))
  index<-data.frame(asset=allfig,file=paste0('figures/',allfig,'.pdf'),citation=c(paste('Fig.',2:6),paste('Supplementary Fig.',paste0('S',1:8))))
  data.table::fwrite(index,file.path(OUT,'text','FIGURE_REFERENCE_INDEX.tsv'),sep='\t')
}
