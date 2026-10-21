figure_discovery <- function() {
  p <- primary(); counts <- list()
  for(t in c("030","025","020")) for(a in c("unadjusted","CPE","LUMP")) {
    folder <- if(t=="030" && a!="LUMP") "primary" else "sensitivity"
    d <- read_run(sprintf("tables/%s/triplets_rho%s_%s.tsv",folder,t,a),"qualifies")
    counts[[length(counts)+1]] <- data.frame(threshold=as.numeric(t)/100,
      adjustment=a,n=sum(truth(d$qualifies)))
  }
  counts <- source_table(do.call(rbind,counts),"Figure_2A_counts")
  A <- ggplot2::ggplot(counts,ggplot2::aes(factor(threshold),n,fill=adjustment))+
    ggplot2::geom_col(position="dodge")+ggplot2::scale_y_log10()+
    ggplot2::scale_fill_manual(values=CONFIG$colours)+
    ggplot2::labs(x="Absolute correlation threshold",y="Qualifying pairs (log scale)",fill="Adjustment")
  cancers <- read_run("processed/stage05/eligible_cancers_cpe.tsv","cancer")$cancer
  cells <- expand.grid(pair=pair_label(p),cancer=sort(cancers),stringsAsFactors=FALSE)
  lookup <- setNames(strsplit(ifelse(is.na(p$supporting_cancers),"",p$supporting_cancers),";"),pair_label(p))
  cells$supported <- mapply(function(a,b) b %in% trimws(lookup[[a]]),cells$pair,cells$cancer)
  cells$pair <- factor(cells$pair,levels=rev(pair_label(p)))
  source_table(cells,"Figure_2B_recurrence")
  B <- ggplot2::ggplot(cells,ggplot2::aes(cancer,pair,fill=supported))+
    ggplot2::geom_tile(colour="white",linewidth=.15)+
    ggplot2::scale_fill_manual(values=c("FALSE"="#F0F0F0","TRUE"=unname(CONFIG$colours["CPE"])),
      labels=c("Unsupported","Supported"))+
    ggplot2::labs(x=NULL,y=NULL,fill=NULL)+ggplot2::guides(fill=ggplot2::guide_legend(nrow=2))+
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,size=7),
      axis.text.y=ggplot2::element_text(size=6))
  m <- down("meta_analysis_all_arms",c("mirna_id","gene","arm","adjustment","pooled_rho","ci_lower","ci_upper"))
  m <- m[m$adjustment=="CPE" & m$status=="ESTIMATED",,drop=FALSE]
  unique_key(m,c("mirna_id","gene","arm"),"CPE meta-analysis")
  m$pair <- factor(pair_label(m),levels=rev(pair_label(p)));source_table(m,"Figure_2C_meta")
  C <- ggplot2::ggplot(m,ggplot2::aes(pooled_rho,pair,colour=arm))+
    ggplot2::geom_vline(xintercept=0,colour="grey70")+
    ggplot2::geom_segment(ggplot2::aes(x=ci_lower,xend=ci_upper,yend=pair),linewidth=.35)+
    ggplot2::geom_point(size=.9)+ggplot2::labs(x="Pooled Spearman rho (95% CI)",y=NULL,colour="Association")+
    ggplot2::guides(colour=ggplot2::guide_legend(nrow=3))+
    ggplot2::theme(axis.text.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank())
  save_panels("Figure_2",list(A,B,C),paste(
    "Discovery and recurrent associations. A, qualifying pair counts by threshold and adjustment;",
    "analysis eligibility differs between adjustment settings, so count differences do not isolate purity effects.",
    "B, cancer-specific three-arm support for all primary pairs, ordered by recurrence and identifier.",
    "C, CPE-adjusted pooled correlations and 95% confidence intervals for the same pairs; colours identify arms.",
    "Pooling uses all eligible cancers, not only discovery-supporting cancers. Shared arms are not independent observations."))
}
figure_composition <- function() {
  p <- primary()
  d <- down("immune_adjusted_recurrence_by_population",c("mirna_id","gene","adjustment_population",
    "recurrence","n_evaluable","is_biological_overlap_sensitivity"))
  unique_key(d,c("mirna_id","gene","adjustment_population"),"immune recurrence")
  grid<-merge(p[c("mirna_id","gene")],
    data.frame(adjustment_population=sort(unique(d$adjustment_population))),by=NULL)
  d<-merge(grid,d,by=c("mirna_id","gene","adjustment_population"),all.x=TRUE,sort=FALSE)
  d$pair <- factor(pair_label(d),levels=rev(pair_label(p)))
  d$display <- ifelse(d$n_evaluable>0,d$recurrence,NA_real_)
  d$category <- ifelse(d$adjustment_population %in% c("T cells","CD8 T cells","Cytotoxic lymphocytes","NK cells","Natural killer cells"),
    "CYT-overlapping populations","Other populations")
  source_table(d,"Figure_3A_immune")
  A <- ggplot2::ggplot(d,ggplot2::aes(adjustment_population,pair,fill=display))+
    ggplot2::geom_tile(colour="white",linewidth=.2)+
    ggplot2::geom_text(ggplot2::aes(label=ifelse(is.na(display),"NA",display)),size=2)+
    ggplot2::scale_fill_gradient(low="white",high="#2166AC",na.value="grey70",
      limits=c(0,20),breaks=c(0,5,10,15,20))+
    ggplot2::facet_grid(.~category,scales="free_x",space="free_x")+
    ggplot2::labs(x=NULL,y=NULL,fill="Supporting cancers")+
    ggplot2::guides(fill=ggplot2::guide_colourbar(direction="vertical",barheight=grid::unit(45,"mm")))+
    ggplot2::theme(legend.position="right",axis.text.x=ggplot2::element_text(angle=55,hjust=1,size=7),
      axis.text.y=ggplot2::element_text(size=6))
  n <- audit("matched_null_target_specificity_results",c("mirna_id","gene","empirical_p","fdr"))
  source_table(n,"Figure_3B_matched_null")
  B <- ggplot2::ggplot(n,ggplot2::aes(empirical_p))+
    ggplot2::geom_histogram(binwidth=.1,boundary=0,fill="#636363",colour="white")+
    ggplot2::labs(x="Matched-null empirical P value",y="Primary pairs")
  save_panels("Figure_3",list(A,B),paste(
    "Cell-composition sensitivity and matched-null assessment. A, three-arm recurrence following separate adjustment",
    "for each MCPcounter population. NA indicates no evaluable context, including marker-overlap exclusions, not zero support.",
    "CYT-overlapping populations are displayed separately. The adjustments are individual, not simultaneous.",
    "B, empirical P-value distribution for the matched-null test; multiplicity conclusions use the supplied BH-adjusted values.",
    "Failure to outperform these matched controls does not demonstrate absence of biological regulation."))
}
figure_survival <- function() {
  all <- do.call(rbind,lapply(c("OS","DSS","PFI"),function(e) {
    nm <- if(e=="OS") "survival_models" else paste0("survival_models_",e)
    d <- down(nm,c("feature","cancer","tier","status","HR","fdr","ph_violated"));d$endpoint<-e;d
  }))
  d <- all[all$tier==CONFIG$survival_tier & all$status=="ESTIMATED",,drop=FALSE]
  sig <- d[!is.na(d$fdr)&d$fdr<CONFIG$fdr,,drop=FALSE]
  sig$PH <- ifelse(truth(sig$ph_violated),"PH flagged","No PH flag")
  ct <- as.data.frame(table(sig$endpoint,sig$cancer,sig$PH),stringsAsFactors=FALSE)
  names(ct)<-c("endpoint","cancer","PH","n");source_table(ct,"Figure_4A_counts")
  A <- ggplot2::ggplot(ct,ggplot2::aes(cancer,n,fill=PH))+ggplot2::geom_col()+
    ggplot2::facet_grid(endpoint~.)+ggplot2::scale_fill_manual(values=c("#2166AC","#D95F02"))+
    ggplot2::labs(x=NULL,y="Associations with BH-FDR < 0.05",fill=NULL)+
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1))
  f <- d[d$endpoint=="OS" & !is.na(d$fdr)& d$fdr<CONFIG$fdr & !truth(d$ph_violated),,drop=FALSE]
  if(!is.null(CONFIG$survival_features)) f<-f[f$feature %in% CONFIG$survival_features,,drop=FALSE]
  f<-head(f[order(f$fdr,f$feature,f$cancer),,drop=FALSE],CONFIG$survival_max_rows)
  if(!nrow(f)) stop("No survival models meet the configured display criteria")
  f$label<-factor(paste(f$feature,f$cancer,sep=" / "),levels=rev(paste(f$feature,f$cancer,sep=" / ")))
  source_table(f,"Figure_4B_forest")
  B<-ggplot2::ggplot(f,ggplot2::aes(HR,label))+ggplot2::geom_vline(xintercept=1,colour="grey60")+
    ggplot2::geom_segment(ggplot2::aes(x=ci_lower,xend=ci_upper,yend=label))+
    ggplot2::geom_point(size=1.5)+ggplot2::scale_x_log10()+
    ggplot2::labs(x="Hazard ratio per feature SD (95% CI)",y=NULL)+
    ggplot2::theme(axis.text.y=ggplot2::element_text(size=6))
  inc<-down("cyt_incremental_models",c("feature","cancer","status","lrt_fdr","concordance_delta"))
  inc<-inc[inc$status=="ESTIMATED",,drop=FALSE]
  inc$significant<-factor(ifelse(!is.na(inc$lrt_fdr)&inc$lrt_fdr<CONFIG$fdr,
    "LRT BH-FDR < 0.05","LRT BH-FDR >= 0.05"),levels=c("LRT BH-FDR >= 0.05","LRT BH-FDR < 0.05"))
  source_table(inc,"Figure_4C_incremental")
  C<-ggplot2::ggplot(inc,ggplot2::aes(concordance_delta,fill=significant))+
    ggplot2::geom_histogram(bins=35,position="stack",colour="white")+
    ggplot2::scale_fill_manual(values=c("#F8766D","#00BFC4"))+
    ggplot2::guides(fill=ggplot2::guide_legend(ncol=1))+
    ggplot2::labs(x="In-sample change in Harrell's C",y="Estimable comparisons",fill=NULL)
  save_panels("Figure_4",list(A,B,C),paste(
    "Cancer-specific survival associations. A, significant associations in",CONFIG$survival_tier,
    "by endpoint, with global proportional-hazards flags shown separately.",
    "B, OS models with FDR < 0.05 and no global PH flag, displayed in ascending FDR order, up to",
    CONFIG$survival_max_rows,"rows; this is a descriptive display selection, not a new testing family.",
    "C, in-sample concordance changes for all estimable Tier 3 comparisons against CYT plus age and sex.",
    "Significance uses the original LRT correction. No independent prognostic validation is implied.",
    "OS, overall survival; DSS, disease-specific survival; PFI, progression-free interval."))
}
figure_cell_context <- function() {
  f<-down("stage24_primary_evidence_funnel",c("n_rows","n_triplets"))
  if(nrow(f)!=5) stop("Expected five sequential primary funnel steps")
  f$display_step<-factor(c("Conjunction BH","LUMP agreement","Nonoverlap","Collinearity check","Permutation"),
    levels=c("Conjunction BH","LUMP agreement","Nonoverlap","Collinearity check","Permutation"))
  if(any(diff(f$n_rows)>0)) stop("Primary funnel is not sequential")
  source_table(f,"Figure_5A_funnel")
  A<-ggplot2::ggplot(f,ggplot2::aes(display_step,n_rows,group=1))+
    ggplot2::geom_line(colour="#2166AC")+ggplot2::geom_point(size=2)+
    ggplot2::geom_text(ggplot2::aes(label=n_rows),vjust=-.7,size=3)+
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(.03,.15)))+
    ggplot2::labs(x=NULL,y="Triplet-cancer-population combinations")
  d<-down("stage24_final_robust_nonoverlapping",c("mirna_id","gene","cancer","population",
    "conjunction_FDR_BH","conjunction_FDR_BY","nonlinearity_flagged","collinearity_assessment"))
  unique_key(d,c("mirna_id","gene","cancer","population"),"final cellular-context table")
  if(nrow(d)!=CONFIG$expected_final_contexts || length(unique(key(d)))!=CONFIG$expected_final_pairs)
    stop("Final cellular-context counts differ from configured expectations")
  d$pair<-pair_label(d);d$context<-paste(d$cancer,d$population,sep=" / ")
  d$evidence<- -log10(pmax(d$conjunction_FDR_BH,.Machine$double.xmin))
  d$sensitivity<-ifelse(truth(d$nonlinearity_flagged),"Flagged","Not flagged")
  d$BY<-factor(ifelse(d$conjunction_FDR_BY<CONFIG$fdr,"BY-FDR < 0.05","BY-FDR >= 0.05"),
    levels=c("BY-FDR >= 0.05","BY-FDR < 0.05"))
  source_table(d,"Figure_5B_final_contexts")
  B<-ggplot2::ggplot(d,ggplot2::aes(context,pair,size=evidence,shape=sensitivity,
    colour=BY))+ggplot2::geom_point()+
    ggplot2::scale_size_continuous(range=c(1.5,5))+
    ggplot2::scale_colour_manual(values=c("#969696","#2166AC"),na.value="black")+
    ggplot2::guides(size=ggplot2::guide_legend(nrow=1,order=3),
      colour=ggplot2::guide_legend(nrow=1,order=1),shape=ggplot2::guide_legend(nrow=1,order=2))+
    ggplot2::labs(x=NULL,y=NULL,size="-log10(BH-FDR)",colour="BY sensitivity",shape="Nonlinearity")+
    ggplot2::theme(legend.position="right",axis.text.x=ggplot2::element_text(angle=55,hjust=1,size=7))
  save_panels("Figure_5",list(A,B),paste(
    "Cancer-specific cellular-context associations. A, the primary sequential funnel, beginning with all",
    "conjunction-BH-supported combinations; it is distinct from the historical reconciliation.",
    "B, final combinations; dot size represents conjunction BH-FDR, colour BY sensitivity and shape exploratory nonlinearity.",
    "Blank cells are not final-supported combinations and do not establish absence of an effect.",
    "Collinearity checks with no correlated population are not evidence of demonstrated stability.",
    "These findings do not establish cell of origin or causal regulation."))
}
