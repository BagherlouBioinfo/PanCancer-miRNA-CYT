figure_meta_sensitivity <- function() {
  d<-down("meta_analysis_all_arms",c("pooled_rho","I2","arm","adjustment","pi_lower","pi_upper"))
  d<-d[d$status=="ESTIMATED",,drop=FALSE]
  d$setting<-ifelse(d$adjustment=="immune",paste("Immune:",d$population),d$adjustment)
  source_table(d,"Figure_S1_meta")
  A<-ggplot2::ggplot(d,ggplot2::aes(setting,I2,fill=arm))+ggplot2::geom_boxplot(outlier.size=.3)+
    ggplot2::labs(x=NULL,y="I-squared (%)",fill="Association")+
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=60,hjust=1,size=6))
  cp<-d[d$adjustment=="CPE",];lp<-d[d$adjustment=="LUMP",]
  unique_key(cp,c("mirna_id","gene","arm"),"CPE meta")
  unique_key(lp,c("mirna_id","gene","arm"),"LUMP meta")
  m<-merge(cp,lp,by=c("mirna_id","gene","arm"),suffixes=c("_CPE","_LUMP"))
  source_table(m,"Figure_S1_concordance")
  B<-ggplot2::ggplot(m,ggplot2::aes(pooled_rho_CPE,pooled_rho_LUMP,colour=arm))+
    ggplot2::geom_abline(slope=1,intercept=0,colour="grey70")+ggplot2::geom_point(size=1.2)+
    ggplot2::labs(x="CPE pooled rho",y="LUMP pooled rho",colour="Association")
  save_panels("Figure_S1",list(A,B),paste(
    "Meta-analysis sensitivity. A, heterogeneity distributions across arms and adjustment settings.",
    "B, CPE and LUMP pooled estimates matched by pair and arm; eligible cancer sets can differ.",
    "Repeated shared association arms are retained as reported and are not independent replications."))
}
figure_mediation <- function() {
  d<-down("mediation_all_contexts",c("point_indirect","point_total","prop_mediated",
    "prop_mediated_meaningful","bh_fdr","primary_support","cancer"))
  d$subset<-ifelse(truth(d$primary_support),"Primary-support contexts","Other eligible contexts")
  d$significant<-factor(ifelse(!is.na(d$bh_fdr)&d$bh_fdr<CONFIG$fdr,
    "BH-FDR < 0.05","BH-FDR >= 0.05"),levels=c("BH-FDR >= 0.05","BH-FDR < 0.05"))
  source_table(d,"Figure_S2_mediation")
  A<-ggplot2::ggplot(d,ggplot2::aes(point_indirect,fill=significant))+
    ggplot2::geom_histogram(bins=40,colour="white")+ggplot2::facet_grid(subset~.)+
    ggplot2::scale_fill_manual(values=c("#BDBDBD","#2166AC"))+
    ggplot2::guides(fill=ggplot2::guide_legend(nrow=1))+
    ggplot2::labs(x="Estimated indirect effect",y="Contexts",fill=NULL)
  m<-d[truth(d$prop_mediated_meaningful)&is.finite(d$prop_mediated),,drop=FALSE]
  B<-ggplot2::ggplot(m,ggplot2::aes(point_total,prop_mediated,colour=subset))+
    ggplot2::geom_hline(yintercept=1,linetype="dotted",colour="grey50")+
    ggplot2::geom_point(size=.9,alpha=.6)+
    ggplot2::guides(colour=ggplot2::guide_legend(ncol=1))+
    ggplot2::labs(x="Estimated total effect",y="Indirect-to-total effect ratio (ab / c)",colour=NULL)
  save_panels("Figure_S2",list(A,B),paste(
    "Exploratory mediation across all eligible contexts, separating primary-support contexts.",
    "A, indirect-effect distributions using the original BH family across all contexts.",
    "B, indirect-to-total effect ratios for contexts satisfying the stored sign-based reporting criterion.",
    "Ratios above one are retained and are not fractions or percentages mediated; the dotted line marks one.",
    "Cross-sectional data, preselection and potential composition confounding prevent causal mediation claims."))
}
figure_enrichment <- function() {
  robust<-down("enrichment_terms_robust_to_leave_one_out",c("term_id","family"))
  plots<-list()
  for(fam in c("GO_BP","Reactome")) {
    file<-if(fam=="GO_BP") "enrichment_go_bp_results" else "enrichment_reactome_results"
    d<-down(file,c("term_id","term_name","fold_enrichment","fdr","n_selected_in_term"))
    d$robust<-d$term_id %in% robust$term_id[robust$family==fam]
    d<-d[!is.na(d$fdr)&d$fdr<CONFIG$fdr,,drop=FALSE]
    d<-head(d[order(d$fdr,d$term_id),,drop=FALSE],CONFIG$enrichment_top_n)
    d$label<-factor(wrap_text(d$term_name,48),levels=rev(wrap_text(d$term_name,48)))
    source_table(d,paste0("Figure_S3_",fam))
    plots[[length(plots)+1]]<-ggplot2::ggplot(d,ggplot2::aes(fold_enrichment,label,
      colour=-log10(pmax(fdr,.Machine$double.xmin)),size=n_selected_in_term,shape=robust))+
      ggplot2::geom_point()+ggplot2::scale_colour_gradient(low="#9ECAE1",high="#08306B")+
      ggplot2::scale_size_continuous(breaks=c(5,10,15,20,25),range=c(1.5,6))+
      ggplot2::scale_shape_manual(values=c("FALSE"=17,"TRUE"=16),
        labels=c("FALSE"="Not stable","TRUE"="Stable"),drop=FALSE)+
      ggplot2::guides(colour=ggplot2::guide_colourbar(barwidth=grid::unit(40,"mm"),order=1),
        size=ggplot2::guide_legend(nrow=1,order=2),shape=ggplot2::guide_legend(nrow=1,order=3))+
      ggplot2::labs(x="Fold enrichment",y=NULL,
        colour="-log10(BH-FDR)",size="Selected genes",shape="Leave-one-gene-out")+
      ggplot2::theme(axis.text.y=ggplot2::element_text(size=7))
  }
  save_panels("Figure_S3",plots,paste(
    "Functional annotation of the selected genes. A, GO biological processes; B, Reactome pathways.",
    "Terms are ordered by their existing FDR and truncated at the configured display count.",
    "Circles mark stability under every leave-one-gene-out removal; triangles mark loss under at least one removal.",
    "Enrichment is conditional on the study-specific eligible gene background."))
}
figure_clustering <- function() {
  v<-down("cluster_validity")
  k<-pick(v,c("k","n_clusters"),"cluster validity")
  sil<-pick(v,c("avg_silhouette","mean_silhouette","silhouette"),"cluster silhouette")
  v$plot_k<-v[[k]];v$plot_silhouette<-v[[sil]];source_table(v,"Figure_S4_validity")
  A<-ggplot2::ggplot(v,ggplot2::aes(plot_k,plot_silhouette))+
    ggplot2::geom_line()+ggplot2::geom_point()+
    ggplot2::labs(x="Number of clusters",y="Mean silhouette")
  d<-down("cluster_pca_coordinates")
  pc1<-pick(d,c("PC1","pc1","Dim.1"),"PCA axis 1")
  pc2<-pick(d,c("PC2","pc2","Dim.2"),"PCA axis 2")
  if(!any(c("cluster","cluster_id","cluster_label") %in% names(d))) {
    mem<-down("cluster_membership")
    id<-intersect(c("sample_id","sample","Sample"),intersect(names(d),names(mem)))
    if(!length(id)) stop("PCA/membership sample join key not documented in these tables")
    id<-id[1];unique_key(d,id,"PCA");unique_key(mem,id,"cluster membership")
    old_n<-nrow(d);d<-merge(d,mem,by=id,all.x=TRUE)
    if(nrow(d)!=old_n) stop("PCA join changed row count")
  }
  cl<-pick(d,c("cluster","cluster_id","cluster_label"),"cluster assignment")
  d$axis1<-d[[pc1]];d$axis2<-d[[pc2]];d$plot_cluster<-factor(d[[cl]])
  if(anyNA(d$plot_cluster)) stop("Missing cluster membership after join")
  # Source data stay local; TCGA sample identifiers must not be posted publicly without review.
  source_table(d,"Figure_S4_PCA")
  B<-ggplot2::ggplot(d,ggplot2::aes(axis1,axis2,colour=plot_cluster))+
    ggplot2::geom_point(size=.35,alpha=.3)+ggplot2::scale_colour_manual(values=CONFIG$clusters)+
    ggplot2::labs(x="PC1 score",y="PC2 score",colour="Cluster")
  save_panels("Figure_S4",list(A,B),paste(
    "Internal clustering characterization. A, observed silhouette across candidate cluster numbers.",
    "B, PCA scores from the existing run, coloured by its recorded cluster membership.",
    "Features were selected through CYT-associated discovery; separation does not constitute independent validation.",
    "Permutation-null convergence diagnostics require separate review; no new null runs are performed here."))
}
figure_specificity <- function() {
  d<-audit("matched_null_target_specificity_results",c("mirna_id","gene","nominal_recurrence",
    "null_recurrence_mean","null_recurrence_min","null_recurrence_max","fdr"))
  d<-d[order(d$mirna_id,d$gene),];d$pair<-factor(pair_label(d),levels=rev(pair_label(d)))
  source_table(d,"Figure_S5_matched_controls")
  A<-ggplot2::ggplot(d,ggplot2::aes(null_recurrence_mean,pair))+
    ggplot2::geom_segment(ggplot2::aes(x=null_recurrence_min,xend=null_recurrence_max,yend=pair),colour="grey70")+
    ggplot2::geom_point(colour="grey40",size=1)+
    ggplot2::geom_point(ggplot2::aes(x=nominal_recurrence),colour="#2166AC",shape=4,size=2)+
    ggplot2::labs(x="Nominal recurrence (cross: observed; dot/range: matched controls)",y=NULL)+
    ggplot2::theme(axis.text.y=ggplot2::element_text(size=6))
  save_panels("Figure_S5",list(A),paste(
    "Observed nominal recurrence and matched-control distributions for every primary pair.",
    "Crosses show observed nominal recurrence; grey points and ranges show control means and minima/maxima.",
    "Nominal recurrence differs from the primary FDR-qualified discovery rule.",
    "These controls test the specific matched-null comparison and are not random unselected genes."))
}
