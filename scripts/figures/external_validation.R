# Publication-only use of stored Stage 25 results. No correlations, bootstrap,
# cohort selection or inferential FDR values are refitted or replaced here.
external_files <- c("coverage_summary", "ensembl_symbol_mapping", "external_pair_arm_comparison",
  "external_pair_validation", "external_sample_annotation", "external_unique_arm_tests",
  "gene_mapping_audit", "mirna_mapping_audit", "testing_families", "validation_summary")
external_dir <- function() {
  if(is.null(CONFIG$external_root)) return(file.path(project_root,"reference/reported_results/external_data"))
  root<-normalizePath(CONFIG$external_root,winslash="/",mustWork=TRUE)
  if(dir.exists(file.path(root,"tables"))) file.path(root,"tables") else root
}
external_bundle <- function() {
  if(exists("EXTERNAL_DATA",envir=environment(external_bundle),inherits=FALSE))
    return(get("EXTERNAL_DATA",envir=environment(external_bundle),inherits=FALSE))
  root<-external_dir();result<-list();manifest<-list()
  for(nm in external_files) {
    path<-file.path(root,paste0(nm,".tsv"))
    if(!file.exists(path)) stop("Missing external input: ",path)
    d<-data.table::fread(path,na.strings=c("NA","NaN",""),check.names=FALSE,data.table=FALSE)
    result[[nm]]<-d
    sha<-digest::digest(file=path,algo="sha256")
    manifest[[nm]]<-data.frame(file=basename(path),path=normalizePath(path,winslash="/"),rows=nrow(d),sha256=sha)
    sources[[paste0("external/",nm)]]<<-data.frame(path=paste0("external/",basename(path)),rows=nrow(d),columns=ncol(d))
  }
  manifest<-do.call(rbind,manifest);rownames(manifest)<-NULL
  if(is.null(CONFIG$external_root)) {
    expected<-data.table::fread(file.path(root,"INPUT_SHA256.tsv"),data.table=FALSE)
    if(!setequal(expected$file,manifest$file) || any(expected$sha256[match(manifest$file,expected$file)]!=manifest$sha256))
      stop("Bundled external input checksum mismatch. Restore the original package.")
  }
  data.table::fwrite(manifest,file.path(OUT,"logs","EXTERNAL_INPUTS.tsv"),sep="\t")
  assign("EXTERNAL_DATA",result,envir=environment(external_bundle))
  result
}
external_unique_comparison <- function(adjustment) {
  d<-external_bundle()$external_pair_arm_comparison
  d<-d[d$adjustment==adjustment,,drop=FALSE]
  d$hypothesis<-ifelse(d$arm=="mirna_cyt",d$mirna_id,
    ifelse(d$arm=="gene_cyt",d$gene,paste(d$mirna_id,d$gene,sep=" / ")))
  cols<-c("arm","hypothesis","rho","TCGA_BRCA_CPE_rho","TCGA_BRCA_unadjusted_rho")
  d<-unique(d[c(cols,"mirna_id","gene")])
  # Same miRNA-CYT/gene-CYT test is repeated under several pair labels upstream.
  # Its numeric values must agree before it is represented by one plotted point.
  numbers<-unique(d[cols])
  unique_key(numbers,c("arm","hypothesis"),"External shared-arm numeric effects")
  numbers
}
validate_external_outputs <- function() {
  ext<-external_bundle();r<-ext$external_pair_validation;a<-ext$external_pair_arm_comparison
  u<-ext$external_unique_arm_tests;s<-ext$external_sample_annotation;p<-primary()
  checks<-list();add<-function(name,ok) checks[[length(checks)+1L]]<<-data.frame(check=name,status=if(isTRUE(ok))"PASS" else "FAIL")
  need(r,c("mirna_id","gene","adjustment","evaluable","p_conjunction","minimum_abs_rho","conjunction_fdr_BH","conjunction_fdr_BY","association_support_BH","all_directions_agree_CPE","all_directions_agree_unadjusted"),"external pair results")
  need(a,c("mirna_id","gene","arm","adjustment","rho","p_raw","fdr_BH","fdr_BY","TCGA_BRCA_CPE_rho","TCGA_BRCA_unadjusted_rho","direction_agrees_CPE","direction_agrees_unadjusted"),"external arm comparisons")
  need(u,c("arm","mirna_id","gene","adjustment","p_raw","fdr_BH","fdr_BY"),"external unique arms")
  need(s,c("sample_id","CYT","GZMA","PRF1","WES_purity"),"external cohort")
  unique_key(r,c("mirna_id","gene","adjustment"),"External pair results")
  unique_key(a,c("mirna_id","gene","arm","adjustment"),"External arm comparisons")
  unique_key(s,"sample_id","External samples")
  add("External settings are explicit",setequal(unique(r$adjustment),c("unadjusted","WES_purity")))
  add("External CYT matches stored gene expression",isTRUE(all.equal(s$CYT,(s$GZMA+s$PRF1)/2,tolerance=1e-10)))
  add("External WES purity is a fraction",all(is.na(s$WES_purity)|(s$WES_purity>=0&s$WES_purity<=1)))
  same_num<-function(x,y) isTRUE(all.equal(as.numeric(x),as.numeric(y),tolerance=1e-9,check.attributes=FALSE))
  ref<-down("meta_analysis_input_all_cancers",c("mirna_id","gene","arm","cancer","adjustment","rho","eligible"))
  arm_names<-c(mirna_cyt="miRNA-CYT",mirna_gene="miRNA-gene",gene_cyt="gene-CYT")
  for(adj in c("unadjusted","WES_purity")) {
    q<-r[r$adjustment==adj,,drop=FALSE];aa<-a[a$adjustment==adj,,drop=FALSE]
    add(paste(adj,"contains the fixed discovery set"),nrow(q)==nrow(p)&&setequal(key(q),key(p)))
    add(paste(adj,"three stored arms per pair"),nrow(aa)==3L*nrow(p)&&all(table(key(aa))==3L))
    pos<-match(key(q),key(p))
    add(paste(adj,"recurrence equals discovery source"),all(q$discovery_recurrence==p$strict_recurrence[pos]))
    computed_p<-vapply(key(q),function(k) {v<-aa$p_raw[key(aa)==k];if(anyNA(v))NA_real_ else max(v)},numeric(1))
    add(paste(adj,"conjunction is maximum raw arm p"),same_num(computed_p,q$p_conjunction))
    computed_effect<-vapply(key(q),function(k) {v<-aa$rho[key(aa)==k];if(anyNA(v))NA_real_ else min(abs(v))},numeric(1))
    add(paste(adj,"minimum absolute effect summary"),same_num(computed_effect,q$minimum_abs_rho))
    for(method in c("BH","BY")) {
      ok<-is.finite(q$p_conjunction);expected<-rep(NA_real_,nrow(q))
      expected[ok]<-stats::p.adjust(q$p_conjunction[ok],method=method,n=nrow(p))
      add(paste(adj,method,"pair FDR agrees with saved values"),same_num(expected,q[[paste0("conjunction_fdr_",method)]]))
    }
    expected_support<-ifelse(truth(q$evaluable),q$conjunction_fdr_BH<CONFIG$fdr&q$minimum_abs_rho>=CONFIG$primary_threshold,NA)
    add(paste(adj,"association support flag"),identical(as.logical(q$association_support_BH),expected_support))
    for(ra in c("CPE","unadjusted")) {
      ref_rows<-ref[ref$cancer=="BRCA"&ref$adjustment==ra,,drop=FALSE]
      ref_rows$arm<-names(arm_names)[match(ref_rows$arm,unname(arm_names))]
      unique_key(ref_rows,c("mirna_id","gene","arm"),paste("BRCA reference",ra))
      idx<-match(key(aa,c("mirna_id","gene","arm")),key(ref_rows,c("mirna_id","gene","arm")))
      expected<-ref_rows$rho[idx];expected[!truth(ref_rows$eligible[idx])]<-NA_real_
      add(paste(adj,"reference matches original run",ra),same_num(expected,aa[[paste0("TCGA_BRCA_",ra,"_rho")]]))
      expected_dir<-ifelse(is.finite(expected)&expected!=0&is.finite(aa$rho)&aa$rho!=0,sign(expected)==sign(aa$rho),NA)
      add(paste(adj,"stored direction agrees with coefficients",ra),identical(as.logical(aa[[paste0("direction_agrees_",ra)]]),as.logical(expected_dir)))
      pair_dir<-vapply(key(q),function(k) {
        v<-aa[[paste0("direction_agrees_",ra)]][key(aa)==k]
        if(anyNA(v)) NA else all(as.logical(v))
      },logical(1))
      add(paste(adj,"pair direction agreement",ra),identical(unname(pair_dir),as.logical(q[[paste0("all_directions_agree_",ra)]])))
      expected_rep<-expected_support & pair_dir
      add(paste(adj,"pair replication flag",ra),identical(unname(expected_rep),as.logical(q[[paste0("directional_replication_vs_",if(ra=="CPE")"CPE" else "unadjusted","_BRCA")]])))
    }
    for(arm in names(arm_names)) {
      b<-u[u$adjustment==adj&u$arm==arm,,drop=FALSE]
      planned<-switch(arm,mirna_cyt=length(unique(p$mirna_id)),gene_cyt=length(unique(p$gene)),mirna_gene=nrow(p))
      add(paste(adj,arm,"unique hypothesis count"),nrow(b)==planned)
      arm_a<-aa[aa$arm==arm,,drop=FALSE]
      arm_key<-function(z) if(arm=="mirna_cyt")z$mirna_id else if(arm=="gene_cyt")z$gene else key(z)
      ii<-match(arm_key(arm_a),arm_key(b))
      add(paste(adj,arm,"repeated-arm estimates equal unique source"),
        same_num(arm_a$rho,b$rho[ii])&&same_num(arm_a$p_raw,b$p_raw[ii])&&same_num(arm_a$fdr_BH,b$fdr_BH[ii]))
      for(method in c("BH","BY")) {
        ok<-is.finite(b$p_raw);expected<-rep(NA_real_,nrow(b))
        expected[ok]<-stats::p.adjust(b$p_raw[ok],method=method,n=planned)
        add(paste(adj,arm,method,"unique-arm FDR"),same_num(expected,b[[paste0("fdr_",method)]]))
      }
    }
    external_unique_comparison(adj)
    reported<-ext$validation_summary[ext$validation_summary$adjustment==adj,,drop=FALSE]
    add(paste(adj,"reported summary matches detail"),nrow(reported)==1L&&
      reported$fixed_pairs==nrow(q)&&reported$evaluable==sum(truth(q$evaluable))&&
      reported$conjunction_BH_and_effect==sum(truth(q$association_support_BH))&&
      reported$directional_vs_CPE_BRCA==sum(truth(q$directional_replication_vs_CPE_BRCA))&&
      reported$directional_vs_unadjusted_BRCA==sum(truth(q$directional_replication_vs_unadjusted_BRCA)))
  }
  z<-do.call(rbind,checks)
  data.table::fwrite(z,file.path(OUT,"logs","EXTERNAL_VALIDATION.tsv"),sep="\t")
  if(any(z$status=="FAIL")) stop("External output consistency check failed; see EXTERNAL_VALIDATION.tsv. No scientific result was replaced.")
  message(nrow(z)," external publication checks passed.")
}
external_arm_scale <- function() ggplot2::scale_colour_manual(
  values=c(mirna_cyt="#2166AC",mirna_gene="#B35806",gene_cyt="#1B7837"),
  labels=c(mirna_cyt="miRNA-CYT",mirna_gene="miRNA-gene",gene_cyt="Gene-CYT"),name="Association")
external_heatmap <- function(adjustment) {
  a<-external_bundle()$external_pair_arm_comparison
  d<-a[a$adjustment==adjustment,,drop=FALSE];p<-primary()
  d$pair<-factor(pair_label(d),levels=rev(pair_label(p)))
  d$arm_label<-factor(d$arm,levels=c("mirna_cyt","mirna_gene","gene_cyt"),labels=c("miRNA-CYT","miRNA-gene","Gene-CYT"))
  d$label<-ifelse(is.finite(d$rho),sprintf("%.2f",d$rho),"NA")
  source_table(d,paste0("External_heatmap_",adjustment))
  ggplot2::ggplot(d,ggplot2::aes(arm_label,pair,fill=rho))+
    ggplot2::geom_tile(colour="white",linewidth=.15)+
    ggplot2::geom_text(ggplot2::aes(label=label),size=2.15,colour="white")+
    ggplot2::scale_fill_gradient(low="white",high="#2166AC",limits=c(0,1),na.value="grey85",name="Spearman rho")+
    ggplot2::guides(fill=ggplot2::guide_colourbar(barwidth=grid::unit(45,"mm"),barheight=grid::unit(3,"mm")))+
    ggplot2::labs(x=NULL,y=NULL)+ggplot2::theme(axis.text.y=ggplot2::element_text(size=6.5),axis.ticks=ggplot2::element_blank())
}
figure_external <- function() {
  d<-external_unique_comparison("unadjusted");source_table(d,"Figure_6_external_comparison")
  s<-external_summary_table();s<-s[s$Analysis=="Unadjusted",,drop=FALSE]
  p<-primary();arm_n<-table(d$arm)
  B<-ggplot2::ggplot(d,ggplot2::aes(TCGA_BRCA_unadjusted_rho,rho,colour=arm))+
    ggplot2::geom_abline(slope=1,intercept=0,colour="grey65",linetype=2)+
    ggplot2::geom_point(size=1.6,alpha=.75)+external_arm_scale()+
    ggplot2::coord_fixed(xlim=c(0,1),ylim=c(0,1))+
    ggplot2::guides(colour=ggplot2::guide_legend(ncol=1))+
    ggplot2::labs(x="TCGA-BRCA Spearman rho",y="CPTAC-BRCA Spearman rho")
  save_panels("Figure_6",list(external_heatmap("unadjusted"),B),paste(
    "External replication of the fixed discovery pairs in CPTAC breast cancer.",
    sprintf("A, unadjusted Spearman correlations for all %d fixed pairs across %d matched tumours, ordered as in the primary discovery display.",nrow(p),s$Samples),
    "Cell labels show effect sizes, not significance; repeated miRNA-CYT and gene-CYT values represent shared hypotheses.",
    sprintf("B, external versus TCGA-BRCA unadjusted correlations for %d unique hypotheses (%d miRNA-CYT, %d miRNA-gene and %d gene-CYT tests).",nrow(d),arm_n['mirna_cyt'],arm_n['mirna_gene'],arm_n['gene_cyt']),
    "The dashed line denotes equality. The reference discovery set was not reselected using external results.",
    sprintf("%d pairs met the stored three-arm conjunction BH and effect-size criteria; complete direction agreement was observed for %d pairs against the unadjusted TCGA-BRCA reference and %d against the CPE-adjusted reference.",s$Three_arm_BH_and_effect,s$Direction_agreement_TCGA_unadjusted,s$Direction_agreement_TCGA_CPE),
    "These findings concern association replication in breast cancer and do not establish direct targeting, causal mediation, survival prediction or generalisation to all cancers (Table 2; Supplementary Table S11)."))
}
figure_external_sensitivity <- function() {
  a<-external_unique_comparison("unadjusted");b<-external_unique_comparison("WES_purity")
  s<-external_summary_table();unadj_s<-s[s$Analysis=="Unadjusted",,drop=FALSE];wes_s<-s[s$Analysis!="Unadjusted",,drop=FALSE]
  d<-merge(a[c("arm","hypothesis","rho")],b,by=c("arm","hypothesis"),suffixes=c("_unadjusted","_WES"))
  source_table(d,"Figure_S6_external_sensitivity")
  A<-ggplot2::ggplot(d,ggplot2::aes(rho_unadjusted,rho_WES,colour=arm))+
    ggplot2::geom_abline(slope=1,intercept=0,colour="grey65",linetype=2)+ggplot2::geom_point(size=1.7,alpha=.75)+external_arm_scale()+
    ggplot2::coord_fixed(xlim=c(0,1),ylim=c(0,1))+
    ggplot2::guides(colour=ggplot2::guide_legend(ncol=1))+
    ggplot2::labs(x="CPTAC unadjusted Spearman rho",y="CPTAC WES-purity partial Spearman rho")
  d$discordant<-is.finite(d$TCGA_BRCA_CPE_rho)&is.finite(d$rho_WES)&sign(d$TCGA_BRCA_CPE_rho)!=sign(d$rho_WES)
  B<-ggplot2::ggplot(d,ggplot2::aes(TCGA_BRCA_CPE_rho,rho_WES,colour=arm))+
    ggplot2::geom_hline(yintercept=0,colour="grey80")+ggplot2::geom_vline(xintercept=0,colour="grey80")+
    ggplot2::geom_abline(slope=1,intercept=0,colour="grey65",linetype=2)+ggplot2::geom_point(size=1.7,alpha=.75)+external_arm_scale()+
    ggplot2::geom_point(data=d[d$discordant,,drop=FALSE],shape=1,colour="black",size=3.2,stroke=.7)+
    ggplot2::coord_fixed(xlim=c(-.25,1),ylim=c(0,1))+
    ggplot2::guides(colour=ggplot2::guide_legend(ncol=1))+
    ggplot2::labs(x="TCGA-BRCA CPE-adjusted rho",y="CPTAC WES-purity partial Spearman rho")
  save_panels("Figure_S6",list(A,B),paste(
    "Purity sensitivity of external breast-cancer associations. A, external unadjusted versus WES-purity-adjusted coefficients.",
    "B, TCGA-BRCA CPE-adjusted coefficients versus external WES-purity-adjusted coefficients; the outlined point identifies a direction-discordant hypothesis.",
    sprintf("Each panel shows unique hypotheses rather than repeating shared arms for each pair. Unadjusted estimates use %d tumours; WES sensitivity uses %d.",unadj_s$Samples,wes_s$Samples),
    sprintf("%d of %d fixed pairs met the stored conjunction and effect-size criteria after WES adjustment; %d had complete CPE-reference direction agreement. Discordant hypotheses and their reference coefficients are retained in Supplementary Table S11.",wes_s$Three_arm_BH_and_effect,wes_s$Fixed_pairs,wes_s$Direction_agreement_TCGA_CPE),
    "CPE and WES purity differ, and the sensitivity uses residualised ranks; this is not an identical CPE-analysis replication (Supplementary Table S11)."))
}
external_summary_table <- function() {
  ext<-external_bundle();r<-ext$external_pair_validation;s<-ext$external_sample_annotation
  do.call(rbind,lapply(c("unadjusted","WES_purity"),function(adj) {
    d<-r[r$adjustment==adj,,drop=FALSE]
    data.frame(Analysis=if(adj=="unadjusted")"Unadjusted" else "WES purity sensitivity",
      Samples=if(adj=="unadjusted")nrow(s) else sum(is.finite(s$WES_purity)),
      Fixed_pairs=nrow(d),Evaluable_pairs=sum(truth(d$evaluable)),
      Three_arm_BH_and_effect=sum(truth(d$association_support_BH)),
      Three_arm_BY_and_effect=sum(d$conjunction_fdr_BY<CONFIG$fdr & d$minimum_abs_rho>=CONFIG$primary_threshold,na.rm=TRUE),
      Direction_agreement_TCGA_unadjusted=sum(truth(d$directional_replication_vs_unadjusted_BRCA)),
      Direction_agreement_TCGA_CPE=sum(truth(d$directional_replication_vs_CPE_BRCA)))
  }))
}
write_main_tables_workbook <- function(mt) {
  main_wb<-openxlsx::createWorkbook()
  main_tabs<-list(Table_1=mt,Table_2=external_summary_table())
  names(main_tabs$Table_1)<-c("miRNA","Pairs","Genes","Median recurrence","Minimum recurrence","Maximum recurrence")
  names(main_tabs$Table_2)<-c("Analysis","Samples","Fixed pairs","Evaluable pairs","Three-arm support (BH)","Three-arm support (BY)","Direction concordance vs TCGA without adjustment","Direction concordance vs TCGA with CPE")
  for(nm in names(main_tabs)) {
    openxlsx::addWorksheet(main_wb,nm)
    openxlsx::writeData(main_wb,nm,main_tabs[[nm]],keepNA=TRUE,na.string="NA")
    openxlsx::addStyle(main_wb,nm,openxlsx::createStyle(textDecoration="bold",fgFill="#D9E5F2",wrapText=TRUE,valign="center"),rows=1,cols=seq_len(ncol(main_tabs[[nm]])),gridExpand=TRUE)
    openxlsx::setColWidths(main_wb,nm,cols=seq_len(ncol(main_tabs[[nm]])),widths=22)
    openxlsx::setColWidths(main_wb,nm,cols=1,widths=30)
    openxlsx::setRowHeights(main_wb,nm,rows=1,heights=64)
    openxlsx::freezePane(main_wb,nm,firstRow=TRUE)
  }
  openxlsx::saveWorkbook(main_wb,file.path(OUT,"tables","Main_Tables_1_2.xlsx"),overwrite=TRUE)
}
write_external_tables <- function() {
  ext<-external_bundle();summary<-external_summary_table()
  data.table::fwrite(summary,file.path(OUT,"tables","Table_2_External_validation_summary.tsv"),sep="\t")
  for(nm in external_files) source_table(ext[[nm]],paste0("External_",nm))
  for(nm in c("external_pair_validation","external_pair_arm_comparison","external_unique_arm_tests"))
    data.table::fwrite(ext[[nm]],file.path(OUT,"tables",paste0("S11_",nm,".tsv")),sep="\t",na="NA")
}
write_external_handoff <- function() {
  ext<-external_bundle();r<-ext$external_pair_validation
  summary<-external_summary_table()
  writeLines(c("# External validation output map", "",
    "Figure 6: external primary effect-size heatmap and comparison with unadjusted TCGA-BRCA effects.",
    "Supplementary Figure S6: purity sensitivity and CPE-reference comparison.",
    "Table 2: external validation counts. Supplementary Table S11: complete external results, coverage and mapping.",
    "Original supplementary tables S1-S10 retain their numbering.","",capture.output(print(summary,row.names=FALSE)),"",
    "For the subsequent manuscript edit: add the independent cohort and locked 47-pair validation design to Methods; add external results; revise Discussion, Limitations and Conclusion; replace the blanket statement that no external validation was performed.",
    "Retain limitations: BRCA only; WES is not CPE; no independent survival, targeting, mediation or cellular-context interaction validation.",
    "Report direction agreement against both unadjusted and CPE-adjusted TCGA-BRCA references using the counts above. Retain all discordant hypotheses and their coefficients in the report.",
    "The existing manuscript Word file is not edited by this publication package.",
    "Reference: Krug K et al. Cell. 2020;183:1436-1456.e31. doi:10.1016/j.cell.2020.10.036.",
    "Source: https://www.linkedomics.org/data_download/CPTAC-pancan-BRCA/"),
    file.path(OUT,"text","External_validation_handoff.md"),useBytes=TRUE)
  writeLines(c("# Table legends","",
    "Table 1. Distribution of the primary miRNA-gene pairs across miRNAs. Recurrence is the number of cancer types meeting the three-arm discovery criteria.","",
    "Table 2. External validation of fixed discovery pairs in CPTAC breast cancer. BH and BY support require conjunction FDR<0.05 and absolute correlation>=0.30 for all three arms. Direction counts additionally require agreement of all three coefficients with the stated TCGA-BRCA reference. WES-purity sensitivity is distinct from CPE adjustment; missing or non-estimable pairs remain in the fixed denominator.","",
    "Supplementary Table S11. Complete external-validation results, unique-arm tests, coverage, mapping and testing-family definitions. FDR values are the stored Stage 25 estimates. WES sensitivity uses partial Spearman coefficients after winsorisation and rank residualisation. Shared arms are not independent repetitions."),
    file.path(OUT,"text","Table_legends.md"),useBytes=TRUE)
}
