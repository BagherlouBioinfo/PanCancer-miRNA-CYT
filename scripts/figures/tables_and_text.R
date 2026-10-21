write_tables <- function() {
  technical<-file.path(OUT,"source_data","technical_tables")
  dir.create(technical,recursive=TRUE,showWarnings=FALSE)
  specs<-list(
    S01_Primary=list("tables/primary/triplets_rho030_CPE.tsv"),
    S02_Immune=list("tables/downstream/immune_adjusted_recurrence_by_population.tsv"),
    S03_Meta=list("tables/downstream/meta_analysis_all_arms.tsv"),
    S04_OS=list("tables/downstream/survival_models.tsv"),
    S05_DSS=list("tables/downstream/survival_models_DSS.tsv"),
    S06_PFI=list("tables/downstream/survival_models_PFI.tsv"),
    S07_Incremental=list("tables/downstream/cyt_incremental_models.tsv"),
    S08_Mediation=list("tables/downstream/mediation_all_contexts.tsv"),
    S09_CellContext=list("tables/downstream/stage24_final_robust_nonoverlapping.tsv"),
    S10_CellFunnel=list("tables/downstream/stage24_primary_evidence_funnel.tsv"),
    S11_MatchedNull=list("tables/audit/matched_null_target_specificity_results.tsv"),
    S12_Redundancy=list("tables/audit/mirna_redundancy_summary.tsv"),
    S13_GO=list("tables/downstream/enrichment_go_bp_results.tsv"),
    S14_Reactome=list("tables/downstream/enrichment_reactome_results.tsv"),
    S15_LOOStable=list("tables/downstream/enrichment_terms_robust_to_leave_one_out.tsv"),
    S16_DrugAnnotations=list("tables/downstream/drug_gene_interactions.tsv"),
    S17_DrugResolution=list("tables/downstream/drug_gene_id_resolution.tsv"),
    S18_ClusterValidity=list("tables/downstream/cluster_validity.tsv"),
    S19_ClusterStability=list("tables/downstream/cluster_stability.tsv"),
    S20_ClusterNull=list("tables/downstream/cluster_null_calibration.tsv"),
    S21_ClusterCharacteristics=list("tables/downstream/cluster_characteristics.tsv"),
    S22_CellDiagnostics=list("tables/downstream/stage24_lump_and_stage_sensitivity.tsv"),
    S23_CellNonlinearity=list("tables/downstream/stage24_nonlinearity_sensitivity.tsv"),
    S24_CellCollinearity=list("tables/downstream/stage24_collinearity_sensitivity.tsv"),
    S25_CellPermutation=list("tables/downstream/stage24_negative_controls.tsv"),
    S26_CellConjunction=list("tables/downstream/stage24_corrected_conjunction.tsv"),
    S28_ClinicalCompleteness=list("tables/downstream/clinical_data_completeness_by_cancer.tsv"),
    S29_EndpointCompleteness=list("tables/downstream/survival_endpoint_completeness.tsv"),
    S30_PurityExclusion=list("tables/audit/purity_sensitivity_option1_vs_option2_comparison.tsv"),
    S31_SeedAssignment=list("tables/audit/mirna_redundancy_seed_family_assignment.tsv"),
    S32_CancerContext=list("tables/downstream/stage24_power_by_cancer.tsv"),
    S33_EnrichmentLOO=list("tables/downstream/enrichment_sensitivity_leave_one_out.tsv"),
    S34_EnrichmentFragile=list("tables/downstream/enrichment_terms_single_gene_driven.tsv"),
    S35_HistoricalComparison=list("tables/downstream/stage24_historical_reconciliation.tsv"))
  tabs<-list();index<-list();dictionary<-list()
  for(nm in names(specs)) {
    rel<-specs[[nm]][[1]]
    if(!file.exists(file.path(CONFIG$run_root,rel))) {
      # Historical reconciliation is not used in any publication table/figure.
      if(nm=="S35_HistoricalComparison") next
      table_issues[[nm]] <<- data.frame(table=nm,detail=paste("Missing input:",rel));next
    }
    d<-read_run(rel)
    note<-"All source rows retained; original testing families and missing estimates preserved."
    if(nm=="S01_Primary") {
      d<-d[truth(d$qualifies),,drop=FALSE];d<-d[order(-d$strict_recurrence,d$mirna_id,d$gene),]
      note<-"Qualifying primary pairs only; qualifies == TRUE."
    }
    data.table::fwrite(d,file.path(technical,paste0(nm,".tsv")),sep="\t",na="NA")
    tabs[[nm]]<-d
    index[[nm]]<-data.frame(table=nm,source=rel,rows=nrow(d),columns=ncol(d),note=note)
    dictionary[[nm]]<-data.frame(table=nm,column=names(d),
      parsed_type=vapply(d,function(x)paste(class(x),collapse="/"),character(1)),
      missing=vapply(d,function(x)sum(is.na(x)),integer(1)))
  }
  # Cancer-level availability summary; individual identifiers are not exported here.
  ann<-read_run("processed/stage05/sample_annotation_with_purity.tsv",c("cancer","CPE","LUMP"))
  cohort<-do.call(rbind,lapply(sort(unique(ann$cancer)),function(ca) {
    x<-ann[ann$cancer==ca,,drop=FALSE]
    data.frame(cancer=ca,n_matched=nrow(x),n_CPE=sum(is.finite(x$CPE)),
      n_LUMP=sum(is.finite(x$LUMP)))
  }))
  for(a in c("unadjusted","cpe","lump")) {
    e<-read_run(paste0("processed/stage05/eligible_cancers_",a,".tsv"),"cancer")
    cohort[[paste0("eligible_",a)]]<-cohort$cancer %in% e$cancer
  }
  tabs$S27_Cohort<-cohort
  data.table::fwrite(cohort,file.path(technical,"S27_Cohort.tsv"),sep="\t")
  index$S27_Cohort<-data.frame(table="S27_Cohort",source="stage05 annotation and eligibility lists",
    rows=nrow(cohort),columns=ncol(cohort),note="Availability and eligibility are separate quantities.")
  # Full sensitivity listings with original qualifies flags, not significant rows only.
  for(t in c("030","025","020")) for(a in c("unadjusted","CPE","LUMP")) {
    folder<-if(t=="030" && a!="LUMP") "primary" else "sensitivity"
    rel<-sprintf("tables/%s/triplets_rho%s_%s.tsv",folder,t,a)
    d<-read_run(rel);nm<-paste0("Sensitivity_",t,"_",a)
    data.table::fwrite(d,file.path(technical,paste0(nm,".tsv")),sep="\t",na="NA")
    tabs[[nm]]<-d[truth(d$qualifies),,drop=FALSE]
    index[[nm]]<-data.frame(table=nm,source=rel,rows=nrow(d),columns=ncol(d),
      note="All tested pairs retained, including qualifies == FALSE.")
  }
  # Six-row main-text summary; the full pair list is in Supplementary Table S2.
  p<-primary()
  mt<-do.call(rbind,lapply(sort(unique(p$mirna_id)),function(mi) {
    d<-p[p$mirna_id==mi,,drop=FALSE]
    data.frame(miRNA=mi,n_pairs=nrow(d),n_genes=length(unique(d$gene)),
      median_recurrence=median(d$strict_recurrence),min_recurrence=min(d$strict_recurrence),
      max_recurrence=max(d$strict_recurrence))
  }))
  data.table::fwrite(mt,file.path(OUT,"tables","Table_1_Primary_summary.tsv"),sep="\t")
  external<-external_bundle()
  external_sections<-c("validation_summary","coverage_summary","external_pair_validation",
    "external_pair_arm_comparison","external_unique_arm_tests","gene_mapping_audit",
    "mirna_mapping_audit","testing_families")
  for(sec in external_sections) {
    nm<-paste0("External_",sec);d<-external[[sec]];tabs[[nm]]<-d
    data.table::fwrite(d,file.path(technical,paste0(nm,".tsv")),sep="\t",na="NA")
    index[[nm]]<-data.frame(table=nm,source=paste0("external/",sec,".tsv"),rows=nrow(d),columns=ncol(d),note="Stored external validation; all fixed pairs retained.")
    dictionary[[nm]]<-data.frame(table=nm,column=names(d),parsed_type=vapply(d,function(x)paste(class(x),collapse="/"),character(1)),missing=vapply(d,function(x)sum(is.na(x)),integer(1)))
  }
  added <- review_table_sections()
  for(nm in names(added)) {
    d<-added[[nm]];tabs[[nm]]<-d
    data.table::fwrite(d,file.path(technical,paste0(nm,".tsv")),sep="\t",na="NA")
    index[[nm]]<-data.frame(table=nm,source="verified post-review add-on",rows=nrow(d),columns=ncol(d),note="Stored inference retained; no model refit.")
    dictionary[[nm]]<-data.frame(table=nm,column=names(d),parsed_type=vapply(d,function(x)paste(class(x),collapse="/"),character(1)),missing=vapply(d,function(x)sum(is.na(x)),integer(1)))
  }
  data.table::fwrite(do.call(rbind,index),file.path(OUT,"logs","SOURCE_TABLE_INDEX.tsv"),sep="\t")
  data.table::fwrite(do.call(rbind,dictionary),file.path(OUT,"logs","DATA_DICTIONARY.tsv"),sep="\t")
  groups<-list(
    S1_Cohort=c("S27_Cohort","S28_ClinicalCompleteness","S29_EndpointCompleteness"),
    S2_Discovery=c("S01_Primary",setdiff(grep("^Sensitivity_",names(tabs),value=TRUE),"Sensitivity_030_CPE"),"S30_PurityExclusion"),
    S3_Immune="S02_Immune",
    S4_Meta="S03_Meta",
    S5_Survival=c("S04_OS","S05_DSS","S06_PFI","S07_Incremental","Cluster_Summary","Cluster_Models","Cluster_PH","Cluster_Time_HR","Cluster_SampleFlow"),
    S6_Mediation="S08_Mediation",
    S7_Clustering=c("S18_ClusterValidity","S19_ClusterStability","S20_ClusterNull","S21_ClusterCharacteristics"),
    S8_Specificity=c("S11_MatchedNull","S12_Redundancy","S31_SeedAssignment","Mir146_Summary","Mir146_Conditional","Mir146_Both","Mir146_Expression"),
    S9_Annotation=c("S13_GO","S14_Reactome","S15_LOOStable","S33_EnrichmentLOO",
      "S34_EnrichmentFragile","S16_DrugAnnotations","S17_DrugResolution"),
    S10_CellContext=c("S10_CellFunnel","S09_CellContext","S22_CellDiagnostics",
      "S23_CellNonlinearity","S24_CellCollinearity","S25_CellPermutation","S26_CellConjunction","S32_CancerContext"),
    S11_External=paste0("External_",external_sections))
  titles<-c(S1_Cohort="Cohort and analysis eligibility",S2_Discovery="Discovery and purity sensitivity",
    S3_Immune="Immune-composition sensitivity",S4_Meta="Cross-cancer meta-analysis",
    S5_Survival="Feature and frozen-cluster survival associations",S6_Mediation="Exploratory mediation",
    S7_Clustering="Clustering diagnostics",S8_Specificity="Matched-null, miRNA-family and conditional miR-146 assessments",
    S9_Annotation="Functional and drug annotations",S10_CellContext="Corrected cellular-context analysis",
    S11_External="External validation in CPTAC breast cancer")
  esc<-function(x) {
    x<-as.character(x);x[is.na(x)]<-"NA"
    x<-gsub("&","&amp;",x,fixed=TRUE);x<-gsub("<","&lt;",x,fixed=TRUE)
    gsub(">","&gt;",x,fixed=TRUE)
  }
  group_index<-list()
  wb<-if(CONFIG$write_xlsx_if_available && requireNamespace("openxlsx",quietly=TRUE)) openxlsx::createWorkbook() else NULL
  for(g in names(groups)) {
    sections<-groups[[g]]
    missing_sections<-setdiff(sections,names(tabs))
    sections<-intersect(sections,names(tabs))
    if(!is.null(wb)) openxlsx::addWorksheet(wb,g)
    path<-file.path(OUT,"tables",paste0(g,".html"))
    con<-file(path,open="wt",encoding="UTF-8")
    writeLines(c('<!doctype html><meta charset="utf-8"><style>body{font:12px Arial;margin:24px}table{border-collapse:collapse;margin-bottom:24px}th,td{border:1px solid #ccc;padding:5px;vertical-align:top}th{background:#d9e5f2}h2{margin-top:30px}</style>',
      paste0("<h1>",esc(sub("_.*","",g)),". ",esc(titles[g]),"</h1>"),
      "<p>Sections retain their original schemas. NA means missing or not estimable, not zero. Original FDR values are preserved.</p>"),con)
    if(length(missing_sections)) writeLines(paste0("<p>INCOMPLETE: missing sections: ",
      esc(paste(missing_sections,collapse=", ")),"</p>"),con)
    row<-4L
    maxcols<-max(6L,vapply(tabs[sections],ncol,integer(1)))
    if(!is.null(wb)) {
      openxlsx::writeData(wb,g,paste(sub("_.*","",g),titles[g]),startRow=1)
      note<-if(g=="S11_External") "Fixed discovery pairs; stored external FDR values. WES purity is not CPE. Direction agreement is reported against both TCGA-BRCA references. NA is missing, not zero." else if(g=="S5_Survival") "Feature and cluster models are separate analyses. Cluster FDR families are endpoint x tier across all attempted cancers. PH flags and time-varying OS sensitivities must accompany interpretation; concordance changes are apparent, not validated." else if(g=="S8_Specificity") "Matched-null and conditional miR-146 tests address different questions. Conditional families pool both miRNAs by purity setting x outcome class. Stored HC3, BH and BY results retained; no causal or biological independence claim." else "NA is missing, not zero. Original FDR families and model statuses are retained."
      openxlsx::writeData(wb,g,note,startRow=2)
      openxlsx::mergeCells(wb,g,cols=seq_len(maxcols),rows=1)
      openxlsx::mergeCells(wb,g,cols=seq_len(maxcols),rows=2)
      openxlsx::addStyle(wb,g,openxlsx::createStyle(wrapText=TRUE,valign="center"),rows=2,cols=1)
      openxlsx::setRowHeights(wb,g,rows=c(1,2),heights=c(24,36))
    }
    for(sec in sections) {
      d<-tabs[[sec]]
      # Keep provenance in the technical archive rather than article-facing sections.
      d<-d[setdiff(names(d),c("run_id","generating_script","in_historical_stage23"))]
      label<-gsub("_"," ",sub("^S[0-9]+_","",sec))
      if(sec %in% names(review_section_labels)) label<-unname(review_section_labels[sec])
      writeLines(paste0("<h2>",esc(label)," (",nrow(d)," rows)</h2><table><thead><tr>",
        paste0("<th>",esc(names(d)),"</th>",collapse=""),"</tr></thead><tbody>"),con)
      if(nrow(d)) for(i in seq_len(nrow(d))) writeLines(paste0("<tr>",
        paste0("<td>",esc(unlist(d[i,,drop=FALSE],use.names=FALSE)),"</td>",collapse=""),"</tr>"),con)
      writeLines("</tbody></table>",con)
      if(!is.null(wb)) {
        openxlsx::writeData(wb,g,label,startRow=row)
        openxlsx::writeData(wb,g,d,startRow=row+1L,withFilter=FALSE,keepNA=TRUE,na.string="NA")
        openxlsx::addStyle(wb,g,openxlsx::createStyle(textDecoration="bold",fgFill="#D9E5F2",wrapText=TRUE,valign="center"),
          rows=row+1L,cols=seq_len(ncol(d)),gridExpand=TRUE)
        openxlsx::setRowHeights(wb,g,rows=row+1L,heights=42)
        if(nrow(d)) {
          pcols<-which(grepl("(^fdr$|fdr_|_fdr$|^bh_fdr$|^raw_p$|^raw_bootstrap_p$|^p_value$|^p_|_p$|^QEp$)",names(d),ignore.case=TRUE))
          if(length(pcols)) openxlsx::addStyle(wb,g,openxlsx::createStyle(numFmt="0.000E+00"),
            rows=seq.int(row+2L,length.out=nrow(d)),cols=pcols,gridExpand=TRUE,stack=TRUE)
        }
      }
      group_index[[length(group_index)+1L]]<-data.frame(table=g,section=label,
        rows=nrow(d),xlsx_title_row=row,note=if(g=="S2_Discovery")
          "Sensitivity listings show qualifying pairs; all attempted pairs are archived in source_data/technical_tables." else
          "All available source rows retained, including not-estimable models.")
      row<-row+nrow(d)+5L
    }
    close(con)
    if(!is.null(wb)) {
      openxlsx::freezePane(wb,g,firstActiveRow=6)
      openxlsx::setColWidths(wb,g,cols=seq_len(maxcols),widths=20)
      openxlsx::setColWidths(wb,g,cols=1,widths=27)
    }
  }
  data.table::fwrite(do.call(rbind,group_index),file.path(OUT,"tables","Supplementary_INDEX.tsv"),sep="\t")
  if(!is.null(wb)) openxlsx::saveWorkbook(wb,file.path(OUT,"tables","Supplementary_Tables_S1_S11.xlsx"),overwrite=TRUE)
  else {
    table_issues[["Excel_workbook"]]<<-data.frame(table="Excel_workbook",detail="openxlsx is not installed; HTML tables were created, XLSX was not.")
    message("Excel not created: openxlsx is not installed. Eleven HTML tables are available; install openxlsx for the grouped workbook.")
  }
  if(!is.null(wb)) write_main_tables_workbook(mt)
  invisible(NULL)
}
write_results_draft <- function() {
  p<-primary()
  med<-down("mediation_all_contexts",c("bh_fdr","primary_support"))
  null<-audit("matched_null_target_specificity_results",c("fdr"))
  s<-down("survival_models",c("tier","status","fdr","ph_violated"))
  final<-down("stage24_final_robust_nonoverlapping",c("mirna_id","gene","cancer","population"))
  drug<-down("drug_gene_interactions",c("record_status"))
  estimated<-s[s$status=="ESTIMATED",,drop=FALSE]
  tx<-c("# Results draft: verified output summaries",
    "", "Draft only: complete interpretation, citations to supplementary items and journal formatting require author review.","",
    sprintf("The primary CPE-adjusted analysis identified %d recurrent miRNA-gene pairs involving %d miRNAs and %d genes (Fig. 2; Table 1). The pairs contributed %d qualifying pair-cancer records.",
      nrow(p),length(unique(p$mirna_id)),length(unique(p$gene)),sum(p$strict_recurrence)),"",
    sprintf("In the matched-null comparison, %d of %d pairs had BH-FDR < 0.05 (Supplementary Fig. S5; Supplementary Table S8). This comparison provided no evidence of target specificity under the implemented matched-control design when none passed; it does not establish absence of regulation.",
      sum(null$fdr<CONFIG$fdr,na.rm=TRUE),nrow(null)),"",
    sprintf("Across all four OS model tiers, %d of %d attempted models were estimable; %d estimated models had a global proportional-hazards flag. Associations are reported by tier and cancer rather than as independent replications (Fig. 4; Supplementary Table S5).",
      nrow(estimated),nrow(s),sum(truth(estimated$ph_violated))),"",
    sprintf("Exploratory mediation met BH-FDR < 0.05 in %d of %d eligible contexts, including %d of %d primary-support contexts (Supplementary Fig. S2; Supplementary Table S6). Preselection and cross-sectional sampling limit causal interpretation.",
      sum(med$bh_fdr<CONFIG$fdr,na.rm=TRUE),nrow(med),
      sum(truth(med$primary_support)&med$bh_fdr<CONFIG$fdr,na.rm=TRUE),sum(truth(med$primary_support))),"",
    sprintf("The corrected cellular-context funnel retained %d triplet-cancer-population combinations representing %d unique miRNA-gene pairs across %d cancers and %d populations (Fig. 5; Supplementary Table S10). These results describe context-dependent associations.",
      nrow(final),length(unique(key(final))),length(unique(final$cancer)),length(unique(final$population))),"",
    sprintf("The drug annotation table contained %d database-listed interaction records and %d other status rows; status rows were not counted as interactions (Supplementary Table S9). These annotations do not establish therapeutic efficacy.",
      sum(drug$record_status=="database_listed_interaction",na.rm=TRUE),
      sum(drug$record_status!="database_listed_interaction",na.rm=TRUE)),"",
    sprintf("The indirect-to-total effect ratio exceeded one in %d contexts, including %d primary-support contexts. These estimates are retained as effect ratios, not fractions or percentages mediated.",
      sum(med$prop_mediated>1,na.rm=TRUE),sum(truth(med$primary_support)&med$prop_mediated>1,na.rm=TRUE)),"",
    "The complete inferential narrative must additionally consider immune-adjustment exclusions, prediction intervals, clustering diagnostics and seed-family mapping limitations. Automated counts alone do not settle these interpretations.")
  writeLines(tx,file.path(OUT,"text","Results_counts_summary.md"),useBytes=TRUE)
  # The Word manuscript will be updated separately after the publication outputs
  # are reviewed. Do not copy the obsolete pre-external-validation narrative.
}
