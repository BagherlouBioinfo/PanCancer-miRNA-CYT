validate_outputs <- function() {
  checks<-list()
  add<-function(label,ok,detail="") {
    checks[[length(checks)+1L]]<<-data.frame(check=label,status=if(isTRUE(ok)) "PASS" else "FAIL",detail=as.character(detail))
  }
  p<-primary()
  add("Primary pair keys are unique",!anyDuplicated(key(p)))
  n_support<-vapply(strsplit(p$supporting_cancers,";"),function(x)length(unique(trimws(x[nzchar(trimws(x))]))),integer(1))
  add("Recurrence equals supporting cancer count",all(n_support==p$strict_recurrence))
  add("Primary recurrence threshold",all(p$strict_recurrence>=3))
  cohort<-read_run("processed/stage05/sample_annotation_with_purity.tsv",c("sample_id","cancer","CPE","LUMP"))
  add("Matched sample IDs unique",!anyDuplicated(cohort$sample_id))
  add("Observed purity fractions within zero and one",all(is.na(cohort$CPE)|(cohort$CPE>=0&cohort$CPE<=1)) &&
    all(is.na(cohort$LUMP)|(cohort$LUMP>=0&cohort$LUMP<=1)))
  meta<-down("meta_analysis_all_arms",c("arm","adjustment","population","pooled_rho","ci_lower","ci_upper","status"))
  add("Meta model keys unique",!anyDuplicated(key(meta,c("mirna_id","gene","arm","adjustment","population"))))
  est<-meta[meta$status=="ESTIMATED",,drop=FALSE]
  add("Meta estimates inside reported intervals",all(est$ci_lower<=est$pooled_rho & est$pooled_rho<=est$ci_upper & abs(est$pooled_rho)<=1))
  for(endpoint in c("OS","DSS","PFI")) {
    d<-down(if(endpoint=="OS") "survival_models" else paste0("survival_models_",endpoint))
    need(d,c("feature","cancer","tier","status","HR","ci_lower","ci_upper"),endpoint)
    add(paste(endpoint,"model keys unique"),!anyDuplicated(key(d,c("feature","cancer","tier"))))
    e<-d[d$status=="ESTIMATED",,drop=FALSE]
    add(paste(endpoint,"HR intervals ordered"),all(e$ci_lower>0 & e$ci_lower<=e$HR & e$HR<=e$ci_upper))
  }
  med<-down("mediation_all_contexts",c("point_a","point_b","point_indirect","point_c_prime","point_total",
    "prop_mediated","prop_mediated_meaningful","primary_support","bh_fdr"))
  add("Mediation indirect effect equals a times b",isTRUE(all.equal(med$point_indirect,med$point_a*med$point_b,tolerance=1e-8)))
  add("Mediation effect decomposition",isTRUE(all.equal(med$point_total,med$point_c_prime+med$point_indirect,tolerance=1e-8)))
  add("Mediation context keys unique",!anyDuplicated(key(med,c("mirna_id","gene","cancer"))))
  k<-c("mirna_id","gene","cancer","population")
  con<-down("stage24_corrected_conjunction",c(k,"all_three_arms_supported","overlap_category"))
  lump<-down("stage24_lump_and_stage_sensitivity",c(k,"lump_sensitivity_pass"))
  coll<-down("stage24_collinearity_sensitivity",c(k,"coefficient_stable","collinearity_assessment"))
  perm<-down("stage24_negative_controls",c(k,"passes_permutation_criterion"))
  final<-down("stage24_final_robust_nonoverlapping",k)
  funnel<-down("stage24_primary_evidence_funnel",c("n_rows","n_triplets"))
  for(nm in c("con","lump","coll","perm","final")) add(paste("Unique cellular-context keys",nm),!anyDuplicated(key(get(nm),k)))
  q<-con[truth(con$all_three_arms_supported),,drop=FALSE];sizes<-nrow(q)
  q<-q[key(q,k) %in% key(lump[truth(lump$lump_sensitivity_pass),,drop=FALSE],k),,drop=FALSE];sizes<-c(sizes,nrow(q))
  q<-q[q$overlap_category=="nonoverlapping_population",,drop=FALSE];sizes<-c(sizes,nrow(q))
  keep<-coll[truth(coll$coefficient_stable)|coll$collinearity_assessment=="not_applicable_no_correlated_population",,drop=FALSE]
  q<-q[key(q,k) %in% key(keep,k),,drop=FALSE];sizes<-c(sizes,nrow(q))
  q<-q[key(q,k) %in% key(perm[truth(perm$passes_permutation_criterion),,drop=FALSE],k),,drop=FALSE];sizes<-c(sizes,nrow(q))
  add("Reconstructed sequential funnel",identical(as.integer(sizes),as.integer(funnel$n_rows)),paste(sizes,collapse=" -> "))
  add("Final rows match reconstructed final keys",setequal(key(q,k),key(final,k)))
  add("Final context count",nrow(final)==CONFIG$expected_final_contexts,nrow(final))
  add("Final unique pair count",length(unique(key(final)))==CONFIG$expected_final_pairs,length(unique(key(final))))
  drug<-down("drug_gene_interactions",c("record_status","gene_symbol"))
  status_counts<-as.data.frame(table(drug$record_status),stringsAsFactors=FALSE)
  names(status_counts)<-c("record_status","n_rows")
  source_table(status_counts,"drug_record_accounting")
  add("Drug status coverage",sum(status_counts$n_rows)==nrow(drug) && !anyNA(drug$record_status))
  diag<-data.frame(metric=c("n_mediation_contexts","n_primary_support_contexts","n_effect_ratios_gt_one",
    "n_primary_support_effect_ratios_gt_one"),value=c(nrow(med),sum(truth(med$primary_support)),
      sum(med$prop_mediated>1,na.rm=TRUE),sum(truth(med$primary_support)&med$prop_mediated>1,na.rm=TRUE)))
  data.table::fwrite(diag,file.path(OUT,"logs","INTERPRETATION_DIAGNOSTICS.tsv"),sep="\t")
  result<-do.call(rbind,checks)
  data.table::fwrite(result,file.path(OUT,"logs","VALIDATION.tsv"),sep="\t")
  if(any(result$status=="FAIL")) stop("Output validation failed. Read logs/VALIDATION.tsv; scientific inputs were not modified.")
  message(sum(result$status=="PASS")," data consistency checks passed.")
}
