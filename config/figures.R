CONFIG <- list(
  write_workflow = TRUE,
  python = "python3",
  run_root = "results",
  output_root = "results/manuscript",
  version = "3.2.0",
  # Stage 27 resolves these paths to current completed analysis outputs.
  external_root = NULL,
  # Archived tables are not substituted by the main runner.
  cluster_survival_root = NULL, mir146_root = NULL,
  cluster_time_example_cancers = c("LGG","LUAD"),
  run_id = "analysis",
  font = "sans", base_size = 9, dpi = 600,
  width_mm = 220, height_mm = 250,
  figure_sizes = list(Figure_2=c(220,260), Figure_3=c(220,259),
    Figure_4=c(220,230), Figure_5=c(240,270),
    Figure_S2=c(220,260),Figure_S3=c(240,270),
    Figure_6=c(220,249),Figure_S6=c(260,160),
    Figure_S7=c(285,310),Figure_S8=c(285,200)),
  panel_letter_size = 12, legend_size = 8,
  align_panels = list(Figure_2=c(2,3),Figure_S7=c(1,2)),
  # Each row: panel index, grid row, grid column, row span, column span.
  layouts = list(
    Figure_2=rbind(c(1,1,1,1,2),c(2,2,1,2,1),c(3,2,2,2,1)),
    Figure_3=rbind(c(1,1,1,4,2),c(2,5,1,1,2)),
    Figure_4=rbind(c(1,1,1,1,2),c(2,2,1,1,1),c(3,2,2,1,1)),
    Figure_5=rbind(c(1,1,1,1,2),c(2,2,1,2,2)),
    Figure_6=rbind(c(1,1,1,1,2),c(2,1,3,1,1)),
    Figure_S6=rbind(c(1,1,1,1,1),c(2,1,2,1,1)),
    Figure_S7=rbind(c(1,1,1,2,1),c(2,1,2,2,1),c(3,3,1,1,2)),
    Figure_S8=rbind(c(1,1,1,1,1),c(2,1,2,1,1),c(3,2,1,1,1),c(4,2,2,1,1))),
  colours = c(CPE="#2166AC", LUMP="#1B7837", unadjusted="#969696"),
  clusters = c("1"="#2166AC", "2"="#B35806"),
  primary_threshold = 0.30, fdr = 0.05,
  # Descriptive display selection only; underlying hypothesis families stay unchanged.
  survival_tier = "Tier3", survival_max_rows = 18,
  survival_features = NULL,  # e.g. c("hsa-miR-150-5p", "CARD11")
  enrichment_top_n = 12,
  write_xlsx_if_available = TRUE,
  # Publication-only regeneration; scientific results are never refitted.
  modules = c("Figure_2","Figure_3","Figure_4","Figure_5",
    "Figure_6","Figure_S1","Figure_S2","Figure_S3","Figure_S4","Figure_S5","Figure_S6","Figure_S7","Figure_S8",
    "External_tables","Supplementary_tables","Results_draft"),
  expected_primary_pairs = 47L, expected_final_contexts = 31L,
  expected_final_pairs = 19L
)
