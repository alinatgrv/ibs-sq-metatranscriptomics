suppressPackageStartupMessages({
  library(Maaslin2)
})

input_data <- "/home/alina_tgrv/beegfs/IBS_SQ/results/humann_metatranscriptome/joined_all/pathabundance_311_unstratified_for_maaslin2.tsv"
input_metadata <- "/home/alina_tgrv/beegfs/IBS_SQ/metadata/metadata_311_for_maaslin2_pathways.tsv"
output_dir <- "/home/alina_tgrv/beegfs/IBS_SQ/results/humann_metatranscriptome/maaslin2_pathways_311_tss_log"

fit_data <- Maaslin2(
  input_data = input_data,
  input_metadata = input_metadata,
  output = output_dir,
  fixed_effects = c("ibs_status", "age", "sex", "bmi"),
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  min_prevalence = 0.1,
  reference = c("ibs_status,Control"),
  standardize = FALSE,
  plot_heatmap = TRUE,
  plot_scatter = TRUE
)
