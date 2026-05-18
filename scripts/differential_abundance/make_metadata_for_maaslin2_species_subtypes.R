#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

base_dir <- "/home/alina_tgrv/beegfs/IBS_SQ"

metadata_file <- file.path(
  base_dir,
  "metadata/metadata_326_for_maaslin2_species.tsv"
)

species_file <- file.path(
  base_dir,
  "results/metaphlan_metatranscriptome/joined/species_326_for_maaslin2.tsv"
)

out_metadata_file <- file.path(
  base_dir,
  "metadata/metadata_305_for_maaslin2_species_subtypes.tsv"
)

out_species_file <- file.path(
  base_dir,
  "results/metaphlan_metatranscriptome/joined/species_305_for_maaslin2_subtypes.tsv"
)

# Need original metadata because bowel_habit is not in metadata_326_for_maaslin2_species.tsv
raw_metadata_file <- file.path(
  base_dir,
  "metadata/metadata_326_clean_v2.tsv"
)

meta <- read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
raw_meta <- read.delim(raw_metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
species <- read.delim(species_file, check.names = FALSE, stringsAsFactors = FALSE)

cat("Input metadata dimensions:\n")
print(dim(meta))

cat("\nInput raw metadata dimensions:\n")
print(dim(raw_meta))

cat("\nInput species dimensions:\n")
print(dim(species))

# Add bowel_habit back
meta2 <- meta %>%
  left_join(raw_meta[, c("sample", "bowel_habit")], by = "sample") %>%
  mutate(
    phenotype_group = case_when(
      ibs_status == "Control" & bowel_habit == "Normal" ~ "Control",
      ibs_status == "IBS" & bowel_habit == "Constipation" ~ "IBS_C",
      ibs_status == "IBS" & bowel_habit == "Diarrhea" ~ "IBS_D",
      ibs_status == "IBS" & bowel_habit == "Mixed" ~ "IBS_M",
      ibs_status == "IBS" & bowel_habit == "Unspecified" ~ "IBS_U",
      TRUE ~ NA_character_
    )
  )

cat("\nPhenotype group counts before filtering:\n")
print(table(meta2$phenotype_group, useNA = "ifany"))

# Main analysis: exclude IBS_U / Unspecified
meta_out <- meta2 %>%
  filter(phenotype_group %in% c("Control", "IBS_C", "IBS_D", "IBS_M")) %>%
  select(
    sample,
    phenotype_group,
    ibs_status,
    bowel_habit,
    age,
    sex,
    bmi,
    diet_original,
    diet_grouped,
    ethnicity_grouped
  )

species_out <- species %>%
  filter(sample %in% meta_out$sample)

# Same order
species_out <- species_out[match(meta_out$sample, species_out$sample), ]

cat("\nFinal phenotype group counts:\n")
print(table(meta_out$phenotype_group))

cat("\nFinal metadata dimensions:\n")
print(dim(meta_out))

cat("\nFinal species dimensions:\n")
print(dim(species_out))

write.table(
  meta_out,
  out_metadata_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  species_out,
  out_species_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nSaved metadata:\n")
cat(out_metadata_file, "\n")

cat("\nSaved species table:\n")
cat(out_species_file, "\n")
