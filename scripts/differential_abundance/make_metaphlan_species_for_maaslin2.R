#!/usr/bin/env Rscript

# ============================================================
# Prepare MetaPhlAn species-level data for MaAsLin2
# IBS vs Control with diet covariate
#
# Input:
#   results/metaphlan_metatranscriptome/joined/species_sample_matrix_with_groups.tsv
#   metadata/metadata_326_clean_v2.tsv
#
# Output:
#   results/metaphlan_metatranscriptome/joined/species_326_for_maaslin2.tsv
#   metadata/metadata_326_for_maaslin2_species.tsv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

base_dir <- "/home/alina_tgrv/beegfs/IBS_SQ"

species_file <- file.path(
  base_dir,
  "results/metaphlan_metatranscriptome/joined/species_sample_matrix_with_groups.tsv"
)

metadata_file <- file.path(
  base_dir,
  "metadata/metadata_326_clean_v2.tsv"
)

out_species_file <- file.path(
  base_dir,
  "results/metaphlan_metatranscriptome/joined/species_326_for_maaslin2.tsv"
)

out_metadata_file <- file.path(
  base_dir,
  "metadata/metadata_326_for_maaslin2_species.tsv"
)

species_df <- read.delim(species_file, check.names = FALSE, stringsAsFactors = FALSE)
metadata_df <- read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)

cat("Species table dimensions:\n")
print(dim(species_df))

cat("Metadata dimensions:\n")
print(dim(metadata_df))

cat("\nSpecies table columns preview:\n")
print(head(colnames(species_df), 10))

cat("\nMetadata columns:\n")
print(colnames(metadata_df))

# -----------------------------
# Basic checks
# -----------------------------

if (!"sample" %in% colnames(species_df)) {
  stop("Column 'sample' not found in species table.")
}

if (!"_group" %in% colnames(species_df)) {
  stop("Column '_group' not found in species table.")
}

if (!"sample" %in% colnames(metadata_df)) {
  stop("Column 'sample' not found in metadata.")
}

required_metadata_cols <- c("sample", "ibs_status", "age", "sex", "bmi", "diet", "ethnicity")
missing_cols <- setdiff(required_metadata_cols, colnames(metadata_df))

if (length(missing_cols) > 0) {
  stop(paste("Missing metadata columns:", paste(missing_cols, collapse = ", ")))
}

# -----------------------------
# Keep only shared samples
# -----------------------------

shared_samples <- intersect(species_df$sample, metadata_df$sample)

cat("\nShared samples:\n")
print(length(shared_samples))

species_df <- species_df %>%
  filter(sample %in% shared_samples)

metadata_df <- metadata_df %>%
  filter(sample %in% shared_samples)

# Order both tables identically
species_df <- species_df[match(shared_samples, species_df$sample), ]
metadata_df <- metadata_df[match(shared_samples, metadata_df$sample), ]

# Check group consistency
group_check <- data.frame(
  sample = species_df$sample,
  group_from_species_table = species_df$`_group`,
  ibs_status_from_metadata = metadata_df$ibs_status
)

mismatch <- group_check[group_check$group_from_species_table != group_check$ibs_status_from_metadata, ]

cat("\nGroup mismatch count:\n")
print(nrow(mismatch))

if (nrow(mismatch) > 0) {
  cat("\nGroup mismatches preview:\n")
  print(head(mismatch, 20))
  stop("Group labels in species table and metadata do not match.")
}

# -----------------------------
# Prepare metadata
# -----------------------------

clean_factor <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA"] <- "Unknown"
  x <- gsub("[ /-]+", "_", x)
  x <- gsub("[^A-Za-z0-9_]", "", x)
  x
}

metadata_out <- metadata_df %>%
  transmute(
    sample = sample,
    ibs_status = ibs_status,
    age = suppressWarnings(as.numeric(age)),
    sex = clean_factor(sex),
    bmi = suppressWarnings(as.numeric(bmi)),
    diet_original = diet,
    diet_grouped = case_when(
      diet %in% c("Standard_American", "Modified_American") ~ "Standard",
      diet %in% c("FODMAPS", "Gluten_free", "Lactose_free") ~ "Restrictive",
      diet %in% c("Mediterranean", "Paleo", "Pescetarian", "Vegan", "Vegetarian") ~ "Other",
      diet %in% c("Unknown", "", NA) ~ "Unknown",
      TRUE ~ "Other"
    ),
    ethnicity_grouped = clean_factor(ethnicity)
  )

# Make sure IBS status has expected labels
metadata_out$ibs_status <- as.character(metadata_out$ibs_status)
metadata_out$ibs_status[metadata_out$ibs_status %in% c("HC", "Healthy", "healthy", "control")] <- "Control"

# Remove samples with missing numeric covariates
before_n <- nrow(metadata_out)

keep <- complete.cases(metadata_out[, c("ibs_status", "age", "sex", "bmi", "diet_grouped", "ethnicity_grouped")])

metadata_out <- metadata_out[keep, ]

after_n <- nrow(metadata_out)

cat("\nSamples before removing missing covariates:\n")
print(before_n)

cat("\nSamples after removing missing covariates:\n")
print(after_n)

cat("\nDropped samples due to missing covariates:\n")
print(before_n - after_n)

# -----------------------------
# Prepare species feature table
# -----------------------------

species_cols <- colnames(species_df)[grepl("^s__", colnames(species_df))]

cat("\nSpecies columns detected:\n")
print(length(species_cols))

species_out <- species_df %>%
  filter(sample %in% metadata_out$sample) %>%
  select(sample, all_of(species_cols))

# Order again
species_out <- species_out[match(metadata_out$sample, species_out$sample), ]

# Remove species that are zero in all remaining samples
nonzero_species <- species_cols[colSums(species_out[, species_cols, drop = FALSE] > 0, na.rm = TRUE) > 0]

cat("\nNon-zero species retained:\n")
print(length(nonzero_species))

species_out <- species_out %>%
  select(sample, all_of(nonzero_species))

# -----------------------------
# Write outputs
# -----------------------------

write.table(
  species_out,
  out_species_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  metadata_out,
  out_metadata_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nFinal metadata dimensions:\n")
print(dim(metadata_out))

cat("\nFinal species table dimensions:\n")
print(dim(species_out))

cat("\nIBS status counts:\n")
print(table(metadata_out$ibs_status))

cat("\nSex counts:\n")
print(table(metadata_out$sex))

cat("\nDiet grouped counts:\n")
print(table(metadata_out$diet_grouped))

cat("\nEthnicity grouped counts:\n")
print(table(metadata_out$ethnicity_grouped))

cat("\nSaved species table:\n")
cat(out_species_file, "\n")

cat("\nSaved metadata table:\n")
cat(out_metadata_file, "\n")
