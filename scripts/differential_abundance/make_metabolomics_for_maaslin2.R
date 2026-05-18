#!/usr/bin/env Rscript

# ============================================================
# Prepare metabolomics data for MaAsLin2
# IBS vs Control
#
# Input:
#   data/metabolomics/raw/metabolite_abundance_368.csv
#   data/metabolomics/raw/metadata_metabolomics_368.csv
#
# Output:
#   data/metabolomics/processed/metabolites_368_for_maaslin2.tsv
#   data/metabolomics/processed/metadata_368_for_maaslin2.tsv
#   data/metabolomics/processed/metabolite_name_mapping.tsv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

base_dir <- "/home/alina_tgrv/beegfs/IBS_SQ"

metabolite_file <- file.path(
  base_dir,
  "data/metabolomics/raw/metabolite_abundance_368.csv"
)

metadata_file <- file.path(
  base_dir,
  "data/metabolomics/raw/metadata_metabolomics_368.csv"
)

out_dir <- file.path(
  base_dir,
  "data/metabolomics/processed"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_metabolites <- file.path(
  out_dir,
  "metabolites_368_for_maaslin2.tsv"
)

out_metadata <- file.path(
  out_dir,
  "metadata_368_for_maaslin2.tsv"
)

out_mapping <- file.path(
  out_dir,
  "metabolite_name_mapping.tsv"
)

# -----------------------------
# Read input
# -----------------------------

metab_raw <- read.csv(
  metabolite_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

metadata_raw <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Raw metabolite table dimensions:\n")
print(dim(metab_raw))

cat("\nRaw metadata dimensions:\n")
print(dim(metadata_raw))

if (!"Metabolite" %in% colnames(metab_raw)) {
  stop("Column 'Metabolite' not found in metabolite table.")
}

if (!"Patient" %in% colnames(metadata_raw)) {
  stop("Column 'Patient' not found in metadata.")
}

# -----------------------------
# Clean metabolite names
# -----------------------------

original_metabolite_names <- metab_raw$Metabolite

clean_feature_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x <- paste0("met_", x)
  make.unique(x)
}

clean_metabolite_names <- clean_feature_name(original_metabolite_names)

mapping <- data.frame(
  original_metabolite_name = original_metabolite_names,
  maaslin2_feature_name = clean_metabolite_names,
  stringsAsFactors = FALSE
)

metab_raw$Metabolite <- clean_metabolite_names

# -----------------------------
# Transpose metabolite matrix
# -----------------------------

metab_matrix <- metab_raw
rownames(metab_matrix) <- metab_matrix$Metabolite
metab_matrix$Metabolite <- NULL

# rows = samples, columns = metabolites
metab_t <- as.data.frame(t(metab_matrix), check.names = FALSE)
metab_t$Patient <- rownames(metab_t)

# Convert all metabolite columns to numeric
metabolite_cols <- setdiff(colnames(metab_t), "Patient")

for (col in metabolite_cols) {
  metab_t[[col]] <- suppressWarnings(as.numeric(metab_t[[col]]))
}

# -----------------------------
# Match with metadata
# -----------------------------

shared_patients <- intersect(metab_t$Patient, metadata_raw$Patient)

cat("\nShared patients:\n")
print(length(shared_patients))

metab_t <- metab_t %>%
  filter(Patient %in% shared_patients)

metadata_raw <- metadata_raw %>%
  filter(Patient %in% shared_patients)

metab_t <- metab_t[match(shared_patients, metab_t$Patient), ]
metadata_raw <- metadata_raw[match(shared_patients, metadata_raw$Patient), ]

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

metadata_out <- metadata_raw %>%
  transmute(
    Patient = Patient,
    Group = clean_factor(Group),
    Age = suppressWarnings(as.numeric(Age)),
    Sex = clean_factor(Sex),
    BMI = suppressWarnings(as.numeric(BMI)),
    Race = clean_factor(Race),
    Diet_Category = clean_factor(Diet_Category),
    Diet_Pattern = clean_factor(Diet_Pattern),
    Batch_metabolomics = clean_factor(Batch_metabolomics),
    BH = clean_factor(BH)
  )

# Remove samples with missing required covariates
before_n <- nrow(metadata_out)

keep <- complete.cases(
  metadata_out[, c(
    "Group",
    "Age",
    "Sex",
    "BMI",
    "Race",
    "Diet_Category",
    "Batch_metabolomics"
  )]
)

metadata_out <- metadata_out[keep, ]

after_n <- nrow(metadata_out)

cat("\nSamples before removing missing covariates:\n")
print(before_n)

cat("\nSamples after removing missing covariates:\n")
print(after_n)

cat("\nDropped samples due to missing covariates:\n")
print(before_n - after_n)

# Keep same samples in metabolite table
metab_out <- metab_t %>%
  filter(Patient %in% metadata_out$Patient)

metab_out <- metab_out[match(metadata_out$Patient, metab_out$Patient), ]

# Remove metabolites that are zero or NA in all samples
metabolite_cols <- setdiff(colnames(metab_out), "Patient")

nonzero_metabolites <- metabolite_cols[
  colSums(metab_out[, metabolite_cols, drop = FALSE] > 0, na.rm = TRUE) > 0
]

cat("\nMetabolites before nonzero filtering:\n")
print(length(metabolite_cols))

cat("\nNon-zero metabolites retained:\n")
print(length(nonzero_metabolites))

metab_out <- metab_out %>%
  select(Patient, all_of(nonzero_metabolites))

# -----------------------------
# Write outputs
# -----------------------------

write.table(
  metab_out,
  out_metabolites,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  metadata_out,
  out_metadata,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  mapping,
  out_mapping,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nFinal metabolite table dimensions:\n")
print(dim(metab_out))

cat("\nFinal metadata dimensions:\n")
print(dim(metadata_out))

cat("\nGroup counts:\n")
print(table(metadata_out$Group))

cat("\nBH counts:\n")
print(table(metadata_out$BH))

cat("\nSex counts:\n")
print(table(metadata_out$Sex))

cat("\nDiet category counts:\n")
print(table(metadata_out$Diet_Category))

cat("\nRace counts:\n")
print(table(metadata_out$Race))

cat("\nBatch counts:\n")
print(table(metadata_out$Batch_metabolomics))

cat("\nSaved metabolite table:\n")
cat(out_metabolites, "\n")

cat("\nSaved metadata:\n")
cat(out_metadata, "\n")

cat("\nSaved metabolite name mapping:\n")
cat(out_mapping, "\n")
