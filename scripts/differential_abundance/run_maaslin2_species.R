#!/usr/bin/env Rscript

# ============================================================
# Final MaAsLin2 run for MetaPhlAn species-level taxonomy
#
# Question:
#   Which species-level taxa are associated with IBS status
#   after adjustment for age, sex, BMI, diet, and ethnicity?
#
# Model:
#   species_abundance ~ ibs_status + age + sex + bmi + diet_grouped + ethnicity_grouped
#
# Input:
#   results/metaphlan_metatranscriptome/joined/species_326_for_maaslin2.tsv
#   metadata/metadata_326_for_maaslin2_species.tsv
#
# Output:
#   results/metaphlan_metatranscriptome/maaslin2_species_326_tss_log_diet_ethnicity_refwhite/
# ============================================================

suppressPackageStartupMessages({
  library(Maaslin2)
})

base_dir <- "/home/alina_tgrv/beegfs/IBS_SQ"

input_data <- file.path(
  base_dir,
  "results/metaphlan_metatranscriptome/joined/species_326_for_maaslin2.tsv"
)

input_metadata <- file.path(
  base_dir,
  "metadata/metadata_326_for_maaslin2_species.tsv"
)

output_dir <- file.path(
  base_dir,
  "results/metaphlan_metatranscriptome/maaslin2_species_326_tss_log_diet_ethnicity_refwhite"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Read input tables
# -----------------------------

data <- read.delim(input_data, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.delim(input_metadata, check.names = FALSE, stringsAsFactors = FALSE)

cat("Input data dimensions before rownames:\n")
print(dim(data))

cat("\nInput metadata dimensions before rownames:\n")
print(dim(metadata))

# -----------------------------
# Move sample IDs to rownames
# -----------------------------

if (!"sample" %in% colnames(data)) {
  stop("Column 'sample' not found in input data.")
}

if (!"sample" %in% colnames(metadata)) {
  stop("Column 'sample' not found in metadata.")
}

rownames(data) <- data$sample
data$sample <- NULL

rownames(metadata) <- metadata$sample
metadata$sample <- NULL

# -----------------------------
# Match samples and order tables
# -----------------------------

shared_samples <- intersect(rownames(data), rownames(metadata))

data <- data[shared_samples, , drop = FALSE]
metadata <- metadata[shared_samples, , drop = FALSE]

cat("\nInput data dimensions after matching:\n")
print(dim(data))

cat("\nInput metadata dimensions after matching:\n")
print(dim(metadata))

# -----------------------------
# Convert metadata types
# -----------------------------

metadata$age <- as.numeric(metadata$age)
metadata$bmi <- as.numeric(metadata$bmi)

metadata$ibs_status <- factor(metadata$ibs_status)
metadata$sex <- factor(metadata$sex)
metadata$diet_grouped <- factor(metadata$diet_grouped)
metadata$ethnicity_grouped <- factor(metadata$ethnicity_grouped)

# -----------------------------
# Set reference levels
# -----------------------------

if ("Control" %in% levels(metadata$ibs_status)) {
  metadata$ibs_status <- relevel(metadata$ibs_status, ref = "Control")
} else {
  stop("Control is not present in ibs_status levels.")
}

if ("Female" %in% levels(metadata$sex)) {
  metadata$sex <- relevel(metadata$sex, ref = "Female")
}

if ("Standard" %in% levels(metadata$diet_grouped)) {
  metadata$diet_grouped <- relevel(metadata$diet_grouped, ref = "Standard")
}

if ("Non_Hispanic_White" %in% levels(metadata$ethnicity_grouped)) {
  metadata$ethnicity_grouped <- relevel(metadata$ethnicity_grouped, ref = "Non_Hispanic_White")
} else {
  stop("Non_Hispanic_White is not present in ethnicity_grouped levels.")
}

# -----------------------------
# Print analysis design
# -----------------------------

cat("\nIBS status counts:\n")
print(table(metadata$ibs_status))

cat("\nSex counts:\n")
print(table(metadata$sex))

cat("\nDiet grouped counts:\n")
print(table(metadata$diet_grouped))

cat("\nEthnicity grouped counts:\n")
print(table(metadata$ethnicity_grouped))

cat("\nFinal factor levels:\n")
cat("ibs_status:\n")
print(levels(metadata$ibs_status))
cat("sex:\n")
print(levels(metadata$sex))
cat("diet_grouped:\n")
print(levels(metadata$diet_grouped))
cat("ethnicity_grouped:\n")
print(levels(metadata$ethnicity_grouped))

reference_vec <- c(
  "ibs_status,Control",
  "sex,Female",
  "diet_grouped,Standard",
  "ethnicity_grouped,Non_Hispanic_White"
)

cat("\nReference vector:\n")
print(reference_vec)

# -----------------------------
# Run MaAsLin2
# -----------------------------

fit_data <- Maaslin2(
  input_data = data,
  input_metadata = metadata,
  output = output_dir,

  fixed_effects = c(
    "ibs_status",
    "age",
    "sex",
    "bmi",
    "diet_grouped",
    "ethnicity_grouped"
  ),

  reference = reference_vec,

  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",

  min_prevalence = 0.1,
  min_abundance = 0.0,

  correction = "BH",
  standardize = TRUE,

  plot_heatmap = TRUE,
  plot_scatter = TRUE
)

cat("\nMaAsLin2 finished.\n")
cat("Main results saved to:\n")
cat(output_dir, "\n")

# -----------------------------
# Post-processing summary tables
# -----------------------------

all_results_file <- file.path(output_dir, "all_results.tsv")

if (!file.exists(all_results_file)) {
  stop("all_results.tsv was not created.")
}

res <- read.delim(all_results_file, check.names = FALSE, stringsAsFactors = FALSE)

# 1. All IBS-status results
ibs <- res[res$metadata == "ibs_status", ]
ibs <- ibs[order(ibs$qval, ibs$pval), ]

ibs_file <- file.path(output_dir, "ibs_status_all_results.tsv")
write.table(
  ibs,
  ibs_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 2. Nominal IBS results p < 0.05
ibs_nominal <- ibs[ibs$pval < 0.05, ]

ibs_nominal_file <- file.path(output_dir, "ibs_status_nominal_p005_results.tsv")
write.table(
  ibs_nominal,
  ibs_nominal_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 3. Significant IBS results q <= 0.25
ibs_q025 <- ibs[ibs$qval <= 0.25, ]

ibs_q025_file <- file.path(output_dir, "ibs_status_q025_results.tsv")
write.table(
  ibs_q025,
  ibs_q025_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 4. Significant associations by metadata variable
sig_file <- file.path(output_dir, "significant_results.tsv")

if (file.exists(sig_file)) {
  sig <- read.delim(sig_file, check.names = FALSE, stringsAsFactors = FALSE)

  if (nrow(sig) > 0) {
    sig_counts <- as.data.frame(table(sig$metadata))
    colnames(sig_counts) <- c("metadata", "n_significant_associations")

    sig_counts_file <- file.path(output_dir, "significant_associations_by_metadata.tsv")
    write.table(
      sig_counts,
      sig_counts_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
}

cat("\nPost-processing completed.\n")

cat("\nIBS-status results:\n")
cat("Total IBS-tested species:", nrow(ibs), "\n")
cat("Nominal IBS results p < 0.05:", nrow(ibs_nominal), "\n")
cat("FDR-significant IBS results q <= 0.25:", nrow(ibs_q025), "\n")

cat("\nSaved summary files:\n")
cat(ibs_file, "\n")
cat(ibs_nominal_file, "\n")
cat(ibs_q025_file, "\n")

cat("\nTop 20 IBS-status results:\n")
cols_to_print <- intersect(c("feature", "metadata", "value", "coef", "stderr", "pval", "qval"), colnames(ibs))
print(head(ibs[, cols_to_print], 20))
