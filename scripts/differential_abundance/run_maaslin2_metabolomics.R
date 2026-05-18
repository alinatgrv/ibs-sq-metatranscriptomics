#!/usr/bin/env Rscript

# ============================================================
# Run MaAsLin2 for metabolomics
# Model:
#   metabolite_abundance ~ Group + Age + Sex + BMI + Race + Diet_Category + Batch_metabolomics
#
# Main comparison:
#   IBS vs Control
# ============================================================

suppressPackageStartupMessages({
  library(Maaslin2)
})

base_dir <- "/home/alina_tgrv/beegfs/IBS_SQ"

input_data <- file.path(
  base_dir,
  "data/metabolomics/processed/metabolites_368_for_maaslin2.tsv"
)

input_metadata <- file.path(
  base_dir,
  "data/metabolomics/processed/metadata_368_for_maaslin2.tsv"
)

output_dir <- file.path(
  base_dir,
  "results/metabolomics/maaslin2_metabolites_368_group_tss_log_diet_batch"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- read.delim(input_data, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.delim(input_metadata, check.names = FALSE, stringsAsFactors = FALSE)

cat("Input data dimensions before rownames:\n")
print(dim(data))

cat("\nInput metadata dimensions before rownames:\n")
print(dim(metadata))

rownames(data) <- data$Patient
data$Patient <- NULL

rownames(metadata) <- metadata$Patient
metadata$Patient <- NULL

# Match samples
shared_samples <- intersect(rownames(data), rownames(metadata))

data <- data[shared_samples, , drop = FALSE]
metadata <- metadata[shared_samples, , drop = FALSE]

cat("\nInput data dimensions after matching:\n")
print(dim(data))

cat("\nInput metadata dimensions after matching:\n")
print(dim(metadata))

cat("\nGroup counts:\n")
print(table(metadata$Group))

cat("\nSex counts:\n")
print(table(metadata$Sex))

cat("\nDiet category counts:\n")
print(table(metadata$Diet_Category))

cat("\nRace counts:\n")
print(table(metadata$Race))

cat("\nBatch counts:\n")
print(table(metadata$Batch_metabolomics))

# Convert categorical variables to factors
metadata$Group <- factor(metadata$Group)
metadata$Sex <- factor(metadata$Sex)
metadata$Race <- factor(metadata$Race)
metadata$Diet_Category <- factor(metadata$Diet_Category)
metadata$Batch_metabolomics <- factor(metadata$Batch_metabolomics)

# Set references
if ("Control" %in% levels(metadata$Group)) {
  metadata$Group <- relevel(metadata$Group, ref = "Control")
}

if ("Female" %in% levels(metadata$Sex)) {
  metadata$Sex <- relevel(metadata$Sex, ref = "Female")
}

if ("Non_Hispanic_White" %in% levels(metadata$Race)) {
  metadata$Race <- relevel(metadata$Race, ref = "Non_Hispanic_White")
}

if ("Standard" %in% levels(metadata$Diet_Category)) {
  metadata$Diet_Category <- relevel(metadata$Diet_Category, ref = "Standard")
}

if ("One" %in% levels(metadata$Batch_metabolomics)) {
  metadata$Batch_metabolomics <- relevel(metadata$Batch_metabolomics, ref = "One")
}

cat("\nFinal factor levels:\n")

cat("Group:\n")
print(levels(metadata$Group))

cat("Sex:\n")
print(levels(metadata$Sex))

cat("Race:\n")
print(levels(metadata$Race))

cat("Diet_Category:\n")
print(levels(metadata$Diet_Category))

cat("Batch_metabolomics:\n")
print(levels(metadata$Batch_metabolomics))

reference_vec <- c(
  paste0("Group,", levels(metadata$Group)[1]),
  paste0("Sex,", levels(metadata$Sex)[1]),
  paste0("Race,", levels(metadata$Race)[1]),
  paste0("Diet_Category,", levels(metadata$Diet_Category)[1]),
  paste0("Batch_metabolomics,", levels(metadata$Batch_metabolomics)[1])
)

cat("\nReference vector:\n")
print(reference_vec)

fit_data <- Maaslin2(
  input_data = data,
  input_metadata = metadata,
  output = output_dir,

  fixed_effects = c(
    "Group",
    "Age",
    "Sex",
    "BMI",
    "Race",
    "Diet_Category",
    "Batch_metabolomics"
  ),

  reference = reference_vec,

  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM",

  min_prevalence = 0.1,
  min_abundance = 0.0,

  correction = "BH",
  standardize = TRUE,

  plot_heatmap = TRUE,
  plot_scatter = TRUE
)

cat("\nMaAsLin2 metabolomics finished.\n")
cat("Main results saved to:\n")
cat(output_dir, "\n")

# -----------------------------
# Post-processing: Group IBS results
# -----------------------------

all_results_file <- file.path(output_dir, "all_results.tsv")

if (file.exists(all_results_file)) {
  res <- read.delim(all_results_file, check.names = FALSE, stringsAsFactors = FALSE)

  group_res <- res[res$metadata == "Group", ]
  group_res <- group_res[order(group_res$qval, group_res$pval), ]

  group_all_file <- file.path(output_dir, "group_IBS_all_results.tsv")
  group_nominal_file <- file.path(output_dir, "group_IBS_nominal_p005_results.tsv")
  group_q025_file <- file.path(output_dir, "group_IBS_q025_results.tsv")

  write.table(
    group_res,
    group_all_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    group_res[group_res$pval < 0.05, ],
    group_nominal_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    group_res[group_res$qval <= 0.25, ],
    group_q025_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("\nPost-processing completed.\n")
  cat("\nGroup IBS results:\n")
  cat("Total tested metabolites:", nrow(group_res), "\n")
  cat("Nominal p < 0.05:", sum(group_res$pval < 0.05, na.rm = TRUE), "\n")
  cat("FDR q <= 0.25:", sum(group_res$qval <= 0.25, na.rm = TRUE), "\n")

  cat("\nTop 20 Group IBS results:\n")
  print(head(group_res[, c("feature", "metadata", "value", "coef", "stderr", "N", "N.not.0", "pval", "qval")], 20))

  cat("\nSaved summary files:\n")
  cat(group_all_file, "\n")
  cat(group_nominal_file, "\n")
  cat(group_q025_file, "\n")
}
