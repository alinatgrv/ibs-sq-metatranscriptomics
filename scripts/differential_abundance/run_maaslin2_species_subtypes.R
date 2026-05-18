#!/usr/bin/env Rscript

# ============================================================
# MaAsLin2 species-level subtype analysis
#
# Question:
#   Are species-level taxa associated with IBS subtype
#   compared with Control?
#
# Model:
#   species_abundance ~ phenotype_group + age + sex + bmi + diet_grouped + ethnicity_grouped
#
# phenotype_group:
#   Control
#   IBS_C
#   IBS_D
#   IBS_M
#
# Input:
#   results/metaphlan_metatranscriptome/joined/species_305_for_maaslin2_subtypes.tsv
#   metadata/metadata_305_for_maaslin2_species_subtypes.tsv
#
# Output:
#   results/metaphlan_metatranscriptome/maaslin2_species_305_subtypes_tss_log_diet_ethnicity_refwhite/
# ============================================================

suppressPackageStartupMessages({
  library(Maaslin2)
})

base_dir <- "/home/alina_tgrv/beegfs/IBS_SQ"

input_data <- file.path(
  base_dir,
  "results/metaphlan_metatranscriptome/joined/species_305_for_maaslin2_subtypes.tsv"
)

input_metadata <- file.path(
  base_dir,
  "metadata/metadata_305_for_maaslin2_species_subtypes.tsv"
)

output_dir <- file.path(
  base_dir,
  "results/metaphlan_metatranscriptome/maaslin2_species_305_subtypes_tss_log_diet_ethnicity_refwhite"
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

metadata$phenotype_group <- factor(metadata$phenotype_group)
metadata$sex <- factor(metadata$sex)
metadata$diet_grouped <- factor(metadata$diet_grouped)
metadata$ethnicity_grouped <- factor(metadata$ethnicity_grouped)

# -----------------------------
# Set reference levels
# -----------------------------

if ("Control" %in% levels(metadata$phenotype_group)) {
  metadata$phenotype_group <- relevel(metadata$phenotype_group, ref = "Control")
} else {
  stop("Control is not present in phenotype_group levels.")
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

cat("\nPhenotype group counts:\n")
print(table(metadata$phenotype_group))

cat("\nSex counts:\n")
print(table(metadata$sex))

cat("\nDiet grouped counts:\n")
print(table(metadata$diet_grouped))

cat("\nEthnicity grouped counts:\n")
print(table(metadata$ethnicity_grouped))

cat("\nFinal factor levels:\n")
cat("phenotype_group:\n")
print(levels(metadata$phenotype_group))
cat("sex:\n")
print(levels(metadata$sex))
cat("diet_grouped:\n")
print(levels(metadata$diet_grouped))
cat("ethnicity_grouped:\n")
print(levels(metadata$ethnicity_grouped))

reference_vec <- c(
  "phenotype_group,Control",
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
    "phenotype_group",
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

cat("\nMaAsLin2 subtype analysis finished.\n")
cat("Main results saved to:\n")
cat(output_dir, "\n")

# -----------------------------
# Post-processing subtype results
# -----------------------------

all_results_file <- file.path(output_dir, "all_results.tsv")

if (!file.exists(all_results_file)) {
  stop("all_results.tsv was not created.")
}

res <- read.delim(all_results_file, check.names = FALSE, stringsAsFactors = FALSE)

# All phenotype_group results
subtype_res <- res[res$metadata == "phenotype_group", ]
subtype_res <- subtype_res[order(subtype_res$qval, subtype_res$pval), ]

subtype_all_file <- file.path(output_dir, "phenotype_group_all_results.tsv")
write.table(
  subtype_res,
  subtype_all_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Split by subtype
for (subtype in c("IBS_C", "IBS_D", "IBS_M")) {
  x <- subtype_res[subtype_res$value == subtype, ]
  x <- x[order(x$qval, x$pval), ]

  all_file <- file.path(output_dir, paste0("phenotype_group_", subtype, "_all_results.tsv"))
  nominal_file <- file.path(output_dir, paste0("phenotype_group_", subtype, "_nominal_p005_results.tsv"))
  q025_file <- file.path(output_dir, paste0("phenotype_group_", subtype, "_q025_results.tsv"))

  write.table(
    x,
    all_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    x[x$pval < 0.05, ],
    nominal_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    x[x$qval <= 0.25, ],
    q025_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# Significant associations count by metadata variable
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

    # Significant phenotype_group only
    sig_subtype <- sig[sig$metadata == "phenotype_group", ]
    sig_subtype_file <- file.path(output_dir, "phenotype_group_significant_results.tsv")
    write.table(
      sig_subtype,
      sig_subtype_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
}

cat("\nPost-processing completed.\n")

cat("\nSubtype results summary:\n")
for (subtype in c("IBS_C", "IBS_D", "IBS_M")) {
  x <- subtype_res[subtype_res$value == subtype, ]

  cat("\n", subtype, " vs Control:\n", sep = "")
  cat("Total tested species:", nrow(x), "\n")
  cat("Nominal p < 0.05:", sum(x$pval < 0.05), "\n")
  cat("FDR q <= 0.25:", sum(x$qval <= 0.25), "\n")

  cat("Top 10 results:\n")
  cols_to_print <- intersect(c("feature", "metadata", "value", "coef", "stderr", "pval", "qval"), colnames(x))
  print(head(x[, cols_to_print], 10))
}

cat("\nSaved subtype summary files to:\n")
cat(output_dir, "\n")
