#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Maaslin2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

BASE <- "/home/alina_tgrv/beegfs/IBS_SQ"

abund_file <- file.path(BASE, "data/metabolomics/raw/metabolite_abundance_368.csv")
metadata_file <- file.path(BASE, "data/metabolomics/raw/metadata_metabolomics_368.csv")

out_dir <- file.path(
  BASE,
  "results/metabolomics/maaslin2_metabolites_368_group_log_diet_batch_HAD"
)

processed_dir <- file.path(BASE, "data/metabolomics/processed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

cat("Reading metabolomics abundance table...\n")
abund_raw <- read.csv(abund_file, check.names = FALSE)

cat("Reading metadata...\n")
metadata <- read.csv(metadata_file, check.names = FALSE)

# -----------------------------
# 1. Prepare metabolite names
# -----------------------------

metabolite_names <- abund_raw$Metabolite

safe_names <- metabolite_names %>%
  str_replace_all("[^A-Za-z0-9]+", "_") %>%
  str_replace_all("^_|_$", "") %>%
  paste0("met_", .)

# make unique if duplicated after cleaning
safe_names <- make.unique(safe_names, sep = "_")

mapping <- data.frame(
  original_metabolite_name = metabolite_names,
  feature = safe_names
)

write.table(
  mapping,
  file.path(processed_dir, "metabolite_name_mapping_HAD_run.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

abund_clean <- abund_raw
abund_clean$Metabolite <- safe_names

# MaAsLin2 wants samples in rows, features in columns
feature_table <- abund_clean
rownames(feature_table) <- feature_table$Metabolite
feature_table$Metabolite <- NULL
feature_table <- as.data.frame(t(feature_table), check.names = FALSE)

feature_table$Patient <- rownames(feature_table)

# -----------------------------
# 2. Prepare metadata
# -----------------------------

metadata_clean <- metadata %>%
  mutate(
    Group = factor(Group, levels = c("Control", "IBS")),
    Sex = factor(Sex),
    Race = factor(Race),
    Diet_Category = factor(Diet_Category),
    Batch_metabolomics = factor(Batch_metabolomics)
  )

# Keep only samples shared between abundance and metadata
shared_samples <- intersect(feature_table$Patient, metadata_clean$Patient)

cat("Shared samples:", length(shared_samples), "\n")

feature_table <- feature_table %>%
  filter(Patient %in% shared_samples) %>%
  arrange(Patient)

rownames(feature_table) <- feature_table$Patient
feature_table$Patient <- NULL

metadata_clean <- metadata_clean %>%
  filter(Patient %in% shared_samples) %>%
  arrange(Patient)

rownames(metadata_clean) <- metadata_clean$Patient
metadata_clean$Patient <- NULL

cat("Group distribution:\n")
print(table(metadata_clean$Group))

cat("Batch distribution:\n")
print(table(metadata_clean$Batch_metabolomics))

cat("HAD_Anxiety summary:\n")
print(summary(metadata_clean$HAD_Anxiety))

# -----------------------------
# 3. Save prepared input
# -----------------------------

write.table(
  feature_table,
  file.path(processed_dir, "metabolite_abundance_368_for_maaslin2_HAD.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

write.table(
  metadata_clean,
  file.path(processed_dir, "metadata_metabolomics_368_for_maaslin2_HAD.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

# -----------------------------
# 4. Run MaAsLin2
# -----------------------------

cat("Running MaAsLin2 metabolomics model with HAD_Anxiety...\n")

fit_data <- Maaslin2(
  input_data = feature_table,
  input_metadata = metadata_clean,
  output = out_dir,
  fixed_effects = c(
    "Group",
    "Age",
    "Sex",
    "BMI",
    "Race",
    "Diet_Category",
    "Batch_metabolomics",
    "HAD_Anxiety"
  ),
  reference = c(
    "Group,Control",
    "Sex,Female",
    "Diet_Category,Standard",
    "Batch_metabolomics,One"
  ),
  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM",
  min_prevalence = 0.10,
  correction = "BH",
  standardize = FALSE,
  plot_heatmap = TRUE,
  plot_scatter = TRUE
)

# -----------------------------
# 5. Extract GroupIBS results
# -----------------------------

all_results <- read.delim(file.path(out_dir, "all_results.tsv"), check.names = FALSE)

group_results <- all_results %>%
  filter(metadata == "Group", value == "IBS") %>%
  left_join(mapping, by = "feature") %>%
  mutate(
    direction_in_IBS = ifelse(coef > 0, "Higher in IBS", "Lower in IBS")
  ) %>%
  select(
    original_metabolite_name,
    feature,
    direction_in_IBS,
    coef,
    stderr,
    pval,
    qval,
    N,
    N.not.0,
    everything()
  ) %>%
  arrange(qval, pval)

write.table(
  group_results,
  file.path(out_dir, "group_IBS_all_results_with_original_names.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

group_nominal <- group_results %>%
  filter(pval < 0.05)

group_q025 <- group_results %>%
  filter(qval <= 0.25)

write.table(
  group_nominal,
  file.path(out_dir, "group_IBS_nominal_p005_results_with_original_names.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  group_q025,
  file.path(out_dir, "group_IBS_q025_results_with_original_names.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Done.\n")
cat("Total GroupIBS results:", nrow(group_results), "\n")
cat("Nominal p < 0.05:", nrow(group_nominal), "\n")
cat("FDR q <= 0.25:", nrow(group_q025), "\n")
cat("Output directory:\n")
cat(out_dir, "\n")
