#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Maaslin2)
  library(dplyr)
  library(stringr)
  library(readr)
})

BASE <- "/home/alina_tgrv/beegfs/IBS_SQ"

species_file <- file.path(
  BASE,
  "results/metaphlan_metatranscriptome/joined/species_sample_matrix_with_groups.tsv"
)

raw_sra_metadata_file <- file.path(
  BASE,
  "metadata/metadata_326_raw.csv"
)

metatranscriptomics_metabolomics_meta_file <- file.path(
  BASE,
  "data/metabolomics/raw/metadata_metatranscriptomics_327.csv"
)

metabolomics_abund_file <- file.path(
  BASE,
  "data/metabolomics/raw/metabolite_abundance_368.csv"
)

out_base <- file.path(
  BASE,
  "results/crossomics/species_vs_sulfur_metabolites_maaslin2"
)

input_dir <- file.path(out_base, "input")

dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

cat("Reading species table...\n")
species_raw <- read.delim(species_file, check.names = FALSE)
species_raw$sample <- as.character(species_raw$sample)

cat("Reading raw SRA metadata...\n")
raw_meta <- read.csv(raw_sra_metadata_file, check.names = FALSE)
raw_meta$Run <- as.character(raw_meta$Run)
raw_meta$`Library Name` <- as.character(raw_meta$`Library Name`)

cat("Reading metatranscriptomics-metabolomics metadata...\n")
mt_meta <- read.csv(metatranscriptomics_metabolomics_meta_file, check.names = FALSE)
mt_meta$SampleID_metatranscriptomics <- as.character(mt_meta$SampleID_metatranscriptomics)
mt_meta$Patient <- as.character(mt_meta$Patient)

cat("Reading metabolomics abundance...\n")
metab_raw <- read.csv(metabolomics_abund_file, check.names = FALSE)
metab_raw$Metabolite <- as.character(metab_raw$Metabolite)

# ------------------------------------------------------------
# 1. Prepare species feature table
# ------------------------------------------------------------

species_cols <- colnames(species_raw)[str_detect(colnames(species_raw), "^s__")]

cat("Detected species-level features:", length(species_cols), "\n")

species_features <- species_raw %>%
  select(sample, all_of(species_cols))

taxa_mapping <- data.frame(
  original_taxon_name = species_cols,
  feature = species_cols %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    make.unique(sep = "_"),
  stringsAsFactors = FALSE
)

colnames(species_features) <- c("sample", taxa_mapping$feature)

write.table(
  taxa_mapping,
  file.path(input_dir, "taxa_name_mapping.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# 2. Correct bridge:
#    SRR sample -> 1UCLA ID -> Patient ID
# ------------------------------------------------------------

bridge <- raw_meta %>%
  transmute(
    sample = as.character(Run),
    SampleID_metatranscriptomics = as.character(`Library Name`)
  ) %>%
  inner_join(mt_meta, by = "SampleID_metatranscriptomics") %>%
  filter(sample %in% species_features$sample) %>%
  mutate(
    sample = as.character(sample),
    Patient = as.character(Patient),
    ibs_status = factor(Group, levels = c("Control", "IBS")),
    sex = factor(Sex),
    age = as.numeric(Age),
    bmi = as.numeric(BMI),
    Batch_metatranscriptomics = factor(Batch_metatranscriptomics),
    HAD_Anxiety = as.numeric(HAD_Anxiety),
    HAD_Depression = as.numeric(HAD_Depression),
    STAI_Tanxiety = as.numeric(STAI_Tanxiety),
    diet_grouped = case_when(
      Diet_Category == "Standard" ~ "Standard",
      Diet_Category == "Restrictive" ~ "Restrictive",
      Diet_Category == "Other" ~ "Other",
      TRUE ~ "Unknown"
    ),
    diet_grouped = factor(
      diet_grouped,
      levels = c("Standard", "Restrictive", "Other", "Unknown")
    ),
    ethnicity_grouped = case_when(
      Race == "Non_Hispanic_White" ~ "Non_Hispanic_White",
      Race == "Asian" ~ "Asian",
      Race == "Hispanic" ~ "Hispanic",
      TRUE ~ "Other"
    ),
    ethnicity_grouped = factor(ethnicity_grouped)
  )

cat("Bridge rows:", nrow(bridge), "\n")
cat("Group distribution in bridge:\n")
print(table(bridge$ibs_status, useNA = "ifany"))

# ------------------------------------------------------------
# 3. Extract selected sulfur/sulfate-related metabolites
# ------------------------------------------------------------

target_metabolites <- c(
  "phenol sulfate",
  "androstenediol (3beta,17beta) disulfate (2)"
)

cat("Checking target metabolites in metabolomics table:\n")
print(target_metabolites)
print(target_metabolites %in% metab_raw$Metabolite)

if (!all(target_metabolites %in% metab_raw$Metabolite)) {
  missing <- target_metabolites[!target_metabolites %in% metab_raw$Metabolite]
  stop(paste("Missing target metabolites:", paste(missing, collapse = "; ")))
}

metab_targets <- metab_raw %>%
  filter(Metabolite %in% target_metabolites)

metab_t <- as.data.frame(t(metab_targets[, -1]), check.names = FALSE)
colnames(metab_t) <- metab_targets$Metabolite
metab_t$Patient <- rownames(metab_t)

metab_t <- metab_t %>%
  rename(
    phenol_sulfate = `phenol sulfate`,
    androstenediol_disulfate = `androstenediol (3beta,17beta) disulfate (2)`
  ) %>%
  mutate(
    Patient = as.character(Patient),
    phenol_sulfate = as.numeric(phenol_sulfate),
    androstenediol_disulfate = as.numeric(androstenediol_disulfate)
  )

# ------------------------------------------------------------
# 4. Merge species features + metadata + metabolites
# ------------------------------------------------------------

metadata_cross <- bridge %>%
  left_join(metab_t, by = "Patient")

cat("Rows after joining metabolomics by Patient:", nrow(metadata_cross), "\n")
cat("Non-missing phenol_sulfate:", sum(!is.na(metadata_cross$phenol_sulfate)), "\n")
cat("Non-missing androstenediol_disulfate:", sum(!is.na(metadata_cross$androstenediol_disulfate)), "\n")

metadata_cross <- metadata_cross %>%
  filter(!is.na(phenol_sulfate), !is.na(androstenediol_disulfate))

species_features2 <- species_features %>%
  filter(sample %in% metadata_cross$sample)

species_features2 <- species_features2 %>%
  arrange(sample)

metadata_cross <- metadata_cross %>%
  arrange(sample)

stopifnot(all(species_features2$sample == metadata_cross$sample))

cat("Final matched cross-omics samples:", nrow(metadata_cross), "\n")
cat("Group distribution:\n")
print(table(metadata_cross$ibs_status, useNA = "ifany"))

cat("Diet_grouped distribution:\n")
print(table(metadata_cross$diet_grouped, useNA = "ifany"))

cat("Ethnicity_grouped distribution:\n")
print(table(metadata_cross$ethnicity_grouped, useNA = "ifany"))

# ------------------------------------------------------------
# 5. Save matched input tables
# ------------------------------------------------------------

feature_table <- species_features2
rownames(feature_table) <- feature_table$sample
feature_table$sample <- NULL

metadata_table <- metadata_cross %>%
  select(
    Patient,
    SampleID_metatranscriptomics,
    ibs_status,
    age,
    sex,
    bmi,
    diet_grouped,
    ethnicity_grouped,
    Batch_metatranscriptomics,
    HAD_Anxiety,
    HAD_Depression,
    STAI_Tanxiety,
    phenol_sulfate,
    androstenediol_disulfate
  )

rownames(metadata_table) <- metadata_cross$sample

write.table(
  feature_table,
  file.path(input_dir, "species_features_matched.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

write.table(
  metadata_table,
  file.path(input_dir, "metadata_crossomics_matched.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

cat("Matched input tables saved to:\n")
cat(input_dir, "\n")

# ------------------------------------------------------------
# 6. Function to run MaAsLin2 for one metabolite
# ------------------------------------------------------------

run_one_metabolite <- function(metabolite_var) {
  
  out_dir <- file.path(out_base, paste0("species_vs_", metabolite_var))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("\n============================================================\n")
  cat("Running MaAsLin2 for metabolite:", metabolite_var, "\n")
  cat("Output:", out_dir, "\n")
  cat("============================================================\n")
  
  fixed_effects <- c(
    metabolite_var,
    "ibs_status",
    "age",
    "sex",
    "bmi",
    "diet_grouped",
    "ethnicity_grouped"
  )
  
  fit <- Maaslin2(
    input_data = feature_table,
    input_metadata = metadata_table,
    output = out_dir,
    fixed_effects = fixed_effects,
    reference = c(
      "ibs_status,Control",
      "sex,Female",
      "diet_grouped,Standard",
      "ethnicity_grouped,Non_Hispanic_White"
    ),
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",
    min_prevalence = 0.10,
    correction = "BH",
    standardize = TRUE,
    plot_heatmap = TRUE,
    plot_scatter = TRUE
  )
  
  all_results_file <- file.path(out_dir, "all_results.tsv")
  
  if (!file.exists(all_results_file)) {
    stop(paste("MaAsLin2 did not produce all_results.tsv for", metabolite_var))
  }
  
  all_results <- read.delim(all_results_file, check.names = FALSE)
  
  metab_results <- all_results %>%
    filter(metadata == metabolite_var) %>%
    left_join(taxa_mapping, by = "feature") %>%
    mutate(
      association_direction = ifelse(coef > 0, "Positive", "Negative")
    ) %>%
    select(
      original_taxon_name,
      feature,
      association_direction,
      coef,
      stderr,
      pval,
      qval,
      N,
      N.not.0,
      everything()
    ) %>%
    arrange(qval, pval)
  
  nominal <- metab_results %>%
    filter(pval < 0.05)
  
  q025 <- metab_results %>%
    filter(qval <= 0.25)
  
  write.table(
    metab_results,
    file.path(out_dir, paste0(metabolite_var, "_all_taxa_results.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  write.table(
    nominal,
    file.path(out_dir, paste0(metabolite_var, "_nominal_p005_taxa_results.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  write.table(
    q025,
    file.path(out_dir, paste0(metabolite_var, "_q025_taxa_results.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  cat("\nResults for", metabolite_var, "\n")
  cat("Total taxa tested:", nrow(metab_results), "\n")
  cat("Nominal p < 0.05:", nrow(nominal), "\n")
  cat("FDR q <= 0.25:", nrow(q025), "\n")
  
  if (nrow(q025) > 0) {
    cat("Top q <= 0.25 taxa:\n")
    print(head(q025[, c("original_taxon_name", "coef", "pval", "qval", "association_direction")], 20))
  } else {
    cat("No FDR-significant taxa for", metabolite_var, "\n")
    cat("Top nominal taxa:\n")
    print(head(nominal[, c("original_taxon_name", "coef", "pval", "qval", "association_direction")], 20))
  }
}

# ------------------------------------------------------------
# 7. Run targeted cross-omics models
# ------------------------------------------------------------

run_one_metabolite("phenol_sulfate")
run_one_metabolite("androstenediol_disulfate")

cat("\nDone. Cross-omics MaAsLin2 results saved to:\n")
cat(out_base, "\n")