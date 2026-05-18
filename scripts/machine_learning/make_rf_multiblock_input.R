#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(42)

outdir <- "results/random_forest/input_blocks"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta_path <- "metadata/metadata_326_clean_v2.tsv"
taxa_path <- "results/metaphlan_metatranscriptome/joined/species_sample_matrix_with_groups.tsv"
sq_path <- "metadata/sq_scores.csv"
metab_path <- "data/metabolomics/processed/metabolite_abundance_368_for_maaslin2_HAD.tsv"

cat("Reading metadata...\n")
meta <- fread(meta_path)
meta[, subject_id := as.character(subject_id)]
meta[, sample := as.character(sample)]

cat("Reading metabolomics abundance...\n")
metab_df <- read.table(
  metab_path,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

metab_ids <- rownames(metab_df)

cat("Building SRR -> metabolomics bridge...\n")
bridge <- copy(meta)

bridge[, metabolomics_id := fifelse(
  subject_id %in% metab_ids,
  subject_id,
  fifelse(paste0("B", subject_id) %in% metab_ids, paste0("B", subject_id), NA_character_)
)]

bridge_mapped <- bridge[!is.na(metabolomics_id)]

cat("Mapped samples:", nrow(bridge_mapped), "\n")
print(table(bridge_mapped$ibs_status))

base <- bridge_mapped[, .(
  sample,
  subject_id,
  metabolomics_id,
  ibs_status
)]

cat("Preparing clinical block...\n")
clinical_block <- bridge_mapped[, .(
  sample,
  clin__age = as.numeric(age),
  clin__bmi = as.numeric(bmi),
  clin__sex_Female = as.integer(sex == "Female"),
  clin__sex_Male = as.integer(sex == "Male")
)]

cat("Reading taxa...\n")
taxa <- fread(taxa_path, check.names = FALSE)
taxa[, sample := as.character(sample)]

taxa_cols <- grep("^s__", names(taxa), value = TRUE)

cat("Taxa before prevalence filtering:", length(taxa_cols), "\n")

taxa_block_raw <- taxa[sample %in% base$sample, c("sample", taxa_cols), with = FALSE]

prev <- sapply(taxa_block_raw[, ..taxa_cols], function(x) {
  mean(as.numeric(x) > 0, na.rm = TRUE)
})

taxa_keep <- names(prev[prev >= 0.10])

cat("Taxa kept at prevalence >= 10%:", length(taxa_keep), "\n")

taxa_block <- taxa_block_raw[, c("sample", taxa_keep), with = FALSE]

taxa_block[, (taxa_keep) := lapply(.SD, function(x) log1p(as.numeric(x))), .SDcols = taxa_keep]

setnames(
  taxa_block,
  taxa_keep,
  paste0("tax__", make.names(taxa_keep, unique = TRUE))
)

cat("Preparing metabolomics block...\n")
metab_block <- as.data.table(metab_df, keep.rownames = "metabolomics_id")

metab_feature_cols <- setdiff(names(metab_block), "metabolomics_id")

metab_block[, (metab_feature_cols) := lapply(.SD, function(x) log1p(as.numeric(x))), .SDcols = metab_feature_cols]

metab_block <- merge(
  bridge_mapped[, .(sample, metabolomics_id)],
  metab_block,
  by = "metabolomics_id",
  all.x = TRUE,
  sort = FALSE
)

metab_block[, metabolomics_id := NULL]

cat("Preparing SQ-score block...\n")
sq <- fread(sq_path)

# remove empty index column from csv if present
if (names(sq)[1] == "" || names(sq)[1] == "V1") {
  sq[, (names(sq)[1]) := NULL]
}

sq[, sample_id := as.character(sample_id)]
sq[, pathway_model := make.names(pathway_model, unique = FALSE)]

sq_value_cols <- c("SQ_score", "coverage", "n_detected_steps", "score_raw")

sq_wide_list <- list()

for (v in sq_value_cols) {
  d <- dcast(
    sq,
    sample_id ~ pathway_model,
    value.var = v,
    fun.aggregate = mean,
    fill = 0
  )
  
  old_names <- setdiff(names(d), "sample_id")
  setnames(d, old_names, paste0("sq__", v, "__", old_names))
  
  sq_wide_list[[v]] <- d
}

sq_summary <- sq[, .(
  sq__SQ_score_max = max(SQ_score, na.rm = TRUE),
  sq__SQ_score_mean = mean(SQ_score, na.rm = TRUE),
  sq__coverage_max = max(coverage, na.rm = TRUE),
  sq__coverage_mean = mean(coverage, na.rm = TRUE),
  sq__n_models_with_detected_steps = sum(n_detected_steps > 0, na.rm = TRUE)
), by = sample_id]

sq_block <- Reduce(
  function(x, y) merge(x, y, by = "sample_id", all = TRUE, sort = FALSE),
  c(list(sq_summary), sq_wide_list)
)

setnames(sq_block, "sample_id", "sample")

cat("SQ features:", ncol(sq_block) - 1, "\n")

write_block <- function(name, blocks) {
  x <- copy(base)
  
  for (b in blocks) {
    x <- merge(x, b, by = "sample", all.x = TRUE, sort = FALSE)
  }
  
  meta_cols <- c("sample", "subject_id", "metabolomics_id", "ibs_status")
  feature_cols <- setdiff(names(x), meta_cols)
  
  for (col in feature_cols) {
    x[[col]] <- as.numeric(x[[col]])
    x[[col]][is.na(x[[col]])] <- 0
  }
  
  setcolorder(x, c(meta_cols, feature_cols))
  
  outpath <- file.path(outdir, paste0(name, ".tsv"))
  fwrite(x, outpath, sep = "\t")
  
  cat(name, ":", nrow(x), "samples x", length(feature_cols), "features\n")
}

cat("Writing RF blocks...\n")

write_block("rf_234_clinical_only", list(clinical_block))

write_block("rf_234_sq_only", list(sq_block))
write_block("rf_234_taxa_only", list(taxa_block))
write_block("rf_234_metabolites_only", list(metab_block))

write_block("rf_234_taxa_sq", list(taxa_block, sq_block))
write_block("rf_234_metabolites_sq", list(metab_block, sq_block))
write_block("rf_234_taxa_metabolites", list(taxa_block, metab_block))

write_block("rf_234_taxa_clinical", list(taxa_block, clinical_block))
write_block("rf_234_metabolites_clinical", list(metab_block, clinical_block))
write_block("rf_234_sq_clinical", list(sq_block, clinical_block))

write_block("rf_234_taxa_sq_clinical", list(taxa_block, sq_block, clinical_block))
write_block("rf_234_metabolites_sq_clinical", list(metab_block, sq_block, clinical_block))
write_block("rf_234_taxa_metabolites_clinical", list(taxa_block, metab_block, clinical_block))

write_block("rf_234_all_taxa_sq_metabolites", list(taxa_block, sq_block, metab_block))
write_block("rf_234_all_taxa_sq_metabolites_clinical", list(taxa_block, sq_block, metab_block, clinical_block))

cat("Done.\n")
