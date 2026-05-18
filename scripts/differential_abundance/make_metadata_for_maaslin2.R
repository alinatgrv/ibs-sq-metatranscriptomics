suppressPackageStartupMessages(library(data.table))

infile <- "/home/alina_tgrv/beegfs/IBS_SQ/results/humann_metatranscriptome/joined_all/sample_metadata_qc_326.tsv"
outfile <- "/home/alina_tgrv/beegfs/IBS_SQ/metadata/metadata_311_for_maaslin2_pathways.tsv"

meta <- fread(infile, sep = "\t", header = TRUE)

# убрать hard_low_pathways
meta <- meta[!grepl("hard_low_pathways", qc_flags)]

# привести sample IDs к формату HUMAnN pathabundance table
meta[, sample := paste0(sample, "_combined_Abundance")]

# оставить нужные колонки
meta <- meta[, .(sample, ibs_status, age, sex, bmi, diet)]

# убрать NA в ключевых ковариатах
meta <- meta[complete.cases(meta)]

fwrite(meta, outfile, sep = "\t", quote = FALSE, na = "NA")

cat("Saved:", outfile, "\n")
cat("Samples:", nrow(meta), "\n")
print(table(meta$ibs_status))
