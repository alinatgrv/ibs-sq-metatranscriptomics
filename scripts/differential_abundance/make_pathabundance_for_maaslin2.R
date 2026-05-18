suppressPackageStartupMessages(library(data.table))

meta_file <- "/home/alina_tgrv/beegfs/IBS_SQ/metadata/metadata_311_for_maaslin2_pathways.tsv"
path_file <- "/home/alina_tgrv/beegfs/IBS_SQ/results/humann_metatranscriptome/joined_all/pathabundance_all.tsv"
out_file <- "/home/alina_tgrv/beegfs/IBS_SQ/results/humann_metatranscriptome/joined_all/pathabundance_311_unstratified_for_maaslin2.tsv"

meta <- fread(meta_file, sep = "\t", header = TRUE)
path <- fread(path_file, sep = "\t", header = TRUE)

setnames(path, 1, "feature")

# убрать special rows
path <- path[!feature %in% c("UNMAPPED", "UNINTEGRATED")]

# оставить только unstratified pathways
path <- path[!grepl("\\|", feature)]

# оставить только samples из metadata
sample_cols <- intersect(meta$sample, colnames(path))
path <- path[, c("feature", sample_cols), with = FALSE]

fwrite(path, out_file, sep = "\t", quote = FALSE, na = "NA")

cat("Saved:", out_file, "\n")
cat("Features:", nrow(path), "\n")
cat("Samples:", ncol(path) - 1, "\n")
