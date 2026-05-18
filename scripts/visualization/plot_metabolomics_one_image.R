#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

# =========================
# Paths
# =========================

project_dir <- "/home/alina_tgrv/beegfs/IBS_SQ"
setwd(project_dir)

res_file <- "results/metabolomics/maaslin2_metabolites_368_group_tss_log_diet_batch/all_results.tsv"
abundance_file <- "data/metabolomics/processed/metabolites_368_for_maaslin2.tsv"
metadata_file <- "data/metabolomics/processed/metadata_368_for_maaslin2.tsv"
mapping_file <- "data/metabolomics/processed/metabolite_name_mapping.tsv"

out_dir <- "results/metabolomics/presentation_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_png <- file.path(out_dir, "fig_metabolomics_top_IBS_associations_with_boxplot_inset.png")

# =========================
# Helper functions
# =========================

read_table_auto <- function(file) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  first_line <- readLines(file, n = 1)
  sep <- ifelse(grepl("\t", first_line), "\t", ",")
  read.table(
    file,
    header = TRUE,
    sep = sep,
    quote = "\"",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

short_metabolite_name <- function(x, width = 34) {
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  vapply(
    x,
    function(z) paste(strwrap(z, width = width), collapse = "\n"),
    character(1)
  )
}

get_label_mapping <- function(mapping_file, features) {
  labels <- setNames(short_metabolite_name(features), features)
  
  if (!file.exists(mapping_file)) {
    return(labels)
  }
  
  mp <- read_table_auto(mapping_file)
  if (ncol(mp) < 2) {
    return(labels)
  }
  
  # Find mapping key column by maximum overlap with MaAsLin2 feature names
  overlaps <- sapply(mp, function(col) sum(as.character(col) %in% features))
  key_col <- names(mp)[which.max(overlaps)]
  
  if (max(overlaps) == 0) {
    return(labels)
  }
  
  candidate_label_cols <- setdiff(names(mp), key_col)
  preferred <- grep("original|metabolite|name|label", candidate_label_cols, ignore.case = TRUE, value = TRUE)
  
  label_col <- if (length(preferred) > 0) preferred[1] else candidate_label_cols[1]
  
  key <- as.character(mp[[key_col]])
  val <- as.character(mp[[label_col]])
  
  mapped <- val[match(features, key)]
  mapped[is.na(mapped) | mapped == ""] <- features[is.na(mapped) | mapped == ""]
  
  setNames(short_metabolite_name(mapped), features)
}

detect_group_column <- function(meta) {
  candidates <- c("Group", "group", "ibs_status", "IBS_status", "gastrointest_disord")
  found <- intersect(candidates, colnames(meta))
  if (length(found) > 0) {
    return(found[1])
  }
  
  idx <- grep("group|ibs|gastro", colnames(meta), ignore.case = TRUE)
  if (length(idx) == 0) {
    stop("Could not detect group column in metadata.")
  }
  colnames(meta)[idx[1]]
}

detect_sample_column <- function(df, sample_values = NULL) {
  if (!is.null(sample_values)) {
    scores <- sapply(df, function(x) sum(as.character(x) %in% sample_values))
    if (max(scores) > 0) {
      return(names(df)[which.max(scores)])
    }
  }
  names(df)[1]
}

# =========================
# Read MaAsLin2 results
# =========================

cat("Reading MaAsLin2 results...\n")
res <- read_table_auto(res_file)

required_cols <- c("feature", "metadata", "value", "coef", "pval", "qval")
missing_cols <- setdiff(required_cols, colnames(res))
if (length(missing_cols) > 0) {
  stop("Missing required columns in MaAsLin2 results: ", paste(missing_cols, collapse = ", "))
}

res$coef <- as.numeric(res$coef)
res$pval <- as.numeric(res$pval)
res$qval <- as.numeric(res$qval)

# Keep only IBS vs Control effect
res_group <- res[
  grepl("group|ibs", res$metadata, ignore.case = TRUE) &
    grepl("IBS", res$value, ignore.case = TRUE),
]

if (nrow(res_group) == 0) {
  stop("No IBS group rows found in MaAsLin2 results.")
}

sig <- res_group[!is.na(res_group$qval) & res_group$qval <= 0.25, ]

if (nrow(sig) == 0) {
  stop("No FDR-significant metabolite associations found at q <= 0.25.")
}

# Select top metabolites for a clean presentation figure
sig <- sig[order(sig$qval, sig$pval, -abs(sig$coef)), ]
top_n <- min(10, nrow(sig))
plot_df <- sig[seq_len(top_n), ]

label_map <- get_label_mapping(mapping_file, plot_df$feature)
plot_df$metabolite_label <- unname(label_map[plot_df$feature])

plot_df$direction <- ifelse(plot_df$coef > 0, "Higher in IBS", "Higher in Control")
plot_df$q_label <- paste0("q=", formatC(plot_df$qval, format = "g", digits = 2))

# Order bars by coefficient
plot_df$metabolite_label <- factor(
  plot_df$metabolite_label,
  levels = plot_df$metabolite_label[order(plot_df$coef)]
)

max_abs_coef <- max(abs(plot_df$coef), na.rm = TRUE)
x_pad <- max_abs_coef * 0.08

plot_df$label_x <- ifelse(plot_df$coef >= 0, plot_df$coef + x_pad, plot_df$coef - x_pad)
plot_df$hjust_q <- ifelse(plot_df$coef >= 0, 0, 1)

x_min <- min(plot_df$coef, plot_df$label_x, na.rm = TRUE) - max_abs_coef * 0.15
x_max <- max(plot_df$coef, plot_df$label_x, na.rm = TRUE) + max_abs_coef * 0.25

top_feature <- plot_df$feature[1]
top_label <- as.character(plot_df$metabolite_label[1])

# =========================
# Read abundance + metadata for inset boxplot
# =========================

cat("Reading abundance and metadata tables...\n")
abund <- read_table_auto(abundance_file)
meta <- read_table_auto(metadata_file)

# Case 1: samples are rows, metabolites are columns
if (top_feature %in% colnames(abund)) {
  sample_col_abund <- names(abund)[1]
  
  box_df <- data.frame(
    sample = as.character(abund[[sample_col_abund]]),
    abundance = as.numeric(abund[[top_feature]]),
    stringsAsFactors = FALSE
  )
  
# Case 2: metabolites are rows, samples are columns
} else if (top_feature %in% as.character(abund[[1]])) {
  row_idx <- which(as.character(abund[[1]]) == top_feature)[1]
  sample_names <- colnames(abund)[-1]
  
  box_df <- data.frame(
    sample = sample_names,
    abundance = as.numeric(abund[row_idx, -1]),
    stringsAsFactors = FALSE
  )
  
} else {
  stop("Top metabolite feature was not found in abundance table: ", top_feature)
}

sample_col_meta <- detect_sample_column(meta, box_df$sample)
group_col <- detect_group_column(meta)

meta_small <- meta[, c(sample_col_meta, group_col)]
colnames(meta_small) <- c("sample", "Group")
meta_small$sample <- as.character(meta_small$sample)

box_df <- merge(box_df, meta_small, by = "sample")
box_df <- box_df[box_df$Group %in% c("Control", "IBS"), ]

box_df$Group <- factor(box_df$Group, levels = c("Control", "IBS"))
box_df$log_abundance <- log10(box_df$abundance + 1)

# =========================
# Plot style
# =========================

group_colors <- c(
  "Higher in Control" = "#4C78A8",
  "Higher in IBS" = "#C44E52"
)

box_colors <- c(
  "Control" = "#4C78A8",
  "IBS" = "#C44E52"
)

base_theme <- theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11, color = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    plot.caption = element_text(size = 9, color = "grey30", hjust = 0)
  )

# =========================
# Main barplot
# =========================

p_bar <- ggplot(plot_df, aes(x = coef, y = metabolite_label, fill = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey40") +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.2) +
  geom_text(
    aes(x = label_x, label = q_label, hjust = hjust_q),
    size = 3.1,
    color = "grey20"
  ) +
  scale_fill_manual(values = group_colors) +
  coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
  labs(
    title = "Metabolomic shifts associated with IBS",
    subtitle = paste0("Top ", top_n, " FDR-significant MaAsLin2 associations"),
    x = "MaAsLin2 coefficient: IBS vs Control",
    y = NULL,
    caption = "Positive coefficient: higher metabolite abundance in IBS. Model adjusted for diet and batch. FDR threshold: q ≤ 0.25."
  ) +
  base_theme +
  theme(
    axis.text.y = element_text(size = 10),
    plot.margin = margin(8, 35, 8, 8)
  )

# =========================
# Inset boxplot
# =========================

p_box <- ggplot(box_df, aes(x = Group, y = log_abundance, fill = Group)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    color = "grey25",
    linewidth = 0.35
  ) +
  geom_jitter(
    width = 0.12,
    size = 0.8,
    alpha = 0.55,
    color = "grey20"
  ) +
  scale_fill_manual(values = box_colors) +
  labs(
    title = "Top association",
    subtitle = gsub("\n", " ", top_label),
    x = NULL,
    y = "log10 abundance + 1"
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 8),
    plot.subtitle = element_text(size = 6.5),
    axis.title.y = element_text(size = 7),
    axis.text = element_text(size = 7, color = "black"),
    legend.position = "none",
    plot.margin = margin(3, 3, 3, 3)
  )

# =========================
# Save one final image
# =========================

cat("Saving final image...\n")

png(out_png, width = 10, height = 6.2, units = "in", res = 450)

print(p_bar)

print(
  p_box,
  vp = viewport(
    x = 0.76,
    y = 0.29,
    width = 0.30,
    height = 0.38
  ),
  newpage = FALSE
)

dev.off()

cat("\nDone.\n")
cat("Saved image:\n")
cat(out_png, "\n")

