suppressPackageStartupMessages({
  library(ggplot2)
})

setwd("/home/alina_tgrv/beegfs/IBS_SQ")

# =========================
# Input files
# =========================

maaslin_file <- "results/metabolomics/maaslin2_metabolites_368_group_tss_log_diet_batch/all_results.tsv"
abundance_file <- "data/metabolomics/processed/metabolites_368_for_maaslin2.tsv"
metadata_file <- "data/metabolomics/processed/metadata_368_for_maaslin2.tsv"

out_dir <- "results/metabolomics/presentation_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Helper functions
# =========================

shorten_name <- function(x, max_len = 55) {
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  ifelse(nchar(x) > max_len, paste0(substr(x, 1, max_len - 3), "..."), x)
}

classify_metabolite <- function(x) {
  y <- tolower(x)

  if (grepl("sulfate|disulfate|sulfon|sulfate", y)) {
    return("Sulfated / sulfur-related")
  }
  if (grepl("andro|estr|pregnen|steroid|cort|testoster|dhea", y)) {
    return("Steroid-related")
  }
  if (grepl("bile|cholate|deoxychol|lithochol|chenodeoxychol", y)) {
    return("Bile acid-related")
  }
  if (grepl("indole|tryptophan|kynuren|seroton", y)) {
    return("Tryptophan / indole")
  }
  if (grepl("phenol|cresol|benzoate|hippurate", y)) {
    return("Aromatic compound")
  }
  if (grepl("carnitine|acyl|oleoyl|linole|palmit|stear|fatty", y)) {
    return("Lipid-related")
  }
  if (grepl("glutamate|glycine|alanine|leucine|valine|arginine|ornithine|amino", y)) {
    return("Amino acid-related")
  }

  return("Other")
}

# =========================
# Read MaAsLin2 results
# =========================

cat("Reading MaAsLin2 results...\n")
res <- read.delim(maaslin_file, check.names = FALSE)

required_cols <- c("feature", "metadata", "value", "coef", "pval", "qval")
missing_cols <- setdiff(required_cols, colnames(res))

if (length(missing_cols) > 0) {
  stop(paste("Missing columns in MaAsLin2 results:", paste(missing_cols, collapse = ", ")))
}

# Keep IBS vs Control metabolite associations
ibs_res <- res[
  grepl("^Group$", res$metadata, ignore.case = TRUE) &
    grepl("IBS", res$value, ignore.case = TRUE),
]

if (nrow(ibs_res) == 0) {
  stop("No IBS-related rows found. Check metadata/value columns in all_results.tsv")
}

ibs_res$abs_coef <- abs(ibs_res$coef)
ibs_res$direction <- ifelse(ibs_res$coef > 0, "Higher in IBS", "Higher in Control")
ibs_res$metabolite_label <- shorten_name(ibs_res$feature)
ibs_res$metabolite_class <- sapply(ibs_res$feature, classify_metabolite)

# Main presentation selection:
# first take nominally significant metabolites, then keep top by effect size
plot_res <- ibs_res[ibs_res$pval < 0.05, ]

if (nrow(plot_res) == 0) {
  cat("No nominal p < 0.05 metabolites found; using top 12 by p-value.\n")
  plot_res <- ibs_res[order(ibs_res$pval), ][1:min(12, nrow(ibs_res)), ]
} else {
  plot_res <- plot_res[order(-plot_res$abs_coef), ]
  plot_res <- plot_res[1:min(12, nrow(plot_res)), ]
}

plot_res$significance <- ifelse(plot_res$qval <= 0.25, "FDR q ≤ 0.25", "nominal p < 0.05")
plot_res$q_label <- ifelse(
  plot_res$qval <= 0.25,
  paste0("q=", signif(plot_res$qval, 2)),
  paste0("p=", signif(plot_res$pval, 2))
)

plot_res <- plot_res[order(plot_res$coef), ]
plot_res$metabolite_label <- factor(plot_res$metabolite_label, levels = plot_res$metabolite_label)

write.table(
  plot_res,
  file = file.path(out_dir, "metabolomics_top_IBS_associations_for_presentation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================
# Barplot: top IBS-associated metabolites
# =========================

p_bar <- ggplot(
  plot_res,
  aes(x = metabolite_label, y = coef, fill = direction)
) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.35) +
  geom_text(
    aes(label = q_label),
    hjust = ifelse(plot_res$coef > 0, -0.08, 1.08),
    size = 3.0
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c(
      "Higher in IBS" = "#B5524A",
      "Higher in Control" = "#4C78A8"
    )
  ) +
  labs(
    title = "IBS-associated metabolomic shifts",
    subtitle = "Top MaAsLin2 associations adjusted for diet and batch",
    x = NULL,
    y = "MaAsLin2 coefficient",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11),
    legend.position = "top",
    legend.text = element_text(size = 10),
    plot.margin = margin(8, 35, 8, 8)
  )

ggsave(
  filename = file.path(out_dir, "fig_metabolomics_top_IBS_barplot.pdf"),
  plot = p_bar,
  width = 8.2,
  height = 5.4,
  device = cairo_pdf
)

ggsave(
  filename = file.path(out_dir, "fig_metabolomics_top_IBS_barplot.png"),
  plot = p_bar,
  width = 8.2,
  height = 5.4,
  dpi = 400
)

# =========================
# Boxplot for strongest metabolite
# Wide slide-friendly version
# =========================

cat("Preparing boxplot...\n")

abund <- read.delim(abundance_file, check.names = FALSE)
meta <- read.delim(metadata_file, check.names = FALSE)

# first column = sample ID
rownames(abund) <- abund[[1]]
abund <- abund[, -1, drop = FALSE]

rownames(meta) <- meta[[1]]

group_col <- grep("^Group$|group|IBS|status", colnames(meta), ignore.case = TRUE, value = TRUE)[1]

if (is.na(group_col)) {
  stop("Could not identify group column in metadata.")
}

# choose strongest FDR-significant metabolite if present, otherwise strongest nominal
fdr_candidates <- plot_res[plot_res$qval <= 0.25, ]

if (nrow(fdr_candidates) > 0) {
  box_feature <- fdr_candidates[order(fdr_candidates$qval, -fdr_candidates$abs_coef), "feature"][1]
} else {
  box_feature <- plot_res[order(plot_res$pval, -plot_res$abs_coef), "feature"][1]
}

if (!(box_feature %in% colnames(abund))) {
  stop(paste("Selected metabolite was not found in abundance table:", box_feature))
}

common_samples <- intersect(rownames(abund), rownames(meta))

box_df <- data.frame(
  sample = common_samples,
  Group = meta[common_samples, group_col],
  abundance = as.numeric(abund[common_samples, box_feature])
)

box_df <- box_df[!is.na(box_df$Group) & !is.na(box_df$abundance), ]

if (all(c("Control", "IBS") %in% unique(box_df$Group))) {
  box_df$Group <- factor(box_df$Group, levels = c("Control", "IBS"))
}

pseudo <- min(box_df$abundance[box_df$abundance > 0], na.rm = TRUE) / 2
box_df$log_abundance <- log10(box_df$abundance + pseudo)

box_stats <- ibs_res[ibs_res$feature == box_feature, ][1, ]

p_box <- ggplot(box_df, aes(x = log_abundance, y = Group, fill = Group)) +
  geom_boxplot(
    width = 0.48,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.45,
    alpha = 0.88
  ) +
  geom_jitter(
    aes(color = Group),
    height = 0.10,
    width = 0,
    size = 1.45,
    alpha = 0.42,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Control" = "#4C78A8",
      "IBS" = "#B5524A"
    )
  ) +
  scale_color_manual(
    values = c(
      "Control" = "#4C78A8",
      "IBS" = "#B5524A"
    )
  ) +
  labs(
    title = shorten_name(box_feature, max_len = 72),
    subtitle = paste0(
      "MaAsLin2: coef = ", signif(box_stats$coef, 3),
      ", p = ", signif(box_stats$pval, 3),
      ", q = ", signif(box_stats$qval, 3)
    ),
    x = "log10(normalized abundance + pseudocount)",
    y = NULL
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 13),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    legend.position = "none",
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(
  filename = file.path(out_dir, "fig_metabolomics_strongest_metabolite_boxplot.pdf"),
  plot = p_box,
  width = 7.2,
  height = 3.2,
  device = cairo_pdf
)

ggsave(
  filename = file.path(out_dir, "fig_metabolomics_strongest_metabolite_boxplot.png"),
  plot = p_box,
  width = 7.2,
  height = 3.2,
  dpi = 400
)