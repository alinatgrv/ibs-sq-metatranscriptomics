suppressPackageStartupMessages({
  library(ggplot2)
})

# =========================
# Paths
# =========================
infile <- "/home/alina_tgrv/beegfs/IBS_SQ/results/metaphlan_metatranscriptome/maaslin2_species_326_tss_log_diet_ethnicity_refwhite/all_results.tsv"
outdir <- "results/metaphlan_metatranscriptome/presentation_figures"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(outdir, "fig_maaslin2_species_covariates_heatmap.pdf")
out_png <- file.path(outdir, "fig_maaslin2_species_covariates_heatmap.png")
out_tsv <- file.path(outdir, "fig_maaslin2_species_covariates_heatmap_plotted_data.tsv")

if (!file.exists(infile)) {
  stop("Input file not found: ", infile)
}

# =========================
# Read data
# =========================
res <- read.delim(infile, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("feature", "metadata", "value", "coef", "qval")
missing_cols <- setdiff(required_cols, colnames(res))
if (length(missing_cols) > 0) {
  stop("Missing columns in input: ", paste(missing_cols, collapse = ", "))
}

res$coef <- as.numeric(res$coef)
res$qval <- as.numeric(res$qval)

# =========================
# Keep only FDR-significant non-IBS associations
# =========================
allowed_metadata <- c("diet_grouped", "sex", "ethnicity_grouped")

sig <- res[
  !is.na(res$qval) &
    res$qval <= 0.25 &
    res$metadata %in% allowed_metadata,
]

if (nrow(sig) == 0) {
  stop("No FDR-significant associations found for diet_grouped / sex / ethnicity_grouped.")
}

# =========================
# Helper functions
# =========================
clean_value <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("Non_Hispanic_White", "Non-Hispanic White", x, fixed = TRUE)
  x <- gsub("Non Hispanic White", "Non-Hispanic White", x, fixed = TRUE)
  x
}

make_term_label <- function(metadata, value) {
  value <- clean_value(value)
  if (metadata == "diet_grouped") {
    paste0("Diet\n", value)
  } else if (metadata == "sex") {
    paste0("Sex\n", value)
  } else if (metadata == "ethnicity_grouped") {
    paste0("Ethnicity\n", value)
  } else {
    paste0(metadata, "\n", value)
  }
}

abbr_species <- function(feature) {
  x <- sub("^s__", "", feature)
  parts <- strsplit(x, "_")[[1]]
  
  if (length(parts) >= 2) {
    first <- paste0(substr(parts[1], 1, 1), ".")
    rest <- paste(parts[-1], collapse = " ")
    lab <- paste(first, rest)
  } else {
    lab <- gsub("_", " ", x)
  }
  
  lab
}

to_plotmath_italic <- function(label) {
  label <- gsub("'", "\\\\'", label)
  paste0("italic('", label, "')")
}

# =========================
# Build labels
# =========================
sig$term_label <- mapply(make_term_label, sig$metadata, sig$value, USE.NAMES = FALSE)
sig$species_short <- vapply(sig$feature, abbr_species, character(1))
sig$signed_score <- -log10(sig$qval) * sign(sig$coef)

# Order x-axis: Diet -> Sex -> Ethnicity
metadata_rank_map <- c("diet_grouped" = 1, "sex" = 2, "ethnicity_grouped" = 3)

term_info <- unique(sig[, c("term_label", "metadata", "value")])
term_info$metadata_rank <- metadata_rank_map[term_info$metadata]
term_info$value_clean <- vapply(term_info$value, clean_value, character(1))
term_info <- term_info[order(term_info$metadata_rank, term_info$value_clean), ]
term_levels <- term_info$term_label

# Order species by strongest signal
species_best_q <- tapply(sig$qval, sig$feature, min, na.rm = TRUE)
species_levels <- names(sort(species_best_q, decreasing = TRUE))   # weakest -> strongest
species_levels <- rev(species_levels)                              # strongest on top

# Full matrix-like grid (white tiles where no association)
grid <- expand.grid(
  feature = species_levels,
  term_label = term_levels,
  stringsAsFactors = FALSE
)

plotdf <- merge(
  grid,
  sig[, c("feature", "term_label", "coef", "qval", "signed_score", "species_short")],
  by = c("feature", "term_label"),
  all.x = TRUE
)

# fill missing tiles
short_map <- tapply(sig$species_short, sig$feature, function(x) x[1])
plotdf$species_short[is.na(plotdf$species_short)] <- short_map[plotdf$feature[is.na(plotdf$species_short)]]
plotdf$signed_score[is.na(plotdf$signed_score)] <- 0

# factor order
plotdf$feature <- factor(plotdf$feature, levels = species_levels)
plotdf$term_label <- factor(plotdf$term_label, levels = term_levels)

# y-axis labels as italic expressions
species_label_map <- setNames(
  vapply(species_levels, function(f) to_plotmath_italic(short_map[[f]]), character(1)),
  species_levels
)

# save plotted table
write.table(plotdf, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

# =========================
# Plot
# =========================
p <- ggplot(plotdf, aes(x = term_label, y = feature, fill = signed_score)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    limits = c(-max(abs(plotdf$signed_score)), max(abs(plotdf$signed_score))),
    name = expression(sign(beta) %.% -log[10](q))
  ) +
  scale_y_discrete(labels = function(x) parse(text = species_label_map[x])) +
  labs(
    title = "FDR-significant species associations detected by MaAsLin2",
    subtitle = "Significant associations were linked to diet, sex, and ethnicity, not IBS status",
    x = "Covariate level",
    y = "Bacterial species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 14, 10, 10)
  )

n_terms <- length(term_levels)
n_species <- length(species_levels)

fig_width <- max(6.5, 1.1 * n_terms + 3.0)
fig_height <- max(4.0, 0.5 * n_species + 2.0)

ggsave(out_pdf, p, width = fig_width, height = fig_height, units = "in")
ggsave(out_png, p, width = fig_width, height = fig_height, units = "in", dpi = 300)

cat("Done.\n")
cat("PDF: ", out_pdf, "\n", sep = "")
cat("PNG: ", out_png, "\n", sep = "")
cat("TSV: ", out_tsv, "\n", sep = "")
cat("Significant associations used: ", nrow(sig), "\n", sep = "")
cat("Species shown: ", n_species, "\n", sep = "")
cat("Columns shown: ", n_terms, "\n", sep = "")
