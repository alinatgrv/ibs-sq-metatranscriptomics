suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

# =========================
# Paths
# =========================
input_file <- "/home/alina_tgrv/beegfs/IBS_SQ/results/metaphlan_metatranscriptome/maaslin2_species_326_tss_log_diet_ethnicity_refwhite/all_results.tsv"
outdir <- "/home/alina_tgrv/beegfs/IBS_SQ/results/metaphlan_metatranscriptome/maaslin2_species_326_tss_log_diet_ethnicity_refwhite/figures_ready"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Helper: shorten taxa names
# =========================
short_taxon <- function(x) {
  x <- gsub("^s__", "", x)
  x <- gsub("_", " ", x)

  # Common shortening for readability
  x <- gsub("Lachnospiraceae", "Lachnosp.", x)
  x <- gsub("Oscillospiraceae", "Oscillosp.", x)
  x <- gsub("Clostridiaceae", "Clostrid.", x)
  x <- gsub("Mediterraneibacter", "Mediterr.", x)
  x <- gsub("Anaerobutyricum", "Anaerobut.", x)
  x <- gsub("bacterium", "bact.", x)

  # Abbreviate genus if simple Genus species pattern
  x <- ifelse(grepl("^[A-Z][a-z]+\\s+[a-z]", x),
              sub("^([A-Z])[a-z]+\\s+", "\\1. ", x),
              x)

  x
}

# =========================
# Read data
# =========================
cat("Reading:", input_file, "\n")
res <- read_tsv(input_file, show_col_types = FALSE)

cat("Columns detected:\n")
print(colnames(res))

# =========================
# Filter:
# only FDR-significant
# only non-IBS associations of interest
# =========================
plot_df <- res %>%
  filter(
    qval <= 0.25,
    metadata %in% c("sex", "diet_grouped", "ethnicity_grouped")
  ) %>%
  mutate(
    metadata_group = case_when(
      metadata == "sex" ~ "Sex",
      metadata == "diet_grouped" ~ "Diet",
      metadata == "ethnicity_grouped" ~ "Ethnicity",
      TRUE ~ metadata
    ),
    value_clean = value %>%
      str_replace_all("_", " ") %>%
      str_replace_all("Non Hispanic White", "Non-Hispanic White"),
    contrast_label = case_when(
      metadata_group == "Sex" ~ paste0("Sex: ", value_clean),
      metadata_group == "Diet" ~ paste0("Diet: ", value_clean),
      metadata_group == "Ethnicity" ~ paste0("Ethnicity: ", value_clean),
      TRUE ~ paste(metadata_group, value_clean, sep = ": ")
    ),
    species_label = short_taxon(feature),
    logq = -log10(qval),
    direction = ifelse(coef > 0, "Positive", "Negative")
  )

cat("Number of rows after filtering:", nrow(plot_df), "\n")
print(plot_df %>% select(feature, species_label, metadata, value, coef, qval))

if (nrow(plot_df) == 0) {
  stop("No FDR-significant covariate associations found.")
}

# =========================
# Order x-axis
# =========================
x_order <- plot_df %>%
  distinct(metadata, contrast_label) %>%
  mutate(
    meta_order = case_when(
      metadata == "sex" ~ 1,
      metadata == "diet_grouped" ~ 2,
      metadata == "ethnicity_grouped" ~ 3,
      TRUE ~ 99
    )
  ) %>%
  arrange(meta_order, contrast_label) %>%
  pull(contrast_label)

# =========================
# Order y-axis
# by max significance then alphabetically
# =========================
y_order <- plot_df %>%
  group_by(species_label) %>%
  summarise(
    n_assoc = n(),
    max_logq = max(logq, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(max_logq), desc(n_assoc), species_label) %>%
  pull(species_label)

plot_df <- plot_df %>%
  mutate(
    contrast_label = factor(contrast_label, levels = x_order),
    species_label = factor(species_label, levels = rev(y_order))
  )

# =========================
# Build plot
# =========================
p <- ggplot(plot_df, aes(x = contrast_label, y = species_label)) +
  geom_point(
    aes(size = logq, fill = coef),
    shape = 21, color = "black", stroke = 0.35
  ) +
  geom_text(
    aes(label = ifelse(coef > 0, "+", "−")),
    size = 4, fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = "#19B4B8",   # teal
    mid = "white",
    high = "#F17C6B",  # coral
    midpoint = 0,
    name = "MaAsLin2\ncoefficient"
  ) +
  scale_size_continuous(
    range = c(5, 12),
    name = expression(-log[10](q))
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 16)) +
  labs(
    title = "FDR-significant species associations are driven by covariates",
    subtitle = "Species-level MaAsLin2 results (q ≤ 0.25);",
    x = "Covariate contrast",
    y = "Bacterial species"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    plot.subtitle = element_text(size = 12, hjust = 0),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# =========================
# Save
# =========================
pdf_file <- file.path(outdir, "maaslin2_species_covariates_fdr_bubbleplot.pdf")
png_file <- file.path(outdir, "maaslin2_species_covariates_fdr_bubbleplot.png")

ggsave(pdf_file, p, width = 10, height = 5.8, dpi = 300)
ggsave(png_file, p, width = 10, height = 5.8, dpi = 300)

cat("\nSaved files:\n")
cat(pdf_file, "\n")
cat(png_file, "\n")