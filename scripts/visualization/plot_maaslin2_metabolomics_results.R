#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(readr)
})

BASE <- "/home/alina_tgrv/beegfs/IBS_SQ"

raw_abund_file <- file.path(BASE, "data/metabolomics/raw/metabolite_abundance_368.csv")
metadata_file  <- file.path(BASE, "data/metabolomics/raw/metadata_metabolomics_368.csv")

maaslin_dir <- file.path(
  BASE,
  "results/metabolomics/maaslin2_metabolites_368_group_tss_log_diet_batch"
)

sig_file <- file.path(maaslin_dir, "group_IBS_q025_results_with_original_names.tsv")

out_dir <- file.path(BASE, "results/metabolomics/figures_maaslin2_metabolites")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. Read input
# -----------------------------

abund_raw <- read.csv(raw_abund_file, check.names = FALSE)
metadata <- read.csv(metadata_file, check.names = FALSE)
sig <- read.delim(sig_file, check.names = FALSE)

cat("Abundance columns:\n")
print(colnames(abund_raw)[1:min(10, ncol(abund_raw))])

cat("Metadata columns:\n")
print(colnames(metadata))

cat("Significant results columns:\n")
print(colnames(sig))

# -----------------------------
# 2. Detect columns robustly
# -----------------------------

# metabolite abundance table: first column should be metabolite name
met_name_col <- colnames(abund_raw)[1]

# metadata sample column
sample_col <- intersect(
  c("sample", "Sample", "sample_id", "SampleID", "ID", "subject_id", "Patient"),
  colnames(metadata)
)[1]

if (is.na(sample_col)) {
  stop("Could not detect sample column in metadata. Check metadata column names.")
}

# group column
group_col <- intersect(
  c("Group", "group", "ibs_status", "IBS_status", "Diagnosis", "diagnosis"),
  colnames(metadata)
)[1]

if (is.na(group_col)) {
  stop("Could not detect group column in metadata. Check metadata column names.")
}

# original metabolite name column in MaAsLin2 result
orig_name_col <- intersect(
  c("original_metabolite_name", "original_name", "Original_name", "metabolite_original_name", "metabolite", "name"),
  colnames(sig)
)[1]

if (is.na(orig_name_col)) {
  message("No original metabolite name column detected. Will use feature names.")
  orig_name_col <- "feature"
}

# coefficient / p / q columns
coef_col <- intersect(c("coef", "coefficient", "Coefficient"), colnames(sig))[1]
p_col    <- intersect(c("pval", "p.value", "p_value", "pvalue"), colnames(sig))[1]
q_col    <- intersect(c("qval", "q.value", "q_value", "qvalue"), colnames(sig))[1]

if (is.na(coef_col) | is.na(p_col) | is.na(q_col)) {
  stop("Could not detect coef/p/q columns in MaAsLin2 results.")
}

# -----------------------------
# 3. Prepare significant metabolite table
# -----------------------------

sig2 <- sig %>%
  mutate(
    metabolite_name = .data[[orig_name_col]],
    coef = .data[[coef_col]],
    pval = .data[[p_col]],
    qval = .data[[q_col]],
    direction = if_else(coef > 0, "Higher in IBS", "Lower in IBS")
  ) %>%
  arrange(coef)

write.table(
  sig2,
  file.path(out_dir, "significant_metabolites_for_plotting.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Number of significant Group-associated metabolites:", nrow(sig2), "\n")

# -----------------------------
# 4. Prepare abundance matrix
# -----------------------------

abund_long <- abund_raw %>%
  rename(metabolite_name = all_of(met_name_col)) %>%
  pivot_longer(
    cols = -metabolite_name,
    names_to = "sample",
    values_to = "abundance"
  )

metadata2 <- metadata %>%
  rename(
    sample = all_of(sample_col),
    Group = all_of(group_col)
  ) %>%
  mutate(Group = as.character(Group))

# harmonize group labels
metadata2 <- metadata2 %>%
  mutate(
    Group = case_when(
      Group %in% c("HC", "Healthy Control", "Healthy_Control", "Control", "control") ~ "Control",
      Group %in% c("IBS", "ibs") ~ "IBS",
      TRUE ~ Group
    )
  )

plot_df <- abund_long %>%
  inner_join(metadata2, by = "sample") %>%
  semi_join(sig2 %>% select(metabolite_name), by = "metabolite_name") %>%
  mutate(log_abundance = log2(abundance + 1))

# -----------------------------
# 5. Heatmap: mean log abundance by group
# -----------------------------

mean_df <- plot_df %>%
  group_by(metabolite_name, Group) %>%
  summarise(mean_log_abundance = mean(log_abundance, na.rm = TRUE), .groups = "drop") %>%
  left_join(sig2 %>% select(metabolite_name, coef, qval, direction), by = "metabolite_name") %>%
  group_by(metabolite_name) %>%
  mutate(
    scaled_mean = as.numeric(scale(mean_log_abundance))
  ) %>%
  ungroup() %>%
  mutate(
    metabolite_name = fct_reorder(metabolite_name, coef)
  )

p_heatmap <- ggplot(mean_df, aes(x = Group, y = metabolite_name, fill = scaled_mean)) +
  geom_tile(color = "white", linewidth = 0.2) +
  labs(
    title = "FDR-significant metabolites associated with IBS",
    subtitle = "Mean log2 abundance by group; metabolites ordered by MaAsLin2 coefficient",
    x = NULL,
    y = NULL,
    fill = "Scaled\nmean"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 7),
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )

ggsave(
  file.path(out_dir, "IBS_vs_Control_significant_metabolites_heatmap.pdf"),
  p_heatmap,
  width = 7,
  height = 10
)

ggsave(
  file.path(out_dir, "IBS_vs_Control_significant_metabolites_heatmap.png"),
  p_heatmap,
  width = 7,
  height = 10,
  dpi = 300
)

# -----------------------------
# 6. Coefficient barplot
# -----------------------------

coef_df <- sig2 %>%
  mutate(
    metabolite_name = fct_reorder(metabolite_name, coef),
    q_label = paste0("q=", signif(qval, 2))
  )

p_coef <- ggplot(coef_df, aes(x = coef, y = metabolite_name, fill = direction)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Direction and effect size of significant metabolite associations",
    subtitle = "MaAsLin2 coefficient for Group IBS relative to Control",
    x = "MaAsLin2 coefficient: IBS vs Control",
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 7),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "IBS_vs_Control_significant_metabolites_coefficients.pdf"),
  p_coef,
  width = 8,
  height = 10
)

ggsave(
  file.path(out_dir, "IBS_vs_Control_significant_metabolites_coefficients.png"),
  p_coef,
  width = 8,
  height = 10,
  dpi = 300
)

# -----------------------------
# 7. SQ-related metabolites
# -----------------------------
# Edit this list if names in your table differ.
sq_candidates <- c(
  "sulfoquinovose",
  "2,3-dihydroxypropane-1-sulfonate",
  "dihydroxypropanesulfonate",
  "DHPS",
  "sulfolactate",
  "3-sulfolactate",
  "isethionate",
  "sulfoacetate"
)

sq_found <- unique(plot_df$metabolite_name[
  str_detect(
    str_to_lower(plot_df$metabolite_name),
    str_c(str_to_lower(sq_candidates), collapse = "|")
  )
])

cat("SQ-related metabolites detected among significant metabolites:\n")
print(sq_found)

# If automatic search finds more than two, take first two.
# If it finds zero, stop but keep the other figures.
if (length(sq_found) >= 1) {

  sq_to_plot <- sq_found[1:min(2, length(sq_found))]

  sq_df <- plot_df %>%
    filter(metabolite_name %in% sq_to_plot) %>%
    left_join(sig2 %>% select(metabolite_name, coef, pval, qval), by = "metabolite_name") %>%
    mutate(
      label = paste0(
        metabolite_name,
        "\ncoef=", signif(coef, 3),
        ", p=", signif(pval, 3),
        ", q=", signif(qval, 3)
      )
    )

  p_sq <- ggplot(sq_df, aes(x = Group, y = log_abundance)) +
    geom_violin(trim = FALSE, alpha = 0.35) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.65) +
    geom_jitter(width = 0.12, alpha = 0.45, size = 1) +
    facet_wrap(~ label, scales = "free_y") +
    labs(
      title = "SQ-related metabolites associated with IBS",
      subtitle = "Log2-transformed metabolite abundance in Control and IBS",
      x = NULL,
      y = "log2(abundance + 1)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )

  ggsave(
    file.path(out_dir, "IBS_vs_Control_SQ_related_metabolites_boxplots.pdf"),
    p_sq,
    width = 9,
    height = 5
  )

  ggsave(
    file.path(out_dir, "IBS_vs_Control_SQ_related_metabolites_boxplots.png"),
    p_sq,
    width = 9,
    height = 5,
    dpi = 300
  )

} else {
  warning("No SQ-related metabolites were found among significant metabolites by automatic name search.")
}

cat("Done. Figures saved to:\n")
cat(out_dir, "\n")
