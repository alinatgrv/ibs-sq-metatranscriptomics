#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
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

abund <- read.csv(raw_abund_file, check.names = FALSE)
metadata <- read.csv(metadata_file, check.names = FALSE)
sig <- read.delim(sig_file, check.names = FALSE)

# Two sulfur/sulfate-related candidates from significant MaAsLin2 results
target_pattern <- "phenol sulfate|androstenediol.*disulfate"

targets <- sig %>%
  filter(str_detect(str_to_lower(original_metabolite_name), target_pattern)) %>%
  mutate(
    label = paste0(
      original_metabolite_name,
      "\ncoef=", signif(coef, 3),
      ", q=", signif(qval, 3)
    )
  )

cat("Detected sulfur/sulfate-related significant metabolites:\n")
print(targets[, c("original_metabolite_name", "coef", "pval", "qval", "direction_in_IBS")])

if (nrow(targets) == 0) {
  stop("No sulfur/sulfate-related target metabolites were found in significant results.")
}

plot_df <- abund %>%
  rename(metabolite_name = Metabolite) %>%
  filter(metabolite_name %in% targets$original_metabolite_name) %>%
  pivot_longer(
    cols = -metabolite_name,
    names_to = "Patient",
    values_to = "abundance"
  ) %>%
  left_join(metadata %>% select(Patient, Group), by = "Patient") %>%
  left_join(
    targets %>% select(
      metabolite_name = original_metabolite_name,
      coef,
      pval,
      qval,
      direction_in_IBS,
      label
    ),
    by = "metabolite_name"
  ) %>%
  mutate(
    Group = case_when(
      Group %in% c("HC", "Healthy Control", "Healthy_Control", "Control", "control") ~ "Control",
      Group %in% c("IBS", "ibs") ~ "IBS",
      TRUE ~ as.character(Group)
    ),
    Group = factor(Group, levels = c("Control", "IBS")),
    log_abundance = log2(abundance + 1)
  )

p <- ggplot(plot_df, aes(x = Group, y = log_abundance)) +
  geom_violin(trim = FALSE, alpha = 0.35) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, alpha = 0.45, size = 1) +
  facet_wrap(~ label, scales = "free_y") +
  labs(
    title = "Sulfur/sulfate-related metabolites associated with IBS",
    subtitle = "FDR-significant metabolites from MaAsLin2 IBS vs Control analysis",
    x = NULL,
    y = "log2(abundance + 1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "IBS_vs_Control_sulfur_sulfate_related_metabolites_boxplots.pdf"),
  p,
  width = 10,
  height = 5
)

ggsave(
  file.path(out_dir, "IBS_vs_Control_sulfur_sulfate_related_metabolites_boxplots.png"),
  p,
  width = 10,
  height = 5,
  dpi = 300
)

write.table(
  targets,
  file.path(out_dir, "sulfur_sulfate_related_significant_metabolites.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Done. Figure saved to:\n")
cat(file.path(out_dir, "IBS_vs_Control_sulfur_sulfate_related_metabolites_boxplots.png"), "\n")
