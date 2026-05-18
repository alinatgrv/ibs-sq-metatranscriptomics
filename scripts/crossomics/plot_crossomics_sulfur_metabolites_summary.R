#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

base_dir <- "/home/alina_tgrv/beegfs/IBS_SQ/results/crossomics/species_vs_sulfur_metabolites_maaslin2"

phenol_dir <- file.path(base_dir, "species_vs_phenol_sulfate")
andro_dir  <- file.path(base_dir, "species_vs_androstenediol_disulfate")

out_dir <- file.path(base_dir, "summary_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_results <- function(path, metabolite_label) {
  if (!file.exists(path)) {
    warning("File not found: ", path)
    return(data.frame())
  }

  df <- read.delim(
    path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )

  if (nrow(df) == 0) {
    return(data.frame())
  }

  df$metabolite <- metabolite_label

  numeric_cols <- c("coef", "stderr", "pval", "qval", "N", "N.not.0")
  for (col in numeric_cols) {
    if (col %in% colnames(df)) {
      df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
    }
  }

  if (!"original_taxon_name" %in% colnames(df) && "feature" %in% colnames(df)) {
    df$original_taxon_name <- df$feature
  }

  df$taxon_clean <- gsub("^s__", "", df$original_taxon_name)
  df$taxon_clean <- gsub("^s_", "", df$taxon_clean)
  df$taxon_clean <- gsub("_", " ", df$taxon_clean)

  df$association_direction <- ifelse(df$coef >= 0, "Positive", "Negative")

  return(df)
}

phenol_all <- read_results(
  file.path(phenol_dir, "phenol_sulfate_all_taxa_results.tsv"),
  "phenol sulfate"
)

phenol_q025 <- read_results(
  file.path(phenol_dir, "phenol_sulfate_q025_taxa_results.tsv"),
  "phenol sulfate"
)

phenol_nominal <- read_results(
  file.path(phenol_dir, "phenol_sulfate_nominal_p005_taxa_results.tsv"),
  "phenol sulfate"
)

andro_all <- read_results(
  file.path(andro_dir, "androstenediol_disulfate_all_taxa_results.tsv"),
  "androstenediol disulfate"
)

andro_q025 <- read_results(
  file.path(andro_dir, "androstenediol_disulfate_q025_taxa_results.tsv"),
  "androstenediol disulfate"
)

andro_nominal <- read_results(
  file.path(andro_dir, "androstenediol_disulfate_nominal_p005_taxa_results.tsv"),
  "androstenediol disulfate"
)

all_results <- rbind(phenol_all, andro_all)
q025_results <- rbind(phenol_q025, andro_q025)
nominal_results <- rbind(phenol_nominal, andro_nominal)

write.table(
  all_results,
  file = file.path(out_dir, "crossomics_sulfur_all_taxa_results_combined.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  q025_results,
  file = file.path(out_dir, "crossomics_sulfur_q025_results_combined.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  nominal_results,
  file = file.path(out_dir, "crossomics_sulfur_nominal_p005_results_combined.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_df <- data.frame(
  metabolite = c("phenol sulfate", "androstenediol disulfate"),
  n_all_tested = c(nrow(phenol_all), nrow(andro_all)),
  n_nominal_p005 = c(nrow(phenol_nominal), nrow(andro_nominal)),
  n_q025 = c(nrow(phenol_q025), nrow(andro_q025))
)

write.table(
  summary_df,
  file = file.path(out_dir, "crossomics_sulfur_summary_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nSummary counts:\n")
print(summary_df)

plot_bar <- function(df, outfile, title_text, top_n = NULL) {
  if (nrow(df) == 0) {
    cat("No rows for plot:", outfile, "\n")
    return(NULL)
  }

  df <- df[order(df$pval), ]

  if (!is.null(top_n) && nrow(df) > top_n) {
    df <- df[seq_len(top_n), ]
  }

  df$label <- paste0(df$taxon_clean, " | ", df$metabolite)
  df$label <- factor(df$label, levels = rev(df$label))

  p <- ggplot(df, aes(x = label, y = coef, fill = metabolite)) +
    geom_col(width = 0.75) +
    coord_flip() +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    labs(
      title = title_text,
      x = NULL,
      y = "MaAsLin2 coefficient"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    )

  ggsave(
    filename = outfile,
    plot = p,
    width = 10,
    height = max(4, 0.35 * nrow(df) + 2),
    units = "in"
  )

  cat("Saved:", outfile, "\n")
}

plot_bar(
  q025_results,
  file.path(out_dir, "crossomics_sulfur_q025_barplot.pdf"),
  "Species associated with sulfur-related metabolites, q <= 0.25"
)

plot_bar(
  nominal_results,
  file.path(out_dir, "crossomics_sulfur_nominal_p005_top20_barplot.pdf"),
  "Top nominal species-metabolite associations, p < 0.05",
  top_n = 20
)

plot_bar(
  phenol_q025,
  file.path(out_dir, "phenol_sulfate_q025_barplot.pdf"),
  "Species associated with phenol sulfate, q <= 0.25"
)

plot_bar(
  andro_nominal,
  file.path(out_dir, "androstenediol_disulfate_nominal_p005_barplot.pdf"),
  "Species nominally associated with androstenediol disulfate, p < 0.05"
)

cat("\nDone. Summary plots and tables saved to:\n")
cat(out_dir, "\n")
