#!/usr/bin/env Rscript

# ============================================================
# MetaPhlAn species-level overview plots for IBS vs Control
# Input:
#   results/metaphlan_metatranscriptome/joined/species_sample_matrix_with_groups.tsv
# Output:
#   results/metaphlan_metatranscriptome/figures_taxonomy_report/
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

base_dir <- "/home/alina_tgrv/beegfs/IBS_SQ"

input_file <- file.path(
  base_dir,
  "results/metaphlan_metatranscriptome/joined/species_sample_matrix_with_groups.tsv"
)

out_dir <- file.path(
  base_dir,
  "results/metaphlan_metatranscriptome/figures_taxonomy_report"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)

cat("Input dimensions:\n")
print(dim(df))
cat("Column names preview:\n")
print(head(colnames(df), 20))

# -----------------------------
# Detect group and sample columns
# -----------------------------

possible_group_cols <- c("_group", "ibs_status", "IBS_status", "group", "Group", "diagnosis", "status")
group_col <- possible_group_cols[possible_group_cols %in% colnames(df)][1]

if (is.na(group_col)) {
  stop("Could not find group column. Expected one of: ibs_status, group, diagnosis, status")
}

possible_sample_cols <- c("sample", "Sample", "sample_id", "SampleID", "run", "Run")
sample_col <- possible_sample_cols[possible_sample_cols %in% colnames(df)][1]

if (is.na(sample_col)) {
  sample_col <- colnames(df)[1]
  cat("Sample column not detected by name; using first column:", sample_col, "\n")
}

df[[group_col]] <- as.factor(df[[group_col]])

cat("Group column:", group_col, "\n")
cat("Sample column:", sample_col, "\n")
cat("Group counts:\n")
print(table(df[[group_col]]))

# -----------------------------
# Detect species columns
# -----------------------------

metadata_cols <- c(
  sample_col, group_col,
  "subject_id", "age", "sex", "bmi", "diet",
  "bowel_habit", "ethnicity", "qc_flags"
)

candidate_cols <- setdiff(colnames(df), metadata_cols)

numeric_cols <- candidate_cols[sapply(df[candidate_cols], is.numeric)]
species_cols <- numeric_cols[grepl("s__", numeric_cols)]

if (length(species_cols) == 0) {
  species_cols <- numeric_cols
  warning("No columns with 's__' detected; using all numeric non-metadata columns as species.")
}

cat("Detected species columns:", length(species_cols), "\n")

mat <- as.matrix(df[, species_cols, drop = FALSE])
rownames(mat) <- df[[sample_col]]

# Replace NA with 0
mat[is.na(mat)] <- 0

# If MetaPhlAn values look like percentages, convert to proportions for diversity only.
# For plots we keep original scale.
mat_for_div <- mat
if (max(mat_for_div, na.rm = TRUE) > 1.5) {
  mat_for_div <- mat_for_div / 100
}

# -----------------------------
# Helper: clean species names
# -----------------------------

clean_taxon <- function(x) {
  x <- sub("^.*\\|s__", "s__", x)
  x <- gsub("s__", "", x)
  x <- gsub("_", " ", x)
  x
}

# -----------------------------
# Summary statistics by group
# -----------------------------

groups <- levels(df[[group_col]])

if (!all(c("Control", "IBS") %in% groups)) {
  cat("Warning: expected groups Control and IBS. Detected groups are:\n")
  print(groups)
}

mean_by_group <- aggregate(
  df[, species_cols, drop = FALSE],
  by = list(group = df[[group_col]]),
  FUN = mean,
  na.rm = TRUE
)

prev_by_group <- aggregate(
  df[, species_cols, drop = FALSE] > 0,
  by = list(group = df[[group_col]]),
  FUN = mean,
  na.rm = TRUE
)

mean_long <- data.frame()
prev_long <- data.frame()

for (sp in species_cols) {
  tmp_mean <- data.frame(
    species = sp,
    group = mean_by_group$group,
    mean_abundance = as.numeric(mean_by_group[[sp]])
  )
  mean_long <- rbind(mean_long, tmp_mean)

  tmp_prev <- data.frame(
    species = sp,
    group = prev_by_group$group,
    prevalence = as.numeric(prev_by_group[[sp]])
  )
  prev_long <- rbind(prev_long, tmp_prev)
}

# Wide summary for IBS vs Control
get_group_value <- function(long_df, value_col, group_name) {
  out <- long_df[long_df$group == group_name, c("species", value_col)]
  colnames(out)[2] <- group_name
  out
}

mean_control <- get_group_value(mean_long, "mean_abundance", "Control")
mean_ibs <- get_group_value(mean_long, "mean_abundance", "IBS")
prev_control <- get_group_value(prev_long, "prevalence", "Control")
prev_ibs <- get_group_value(prev_long, "prevalence", "IBS")

summary_df <- Reduce(
  function(x, y) merge(x, y, by = "species", all = TRUE),
  list(mean_control, mean_ibs, prev_control, prev_ibs)
)

colnames(summary_df) <- c(
  "species",
  "mean_Control",
  "mean_IBS",
  "prev_Control",
  "prev_IBS"
)

summary_df[is.na(summary_df)] <- 0

eps <- 1e-6
summary_df$log2FC_IBS_vs_Control <- log2((summary_df$mean_IBS + eps) / (summary_df$mean_Control + eps))
summary_df$prevalence_diff_IBS_minus_Control <- summary_df$prev_IBS - summary_df$prev_Control
summary_df$mean_overall <- rowMeans(summary_df[, c("mean_Control", "mean_IBS")])
summary_df$species_clean <- clean_taxon(summary_df$species)

summary_out <- file.path(out_dir, "species_group_summary_for_plots.tsv")
write.table(summary_df, summary_out, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Saved summary table:\n")
cat(summary_out, "\n")

# ============================================================
# 1. Heatmap top 30 mean abundance by group
# ============================================================

top_n <- 30
top_species <- summary_df[order(summary_df$mean_overall, decreasing = TRUE), ][1:min(top_n, nrow(summary_df)), ]

heat_df <- data.frame(
  species = rep(top_species$species_clean, each = 2),
  group = rep(c("Control", "IBS"), times = nrow(top_species)),
  mean_abundance = as.vector(t(top_species[, c("mean_Control", "mean_IBS")]))
)

heat_df$species <- factor(heat_df$species, levels = rev(top_species$species_clean))
heat_df$group <- factor(heat_df$group, levels = c("Control", "IBS"))

p_heat <- ggplot(heat_df, aes(x = group, y = species, fill = log10(mean_abundance + eps))) +
  geom_tile(color = "white") +
  labs(
    title = "Top species by mean abundance",
    x = "",
    y = "",
    fill = "log10(mean abundance)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  file.path(out_dir, "01_top30_species_mean_abundance_heatmap.pdf"),
  p_heat,
  width = 7,
  height = 8
)

# ============================================================
# 2. Top log2FC barplot
# ============================================================

fc_df <- summary_df[
  summary_df$prev_IBS >= 0.05 | summary_df$prev_Control >= 0.05,
]

fc_top_up <- fc_df[order(fc_df$log2FC_IBS_vs_Control, decreasing = TRUE), ][1:min(15, nrow(fc_df)), ]
fc_top_down <- fc_df[order(fc_df$log2FC_IBS_vs_Control, decreasing = FALSE), ][1:min(15, nrow(fc_df)), ]
fc_plot_df <- rbind(fc_top_down, fc_top_up)
fc_plot_df <- fc_plot_df[!duplicated(fc_plot_df$species), ]
fc_plot_df$direction <- ifelse(fc_plot_df$log2FC_IBS_vs_Control > 0, "Higher in IBS", "Higher in Control")
fc_plot_df$species_clean <- factor(fc_plot_df$species_clean, levels = fc_plot_df$species_clean[order(fc_plot_df$log2FC_IBS_vs_Control)])

p_fc <- ggplot(fc_plot_df, aes(x = species_clean, y = log2FC_IBS_vs_Control, fill = direction)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Species with strongest mean abundance differences",
    x = "",
    y = "log2FC IBS vs Control",
    fill = ""
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  file.path(out_dir, "02_top_species_log2fc_barplot.pdf"),
  p_fc,
  width = 8,
  height = 8
)

# ============================================================
# 3. Prevalence difference plot
# ============================================================

prev_df <- summary_df[
  summary_df$prev_IBS >= 0.05 | summary_df$prev_Control >= 0.05,
]

prev_top_up <- prev_df[order(prev_df$prevalence_diff_IBS_minus_Control, decreasing = TRUE), ][1:min(15, nrow(prev_df)), ]
prev_top_down <- prev_df[order(prev_df$prevalence_diff_IBS_minus_Control, decreasing = FALSE), ][1:min(15, nrow(prev_df)), ]
prev_plot_df <- rbind(prev_top_down, prev_top_up)
prev_plot_df <- prev_plot_df[!duplicated(prev_plot_df$species), ]
prev_plot_df$direction <- ifelse(prev_plot_df$prevalence_diff_IBS_minus_Control > 0, "More prevalent in IBS", "More prevalent in Control")
prev_plot_df$species_clean <- factor(prev_plot_df$species_clean, levels = prev_plot_df$species_clean[order(prev_plot_df$prevalence_diff_IBS_minus_Control)])

p_prev <- ggplot(prev_plot_df, aes(x = species_clean, y = prevalence_diff_IBS_minus_Control, fill = direction)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Species with strongest prevalence differences",
    x = "",
    y = "Prevalence difference: IBS - Control",
    fill = ""
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  file.path(out_dir, "03_top_species_prevalence_difference.pdf"),
  p_prev,
  width = 8,
  height = 8
)

# ============================================================
# 4. Alpha diversity: Shannon
# ============================================================

shannon <- function(x) {
  x <- x[x > 0]
  if (length(x) == 0) return(0)
  p <- x / sum(x)
  -sum(p * log(p))
}

alpha_df <- data.frame(
  sample = df[[sample_col]],
  group = df[[group_col]],
  shannon = apply(mat_for_div, 1, shannon),
  observed_species = rowSums(mat > 0)
)

write.table(
  alpha_df,
  file.path(out_dir, "alpha_diversity_species.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_alpha <- ggplot(alpha_df, aes(x = group, y = shannon, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  labs(
    title = "Species-level alpha diversity",
    x = "",
    y = "Shannon diversity"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  file.path(out_dir, "04_alpha_diversity_shannon_boxplot.pdf"),
  p_alpha,
  width = 5,
  height = 5
)

# ============================================================
# 5. Bray-Curtis PCoA
# ============================================================

bray_curtis <- function(x) {
  n <- nrow(x)
  d <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      denom <- sum(x[i, ] + x[j, ])
      if (denom == 0) {
        d[i, j] <- 0
      } else {
        d[i, j] <- sum(abs(x[i, ] - x[j, ])) / denom
      }
    }
  }
  as.dist(d)
}

cat("Calculating Bray-Curtis distance...\n")
bc <- bray_curtis(mat_for_div)

pcoa <- cmdscale(bc, eig = TRUE, k = 2)

pcoa_df <- data.frame(
  sample = df[[sample_col]],
  group = df[[group_col]],
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2]
)

var_explained <- round(100 * pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]), 1)

write.table(
  pcoa_df,
  file.path(out_dir, "pcoa_bray_curtis_coordinates.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_pcoa <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = group)) +
  geom_point(alpha = 0.75, size = 2) +
  labs(
    title = "PCoA of species-level MetaPhlAn profiles",
    x = paste0("PCoA1 (", var_explained[1], "%)"),
    y = paste0("PCoA2 (", var_explained[2], "%)"),
    color = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  file.path(out_dir, "05_pcoa_bray_curtis_ibs_control.pdf"),
  p_pcoa,
  width = 6,
  height = 5
)

# ============================================================
# 6. Candidate species boxplots
# ============================================================

candidate_species <- c(
  "s__Bacteroides_ovatus",
  "s__Blautia_faecis",
  "s__Parabacteroides_merdae",
  "s__Phascolarctobacterium_faecium",
  "s__Parabacteroides_distasonis",
  "s__Bacteroides_caccae",
  "s__Alistipes_onderdonkii",
  "s__Bacteroides_dorei",
  "s__Faecalibacterium_prausnitzii",
  "s__Bilophila_wadsworthia",
  "s__Roseburia_inulinivorans"
)

matched_candidates <- species_cols[
  sapply(candidate_species, function(x) {
    grep(x, species_cols, value = TRUE)[1]
  }) |> unlist() |> is.na() == FALSE
]

# safer matching
matched_candidates <- unique(unlist(lapply(candidate_species, function(x) grep(x, species_cols, value = TRUE))))

if (length(matched_candidates) > 0) {
  box_df <- data.frame()

  for (sp in matched_candidates) {
    tmp <- data.frame(
      sample = df[[sample_col]],
      group = df[[group_col]],
      species = clean_taxon(sp),
      abundance = df[[sp]]
    )
    box_df <- rbind(box_df, tmp)
  }

  box_df$species <- factor(box_df$species, levels = unique(box_df$species))

  write.table(
    box_df,
    file.path(out_dir, "candidate_species_abundance_long.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  p_box <- ggplot(box_df, aes(x = group, y = abundance, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.35, size = 0.7) +
    facet_wrap(~ species, scales = "free_y", ncol = 3) +
    labs(
      title = "Candidate species abundance by IBS status",
      x = "",
      y = "MetaPhlAn relative abundance"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 8),
      plot.title = element_text(hjust = 0.5)
    )

  ggsave(
    file.path(out_dir, "06_candidate_species_boxplots.pdf"),
    p_box,
    width = 10,
    height = 8
  )
} else {
  cat("No candidate species were matched in the matrix.\n")
}

cat("\nDone. Figures saved to:\n")
cat(out_dir, "\n")
