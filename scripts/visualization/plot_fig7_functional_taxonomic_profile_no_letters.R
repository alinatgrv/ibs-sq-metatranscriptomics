suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

base_dir <- "/home/alina_tgrv/beegfs/IBS_SQ"

out_dir <- file.path(base_dir, "results/presentation/main_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)

# -----------------------------
# Helper functions
# -----------------------------

first_existing <- function(paths) {
  full <- file.path(base_dir, paths)
  hit <- full[file.exists(full)]
  if (length(hit) == 0) {
    stop("None of these files exist:\n", paste(full, collapse = "\n"))
  }
  hit[1]
}

clean_sample_id <- function(x) {
  x <- gsub("_combined_Abundance-RPKs$", "", x)
  x <- gsub("_combined_Abundance$", "", x)
  x <- gsub("_combined_Coverage$", "", x)
  x <- gsub("_Abundance-RPKs$", "", x)
  x <- gsub("_Abundance$", "", x)
  x <- gsub("_Coverage$", "", x)
  x <- gsub("-RPKs$", "", x)
  x
}

truncate_label <- function(x, n = 55) {
  x <- as.character(x)
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
}

clean_pathway_label <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\|.*$", "", x)
  truncate_label(x, 58)
}

read_feature_table <- function(path) {
  message("Reading feature table: ", path)

  df <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    comment.char = "",
    quote = "",
    stringsAsFactors = FALSE
  )

  features <- df[[1]]
  df[[1]] <- NULL

  sample_names <- clean_sample_id(colnames(df))

  mat <- do.call(cbind, lapply(df, function(z) suppressWarnings(as.numeric(as.character(z)))))
  rownames(mat) <- features
  colnames(mat) <- sample_names

  if (any(duplicated(colnames(mat)))) {
    message("Duplicated sample names after cleaning were collapsed by summing.")
    mat <- t(rowsum(t(mat), group = colnames(mat), reorder = FALSE))
  }

  mat_t <- t(mat)
  mat_t[is.na(mat_t)] <- 0
  mat_t
}

bray_curtis_dist <- function(mat) {
  n <- nrow(mat)
  d <- matrix(0, n, n)
  rownames(d) <- rownames(mat)
  colnames(d) <- rownames(mat)

  for (i in seq_len(n - 1)) {
    xi <- mat[i, ]
    for (j in (i + 1):n) {
      xj <- mat[j, ]
      denom <- sum(xi + xj, na.rm = TRUE)
      val <- ifelse(denom == 0, 0, sum(abs(xi - xj), na.rm = TRUE) / denom)
      d[i, j] <- val
      d[j, i] <- val
    }
  }

  as.dist(d)
}

theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(size = base_size),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey88"),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

group_colors <- c(
  "Control" = "#F8766D",
  "IBS" = "#00BFC4"
)

# -----------------------------
# Input files
# -----------------------------

metadata_path <- first_existing(c(
  "metadata/metadata_326_clean_v2.tsv"
))

species_path <- first_existing(c(
  "results/metaphlan_metatranscriptome/joined/species_sample_matrix_with_groups.tsv"
))

pathabundance_path <- first_existing(c(
  "results/humann_metatranscriptome/joined_all/split/pathabundance_all_unstratified_cpm.tsv",
  "results/humann_metatranscriptome/joined_all/split/pathabundance_all_unstratified.tsv",
  "results/humann_metatranscriptome/joined_all/pathabundance_all_unstratified_cpm.tsv",
  "results/humann_metatranscriptome/joined_all/pathabundance_all_unstratified.tsv",
  "results/humann_metatranscriptome/joined_all/pathabundance_all.tsv"
))

message("\nSelected input files:")
message("Metadata:      ", metadata_path)
message("Species table: ", species_path)
message("Pathways:      ", pathabundance_path)

# -----------------------------
# Read metadata
# -----------------------------

meta <- read.delim(
  metadata_path,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

meta$sample <- as.character(meta$sample)
meta$ibs_status <- factor(meta$ibs_status, levels = c("Control", "IBS"))

message("\nMetadata samples: ", nrow(meta))
message("Group distribution:")
print(table(meta$ibs_status, useNA = "ifany"))

# -----------------------------
# Read HUMAnN pathway abundance
# -----------------------------

pwy <- read_feature_table(pathabundance_path)

# Remove stratified rows if present, and HUMAnN technical rows
keep_features <- !grepl("\\|", colnames(pwy)) &
  !(colnames(pwy) %in% c("UNMAPPED", "UNINTEGRATED"))

pwy <- pwy[, keep_features, drop = FALSE]

common_pwy_samples <- intersect(rownames(pwy), meta$sample)
pwy <- pwy[common_pwy_samples, , drop = FALSE]
meta_pwy <- meta[match(common_pwy_samples, meta$sample), , drop = FALSE]

message("\nPathway matrix after filtering:")
message("Samples: ", nrow(pwy))
message("Pathways: ", ncol(pwy))

# -----------------------------
# Panel A: detected pathways per sample
# -----------------------------

detected_df <- data.frame(
  sample = rownames(pwy),
  detected_pathways = rowSums(pwy > 0, na.rm = TRUE),
  ibs_status = meta_pwy$ibs_status,
  stringsAsFactors = FALSE
)

pA <- ggplot(detected_df, aes(x = ibs_status, y = detected_pathways, fill = ibs_status)) +
  geom_violin(trim = FALSE, alpha = 0.55, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.85, linewidth = 0.35) +
  geom_jitter(aes(color = ibs_status), width = 0.12, size = 1.25, alpha = 0.55, show.legend = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  labs(
    title = "A. HUMAnN pathway reconstruction",
    subtitle = "Detected MetaCyc pathways per sample",
    x = NULL,
    y = "Detected pathways"
  ) +
  theme_pub(11) +
  theme(legend.position = "none")

# -----------------------------
# Panel B: top MetaCyc pathways heatmap
# -----------------------------

pwy_log <- log10(pwy + 1)
top_n <- 12

top_features <- names(sort(colMeans(pwy_log, na.rm = TRUE), decreasing = TRUE))[seq_len(min(top_n, ncol(pwy_log)))]

long_list <- lapply(top_features, function(feat) {
  data.frame(
    sample = rownames(pwy_log),
    pathway = feat,
    value = as.numeric(pwy_log[, feat]),
    stringsAsFactors = FALSE
  )
})

top_long <- do.call(rbind, long_list)
top_long$ibs_status <- meta_pwy$ibs_status[match(top_long$sample, meta_pwy$sample)]

heat_df <- aggregate(
  value ~ pathway + ibs_status,
  data = top_long,
  FUN = mean,
  na.rm = TRUE
)

path_labels <- data.frame(
  pathway = top_features,
  pathway_label = clean_pathway_label(top_features),
  stringsAsFactors = FALSE
)

heat_df$pathway_label <- path_labels$pathway_label[match(heat_df$pathway, path_labels$pathway)]
heat_df$pathway_label <- factor(heat_df$pathway_label, levels = rev(path_labels$pathway_label))
heat_df$ibs_status <- factor(heat_df$ibs_status, levels = c("Control", "IBS"))

pB <- ggplot(heat_df, aes(x = ibs_status, y = pathway_label, fill = value)) +
  geom_tile(color = "white", linewidth = 0.6) +
  scale_fill_gradient(
    low = "grey95",
    high = "#2166AC",
    name = "Mean\nlog10(CPM+1)"
  ) +
  labs(
    title = "Dominant functional pathways",
    subtitle = paste0("Top ", length(top_features), " MetaCyc pathways"),
    x = NULL,
    y = NULL
  ) +
  theme_pub(10) +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  )

# -----------------------------
# Read MetaPhlAn species table
# -----------------------------

message("\nReading species table: ", species_path)

species_df <- read.delim(
  species_path,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_col <- if ("sample" %in% colnames(species_df)) "sample" else colnames(species_df)[1]
species_df[[sample_col]] <- as.character(species_df[[sample_col]])

species_cols <- grep("^s__", colnames(species_df), value = TRUE)

if (length(species_cols) == 0) {
  stop("No species columns starting with s__ were found in species table.")
}

species_mat <- as.matrix(data.frame(
  lapply(species_df[, species_cols, drop = FALSE], function(z) suppressWarnings(as.numeric(as.character(z)))),
  check.names = FALSE
))

rownames(species_mat) <- species_df[[sample_col]]
colnames(species_mat) <- species_cols
species_mat[is.na(species_mat)] <- 0

common_species_samples <- intersect(rownames(species_mat), meta$sample)
species_mat <- species_mat[common_species_samples, , drop = FALSE]
meta_species <- meta[match(common_species_samples, meta$sample), , drop = FALSE]

message("\nSpecies matrix:")
message("Samples: ", nrow(species_mat))
message("Species: ", ncol(species_mat))

# -----------------------------
# Panel C: PCoA species-level
# -----------------------------

message("Computing Bray-Curtis distance for species-level PCoA...")

bc <- bray_curtis_dist(species_mat)
pcoa <- cmdscale(bc, k = 2, eig = TRUE)

eig <- pcoa$eig
positive_eig <- eig[eig > 0]
var1 <- round(100 * eig[1] / sum(positive_eig), 1)
var2 <- round(100 * eig[2] / sum(positive_eig), 1)

pcoa_df <- data.frame(
  sample = rownames(pcoa$points),
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2],
  ibs_status = meta_species$ibs_status[match(rownames(pcoa$points), meta_species$sample)],
  stringsAsFactors = FALSE
)

pcoa_df$ibs_status <- factor(pcoa_df$ibs_status, levels = c("Control", "IBS"))

pC <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = ibs_status)) +
  geom_point(size = 2.1, alpha = 0.75) +
  stat_ellipse(linewidth = 0.55, alpha = 0.8, show.legend = FALSE) +
  scale_color_manual(values = group_colors, drop = FALSE, name = "Group") +
  labs(
    title = "Species-level taxonomic profile",
    subtitle = "MetaPhlAn profiles, Bray-Curtis PCoA",
    x = paste0("PCoA1 (", var1, "%)"),
    y = paste0("PCoA2 (", var2, "%)")
  ) +
  theme_pub(11)

# -----------------------------
# Panel D: SQ pathway detection summary
# -----------------------------

sq_targets <- data.frame(
  pathway_id = c("PWY-7446", "PWY-7722", "PWY-8213", "PWY-8348", "PWY-8349", "PWY-8350"),
  pathway_name = c(
    "sulfoquinovose degradation I",
    "sulfoquinovose degradation II",
    "SQ-related pathway",
    "SQ-related pathway",
    "SQ-related pathway",
    "SQ-related pathway"
  ),
  stringsAsFactors = FALSE
)

sq_summary <- do.call(rbind, lapply(seq_len(nrow(sq_targets)), function(i) {
  id <- sq_targets$pathway_id[i]
  hits <- grep(id, colnames(pwy), fixed = TRUE, value = TRUE)

  if (length(hits) == 0) {
    n_detected <- 0
    max_abundance <- 0
  } else {
    sub <- pwy[, hits, drop = FALSE]
    n_detected <- sum(rowSums(sub > 0, na.rm = TRUE) > 0)
    max_abundance <- max(sub, na.rm = TRUE)
  }

  data.frame(
    pathway_id = id,
    pathway_name = sq_targets$pathway_name[i],
    n_detected = n_detected,
    total_samples = nrow(pwy),
    max_abundance = max_abundance,
    status = ifelse(n_detected > 0, "Detected", "Not detected"),
    stringsAsFactors = FALSE
  )
}))

sq_summary$label <- paste0(sq_summary$pathway_id, "\n", sq_summary$pathway_name)
sq_summary$label <- factor(sq_summary$label, levels = rev(sq_summary$label))
sq_summary$status <- factor(sq_summary$status, levels = c("Detected", "Not detected"))

write.table(
  sq_summary,
  file = file.path(out_dir, "fig7D_sq_pathway_detection_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pD <- ggplot(sq_summary, aes(x = "HUMAnN\npathabundance", y = label, fill = status)) +
  geom_tile(color = "white", linewidth = 0.7, width = 0.85, height = 0.8) +
  geom_text(
    aes(label = paste0(n_detected, "/", total_samples, "\nsamples")),
    size = 3.1,
    color = "black"
  ) +
  scale_fill_manual(
    values = c("Detected" = "#2A9D8F", "Not detected" = "grey82"),
    name = "Status",
    drop = FALSE
  ) +
  labs(
    title = "SQ-specific pathway check",
    subtitle = "No reconstructed SQ pathway signal in HUMAnN output",
    x = NULL,
    y = NULL
  ) +
  theme_pub(10) +
  theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  )

# -----------------------------
# Save individual figures
# -----------------------------

save_plot <- function(plot, filename, width, height) {
  pdf_path <- file.path(out_dir, paste0(filename, ".pdf"))
  png_path <- file.path(out_dir, paste0(filename, ".png"))

  ggsave(pdf_path, plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 450)

  message("Saved: ", pdf_path)
  message("Saved: ", png_path)
}

save_plot(pA, "fig7A_detected_pathways_by_group", 5.2, 4.2)
save_plot(pB, "fig7B_top_metacyc_pathways_heatmap", 6.2, 4.5)
save_plot(pC, "fig7C_metaphlan_species_pcoa", 5.4, 4.5)
save_plot(pD, "fig7D_sq_pathway_detection_summary", 6.2, 4.1)

# -----------------------------
# Save combined 2x2 figure
# -----------------------------

combined_pdf <- file.path(out_dir, "fig7_combined_functional_taxonomic_profile_no_letters.pdf")
combined_png <- file.path(out_dir, "fig7_combined_functional_taxonomic_profile_no_letters.png")

pdf(combined_pdf, width = 13.33, height = 7.5, family = "sans")
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(pA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(pC, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(pD, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()

png(combined_png, width = 4800, height = 2700, res = 400)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(pA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(pC, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(pD, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()

message("\nDone.")
message("Combined figure:")
message(combined_pdf)
message(combined_png)
message("\nSQ summary:")
print(sq_summary)
