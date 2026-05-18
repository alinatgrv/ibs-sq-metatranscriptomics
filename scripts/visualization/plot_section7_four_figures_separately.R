suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(scales)
})

# ==============================
# Paths
# ==============================

PROJECT_DIR <- "/home/alina_tgrv/beegfs/IBS_SQ"
setwd(PROJECT_DIR)

OUT_DIR <- file.path(
  PROJECT_DIR,
  "results/presentation_figures/section7_functional_taxonomic_profile_separate"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TOP_PATHWAYS_HEATMAP <- 20
TOP_PATHWAYS_BARPLOT <- 12
TOP_SPECIES_HEATMAP <- 20

# ==============================
# Helper functions
# ==============================

pick_file <- function(candidates) {
  for (f in candidates) {
    if (file.exists(f)) return(f)
  }
  stop("None of these files exist:\n", paste(candidates, collapse = "\n"))
}

clean_sample_name <- function(x) {
  x %>%
    as.character() %>%
    str_replace("^#", "") %>%
    str_replace("(_Abundance.*|_Coverage.*)$", "") %>%
    str_replace("_combined$", "") %>%
    str_squish()
}

standardize_group <- function(x) {
  y <- as.character(x)
  y2 <- case_when(
    str_detect(str_to_lower(y), "control|healthy|hc") ~ "Control",
    str_detect(str_to_lower(y), "ibs") ~ "IBS",
    TRUE ~ y
  )
  factor(y2, levels = c("Control", "IBS", sort(setdiff(unique(y2), c("Control", "IBS")))))
}

clean_pathway_label <- function(x, max_width = 58) {
  x %>%
    as.character() %>%
    str_replace("^([A-Z0-9-]+):\\s*", "") %>%
    str_replace_all("_", " ") %>%
    str_replace_all("-", "-") %>%
    str_squish() %>%
    str_to_sentence() %>%
    str_trunc(width = max_width)
}

clean_species_label <- function(x) {
  x <- x %>%
    as.character() %>%
    str_replace("^s__", "") %>%
    str_replace_all("_", " ") %>%
    str_squish()

  vapply(str_split(x, " "), function(parts) {
    if (length(parts) >= 2) {
      paste0(str_sub(parts[1], 1, 1), ". ", paste(parts[-1], collapse = " "))
    } else {
      parts[1]
    }
  }, character(1))
}

theme_pub <- function(base_size = 16) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 3, hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, color = "grey25"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      strip.text = element_text(face = "bold", size = base_size),
      plot.margin = margin(10, 16, 10, 10)
    )
}

save_plot <- function(plot, name, width, height) {
  ggsave(
    filename = file.path(OUT_DIR, paste0(name, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in"
  )

  ggsave(
    filename = file.path(OUT_DIR, paste0(name, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 400,
    bg = "white"
  )
}

# ==============================
# Input files
# ==============================

metadata_file <- pick_file(c(
  "results/humann_metatranscriptome/qc_report_2026-04-13/sample_metadata_qc_326.tsv",
  "results/humann_metatranscriptome/joined_all/sample_metadata_qc_326.tsv",
  "data/metadata_326_clean_v2.tsv",
  "metadata_326_clean_v2.tsv"
))

pathabundance_file <- pick_file(c(
  "results/humann_metatranscriptome/joined_all/split/pathabundance_all_unstratified_cpm.tsv",
  "results/humann_metatranscriptome/joined_all/split/pathabundance_all_unstratified.tsv",
  "results/humann_metatranscriptome/joined_all/pathabundance_all.tsv"
))

species_file <- pick_file(c(
  "results/metaphlan_metatranscriptome/joined/species_sample_matrix_with_groups.tsv"
))

message("Metadata file: ", metadata_file)
message("Pathabundance file: ", pathabundance_file)
message("Species file: ", species_file)

# ==============================
# Metadata
# ==============================

meta <- fread(metadata_file, check.names = FALSE) %>% as_tibble()

sample_col <- intersect(c("sample", "Sample", "sample_id", "SampleID", "run_accession"), names(meta))[1]
if (is.na(sample_col)) sample_col <- names(meta)[1]

group_col <- intersect(c("ibs_status", "Group", "group", "_group"), names(meta))[1]
if (is.na(group_col)) {
  stop("Could not find IBS/Control group column in metadata.")
}

meta2 <- meta %>%
  mutate(
    sample_key = clean_sample_name(.data[[sample_col]]),
    group = standardize_group(.data[[group_col]])
  ) %>%
  distinct(sample_key, .keep_all = TRUE) %>%
  select(sample_key, group)

# ==============================
# HUMAnN pathabundance
# ==============================

path_dt <- fread(pathabundance_file, check.names = FALSE)
feature_col <- names(path_dt)[1]
setnames(path_dt, feature_col, "feature")

path_dt <- path_dt %>%
  as_tibble() %>%
  filter(
    !str_detect(feature, "\\|"),
    !feature %in% c("UNMAPPED", "UNINTEGRATED")
  )

sample_cols <- setdiff(names(path_dt), "feature")

path_long <- path_dt %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample_raw",
    values_to = "abundance"
  ) %>%
  mutate(
    abundance = suppressWarnings(as.numeric(abundance)),
    sample_key = clean_sample_name(sample_raw)
  ) %>%
  left_join(meta2, by = "sample_key") %>%
  filter(!is.na(group), !is.na(abundance))

abundance_axis_label <- ifelse(
  str_detect(str_to_lower(pathabundance_file), "cpm"),
  "Mean HUMAnN pathway abundance, CPM",
  "Mean HUMAnN pathway abundance"
)

# ==============================
# Figure 1: HUMAnN pathway profile
# ==============================

top_pathways <- path_long %>%
  group_by(feature) %>%
  summarise(mean_all = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_all)) %>%
  slice_head(n = TOP_PATHWAYS_HEATMAP) %>%
  mutate(feature_label = make.unique(clean_pathway_label(feature), sep = " "))

fig1_df <- path_long %>%
  filter(feature %in% top_pathways$feature) %>%
  group_by(group, feature) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  left_join(top_pathways, by = "feature") %>%
  mutate(
    feature_label = factor(
      feature_label,
      levels = rev(top_pathways$feature_label[order(top_pathways$mean_all, decreasing = TRUE)])
    ),
    log_mean = log10(mean_abundance + 1)
  )

fig1 <- ggplot(fig1_df, aes(x = group, y = feature_label, fill = log_mean)) +
  geom_tile(color = "white", linewidth = 0.45) +
  scale_fill_viridis_c(
    option = "C",
    name = "log10(mean + 1)"
  ) +
  labs(
    title = "HUMAnN pathway profile",
    subtitle = "Top pathways by mean abundance across all samples",
    x = NULL,
    y = NULL
  ) +
  theme_pub(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 15, face = "bold")
  )

save_plot(fig1, "01_HUMAnN_pathway_profile", width = 8.5, height = 8.2)

# ==============================
# Figure 2: Dominant functional pathways
# ==============================

top_bar <- path_long %>%
  group_by(feature) %>%
  summarise(mean_all = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_all)) %>%
  slice_head(n = TOP_PATHWAYS_BARPLOT) %>%
  mutate(feature_label = make.unique(clean_pathway_label(feature, max_width = 52), sep = " "))

fig2_df <- path_long %>%
  filter(feature %in% top_bar$feature) %>%
  group_by(group, feature) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  left_join(top_bar, by = "feature") %>%
  mutate(
    feature_label = factor(
      feature_label,
      levels = rev(top_bar$feature_label[order(top_bar$mean_all, decreasing = TRUE)])
    )
  )

fig2 <- ggplot(fig2_df, aes(x = mean_abundance, y = feature_label, fill = group)) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "grey20",
    linewidth = 0.2
  ) +
  scale_x_continuous(labels = label_number(big.mark = " ", accuracy = 1)) +
  labs(
    title = "Dominant functional pathways",
    subtitle = "Mean abundance of the most abundant HUMAnN pathways",
    x = abundance_axis_label,
    y = NULL,
    fill = NULL
  ) +
  theme_pub(base_size = 17) +
  theme(
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

save_plot(fig2, "02_Dominant_functional_pathways", width = 9.2, height = 6.4)

# ==============================
# Species-level taxonomic profile
# ==============================

sp_dt <- fread(species_file, check.names = FALSE) %>% as_tibble()

sp_sample_col <- intersect(c("sample", "Sample", "sample_id", "SampleID", "run_accession"), names(sp_dt))[1]
if (is.na(sp_sample_col)) sp_sample_col <- names(sp_dt)[1]

sp_group_col <- intersect(c("_group", "group", "Group", "ibs_status"), names(sp_dt))[1]

species_cols <- names(sp_dt)[str_detect(names(sp_dt), "^s__")]
if (length(species_cols) == 0) {
  exclude_cols <- c(sp_sample_col, sp_group_col)
  numeric_cols <- names(sp_dt)[sapply(sp_dt, is.numeric)]
  species_cols <- setdiff(numeric_cols, exclude_cols)
}

if (length(species_cols) == 0) {
  stop("No species columns were detected in species table.")
}

species_long <- sp_dt %>%
  mutate(sample_key = clean_sample_name(.data[[sp_sample_col]])) %>%
  {
    if (!is.na(sp_group_col)) {
      mutate(., group = standardize_group(.data[[sp_group_col]]))
    } else {
      left_join(., meta2, by = "sample_key")
    }
  } %>%
  pivot_longer(
    cols = all_of(species_cols),
    names_to = "species",
    values_to = "abundance"
  ) %>%
  mutate(abundance = suppressWarnings(as.numeric(abundance))) %>%
  filter(!is.na(group), !is.na(abundance))

top_species <- species_long %>%
  group_by(species) %>%
  summarise(mean_all = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_all)) %>%
  slice_head(n = TOP_SPECIES_HEATMAP) %>%
  mutate(species_label = make.unique(clean_species_label(species), sep = " "))

fig3_df <- species_long %>%
  filter(species %in% top_species$species) %>%
  group_by(group, species) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  left_join(top_species, by = "species") %>%
  mutate(
    species_label = factor(
      species_label,
      levels = rev(top_species$species_label[order(top_species$mean_all, decreasing = TRUE)])
    )
  )

fig3 <- ggplot(fig3_df, aes(x = group, y = species_label, fill = mean_abundance)) +
  geom_tile(color = "white", linewidth = 0.45) +
  scale_fill_viridis_c(
    option = "C",
    name = "Mean relative\nabundance"
  ) +
  labs(
    title = "Species-level taxonomic profile",
    subtitle = "Top species by mean MetaPhlAn abundance",
    x = NULL,
    y = NULL
  ) +
  theme_pub(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 12, face = "italic"),
    axis.text.x = element_text(size = 15, face = "bold")
  )

save_plot(fig3, "03_Species_level_taxonomic_profile", width = 7.8, height = 8.2)

# ==============================
# Figure 4: SQ-specific pathway check
# ==============================

target_sq_pathways <- c(
  "PWY-7446",
  "PWY-7722",
  "PWY-8213",
  "PWY-8348",
  "PWY-8349",
  "PWY-8350"
)

detected_ids <- vapply(
  target_sq_pathways,
  function(id) any(str_detect(path_dt$feature, fixed(id))),
  logical(1)
)

sq_target_df <- tibble(
  layer = "Target MetaCyc SQ pathway IDs in HUMAnN pathabundance",
  item = target_sq_pathways,
  status = if_else(detected_ids, "Detected", "Not detected"),
  label = if_else(detected_ids, "detected", "absent")
)

mapping_files <- c(
  "KO mapping" = "databases/humann/utility_mapping/utility_mapping/map_ko_name.txt.gz",
  "EC mapping" = "databases/humann/utility_mapping/utility_mapping/map_ec_name.txt.gz",
  "GO mapping" = "databases/humann/utility_mapping/utility_mapping/map_go_name.txt.gz",
  "eggNOG mapping" = "databases/humann/utility_mapping/utility_mapping/map_eggnog_name.txt.gz"
)

sq_pattern <- "sulfoquinovose|sulfoquinovos|sulfolactaldehyde|sulfopropanediol"

mapping_df <- imap_dfr(mapping_files, function(f, nm) {
  if (!file.exists(f)) {
    tibble(
      layer = "SQ-related terms in HUMAnN utility mappings",
      item = nm,
      status = "Not detected",
      label = "file missing"
    )
  } else {
    lines <- readLines(gzfile(f), warn = FALSE)
    hits <- lines[str_detect(str_to_lower(lines), sq_pattern)]
    tibble(
      layer = "SQ-related terms in HUMAnN utility mappings",
      item = nm,
      status = if_else(length(hits) > 0, "Detected", "Not detected"),
      label = if_else(length(hits) > 0, paste0(length(hits), " terms"), "0 terms")
    )
  }
})

fig4_df <- bind_rows(sq_target_df, mapping_df) %>%
  mutate(
    item = factor(item, levels = rev(unique(item))),
    status = factor(status, levels = c("Detected", "Not detected"))
  )

fig4 <- ggplot(fig4_df, aes(x = "Search result", y = item, fill = status)) +
  geom_tile(color = "white", linewidth = 0.7, width = 0.8, height = 0.8) +
  geom_text(aes(label = label), size = 5.2, fontface = "bold") +
  facet_grid(layer ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(
    values = c("Detected" = "#2C7BB6", "Not detected" = "#D9D9D9"),
    drop = FALSE
  ) +
  labs(
    title = "SQ-specific pathway check",
    subtitle = "Target SQ MetaCyc pathway IDs were checked directly in HUMAnN pathabundance; SQ-related enzyme terms were checked in utility mappings",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_pub(base_size = 16) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    strip.text.y = element_text(angle = 0, hjust = 0, size = 13)
  )

save_plot(fig4, "04_SQ_specific_pathway_check", width = 9.2, height = 6.4)

# ==============================
# Done
# ==============================

message("\nDone.")
message("Figures saved to: ", OUT_DIR)
message("\nCreated files:")
print(list.files(OUT_DIR, full.names = TRUE))
