#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ==============================
# INPUT / OUTPUT
# ==============================

SUMMARY_FILE = Path("results/random_forest/ibs_multiblock_python/rf_block_summary.tsv")
FIG_DIR = Path("results/random_forest/ibs_multiblock_python/figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)

OUT_PREFIX = FIG_DIR / "03b_rf_metrics_heatmap_selected_blocks"

# ==============================
# READ DATA
# ==============================

df = pd.read_csv(SUMMARY_FILE, sep="\t")

# Оставляем только комбинированные модели, без строк "... only"
selected_blocks = [
    "rf_234_metabolites_clinical",
    "rf_234_metabolites_sq",
    "rf_234_taxa_metabolites_clinical",
    "rf_234_metabolites_sq_clinical",
    "rf_234_all_taxa_sq_metabolites",
    "rf_234_taxa_metabolites",
    "rf_234_taxa_sq_clinical",
    "rf_234_sq_clinical",
]

pretty_names = {
    "rf_234_metabolites_clinical": "metabolites + clinical",
    "rf_234_metabolites_sq": "metabolites + SQ",
    "rf_234_taxa_metabolites_clinical": "taxa + metabolites + clinical",
    "rf_234_metabolites_sq_clinical": "metabolites + SQ + clinical",
    "rf_234_all_taxa_sq_metabolites": "taxa + SQ + metabolites",
    "rf_234_taxa_metabolites": "taxa + metabolites",
    "rf_234_taxa_sq_clinical": "taxa + SQ + clinical",
    "rf_234_sq_clinical": "SQ + clinical",
}

metric_cols = [
    "roc_auc_mean",
    "pr_auc_mean",
    "balanced_accuracy_mean",
    "accuracy_mean",
    "f1_mean",
]

metric_names = [
    "ROC-AUC",
    "PR-AUC",
    "Balanced accuracy",
    "Accuracy",
    "F1",
]

plot_df = df[df["block"].isin(selected_blocks)].copy()

# Сохраняем порядок вручную, а не сортируем автоматически
plot_df["block"] = pd.Categorical(
    plot_df["block"],
    categories=selected_blocks,
    ordered=True
)
plot_df = plot_df.sort_values("block")

plot_df["pretty_block"] = plot_df["block"].map(pretty_names)

mat = plot_df[metric_cols].to_numpy(dtype=float)

# ==============================
# PLOT
# ==============================

fig_height = max(6, 0.7 * len(plot_df) + 2)
fig_width = 11

fig, ax = plt.subplots(figsize=(fig_width, fig_height))

im = ax.imshow(
    mat,
    aspect="auto",
    vmin=0.45,
    vmax=1.0,
    cmap="viridis"
)

# Оси
ax.set_xticks(np.arange(len(metric_names)))
ax.set_xticklabels(metric_names, fontsize=17, rotation=30, ha="right")

ax.set_yticks(np.arange(len(plot_df)))
ax.set_yticklabels(plot_df["pretty_block"], fontsize=17)

# Значения внутри ячеек
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        value = mat[i, j]
        ax.text(
            j,
            i,
            f"{value:.2f}",
            ha="center",
            va="center",
            fontsize=18,
            fontweight="normal",
            color="white"
        )

# Заголовок
ax.set_title(
    "Random Forest IBS classification: selected multiblock models",
    fontsize=24,
    pad=18
)

# Сетка между ячейками
ax.set_xticks(np.arange(-0.5, len(metric_names), 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(plot_df), 1), minor=True)
ax.grid(which="minor", color="white", linestyle="-", linewidth=2)
ax.tick_params(which="minor", bottom=False, left=False)

# Colorbar
cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
cbar.set_label("Metric value", fontsize=18)
cbar.ax.tick_params(labelsize=15)

plt.tight_layout()

# ==============================
# SAVE
# ==============================

fig.savefig(f"{OUT_PREFIX}.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT_PREFIX}.pdf", bbox_inches="tight")

print("Done.")
print(f"Saved:")
print(f"- {OUT_PREFIX}.png")
print(f"- {OUT_PREFIX}.pdf")
