#!/usr/bin/env python3

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from sklearn.ensemble import RandomForestClassifier


# ==============================
# SETTINGS
# ==============================

INPUT_DIR = Path("results/random_forest/input_blocks")
RF_OUTDIR = Path("results/random_forest/ibs_multiblock_python")
SUMMARY_FILE = RF_OUTDIR / "rf_block_summary.tsv"

FIG_DIR = RF_OUTDIR / "figures"
TOP_FEATURE_DIR = FIG_DIR / "top_features_by_block"

FIG_DIR.mkdir(parents=True, exist_ok=True)
TOP_FEATURE_DIR.mkdir(parents=True, exist_ok=True)

N_JOBS = int(os.getenv("RF_N_JOBS", "6"))
N_TREES = int(os.getenv("RF_NUM_TREES_FOR_IMPORTANCE", "1000"))
TOP_N = int(os.getenv("RF_TOP_N_FEATURES", "15"))

META_COLS = {"sample", "subject_id", "metabolomics_id", "ibs_status"}

CATEGORY_COLORS = {
    "metabolites": "#4C78A8",
    "taxa": "#59A14F",
    "SQ-score": "#F28E2B",
    "clinical": "#E15759",
    "other": "#9D9D9D",
}

METRIC_LABELS = {
    "roc_auc_mean": "ROC-AUC",
    "pr_auc_mean": "PR-AUC",
    "balanced_accuracy_mean": "Balanced accuracy",
    "accuracy_mean": "Accuracy",
    "f1_mean": "F1",
}


plt.rcParams.update({
    "figure.dpi": 140,
    "savefig.dpi": 300,
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# ==============================
# HELPER FUNCTIONS
# ==============================

def pretty_block_name(block: str) -> str:
    name = block
    name = re.sub(r"^rf_234_", "", name)
    name = name.replace("all_taxa_sq_metabolites_clinical", "taxa + SQ + metabolites + clinical")
    name = name.replace("all_taxa_sq_metabolites", "taxa + SQ + metabolites")
    name = name.replace("taxa_metabolites_clinical", "taxa + metabolites + clinical")
    name = name.replace("taxa_metabolites", "taxa + metabolites")
    name = name.replace("metabolites_sq_clinical", "metabolites + SQ + clinical")
    name = name.replace("metabolites_sq", "metabolites + SQ")
    name = name.replace("taxa_sq_clinical", "taxa + SQ + clinical")
    name = name.replace("taxa_sq", "taxa + SQ")
    name = name.replace("metabolites_clinical", "metabolites + clinical")
    name = name.replace("taxa_clinical", "taxa + clinical")
    name = name.replace("sq_clinical", "SQ + clinical")
    name = name.replace("clinical_only", "clinical only")
    name = name.replace("metabolites_only", "metabolites only")
    name = name.replace("taxa_only", "taxa only")
    name = name.replace("sq_only", "SQ only")
    return name


def short_feature_name(feature: str, max_len: int = 55) -> str:
    x = feature
    x = re.sub(r"^met_", "met: ", x)
    x = re.sub(r"^tax__", "tax: ", x)
    x = re.sub(r"^sq__", "SQ: ", x)
    x = re.sub(r"^clin__", "clin: ", x)
    x = x.replace("_", " ")
    if len(x) > max_len:
        x = x[:max_len - 3] + "..."
    return x


def feature_category(feature: str) -> str:
    if feature.startswith("met_"):
        return "metabolites"
    if feature.startswith("tax__"):
        return "taxa"
    if feature.startswith("sq__"):
        return "SQ-score"
    if feature.startswith("clin__"):
        return "clinical"
    return "other"


def savefig(fig, path_base: Path):
    fig.savefig(path_base.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(path_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def load_block(block: str):
    path = INPUT_DIR / f"{block}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Missing input block: {path}")

    df = pd.read_csv(path, sep="\t")

    feature_cols = [c for c in df.columns if c not in META_COLS]
    X = df[feature_cols].apply(pd.to_numeric, errors="coerce")

    # impute missing values defensively
    X = X.fillna(X.median(numeric_only=True)).fillna(0)

    # remove constant columns, exactly as in RF script logic
    non_constant = X.nunique(dropna=False) > 1
    X = X.loc[:, non_constant]

    y = df["ibs_status"].map({"Control": 0, "IBS": 1})
    if y.isna().any():
        raise ValueError(f"Unexpected ibs_status values in {path}")

    return df, X, y.astype(int)


# ==============================
# LOAD SUMMARY
# ==============================

if not SUMMARY_FILE.exists():
    raise FileNotFoundError(
        f"Cannot find {SUMMARY_FILE}. First run scripts/run_rf_multiblock_ibs.py"
    )

summary = pd.read_csv(SUMMARY_FILE, sep="\t")
summary["pretty_block"] = summary["block"].apply(pretty_block_name)

numeric_cols = [
    "n_samples", "n_features",
    "roc_auc_mean", "roc_auc_sd",
    "pr_auc_mean", "pr_auc_sd",
    "balanced_accuracy_mean", "balanced_accuracy_sd",
    "accuracy_mean", "accuracy_sd",
    "f1_mean", "f1_sd",
]
for col in numeric_cols:
    if col in summary.columns:
        summary[col] = pd.to_numeric(summary[col], errors="coerce")

summary = summary.sort_values("roc_auc_mean", ascending=False).reset_index(drop=True)

summary.to_csv(FIG_DIR / "rf_block_summary_for_plotting.tsv", sep="\t", index=False)


# ==============================
# 1. ROC-AUC BY BLOCK
# ==============================

fig_h = max(6, 0.45 * len(summary))
fig, ax = plt.subplots(figsize=(10, fig_h))

y_pos = np.arange(len(summary))
ax.barh(
    y_pos,
    summary["roc_auc_mean"],
    xerr=summary["roc_auc_sd"],
    color="#4C78A8",
    alpha=0.9,
    error_kw={"elinewidth": 1, "capsize": 3},
)

ax.axvline(0.5, color="black", linestyle="--", linewidth=1, alpha=0.7)
ax.set_yticks(y_pos)
ax.set_yticklabels(summary["pretty_block"])
ax.invert_yaxis()
ax.set_xlim(0.4, 1.0)
ax.set_xlabel("ROC-AUC, mean ± SD")
ax.set_title("Random Forest IBS classification: comparison of input blocks")

for i, row in summary.iterrows():
    ax.text(
        row["roc_auc_mean"] + 0.01,
        i,
        f'{row["roc_auc_mean"]:.3f}',
        va="center",
        fontsize=8,
    )

savefig(fig, FIG_DIR / "01_rf_roc_auc_by_block")


# ==============================
# 2. BALANCED ACCURACY BY BLOCK
# ==============================

fig, ax = plt.subplots(figsize=(10, fig_h))

ax.barh(
    y_pos,
    summary["balanced_accuracy_mean"],
    xerr=summary["balanced_accuracy_sd"],
    color="#59A14F",
    alpha=0.9,
    error_kw={"elinewidth": 1, "capsize": 3},
)

ax.axvline(0.5, color="black", linestyle="--", linewidth=1, alpha=0.7)
ax.set_yticks(y_pos)
ax.set_yticklabels(summary["pretty_block"])
ax.invert_yaxis()
ax.set_xlim(0.4, 1.0)
ax.set_xlabel("Balanced accuracy, mean ± SD")
ax.set_title("Random Forest IBS classification: balanced accuracy by input block")

for i, row in summary.iterrows():
    ax.text(
        row["balanced_accuracy_mean"] + 0.01,
        i,
        f'{row["balanced_accuracy_mean"]:.3f}',
        va="center",
        fontsize=8,
    )

savefig(fig, FIG_DIR / "02_rf_balanced_accuracy_by_block")


# ==============================
# 3. METRIC HEATMAP
# ==============================

metric_cols = [
    "roc_auc_mean",
    "pr_auc_mean",
    "balanced_accuracy_mean",
    "accuracy_mean",
    "f1_mean",
]
metric_cols = [c for c in metric_cols if c in summary.columns]

heat = summary.set_index("pretty_block")[metric_cols]
heat_labels = [METRIC_LABELS.get(c, c) for c in metric_cols]

fig, ax = plt.subplots(figsize=(9, fig_h))
im = ax.imshow(heat.values, aspect="auto", vmin=0.45, vmax=1.0, cmap="viridis")

ax.set_xticks(np.arange(len(metric_cols)))
ax.set_xticklabels(heat_labels, rotation=30, ha="right")
ax.set_yticks(np.arange(len(heat.index)))
ax.set_yticklabels(heat.index)
ax.set_title("Random Forest performance metrics across input blocks")

for i in range(heat.shape[0]):
    for j in range(heat.shape[1]):
        val = heat.iloc[i, j]
        ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8, color="white")

cbar = fig.colorbar(im, ax=ax)
cbar.set_label("Metric value")

savefig(fig, FIG_DIR / "03_rf_metrics_heatmap")


# ==============================
# 4. FEATURE COMPOSITION BY BLOCK
# ==============================

composition_records = []

for block in summary["block"]:
    try:
        _, X, _ = load_block(block)
    except Exception as e:
        print(f"WARNING: cannot load block {block}: {e}")
        continue

    counts = pd.Series([feature_category(c) for c in X.columns]).value_counts().to_dict()

    record = {
        "block": block,
        "pretty_block": pretty_block_name(block),
        "n_features_after_constant_filter": X.shape[1],
    }
    for cat in CATEGORY_COLORS:
        record[cat] = counts.get(cat, 0)

    composition_records.append(record)

composition = pd.DataFrame(composition_records)
composition = composition.set_index("block").loc[summary["block"]].reset_index()
composition.to_csv(FIG_DIR / "rf_block_feature_composition.tsv", sep="\t", index=False)

composition_prop = composition.copy()
cat_cols = list(CATEGORY_COLORS.keys())
composition_prop[cat_cols] = composition_prop[cat_cols].div(
    composition_prop[cat_cols].sum(axis=1), axis=0
).fillna(0)

fig, ax = plt.subplots(figsize=(10, fig_h))
left = np.zeros(len(composition_prop))

for cat in cat_cols:
    vals = composition_prop[cat].values
    ax.barh(
        np.arange(len(composition_prop)),
        vals,
        left=left,
        label=cat,
        color=CATEGORY_COLORS[cat],
        alpha=0.9,
    )
    left += vals

ax.set_yticks(np.arange(len(composition_prop)))
ax.set_yticklabels(composition_prop["pretty_block"])
ax.invert_yaxis()
ax.set_xlim(0, 1)
ax.set_xlabel("Fraction of non-constant features")
ax.set_title("Feature composition of Random Forest input blocks")
ax.legend(loc="lower right", frameon=False)

savefig(fig, FIG_DIR / "04_rf_feature_block_composition")


# ==============================
# 5. ROC-AUC VS NUMBER OF FEATURES
# ==============================

fig, ax = plt.subplots(figsize=(8, 6))

ax.scatter(
    summary["n_features"],
    summary["roc_auc_mean"],
    s=90,
    color="#4C78A8",
    alpha=0.85,
    edgecolor="black",
    linewidth=0.5,
)

ax.set_xscale("log")
ax.axhline(0.5, color="black", linestyle="--", linewidth=1, alpha=0.7)
ax.set_xlabel("Number of features, log scale")
ax.set_ylabel("ROC-AUC")
ax.set_title("Model performance vs feature-space size")

for _, row in summary.iterrows():
    ax.text(
        row["n_features"] * 1.04,
        row["roc_auc_mean"],
        pretty_block_name(row["block"]),
        fontsize=7,
        va="center",
    )

savefig(fig, FIG_DIR / "05_rf_auc_vs_number_of_features")


# ==============================
# 6. FULL-DATA RF FEATURE IMPORTANCE FOR EACH BLOCK
# ==============================

all_importances = []

print("\nTraining full-data RF models for feature-importance plots...")
print("These importances are descriptive; validation metrics are taken from CV summary.\n")

for block in summary["block"]:
    print(f"Feature importance: {block}")

    _, X, y = load_block(block)

    model = RandomForestClassifier(
        n_estimators=N_TREES,
        max_features="sqrt",
        class_weight="balanced",
        random_state=42,
        n_jobs=N_JOBS,
    )

    model.fit(X, y)

    imp = pd.DataFrame({
        "block": block,
        "pretty_block": pretty_block_name(block),
        "feature": X.columns,
        "feature_pretty": [short_feature_name(c) for c in X.columns],
        "category": [feature_category(c) for c in X.columns],
        "importance": model.feature_importances_,
    }).sort_values("importance", ascending=False)

    imp["rank"] = np.arange(1, len(imp) + 1)

    all_importances.append(imp)

    top = imp.head(TOP_N).iloc[::-1]

    fig, ax = plt.subplots(figsize=(9, max(4.5, 0.35 * TOP_N)))

    colors = [CATEGORY_COLORS.get(cat, "#9D9D9D") for cat in top["category"]]

    ax.barh(
        np.arange(len(top)),
        top["importance"],
        color=colors,
        alpha=0.9,
    )

    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(top["feature_pretty"])
    ax.set_xlabel("Gini feature importance")
    ax.set_title(f"Top {TOP_N} features: {pretty_block_name(block)}")

    handles = [
        plt.Line2D([0], [0], marker="s", color="w",
                   markerfacecolor=CATEGORY_COLORS[c], label=c, markersize=9)
        for c in CATEGORY_COLORS
        if c in set(top["category"])
    ]
    ax.legend(handles=handles, frameon=False, loc="lower right")

    safe_block = re.sub(r"[^A-Za-z0-9_]+", "_", block)
    savefig(fig, TOP_FEATURE_DIR / f"top_features_{safe_block}")


all_imp = pd.concat(all_importances, ignore_index=True)
all_imp.to_csv(FIG_DIR / "rf_full_model_feature_importance_all.tsv", sep="\t", index=False)

top_imp = all_imp.query("rank <= @TOP_N").copy()
top_imp.to_csv(FIG_DIR / "rf_full_model_feature_importance_top.tsv", sep="\t", index=False)


# multipage PDF with all top-feature plots
pdf_path = FIG_DIR / "06_rf_top_features_all_blocks.pdf"

with PdfPages(pdf_path) as pdf:
    for block in summary["block"]:
        imp = top_imp[top_imp["block"] == block].sort_values("importance", ascending=True)

        fig, ax = plt.subplots(figsize=(9, max(4.5, 0.35 * len(imp))))
        colors = [CATEGORY_COLORS.get(cat, "#9D9D9D") for cat in imp["category"]]

        ax.barh(
            np.arange(len(imp)),
            imp["importance"],
            color=colors,
            alpha=0.9,
        )

        ax.set_yticks(np.arange(len(imp)))
        ax.set_yticklabels(imp["feature_pretty"])
        ax.set_xlabel("Gini feature importance")
        ax.set_title(f"Top {TOP_N} features: {pretty_block_name(block)}")

        handles = [
            plt.Line2D([0], [0], marker="s", color="w",
                       markerfacecolor=CATEGORY_COLORS[c], label=c, markersize=9)
            for c in CATEGORY_COLORS
            if c in set(imp["category"])
        ]
        ax.legend(handles=handles, frameon=False, loc="lower right")

        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

print(f"Written multipage PDF: {pdf_path}")


# ==============================
# 7. FINAL PRESENTATION FIGURE
# ==============================

best_block = summary.iloc[0]["block"]
best_name = pretty_block_name(best_block)
best_imp = top_imp[top_imp["block"] == best_block].sort_values("importance", ascending=True)

fig = plt.figure(figsize=(15, 10))
gs = fig.add_gridspec(
    2, 2,
    height_ratios=[1.2, 1],
    width_ratios=[1.25, 1],
    hspace=0.35,
    wspace=0.35,
)

# Panel A: ROC-AUC ranking
ax1 = fig.add_subplot(gs[:, 0])

ax1.barh(
    np.arange(len(summary)),
    summary["roc_auc_mean"],
    xerr=summary["roc_auc_sd"],
    color="#4C78A8",
    alpha=0.9,
    error_kw={"elinewidth": 1, "capsize": 3},
)

ax1.axvline(0.5, color="black", linestyle="--", linewidth=1, alpha=0.7)
ax1.set_yticks(np.arange(len(summary)))
ax1.set_yticklabels(summary["pretty_block"])
ax1.invert_yaxis()
ax1.set_xlim(0.4, 1.0)
ax1.set_xlabel("ROC-AUC, mean ± SD")
ax1.set_title("A. Predictive performance of data blocks")

for i, row in summary.iterrows():
    ax1.text(
        row["roc_auc_mean"] + 0.01,
        i,
        f'{row["roc_auc_mean"]:.3f}',
        va="center",
        fontsize=8,
    )

# Panel B: composition
ax2 = fig.add_subplot(gs[0, 1])

top_blocks_for_comp = composition_prop.copy()
top_blocks_for_comp = top_blocks_for_comp.set_index("block").loc[summary["block"]].reset_index()

left = np.zeros(len(top_blocks_for_comp))
for cat in cat_cols:
    vals = top_blocks_for_comp[cat].values
    ax2.barh(
        np.arange(len(top_blocks_for_comp)),
        vals,
        left=left,
        label=cat,
        color=CATEGORY_COLORS[cat],
        alpha=0.9,
    )
    left += vals

ax2.set_yticks(np.arange(len(top_blocks_for_comp)))
ax2.set_yticklabels(top_blocks_for_comp["pretty_block"], fontsize=7)
ax2.invert_yaxis()
ax2.set_xlim(0, 1)
ax2.set_xlabel("Fraction of features")
ax2.set_title("B. Feature composition")
ax2.legend(frameon=False, loc="lower right")

# Panel C: top features from best model
ax3 = fig.add_subplot(gs[1, 1])

colors = [CATEGORY_COLORS.get(cat, "#9D9D9D") for cat in best_imp["category"]]

ax3.barh(
    np.arange(len(best_imp)),
    best_imp["importance"],
    color=colors,
    alpha=0.9,
)

ax3.set_yticks(np.arange(len(best_imp)))
ax3.set_yticklabels(best_imp["feature_pretty"], fontsize=8)
ax3.set_xlabel("Gini feature importance")
ax3.set_title(f"C. Top features in best block:\n{best_name}")

handles = [
    plt.Line2D([0], [0], marker="s", color="w",
               markerfacecolor=CATEGORY_COLORS[c], label=c, markersize=9)
    for c in CATEGORY_COLORS
    if c in set(best_imp["category"])
]
ax3.legend(handles=handles, frameon=False, loc="lower right")

fig.suptitle(
    "Random Forest classification of IBS status using multi-omics input blocks",
    fontsize=15,
    y=0.98,
)

savefig(fig, FIG_DIR / "00_RF_final_presentation_summary")


# ==============================
# DONE
# ==============================

print("\nDone.")
print(f"Figures written to: {FIG_DIR}")
print("\nMain files:")
print(f"- {FIG_DIR / '00_RF_final_presentation_summary.png'}")
print(f"- {FIG_DIR / '00_RF_final_presentation_summary.pdf'}")
print(f"- {FIG_DIR / '01_rf_roc_auc_by_block.pdf'}")
print(f"- {FIG_DIR / '03_rf_metrics_heatmap.pdf'}")
print(f"- {FIG_DIR / '04_rf_feature_block_composition.pdf'}")
print(f"- {FIG_DIR / '06_rf_top_features_all_blocks.pdf'}")
print(f"- {FIG_DIR / 'rf_full_model_feature_importance_top.tsv'}")
