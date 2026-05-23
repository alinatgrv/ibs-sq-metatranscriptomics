"""
Plot targeted SQ pathway scores for IBS and HC groups.

Input:
sq_scores.tsv

Expected columns:
- sample_id
- condition
- pathway_model
- SQ_score

Optional input:
sq_score_ibs_hc_stats.tsv

If provided, FDR-significant pathway models can be marked on the plot.

Outputs:
1. sq_score_ibs_hc_boxplot.png
2. sq_score_ibs_hc_boxplot.pdf
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def read_sq_scores(sq_scores_path: str | Path) -> pd.DataFrame:
    """
    Read long-format SQ score table.
    """
    sq_scores_path = Path(sq_scores_path)

    if not sq_scores_path.exists():
        raise FileNotFoundError(f"File not found: {sq_scores_path}")

    df = pd.read_csv(sq_scores_path, sep="\t")

    required_columns = {
        "sample_id",
        "condition",
        "pathway_model",
        "SQ_score",
    }

    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise ValueError(
            "SQ score table is missing required columns: "
            f"{sorted(missing_columns)}"
        )

    df["SQ_score"] = pd.to_numeric(df["SQ_score"], errors="coerce")

    df = df.dropna(
        subset=[
            "sample_id",
            "condition",
            "pathway_model",
            "SQ_score",
        ]
    ).copy()

    return df


def read_stats(stats_path: str | Path | None) -> pd.DataFrame | None:
    """
    Read optional IBS vs HC statistics table.
    """
    if stats_path is None:
        return None

    stats_path = Path(stats_path)

    if not stats_path.exists():
        print(f"Warning: statistics file not found: {stats_path}")
        return None

    stats = pd.read_csv(stats_path, sep="\t")

    required_columns = {
        "pathway_model",
        "q_value_FDR",
    }

    missing_columns = required_columns - set(stats.columns)

    if missing_columns:
        print(
            "Warning: statistics table is missing required columns: "
            f"{sorted(missing_columns)}. Significance labels will be skipped."
        )
        return None

    stats["q_value_FDR"] = pd.to_numeric(
        stats["q_value_FDR"],
        errors="coerce",
    )

    return stats


def get_pathway_order(
    sq_scores: pd.DataFrame,
    order_by: str = "median",
) -> list[str]:
    """
    Determine pathway order for plotting.
    """
    if order_by == "alphabetical":
        return sorted(sq_scores["pathway_model"].unique())

    if order_by == "median":
        return (
            sq_scores
            .groupby("pathway_model")["SQ_score"]
            .median()
            .sort_values(ascending=True)
            .index
            .tolist()
        )

    if order_by == "mean":
        return (
            sq_scores
            .groupby("pathway_model")["SQ_score"]
            .mean()
            .sort_values(ascending=True)
            .index
            .tolist()
        )

    raise ValueError(
        "Unknown order_by value. Use one of: median, mean, alphabetical."
    )


def significance_label(q_value: float) -> str:
    """
    Convert FDR-adjusted q-value to significance label.
    """
    if pd.isna(q_value):
        return ""

    if q_value < 0.001:
        return "***"

    if q_value < 0.01:
        return "**"

    if q_value < 0.05:
        return "*"

    return ""


def add_significance_labels(
    ax,
    sq_scores: pd.DataFrame,
    stats: pd.DataFrame,
    order: list[str],
) -> None:
    """
    Add FDR significance labels to the right side of the plot.
    """
    if stats is None:
        return

    stats = stats.copy()
    stats["label"] = stats["q_value_FDR"].map(significance_label)

    label_map = dict(zip(stats["pathway_model"], stats["label"]))

    x_max = sq_scores["SQ_score"].max()

    if pd.isna(x_max) or x_max == 0:
        x_max = 1.0

    x_position = x_max * 1.05

    for y_position, pathway_model in enumerate(order):
        label = label_map.get(pathway_model, "")

        if label:
            ax.text(
                x_position,
                y_position,
                label,
                va="center",
                ha="left",
                fontsize=11,
                fontweight="bold",
            )

    ax.set_xlim(right=x_max * 1.15)


def plot_sq_score_boxplot(
    sq_scores: pd.DataFrame,
    stats: pd.DataFrame | None,
    output_prefix: str | Path,
    case_label: str = "IBS",
    control_label: str = "HC",
    order_by: str = "median",
    width: float = 10,
    height: float = 7,
) -> None:
    """
    Plot SQ score distributions by condition.
    """
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    plot_df = sq_scores.copy()

    hue_order = [control_label, case_label]
    hue_order = [
        condition
        for condition in hue_order
        if condition in plot_df["condition"].unique()
    ]

    order = get_pathway_order(plot_df, order_by=order_by)

    sns.set_theme(style="whitegrid", context="notebook")

    fig, ax = plt.subplots(figsize=(width, height))

    sns.boxplot(
        data=plot_df,
        y="pathway_model",
        x="SQ_score",
        hue="condition",
        order=order,
        hue_order=hue_order,
        showfliers=False,
        ax=ax,
    )

    sns.stripplot(
        data=plot_df,
        y="pathway_model",
        x="SQ_score",
        hue="condition",
        order=order,
        hue_order=hue_order,
        dodge=True,
        alpha=0.35,
        size=2.5,
        jitter=0.2,
        legend=False,
        ax=ax,
    )

    add_significance_labels(
        ax=ax,
        sq_scores=plot_df,
        stats=stats,
        order=order,
    )

    ax.set_xlabel("SQ score")
    ax.set_ylabel("")
    ax.set_title("Targeted SQ pathway scores in IBS and HC samples")

    handles, labels = ax.get_legend_handles_labels()

    if handles:
        ax.legend(
            handles[: len(hue_order)],
            labels[: len(hue_order)],
            title="Condition",
            frameon=True,
            loc="lower right",
        )

    sns.despine(ax=ax, left=False, bottom=False)

    fig.tight_layout()

    png_path = output_prefix.with_suffix(".png")
    pdf_path = output_prefix.with_suffix(".pdf")

    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")

    plt.close(fig)

    print(f"Saved plot to: {png_path}")
    print(f"Saved plot to: {pdf_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot targeted SQ pathway scores for IBS and HC groups."
    )

    parser.add_argument(
        "--sq-scores",
        default="results/tables/targeted_sq/sq_scores.tsv",
        help="Long-format SQ score table.",
    )

    parser.add_argument(
        "--stats",
        default="results/tables/targeted_sq/sq_score_ibs_hc_stats.tsv",
        help=(
            "Optional IBS vs HC statistics table. "
            "Used for FDR significance labels if available."
        ),
    )

    parser.add_argument(
        "--outdir",
        default="results/figures/targeted_sq",
        help="Output directory for figures.",
    )

    parser.add_argument(
        "--case-label",
        default="IBS",
        help="Case group label. Default: IBS.",
    )

    parser.add_argument(
        "--control-label",
        default="HC",
        help="Control group label. Default: HC.",
    )

    parser.add_argument(
        "--order-by",
        default="median",
        choices=["median", "mean", "alphabetical"],
        help="Pathway ordering method.",
    )

    parser.add_argument(
        "--width",
        type=float,
        default=10,
        help="Figure width in inches.",
    )

    parser.add_argument(
        "--height",
        type=float,
        default=7,
        help="Figure height in inches.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("Reading SQ score table...")
    sq_scores = read_sq_scores(args.sq_scores)

    print(f"Rows: {len(sq_scores):,}")
    print(f"Samples: {sq_scores['sample_id'].nunique():,}")
    print(f"Pathway models: {sq_scores['pathway_model'].nunique():,}")

    stats = read_stats(args.stats)

    output_prefix = outdir / "sq_score_ibs_hc_boxplot"

    print("Plotting SQ score boxplot...")

    plot_sq_score_boxplot(
        sq_scores=sq_scores,
        stats=stats,
        output_prefix=output_prefix,
        case_label=args.case_label,
        control_label=args.control_label,
        order_by=args.order_by,
        width=args.width,
        height=args.height,
    )


if __name__ == "__main__":
    main()