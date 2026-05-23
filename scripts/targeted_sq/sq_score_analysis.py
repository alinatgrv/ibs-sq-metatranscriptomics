"""
Compare targeted SQ pathway scores between IBS and HC groups.

Input:
sq_scores.tsv

Expected columns:
- sample_id
- condition
- pathway_model
- SQ_score

Outputs:
1. sq_score_group_summary.tsv
2. sq_score_ibs_hc_stats.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    """
    Apply Benjamini-Hochberg FDR correction.

    Input: p_values : pandas.Series.

    Returns: pandas.Series (FDR-adjusted q-values).
    """
    p_values = p_values.astype(float)
    q_values = pd.Series(np.nan, index=p_values.index, dtype=float)

    valid = p_values.notna()

    if valid.sum() == 0:
        return q_values

    p = p_values[valid]
    n = len(p)

    order = np.argsort(p.values)
    ranked_p = p.values[order]

    adjusted = ranked_p * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    adjusted_original_order = np.empty(n)
    adjusted_original_order[order] = adjusted

    q_values.loc[p.index] = adjusted_original_order

    return q_values


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """
    Calculate Cliff's delta effect size.

    Cliff's delta estimates how often values from group x are larger
    than values from group y, minus the reverse probability.

    Values range from -1 to 1.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]

    if len(x) == 0 or len(y) == 0:
        return np.nan

    greater = 0
    lower = 0

    for value in x:
        greater += np.sum(value > y)
        lower += np.sum(value < y)

    return (greater - lower) / (len(x) * len(y))


def summarize_by_group(
    sq_scores: pd.DataFrame,
    value_col: str = "SQ_score",
) -> pd.DataFrame:
    """
    Calculate descriptive statistics for each pathway and condition.
    """
    summary = (
        sq_scores
        .groupby(["pathway_model", "condition"], as_index=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            mean_score=(value_col, "mean"),
            median_score=(value_col, "median"),
            sd_score=(value_col, "std"),
            q1_score=(value_col, lambda x: x.quantile(0.25)),
            q3_score=(value_col, lambda x: x.quantile(0.75)),
            min_score=(value_col, "min"),
            max_score=(value_col, "max"),
        )
    )

    summary["iqr_score"] = summary["q3_score"] - summary["q1_score"]

    numeric_columns = [
        "mean_score",
        "median_score",
        "sd_score",
        "q1_score",
        "q3_score",
        "iqr_score",
        "min_score",
        "max_score",
    ]

    for column in numeric_columns:
        summary[column] = summary[column].round(6)

    return summary.sort_values(
        ["pathway_model", "condition"]
    ).reset_index(drop=True)


def compare_groups_by_pathway(
    sq_scores: pd.DataFrame,
    case_label: str = "IBS",
    control_label: str = "HC",
    value_col: str = "SQ_score",
) -> pd.DataFrame:
    """
    Compare SQ scores between case and control groups for each pathway model.
    """
    rows = []

    for pathway_model, sub in sq_scores.groupby("pathway_model"):
        case_values = (
            sub.loc[sub["condition"] == case_label, value_col]
            .dropna()
            .astype(float)
            .values
        )

        control_values = (
            sub.loc[sub["condition"] == control_label, value_col]
            .dropna()
            .astype(float)
            .values
        )

        n_case = len(case_values)
        n_control = len(control_values)

        case_median = np.nan
        control_median = np.nan
        median_difference = np.nan
        statistic = np.nan
        p_value = np.nan
        delta = np.nan

        if n_case > 0:
            case_median = float(np.median(case_values))

        if n_control > 0:
            control_median = float(np.median(control_values))

        if n_case > 0 and n_control > 0:
            median_difference = case_median - control_median
            delta = cliffs_delta(case_values, control_values)

        if n_case >= 2 and n_control >= 2:
            test = mannwhitneyu(
                case_values,
                control_values,
                alternative="two-sided",
            )
            statistic = float(test.statistic)
            p_value = float(test.pvalue)

        rows.append(
            {
                "pathway_model": pathway_model,
                "case_group": case_label,
                "control_group": control_label,
                "n_case": n_case,
                "n_control": n_control,
                "case_median": case_median,
                "control_median": control_median,
                "median_difference_case_minus_control": median_difference,
                "mannwhitney_u": statistic,
                "p_value": p_value,
                "cliffs_delta": delta,
            }
        )

    stats = pd.DataFrame(rows)

    stats["q_value_FDR"] = benjamini_hochberg(stats["p_value"])

    numeric_columns = [
        "case_median",
        "control_median",
        "median_difference_case_minus_control",
        "mannwhitney_u",
        "p_value",
        "q_value_FDR",
        "cliffs_delta",
    ]

    for column in numeric_columns:
        stats[column] = stats[column].round(6)

    stats["significant_FDR_0_05"] = stats["q_value_FDR"] < 0.05

    return stats.sort_values(
        ["q_value_FDR", "p_value", "pathway_model"],
        na_position="last",
    ).reset_index(drop=True)


def read_sq_scores(sq_scores_path: str | Path) -> pd.DataFrame:
    """
    Read SQ score table.
    """
    sq_scores_path = Path(sq_scores_path)

    if not sq_scores_path.exists():
        raise FileNotFoundError(f"File not found: {sq_scores_path}")

    sq_scores = pd.read_csv(sq_scores_path, sep="\t")

    required_columns = {
        "sample_id",
        "condition",
        "pathway_model",
        "SQ_score",
    }

    missing_columns = required_columns - set(sq_scores.columns)

    if missing_columns:
        raise ValueError(
            "SQ scores table is missing required columns: "
            f"{sorted(missing_columns)}"
        )

    sq_scores["SQ_score"] = pd.to_numeric(
        sq_scores["SQ_score"],
        errors="coerce",
    )

    sq_scores = sq_scores.dropna(
        subset=[
            "sample_id",
            "condition",
            "pathway_model",
            "SQ_score",
        ]
    ).copy()

    return sq_scores


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare SQ pathway scores between IBS and HC groups."
    )

    parser.add_argument(
        "--sq-scores",
        default="results/tables/targeted_sq/sq_scores.tsv",
        help="Long-format SQ score table.",
    )

    parser.add_argument(
        "--outdir",
        default="results/tables/targeted_sq",
        help="Output directory.",
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

    group_counts = (
        sq_scores[["sample_id", "condition"]]
        .drop_duplicates()
        ["condition"]
        .value_counts()
    )

    print("Samples per condition:")
    print(group_counts.to_string())

    print("Calculating group summaries...")
    group_summary = summarize_by_group(sq_scores)

    print("Running Mann-Whitney U tests...")
    stats = compare_groups_by_pathway(
        sq_scores=sq_scores,
        case_label=args.case_label,
        control_label=args.control_label,
    )

    group_summary_out = outdir / "sq_score_group_summary.tsv"
    stats_out = outdir / "sq_score_ibs_hc_stats.tsv"

    group_summary.to_csv(group_summary_out, sep="\t", index=False)
    stats.to_csv(stats_out, sep="\t", index=False)

    print(f"Saved group summary to: {group_summary_out}")
    print(f"Saved IBS vs HC statistics to: {stats_out}")

    significant = stats["significant_FDR_0_05"].sum()
    print(f"Significant pathway models after FDR < 0.05: {significant}")


if __name__ == "__main__":
    main()