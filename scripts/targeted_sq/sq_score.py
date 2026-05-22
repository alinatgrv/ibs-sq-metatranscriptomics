"""
Compute targeted SQ pathway scores from filtered DIAMOND hits.

Input:
filtered_diamond_hits.tsv

Expected key columns:
- qseqid
- sample_id
- accession
- identity
- bit_score
- aln_len_aa
- subject_cov_pct
- ref_protein_name
- sulfo_pathway
- condition
- length_of_protein

Outputs:
1. sample_gene_abundance.tsv
2. sq_scores.tsv
3. sq_scores_wide.tsv

SQ score calculation:

For each sample and pathway model:
1. RPK_e = reads_e / protein_length_e(kb)
2. TPM_e = 10^6 * RPK_e / sum(RPK_all)
3. score_raw = mean(log1p(TPM_e) for core pathway steps)
4. coverage = number of detected core steps / number of core steps
5. SQ_score = score_raw * coverage
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


CORE_PATHWAYS = {
    "CORE_EMP1": [
        ["YihV"],
        ["YihT"],
        ["YihS"],
        ["YihQ"],
    ],

    "CORE_EMP2": [
        ["SqiA"],
        ["SqiK"],
        ["SqvD"],
        ["SqgA"],
    ],

    "CORE_EMP": [
        ["YihV", "SqiK"],
        ["YihT", "SqiA"],
        ["YihS", "SqvD"],
        ["YihQ", "SqgA"],
        ["YihU", "SlaB"],
    ],

    "CORE_EMP_alt": [
        ["YihV", "SqiK"],
        ["YihT", "SqiA"],
        ["YihS", "SqvD"],
        ["YihQ", "SqgA"],
        ["SqvB", "YihR"],
    ],

    "CORE_TAL": [
        ["SqvA"],
        ["YihS", "SqvD"],
        ["YihQ", "SqgA"],
        ["SqvB", "YihR"],
        ["YihU", "SlaB"],
    ],

    "CORE_TK_var1": [
        ["SqwGH"],
        ["SqwI"],
        ["YihS", "SqvD"],
        ["YihQ", "SqgA"],
    ],

    "CORE_TK_var2": [
        ["SqwGH"],
        ["SqwI"],
        ["YihS", "SqvD"],
        ["YihQ", "SqgA"],
        ["SqwF"],
        ["SqwD"],
    ],

    "CORE_ED_var1": [
        ["SedA"],
        ["SedB"],
        ["SedC"],
        ["SedD"],
        ["YihQ", "SqgA"],
    ],

    "CORE_ED_var2": [
        ["SedA"],
        ["SedB"],
        ["SedC"],
        ["SedD"],
        ["YihQ", "SqgA"],
        ["YihU", "SlaB"],
    ],

    "CORE_ASDO": [
        ["SquD"],
        ["YihQ", "SqgA"],
        ["SquF"],
    ],

    "CORE_ASDO_var2": [
        ["SquD"],
        ["YihQ", "SqgA"],
    ],

    "CORE_ASDO_var3": [
        ["SquD"],
        ["SquF"],
    ],

    "CORE_ASMO": [
        ["SqoD"],
        ["YihQ", "SqgA"],
        ["SquF"],
    ],
}


def read_filtered_hits(filtered_hits_path: str | Path) -> pd.DataFrame:
    """
    Read filtered DIAMOND hits.
    """
    filtered_hits_path = Path(filtered_hits_path)

    if not filtered_hits_path.exists():
        raise FileNotFoundError(f"File not found: {filtered_hits_path}")

    hits = pd.read_csv(filtered_hits_path, sep="\t")

    required_columns = {
        "qseqid",
        "sample_id",
        "accession",
        "identity",
        "bit_score",
        "ref_protein_name",
        "length_of_protein",
    }

    missing_columns = required_columns - set(hits.columns)

    if missing_columns:
        raise ValueError(
            "Filtered hits table is missing required columns: "
            f"{sorted(missing_columns)}"
        )

    numeric_columns = [
        "identity",
        "bit_score",
        "aln_len_aa",
        "subject_cov_pct",
        "length_of_protein",
    ]

    for column in numeric_columns:
        if column in hits.columns:
            hits[column] = pd.to_numeric(hits[column], errors="coerce")

    hits = hits.dropna(
        subset=[
            "qseqid",
            "sample_id",
            "ref_protein_name",
            "length_of_protein",
        ]
    ).copy()

    return hits


def aggregate_sample_gene(filtered_hits: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate filtered hits by sample and reference protein.

    One row corresponds to one sample and one SQ-related reference enzyme/gene.
    """
    group_columns = ["sample_id", "ref_protein_name"]

    agg_dict = {
        "read_count": ("qseqid", "nunique"),
        "total_hit_count": ("qseqid", "size"),
        "n_accessions_detected": ("accession", "nunique"),
        "identity_mean": ("identity", "mean"),
        "identity_median": ("identity", "median"),
        "bit_score_mean": ("bit_score", "mean"),
        "bit_score_median": ("bit_score", "median"),
        "protein_length_median": ("length_of_protein", "median"),
        "protein_length_mean": ("length_of_protein", "mean"),
    }

    if "subject_cov_pct" in filtered_hits.columns:
        agg_dict["subject_cov_mean"] = ("subject_cov_pct", "mean")
        agg_dict["subject_cov_median"] = ("subject_cov_pct", "median")

    sample_gene = (
        filtered_hits
        .groupby(group_columns, as_index=False)
        .agg(**agg_dict)
    )

    if "condition" in filtered_hits.columns:
        condition = (
            filtered_hits[["sample_id", "condition"]]
            .dropna()
            .drop_duplicates()
            .groupby("sample_id", as_index=False)
            .agg(condition=("condition", "first"))
        )

        sample_gene = sample_gene.merge(
            condition,
            on="sample_id",
            how="left",
        )

    if "sulfo_pathway" in filtered_hits.columns:
        pathway = (
            filtered_hits[["ref_protein_name", "sulfo_pathway"]]
            .dropna()
            .drop_duplicates()
            .groupby("ref_protein_name", as_index=False)
            .agg(sulfo_pathway=("sulfo_pathway", "first"))
        )

        sample_gene = sample_gene.merge(
            pathway,
            on="ref_protein_name",
            how="left",
        )

    numeric_round_columns = [
        "identity_mean",
        "identity_median",
        "bit_score_mean",
        "bit_score_median",
        "subject_cov_mean",
        "subject_cov_median",
        "protein_length_median",
        "protein_length_mean",
    ]

    for column in numeric_round_columns:
        if column in sample_gene.columns:
            sample_gene[column] = sample_gene[column].round(2)

    return sample_gene.sort_values(
        ["sample_id", "ref_protein_name"]
    ).reset_index(drop=True)


def add_tpm_like_abundance(sample_gene: pd.DataFrame) -> pd.DataFrame:
    """
    Add RPK and TPM-like normalized abundance.

    Protein length is represented by the median length of detected accessions
    within each reference protein group.
    """
    df = sample_gene.copy()

    df["protein_length_kb"] = df["protein_length_median"] / 1000

    df["RPK"] = np.where(
        df["protein_length_kb"] > 0,
        df["read_count"] / df["protein_length_kb"],
        np.nan,
    )

    df["sum_RPK_all"] = df.groupby("sample_id")["RPK"].transform("sum")

    df["TPM"] = np.where(
        df["sum_RPK_all"] > 0,
        1_000_000 * df["RPK"] / df["sum_RPK_all"],
        0.0,
    )

    df["protein_length_kb"] = df["protein_length_kb"].round(4)
    df["RPK"] = df["RPK"].round(6)
    df["TPM"] = df["TPM"].round(6)

    return df


def compute_sq_scores(
    sample_gene_tpm: pd.DataFrame,
    core_pathways: dict[str, list[list[str]]],
    tpm_presence_threshold: float = 0.0,
) -> pd.DataFrame:
    """
    Compute SQ scores for all samples and pathway models.

    Alternative enzymes within one pathway step are treated using OR logic:
    if several genes can represent the same step, the maximum TPM is used.
    """
    rows = []

    has_condition = "condition" in sample_gene_tpm.columns

    for sample_id, sub in sample_gene_tpm.groupby("sample_id"):
        gene_to_tpm = dict(zip(sub["ref_protein_name"], sub["TPM"]))

        condition = pd.NA

        if has_condition:
            conditions = sub["condition"].dropna().unique()

            if len(conditions) > 0:
                condition = conditions[0]

        for pathway_model, core_steps in core_pathways.items():
            step_tpms = []
            step_detected = []
            detected_step_names = []

            for step in core_steps:
                alternative_tpms = [
                    gene_to_tpm.get(gene, 0.0)
                    for gene in step
                ]

                step_tpm = max(alternative_tpms)
                detected = step_tpm > tpm_presence_threshold

                step_tpms.append(step_tpm)
                step_detected.append(detected)

                if detected:
                    detected_step_names.append("/".join(step))

            score_raw = float(np.mean([np.log1p(value) for value in step_tpms]))
            coverage = float(np.mean(step_detected))
            sq_score = score_raw * coverage

            rows.append(
                {
                    "sample_id": sample_id,
                    "condition": condition,
                    "pathway_model": pathway_model,
                    "n_core_steps": len(core_steps),
                    "n_detected_steps": int(sum(step_detected)),
                    "coverage": coverage,
                    "score_raw": score_raw,
                    "SQ_score": sq_score,
                    "detected_steps": ";".join(detected_step_names),
                }
            )

    sq_scores = pd.DataFrame(rows)

    numeric_columns = ["coverage", "score_raw", "SQ_score"]

    for column in numeric_columns:
        sq_scores[column] = sq_scores[column].round(6)

    return sq_scores.sort_values(
        ["sample_id", "pathway_model"]
    ).reset_index(drop=True)


def make_sq_scores_wide(sq_scores: pd.DataFrame) -> pd.DataFrame:
    """
    Convert long SQ score table to sample x pathway matrix.
    """
    wide = (
        sq_scores
        .pivot_table(
            index="sample_id",
            columns="pathway_model",
            values="SQ_score",
            aggfunc="mean",
        )
        .reset_index()
    )

    if "condition" in sq_scores.columns:
        condition = (
            sq_scores[["sample_id", "condition"]]
            .dropna()
            .drop_duplicates()
        )

        wide = condition.merge(wide, on="sample_id", how="right")

    return wide


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute SQ pathway scores from filtered DIAMOND hits."
    )

    parser.add_argument(
        "--filtered-hits",
        default="results/tables/targeted_sq/filtered_diamond_hits.tsv",
        help="Filtered DIAMOND hits table.",
    )

    parser.add_argument(
        "--outdir",
        default="results/tables/targeted_sq",
        help="Output directory.",
    )

    parser.add_argument(
        "--tpm-presence-threshold",
        type=float,
        default=0.0,
        help=(
            "TPM threshold for considering a core pathway step detected. "
            "Default: 0.0"
        ),
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("Reading filtered DIAMOND hits...")
    filtered_hits = read_filtered_hits(args.filtered_hits)
    print(f"Filtered hits: {len(filtered_hits):,}")
    print(f"Samples: {filtered_hits['sample_id'].nunique():,}")
    print(f"Reference proteins: {filtered_hits['ref_protein_name'].nunique():,}")

    print("Aggregating hits by sample and reference protein...")
    sample_gene = aggregate_sample_gene(filtered_hits)

    print("Computing TPM-like abundance...")
    sample_gene_tpm = add_tpm_like_abundance(sample_gene)

    print("Computing SQ scores...")
    sq_scores = compute_sq_scores(
        sample_gene_tpm=sample_gene_tpm,
        core_pathways=CORE_PATHWAYS,
        tpm_presence_threshold=args.tpm_presence_threshold,
    )

    sq_scores_wide = make_sq_scores_wide(sq_scores)

    sample_gene_out = outdir / "sample_gene_abundance.tsv"
    sq_scores_out = outdir / "sq_scores.tsv"
    sq_scores_wide_out = outdir / "sq_scores_wide.tsv"

    sample_gene_tpm.to_csv(sample_gene_out, sep="\t", index=False)
    sq_scores.to_csv(sq_scores_out, sep="\t", index=False)
    sq_scores_wide.to_csv(sq_scores_wide_out, sep="\t", index=False)

    print(f"Saved sample-gene abundance table to: {sample_gene_out}")
    print(f"Saved SQ scores to: {sq_scores_out}")
    print(f"Saved wide SQ score matrix to: {sq_scores_wide_out}")


if __name__ == "__main__":
    main()