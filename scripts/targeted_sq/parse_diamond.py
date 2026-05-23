"""
Parse, annotate, and filter DIAMOND blastx hits for targeted SQ-pathway analysis.

Inputs:
1. DIAMOND blastx TSV files
2. SQ protein reference FASTA
3. Homolog annotation table
4. Optional metadata table with sample condition

Outputs:
1. annotated_diamond_hits.tsv
2. filtered_diamond_hits.tsv
3. diamond_filtering_qc.tsv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))


from scripts.targeted_sq.fasta_utils import parse_fasta
from scripts.targeted_sq.sq_helpers import (
    extract_accession_from_sseqid,
    extract_sample_id,
    map_enzyme,
    map_pathway,
    normalize_gene_name,
    parse_condition,
)


diamond_col_names = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


num_cols = [
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


def read_diamond_tables(diamond_dir: str | Path) -> pd.DataFrame:
    """
    Read all DIAMOND TSV files from a directory.

    Input: 
    diamond_dir : str or pathlib.Path (directory with DIAMOND .tsv files).

    Returns: pandas.DataFrame (combined DIAMOND table).
    """
    diamond_dir = Path(diamond_dir)
    files = sorted(diamond_dir.glob("*.tsv"))

    if not files:
        raise FileNotFoundError(f"No .tsv files found in {diamond_dir}")

    dfs = []

    for file in files:
        try:
            df = pd.read_csv(file, sep="\t", names=diamond_col_names)
        except pd.errors.EmptyDataError:
            print(f"Warning: empty DIAMOND file skipped: {file}")
            continue

        if df.empty:
            print(f"Warning: empty DIAMOND file skipped: {file}")
            continue

        df["source_file"] = file.name
        dfs.append(df)

    if not dfs:
        raise ValueError(f"All DIAMOND files in {diamond_dir} are empty.")

    diamond = pd.concat(dfs, ignore_index=True)

    for column in num_cols:
        diamond[column] = pd.to_numeric(diamond[column], errors="coerce")

    diamond["sample_id"] = diamond["qseqid"].map(extract_sample_id)
    diamond["accession"] = diamond["sseqid"].map(extract_accession_from_sseqid)

    diamond = diamond.rename(
        columns={
            "pident": "identity",
            "length": "aln_len_aa",
            "bitscore": "bit_score",
        }
    )

    return diamond


def load_homolog_table(homolog_table: str | Path) -> pd.DataFrame:
    """
    Load homolog annotation table and normalize column names.

    Expected original columns may include:
    - ID_protein
    - Homolog_of_gene
    - Identity
    - Similarity
    - Organism
    """
    homolog_table = Path(homolog_table)
    df = pd.read_csv(homolog_table)

    rename_map = {
        "ID_protein": "accession",
        "Homolog_of_gene": "homolog_of_gene",
        "Identity": "homolog_identity",
        "Similarity": "homolog_similarity",
        "Organism": "homolog_organism",
    }

    existing_rename_map = {
        old: new for old, new in rename_map.items() if old in df.columns
    }

    df = df.rename(columns=existing_rename_map)

    required_columns = {"accession", "homolog_of_gene"}
    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise ValueError(
            "Homolog table is missing required columns after renaming: "
            f"{sorted(missing_columns)}. "
            "Expected at least ID_protein/accession and Homolog_of_gene/homolog_of_gene."
        )

    df["ref_protein_name"] = df["homolog_of_gene"].map(normalize_gene_name)
    df["sulfo_pathway"] = df["homolog_of_gene"].map(map_pathway)
    df["enzyme"] = df.apply(
        lambda row: map_enzyme(row["homolog_of_gene"]),
        axis=1,
    )

    return df


def load_metadata(metadata_csv: str | Path) -> pd.DataFrame:
    """
    Load metadata and return sample_id-condition table.

    The function expects:
    - Run column with sample IDs
    - gastrointest_disord column with IBS / healthy / control labels
    """
    metadata_csv = Path(metadata_csv)
    metadata = pd.read_csv(metadata_csv)

    if "Run" not in metadata.columns:
        raise ValueError("Metadata table must contain column 'Run'.")

    if "gastrointest_disord" not in metadata.columns:
        raise ValueError("Metadata table must contain column 'gastrointest_disord'.")

    metadata = metadata.rename(
        columns={
            "Run": "sample_id",
            "gastrointest_disord": "condition_raw",
        }
    )

    metadata["condition"] = metadata["condition_raw"].map(parse_condition)

    return metadata[["sample_id", "condition"]].drop_duplicates()


def build_annotated_hit_table(
    diamond_dir: str | Path,
    fasta_path: str | Path,
    homolog_table: str | Path,
    metadata_csv: str | Path | None = None,
    require_homolog_annotation: bool = True,
) -> pd.DataFrame:
    """
    Build annotated DIAMOND hit table.
    """
    diamond = read_diamond_tables(diamond_dir)
    fasta_annotation = parse_fasta(fasta_path)
    homolog_annotation = load_homolog_table(homolog_table)

    hit_table = (
        diamond
        .merge(fasta_annotation, on="accession", how="left")
        .merge(
            homolog_annotation,
            on="accession",
            how="left",
            suffixes=("", "_homolog"),
        )
    )

    hit_table["subject_cov_pct"] = np.where(
        hit_table["length_of_protein"].notna() & (hit_table["length_of_protein"] > 0),
        100 * hit_table["aln_len_aa"] / hit_table["length_of_protein"],
        np.nan,
    )

    if metadata_csv is not None:
        metadata = load_metadata(metadata_csv)

        hit_table = hit_table.merge(
            metadata,
            on="sample_id",
            how="left",
        )

    if require_homolog_annotation:
        hit_table = hit_table[hit_table["homolog_of_gene"].notna()].copy()

    return hit_table


def filter_diamond_hits(
    hit_table: pd.DataFrame,
    min_identity: float = 35.0,
    min_bitscore: float = 50.0,
    min_aln_len_aa: int = 30,
    min_subject_cov_pct: float = 15.0,
    keep_best_hit_per_read: bool = True,
) -> pd.DataFrame:
    """
    Filter DIAMOND hits and optionally keep the best hit per read.
    """
    mask = (
        (hit_table["identity"] >= min_identity)
        & (hit_table["bit_score"] >= min_bitscore)
        & (hit_table["aln_len_aa"] >= min_aln_len_aa)
        & (hit_table["subject_cov_pct"] >= min_subject_cov_pct)
    )

    filtered = hit_table.loc[mask].copy()

    if keep_best_hit_per_read:
        filtered = (
            filtered
            .sort_values(
                ["qseqid", "bit_score", "identity", "subject_cov_pct"],
                ascending=[True, False, False, False],
            )
            .drop_duplicates(subset=["qseqid"], keep="first")
            .copy()
        )

    return filtered


def make_filtering_qc(
    annotated_hits: pd.DataFrame,
    filtered_hits: pd.DataFrame,
) -> pd.DataFrame:
    """
    Create per-sample QC table for hit filtering.
    """
    raw_qc = (
        annotated_hits
        .groupby("sample_id")
        .agg(
            raw_hits=("qseqid", "size"),
            raw_unique_reads=("qseqid", "nunique"),
        )
        .reset_index()
    )

    filtered_qc = (
        filtered_hits
        .groupby("sample_id")
        .agg(
            filtered_hits=("qseqid", "size"),
            filtered_unique_reads=("qseqid", "nunique"),
        )
        .reset_index()
    )

    qc = raw_qc.merge(filtered_qc, on="sample_id", how="left")

    qc[["filtered_hits", "filtered_unique_reads"]] = (
        qc[["filtered_hits", "filtered_unique_reads"]]
        .fillna(0)
        .astype(int)
    )

    qc["hit_retention_pct"] = (
        100 * qc["filtered_hits"] / qc["raw_hits"]
    )

    qc["read_retention_pct"] = (
        100 * qc["filtered_unique_reads"] / qc["raw_unique_reads"]
    )

    return qc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Parse, annotate and filter DIAMOND hits for targeted SQ analysis."
    )

    parser.add_argument(
        "--diamond-dir",
        required=True,
        help="Directory with DIAMOND .tsv output files.",
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="SQ protein reference FASTA file.",
    )

    parser.add_argument(
        "--homolog-table",
        required=True,
        help="CSV table with SQ homolog annotation.",
    )

    parser.add_argument(
        "--metadata",
        default=None,
        help="Optional metadata CSV file with Run and gastrointest_disord columns.",
    )

    parser.add_argument(
        "--outdir",
        default="results/tables/targeted_sq",
        help="Output directory.",
    )

    parser.add_argument(
        "--min-identity",
        type=float,
        default=35.0,
        help="Minimum percentage identity.",
    )

    parser.add_argument(
        "--min-bitscore",
        type=float,
        default=50.0,
        help="Minimum DIAMOND bit score.",
    )

    parser.add_argument(
        "--min-aln-len-aa",
        type=int,
        default=30,
        help="Minimum alignment length in amino acids.",
    )

    parser.add_argument(
        "--min-subject-cov-pct",
        type=float,
        default=15.0,
        help="Minimum subject coverage percentage.",
    )

    parser.add_argument(
        "--keep-all-hits-per-read",
        action="store_true",
        help="Do not keep only the best hit per read.",
    )

    parser.add_argument(
        "--allow-missing-homolog-annotation",
        action="store_true",
        help="Keep hits without homolog annotation.",
    )

    parser.add_argument(
        "--save-annotated-hits",
        action="store_true",
        help="Save full annotated hit table. This file may be large.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("Building annotated DIAMOND hit table...")

    annotated_hits = build_annotated_hit_table(
        diamond_dir=args.diamond_dir,
        fasta_path=args.fasta,
        homolog_table=args.homolog_table,
        metadata_csv=args.metadata,
        require_homolog_annotation=not args.allow_missing_homolog_annotation,
    )

    print(f"Annotated hits: {len(annotated_hits):,}")

    print("Filtering DIAMOND hits...")

    filtered_hits = filter_diamond_hits(
        hit_table=annotated_hits,
        min_identity=args.min_identity,
        min_bitscore=args.min_bitscore,
        min_aln_len_aa=args.min_aln_len_aa,
        min_subject_cov_pct=args.min_subject_cov_pct,
        keep_best_hit_per_read=not args.keep_all_hits_per_read,
    )

    print(f"Filtered hits: {len(filtered_hits):,}")
    print(f"Filtered unique reads: {filtered_hits['qseqid'].nunique():,}")

    qc = make_filtering_qc(
        annotated_hits=annotated_hits,
        filtered_hits=filtered_hits,
    )

    filtered_out = outdir / "filtered_diamond_hits.tsv"
    qc_out = outdir / "diamond_filtering_qc.tsv"

    filtered_hits.to_csv(filtered_out, sep="\t", index=False)
    qc.to_csv(qc_out, sep="\t", index=False)

    print(f"Saved filtered hits to: {filtered_out}")
    print(f"Saved QC table to: {qc_out}")

    if args.save_annotated_hits:
        annotated_out = outdir / "annotated_diamond_hits.tsv"
        annotated_hits.to_csv(annotated_out, sep="\t", index=False)
        print(f"Saved annotated hits to: {annotated_out}")


if __name__ == "__main__":
    main()