"""
FASTA parsing utilities for targeted SQ-pathway analysis.

This module parses UniProt-like FASTA headers and extracts metadata required
for DIAMOND hit annotation and downstream SQ score calculation.
"""

from __future__ import annotations

import re
import pandas as pd
from pathlib import Path

from scripts.targeted_sq.sq_helpers import (
    extract_accession_from_sseqid,
    normalize_gene_name,
)


def parse_uniprot_like_header(header: str) -> dict:
    """
    Parse a UniProt-like FASTA header.

    Example:
    sp|P32141|SQUT_ECOLI Sulfofructosephosphate aldolase OS=Escherichia coli GN=yihT PE=1 SV=1

    Returns: dict
    Parsed metadata:
     - sseqid
     - accession
     - protein_name
     - gene_from_fasta
     - fasta_header
    """
    header = str(header).strip()

    if header.startswith(">"):
        header = header[1:]

    sseqid = header.split(" ", 1)[0]
    accession = extract_accession_from_sseqid(sseqid)

    protein_name_match = re.search(r"^[^\s]+\s+(.*?)\s+OS=", header)
    protein_name = (
        protein_name_match.group(1).strip()
        if protein_name_match
        else pd.NA
    )

    gene_match = re.search(r"\bGN=([^\s]+)", header)
    gene_name = (
        normalize_gene_name(gene_match.group(1))
        if gene_match
        else pd.NA
    )

    return {
        "sseqid": sseqid,
        "accession": accession,
        "protein_name": protein_name,
        "gene_from_fasta": gene_name,
        "fasta_header": header,
    }


def parse_fasta(fasta_path: str | Path) -> pd.DataFrame:
    """
    Parse a FASTA file and return a table with protein metadata.

    Returns: pandas.DataFrame (table with one row per FASTA record).

    Columns:
        - sseqid
        - accession
        - protein_name
        - gene_from_fasta
        - fasta_header
        - length_of_protein
    """
    fasta_path = Path(fasta_path)
    records = []

    with fasta_path.open("r", encoding="utf-8") as file:
        header = None
        seq_chunks = []

        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    record = parse_uniprot_like_header(header)
                    record["length_of_protein"] = len("".join(seq_chunks))
                    records.append(record)

                header = line
                seq_chunks = []

            else:
                seq_chunks.append(line)

        if header is not None:
            record = parse_uniprot_like_header(header)
            record["length_of_protein"] = len("".join(seq_chunks))
            records.append(record)

    return pd.DataFrame(records)