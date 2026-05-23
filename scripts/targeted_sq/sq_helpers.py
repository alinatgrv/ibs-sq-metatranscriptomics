"""
Helper functions for targeted SQ-pathway analysis.

This module contains small reusable functions used across the DIAMOND parsing,
SQ gene annotation, pathway mapping, and condition assignment steps.
"""

from __future__ import annotations
import pandas as pd


sq_pathway_map= {
    # sulfo-EMP
    "yihs": "sulfo-EMP",
    "yihr": "sulfo-EMP",
    "yihv": "sulfo-EMP",
    "sqia": "sulfo-EMP",
    "sqik": "sulfo-EMP",
    "yihq": "SQ hydrolysis",

    # sulfo-TAL 
    "sqod": "sulfo-TAL",
    "squd": "sulfo-TAL",
    "sqvf": "sulfo-TAL",
    "sqvg": "sulfo-TAL",
    "sqvh": "sulfo-TAL",

    # sulfo-TK 
    "sqwf": "sulfo-TK",
    "sqwg": "sulfo-TK",
    "sqwh": "sulfo-TK",
    "sqwi": "sulfo-TK",

    # sulfo-ED 
    "sedc": "sulfo-ED",
    "sedd": "sulfo-ED",
    "slab": "sulfo-ED",
    "sqve": "sulfo-ED",
}


enzyme_name_map = {
    "yihq": "sulfoquinovosidase",
    "yihr": "sulfoquinovose mutarotase",
    "yihs": "sulfoquinovose isomerase",
    "yihv": "sulfofructose kinase",
    "sqia": "sulfofructosephosphate aldolase",
    "sqik": "sulfofructose kinase",

    "sqod": "sulfolactaldehyde dehydrogenase / related enzyme",
    "squd": "sulfoacetaldehyde acetyltransferase / related enzyme",
    "sqve": "sulfolactate processing enzyme",

    "sqvf": "transaldolase-like sulfo pathway enzyme",
    "sqvg": "sulfo pathway enzyme",
    "sqvh": "sulfo pathway enzyme",

    "sqwf": "transketolase-like sulfo pathway enzyme",
    "sqwg": "sulfo pathway enzyme",
    "sqwh": "sulfo pathway enzyme",
    "sqwi": "sulfo pathway enzyme",

    "sedc": "sulfo-ED pathway enzyme",
    "sedd": "sulfo-ED pathway enzyme",
    "slab": "sulfo-ED pathway enzyme",
}


def extract_sample_id(qseqid: str) -> str:
    """
    Extract sample ID from DIAMOND query sequence ID.

    Example:
    SRR123456.1 to SRR123456
    """
    return str(qseqid).split(".")[0]


def normalize_gene_name(gene) -> str:
    """
    Normalize gene names to a consistent format.

    Examples:
    yihT to YihT
    sqiA to SqiA
    """
    if pd.isna(gene):
        return pd.NA

    gene = str(gene).strip()

    if not gene:
        return pd.NA

    if len(gene) == 1:
        return gene.upper()

    return gene[0].upper() + gene[1:-1].lower() + gene[-1].upper()


def extract_accession_from_sseqid(sseqid: str) -> str:
    """
    Extract UniProt accession from DIAMOND subject ID.

    Examples:
    tr|A0A7W8ULP4|A0A7W8ULP4_9HYPH convert to A0A7W8ULP4
    sp|P32141|SQUT_ECOLI convert to  P32141

    If the subject ID has another format, return it unchanged.
    """
    sseqid = str(sseqid)
    parts = sseqid.split("|")

    if len(parts) >= 3:
        return parts[1]

    return sseqid


def map_pathway(homolog_gene):
    """
    Map a homolog gene name to an SQ-related pathway.
    """
    if pd.isna(homolog_gene):
        return pd.NA

    gene = str(homolog_gene).lower()

    return sq_pathway_map.get(gene, pd.NA)


def map_enzyme(homolog_gene, protein_name=None):
    """
    Map a homolog gene name to a normalized enzyme name.

    If the gene is absent from the dictionary, return the original protein name
    when available.
    """
    if pd.notna(homolog_gene):
        gene = str(homolog_gene).lower()

        if gene in enzyme_name_map :
            return enzyme_name_map [gene]

    return protein_name if pd.notna(protein_name) else pd.NA


def parse_condition(value):
    """
    Convert metadata condition labels to IBS / HC.
    """
    value = str(value).lower()

    if "ibs" in value:
        return "IBS"

    if "healthy" in value or "control" in value:
        return "HC"

    return pd.NA