# Workflow

The analysis workflow consisted of the following main steps:

1. SRA metadata collection and FASTQ download.
2. Quality control, trimming, and host read removal with KneadData.
3. Targeted SQ-related gene search with DIAMOND.
4. DIAMOND hit filtering, enzyme-level aggregation, and SQ score calculation.
5. Taxonomic profiling with MetaPhlAn.
6. Functional profiling with HUMAnN.
7. HUMAnN table joining, normalization, and QC.
8. SQ pathway checks in HUMAnN / MetaCyc output.
9. Species-level and pathway-level association testing with MaAsLin2.
10. Metabolomics association testing with MaAsLin2.
11. Cross-omics association analysis of sulfur-related metabolites and taxa.
12. Multi-block Random Forest classification.

The main scripts used for each step are stored in the `scripts/` directory.
