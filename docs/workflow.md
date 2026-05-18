# Workflow

The analysis workflow consisted of the following main steps:

1. SRA metadata collection and FASTQ download.
2. Quality control, trimming, and host read removal with KneadData.
3. Targeted SQ-related gene search with DIAMOND.
4. Taxonomic profiling with MetaPhlAn.
5. Functional profiling with HUMAnN.
6. HUMAnN table joining, normalization, and QC.
7. SQ pathway checks in HUMAnN / MetaCyc output.
8. Species-level and pathway-level association testing with MaAsLin2.
9. Metabolomics association testing with MaAsLin2.
10. Cross-omics association analysis of sulfur-related metabolites and taxa.
11. Multi-block Random Forest classification.

The main scripts used for each step are stored in the `scripts/` directory.
