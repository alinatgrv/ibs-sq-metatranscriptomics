# Methods

This document provides additional methodological details for the analysis.

## Preprocessing

Raw metatranscriptomic reads were processed with KneadData for read trimming and host read removal. FastQC reports were aggregated with MultiQC.

## Targeted SQ-related gene search

SQ-related genes were searched using a targeted DIAMOND-based approach with a custom SQ protein reference database.

## Taxonomic profiling

Taxonomic profiling was performed with MetaPhlAn. Species-level profiles were merged and used for downstream exploratory analysis and MaAsLin2 association testing.

## Functional profiling

Functional profiling was performed with HUMAnN. HUMAnN output tables were joined, normalized, filtered, and used for pathway-level analysis.

## Differential abundance analysis

MaAsLin2 was used for covariate-adjusted association testing of taxonomic, functional, and metabolomic features.

## Cross-omics analysis

Sulfur-related metabolites were tested for associations with species-level taxonomic profiles.

## Machine learning

Random Forest classification was used to compare clinical, taxonomic, SQ-score, metabolomic, and combined feature blocks.
