# Data

Large raw and intermediate data files are not stored in this repository.

## Raw metatranscriptomic data

Metatranscriptomic sequencing data were obtained from NCBI SRA under BioProject PRJNA812699.

Raw FASTQ/SRA files are not included in the repository because of their size. Download scripts are provided in:

- `scripts/preprocessing/download_PRJNA812699.sbatch`
- `scripts/preprocessing/download_PRJNA812699_fqd.sbatch`

## Metabolomics data

Untargeted fecal metabolomics data were obtained from the supplementary materials of Jacobs et al. (2023), Microbiome.

The full metabolomics abundance matrix is not included in the repository. Selected processed summary tables are available in:

- `results/tables/`

## Intermediate files

Large intermediate outputs from KneadData, MetaPhlAn, HUMAnN, MaAsLin2, and Random Forest preprocessing are not included.

Selected final figures and result tables required for the report are stored in:

- `results/figures/`
- `results/tables/`
