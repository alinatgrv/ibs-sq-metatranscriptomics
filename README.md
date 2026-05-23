# Multi-omics analysis of sulfoglucose metabolism in Irritable Bowel Syndrome

## Project overview

This repository contains the analysis workflow and selected results for a metatranscriptomic and multi-omics study of sulfoquinovose (SQ) metabolism in the gut microbiome of patients with irritable bowel syndrome (IBS).

Sulfoquinovose is a plant-derived sulfosugar that can be degraded by gut bacteria through several microbial pathways. Depending on the microbial community and downstream cross-feeding interactions, SQ degradation may contribute to different metabolic outputs, including short-chain fatty acid production or sulfur-associated metabolites such as hydrogen sulfide. The central motivation of this project was to investigate whether SQ-related microbial metabolism is detectable in public IBS gut microbiome data and whether it is associated with IBS status or IBS subtypes.

The analysis was based on public metatranscriptomic data from NCBI BioProject **PRJNA812699**, together with available metadata and metabolomics data from the same study.

## Aim

The main aim of the project was to analyze gut microbiome metatranscriptomic data in order to evaluate SQ-related microbial activity and its potential association with IBS.

## Objectives

1. Download and organize public metatranscriptomic sequencing data.
2. Perform quality control, read trimming, and host read removal.
3. Generate taxonomic profiles using MetaPhlAn.
4. Generate functional profiles using HUMAnN.
5. Check whether known MetaCyc SQ degradation pathways are reconstructed in HUMAnN output.
6. Perform targeted DIAMOND-based search for SQ-related homologous genes in cleaned metatranscriptomic reads.
7. Aggregate DIAMOND hits at the enzyme level and calculate targeted SQ scores for predefined SQ pathway models.
8. Analyze species-level taxonomic differences between IBS and control samples.
9. Analyze pathway-level functional differences between IBS and control samples.
10. Integrate metabolomics data and test IBS-associated metabolite changes.
11. Explore associations between sulfur-related metabolites and bacterial species.
12. Compare clinical, taxonomic, SQ-score, and metabolomic feature blocks using Random Forest classification.

## Dataset

### Metatranscriptomics

- Source: NCBI Sequence Read Archive
- BioProject: **PRJNA812699**
- Initial number of sequencing runs: **1184**
- Main metatranscriptomic working cohort after preprocessing and metadata matching: **326 samples**
  - IBS: **207 samples**
  - Control: **119 samples**

### Metabolomics

- Metabolomics cohort: **368 samples**
  - IBS: **229 samples**
  - Control: **139 samples**
- Matched metatranscriptomics-metabolomics subset used for multi-block Random Forest analysis: **234 samples**
  - IBS: **128 samples**
  - Control: **106 samples**

Large raw sequencing files and full intermediate tables are not stored in this repository. The repository contains scripts, selected summary tables, selected figures. Full raw and intermediate data were stored on the HPC cluster under the main project directory:

```text
/home/alina_tgrv/beegfs/IBS_SQ
```

## Prerequisites and computational requirements
The workflow was designed for execution on a Linux HPC cluster with SLURM. The full analysis requires substantial disk space because raw SRA files, FASTQ files, KneadData outputs, HUMAnN/MetaPhlAn profiles, DIAMOND outputs, reference databases, and logs are generated as intermediate data.

| Resource | Requirement |
|---|---:|
| Operating system | Linux / HPC environment |
| Workflow manager | SLURM |
| CPU | up to 16 CPU threads per job |
| RAM | 32 GB minimum; 64 GB recommended for memory-intensive steps |
| Storage | at least 5 TB; 8–10 TB recommended for the full workflow |
| Environment management | Conda / Mamba |

Large raw sequencing files, full intermediate outputs, complete DIAMOND hit tables, and reference databases are not stored in this repository.


## Repository structure

This repository is a cleaned reporting version of the analysis project. Large raw data files, complete FASTQ files, full HUMAnN per-sample output directories, and large intermediate matrices are not intended to be stored in GitHub.

The recommended repository structure is:

```text
.
├── README.md
├── LICENSE
├── .gitignore
├── docs/
│   ├── methods.md
│   ├── workflow.md
│   └── report.md
├── envs/
│   └── software_versions.md
├── scripts/
│   ├── preprocessing/
│   │   ├── download_PRJNA812699.sbatch
│   │   ├── download_PRJNA812699_fqd.sbatch
│   │   ├── kneaddata_array.sh
│   │   └── multiqc.sbatch
│   ├── profiling/
│   │   ├── metaphlan_array.sh
│   │   ├── humann_array.sh
│   │   └── prepare_metaphlan_samples.sh
│   ├── differential_abundance/
│   │   ├── make_metadata_for_maaslin2.R
│   │   ├── make_metadata_for_maaslin2_species_subtypes.R
│   │   ├── make_metaphlan_species_for_maaslin2.R
│   │   ├── make_pathabundance_for_maaslin2.R
│   │   ├── make_metabolomics_for_maaslin2.R
│   │   ├── run_maaslin2_pathways.R
│   │   ├── run_maaslin2_species.R
│   │   ├── run_maaslin2_species_subtypes.R
│   │   ├── run_maaslin2_metabolomics.R
│   │   └── run_maaslin2_metabolomics_with_HAD.R
│   ├── crossomics/
│   │   ├── run_maaslin2_species_vs_sulfur_metabolites.R
│   │   └── plot_crossomics_sulfur_metabolites_summary.R
│   ├── machine_learning/
│   │   ├── make_rf_multiblock_input.R
│   │   ├── run_rf_multiblock_ibs.py
│   │   ├── run_rf_multiblock_ibs.slurm
│   │   ├── plot_rf_multiblock_results.py
│   │   └── plot_rf_metrics_heatmap_selected.py
│   └── visualization/
│   │   ├── make_qc_figures.py
│   │   ├── metaphlan_species_stats.py
│   │   ├── plot_metaphlan_species_overview.R
│   │   ├── plot_maaslin2_species_covariates_heatmap.R
│   │   ├── plot_maaslin2_metabolomics_results.R
│   │   ├── plot_sulfur_related_metabolites.R
│   │   ├── plot_metabolomics_presentation_figures.R
│   │   ├── plot_metabolomics_one_image.R
│   │   ├── plot_fig7_functional_taxonomic_profile.R
│   │   ├── plot_fig7_functional_taxonomic_profile_no_letters.R
│   │   ├── plot_section7_four_figures_separately.R
│   │   └── plot_maaslin2_covariates_ready.R
│   └── targeted_sq/
│       ├── SQ_diamond_blastx.sbatch
│       ├── fasta_utils.py
│       ├── parse_diamond.py
│       ├── sq_helpers.py
│       ├── sq_score_analysis.py
│       ├── sq_score_plots.py
│       └── sq_score.py
├── results/
│   ├── figures/
│   └── tables/
└── notebooks/
```

## Workflow

```text
SRA metadata and FASTQ download
        ↓
Quality control and host read removal with KneadData
        ↓
Targeted SQ-related gene search with DIAMOND
        ↓
DIAMOND hit filtering, enzyme-level aggregation, and SQ score calculation
        ↓
Taxonomic profiling with MetaPhlAn
        ↓
Functional profiling with HUMAnN
        ↓
HUMAnN table joining, normalization, and QC
        ↓
SQ pathway checks in HUMAnN / MetaCyc output
        ↓
Species-level and pathway-level association testing with MaAsLin2
        ↓
Metabolomics association testing with MaAsLin2
        ↓
Cross-omics association analysis of sulfur-related metabolites and taxa
        ↓
Random Forest multi-block classification
```

## Methods

### Software versions

The main software versions used in the analysis were:

| Tool | Version |
|---|---:|
| sra-tools (`prefetch`, `fasterq-dump`) | 3.2.1 |
| pigz | 2.8 |
| KneadData | 0.12.4 |
| Trimmomatic | 0.40 |
| FastQC | 0.12.1 |
| MultiQC | 1.35 |
| HUMAnN | 3.9 |
| MetaPhlAn | 4.1.1 |
| Bowtie2 | 2.5.5 |
| DIAMOND | 2.2.0 |
| R | 4.3.3 |
| MaAsLin2 | 1.18.0 |
| Python | 3.11.15 |
| scikit-learn | 1.8.0 |
| pandas | 3.0.3 |
| NumPy | 2.4.4 |
| SciPy | 1.17.1 |
| matplotlib | 3.10.9 |
| seaborn | 0.13.2 |
| scikit-learn | 1.8.0 |

Preprocessing tools were run from the `ibs_env` conda environment. HUMAnN and MetaPhlAn were run from the `humann39_env_fix` environment. MaAsLin2 analyses were run from `maaslin2_env`, and Random Forest models were run from `rf_py_env`. Targeted SQ DIAMOND searches were run from the `diamond_sq_env` environment. Targeted SQ-score parsing, filtering, statistical comparison, and visualization scripts were run using the Python environment listed above.

### Data download

Sequencing run accessions were obtained from NCBI BioProject PRJNA812699. Data were downloaded using `sra-tools`:

- `prefetch` for downloading `.sra` files;
- `fasterq-dump` for conversion to FASTQ;
- `pigz` for parallel compression.

The download and preprocessing steps were executed on an HPC cluster using SLURM jobs.

### Quality control and host read removal

Raw reads were processed using **KneadData**, which internally used Trimmomatic, Bowtie2, and FastQC.

Human reads were removed using the human genome Bowtie2 database downloaded with:

```bash
kneaddata_database --download human_genome bowtie2 databases/kneaddata/
```

The final preprocessing configuration used:

```bash
--sequencer-source none
--trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:74"
```

These parameters were selected after test runs showed no substantial adapter contamination in FastQC reports. Therefore, preset adapter trimming was disabled with `--sequencer-source none`. The minimum read length was set to 74 bp to retain sufficiently informative reads after trimming.

#### Read quality control

FastQC reports were aggregated with MultiQC after preprocessing. The plot below shows the mean per-base sequence quality across 1182 FastQC entries. Most read positions had high Phred quality scores, mostly around 35–40, indicating overall good read quality after preprocessing.

<p align="center">
  <img src="results/figures/qc_fastqc_mean_quality_scores.png" alt="FastQC mean quality scores across read positions" width="900">
</p>

<p align="center">
  <b>Figure 1.</b> MultiQC summary of FastQC mean per-base quality scores across preprocessed paired-end read files.
</p>

The full interactive MultiQC report is available here:

[Open full MultiQC report](docs/qc/SQ_data_multiqc_report.html)

### Taxonomic profiling

Taxonomic profiling was performed with **MetaPhlAn**. Individual sample profiles were merged into a single species-level table. The main downstream taxonomic matrix contained:

- 326 samples;
- 626 unique species-level taxa;
- IBS/control group labels.

Main input table for downstream taxonomic analysis:

```text
results/metaphlan_metatranscriptome/joined/species_sample_matrix_with_groups.tsv
```

### Functional profiling

Functional profiling was performed with **HUMAnN**. The following HUMAnN output tables were generated and merged across samples:

- `genefamilies.tsv`
- `pathabundance.tsv`
- `pathcoverage.tsv`

Merged tables were generated using `humann_join_tables`. Abundance tables were normalized to CPM using `humann_renorm_table`. Stratified and unstratified tables were separated for downstream analysis.

Main HUMAnN output files:

```text
results/humann_metatranscriptome/joined_all/genefamilies_all.tsv
results/humann_metatranscriptome/joined_all/pathabundance_all.tsv
results/humann_metatranscriptome/joined_all/pathcoverage_all.tsv
```

### HUMAnN QC

HUMAnN output quality was evaluated using the number of detected pathways and the abundance of `UNMAPPED` and `UNINTEGRATED` features.

Main QC summary:

```text
results/humann_metatranscriptome/qc_report_2026-04-13/tables/qc_overview.txt
```

QC categories included:

- `ok`: 276 samples;
- `borderline_low_pathways`: 24 samples;
- `hard_low_pathways`: 14 samples;
- additional samples with high `UNMAPPED` or `UNINTEGRATED` values.

For pathway-level MaAsLin2 analysis, the most problematic samples were excluded, resulting in a filtered cohort of 311 samples.

### SQ pathway search in HUMAnN output

Known MetaCyc SQ degradation pathway IDs were checked in HUMAnN outputs:

- `PWY-7446`: sulfoquinovose degradation I;
- `PWY-7722`: sulfoquinovose degradation II;
- `PWY-8213`;
- `PWY-8348`;
- `PWY-8349`;
- `PWY-8350`.

None of these SQ pathways were reconstructed in the final HUMAnN `pathabundance` or `pathcoverage` tables.

Additional checks showed that `PWY-7446` and `PWY-7722` were present in the internal HUMAnN MetaCyc database used in this project, but their supporting UniRef90 gene families were not detected in the merged unstratified HUMAnN gene-family table. The other SQ pathway IDs were not present in the local HUMAnN MetaCyc database version.

Therefore, direct pathway-level evidence for SQ degradation was not detected by HUMAnN in this cohort, and SQ-related analysis was moved toward targeted gene/SQ-score approaches.

### Targeted SQ-related gene search and SQ score calculation

Because known SQ degradation pathways were not reconstructed in HUMAnN output, an additional targeted read-level search was performed against a custom protein reference library of SQ-related homologs.

The custom reference library included homologous proteins associated with several SQ degradation route models, including sulfo-EMP, SQ hydrolysis, sulfo-TAL, sulfo-TK, and sulfo-ED-related pathways. Cleaned metatranscriptomic reads were aligned against this protein reference using DIAMOND `blastx`.

DIAMOND hits were parsed and annotated using protein metadata from the reference FASTA file and a homolog annotation table. Hits were filtered using the following criteria:

- minimum percentage identity;
- minimum bit score;
- minimum alignment length;
- minimum subject coverage.

Subject coverage was calculated as:
```text
subject_coverage = 100 × alignment_length_aa / protein_length_aa
```
When multiple hits were assigned to the same read, the best hit was retained based on bit score, sequence identity, and subject coverage.

Filtered hits were aggregated by sample and reference protein. For each sample and SQ-related enzyme, the number of uniquely assigned reads was counted. Enzyme-level abundance was then normalized using a TPM-like approach based on protein length:
```text
RPK_e = reads_e / protein_length_e(kb)

TPM_e = 10^6 × RPK_e / sum(RPK_all)
```

For each SQ pathway model, a targeted SQ score was calculated using predefined core enzymatic steps. Alternative enzymes within the same step were treated using OR logic, and the maximum TPM among alternative enzymes was used for that step.

The SQ score was calculated as:
```text
score_raw = mean(log1p(TPM_e) for core pathway steps)

coverage = detected_core_steps / total_core_steps

SQ_score = score_raw × coverage
```
This score was used as a targeted read-level proxy for SQ-related transcriptional signal. It should be interpreted as enzyme-level support for SQ-related pathway models rather than as direct evidence of complete pathway reconstruction.

### Differential abundance analysis with MaAsLin2

Differential abundance testing was performed using **MaAsLin2** with multivariable linear models.

For species-level microbiome data, the model was:

```text
species_abundance ~ ibs_status + age + sex + bmi + diet_grouped + ethnicity_grouped
```

For IBS subtype analysis, the model was:

```text
species_abundance ~ phenotype_group + age + sex + bmi + diet_grouped + ethnicity_grouped
```

For pathway-level HUMAnN analysis, the model was:

```text
pathway_abundance ~ ibs_status + age + sex + bmi
```

For metabolomics, the model was:

```text
metabolite_abundance ~ Group + Age + Sex + BMI + Race + Diet_Category + Batch_metabolomics
```

The general settings were:

- microbiome data: TSS normalization + LOG transformation;
- metabolomics data: no TSS normalization + LOG transformation;
- method: linear model;
- multiple testing correction: Benjamini-Hochberg FDR;
- significance threshold: q-value ≤ 0.25;
- prevalence filtering: features present in at least 10% of samples.

## How to reproduce the analysis

### 1. Clone the repository

```bash
git clone https://github.com/<username>/<repository-name>.git
cd <repository-name>
```

### 2. Create environments

Example:

```bash
mamba env create -f envs/maaslin2_env.yml
mamba env create -f envs/random_forest_env.yml
mamba env create -f envs/diamond_sq_env.yml
```

### 3. Run preprocessing and profiling scripts

Example SLURM scripts are stored in:

```text
scripts/preprocessing/
scripts/profiling/
```

### 4. Run targeted SQ-score analysis

Targeted SQ-related gene analysis is implemented in:
```text
scripts/targeted_sq/
```
Before running targeted SQ-score scripts, activate the corresponding conda environment:
```bash
mamba activate targeted_sq_env
```

Example DIAMOND hit parsing and filtering:
```bash
python scripts/targeted_sq/04_parse_filter_diamond_hits.py \
  --diamond-dir output/diamond_out \
  --fasta output/SQ_all_prot.faa \
  --homolog-table output/datasheets/all_homolog_sq_2.csv \
  --metadata data/SraRunTable.csv \
  --outdir results/tables/targeted_sq
```

Example SQ score calculation:
```bash
python scripts/targeted_sq/05_compute_sq_scores.py \
  --filtered-hits results/tables/targeted_sq/filtered_diamond_hits.tsv \
  --outdir results/tables/targeted_sq
```

Example IBS vs control comparison:
```bash
python scripts/targeted_sq/06_compare_sq_scores_ibs_hc.py \
  --sq-scores results/tables/targeted_sq/sq_scores.tsv \
  --outdir results/tables/targeted_sq
```

### 5. Run downstream analyses

Example:

```bash
Rscript scripts/differential_abundance/run_maaslin2_species.R
Rscript scripts/differential_abundance/run_maaslin2_species_subtypes.R
Rscript scripts/differential_abundance/run_maaslin2_metabolomics.R
Rscript scripts/crossomics/run_maaslin2_species_vs_sulfur_metabolites.R
python scripts/machine_learning/run_rf_multiblock_ibs.py
```

## Main results

### 1. HUMAnN pathway-level analysis

After QC filtering, 311 samples and 300 unstratified pathways were prepared for pathway-level MaAsLin2 analysis.

After prevalence filtering, 185 pathways remained in the model. No pathway-level associations with IBS status passed FDR correction at q ≤ 0.25.

This suggests that broad pathway-level HUMAnN profiles did not show robust IBS-associated differences under the tested model.

### 2. Targeted SQ-score analysis

Because HUMAnN did not reconstruct known SQ degradation pathways, a targeted DIAMOND-based SQ-score analysis was performed. Filtered read-level hits were aggregated by SQ-related reference enzymes and normalized using a TPM-like approach. SQ scores were then calculated for predefined pathway models based on core enzymatic step coverage and enzyme-level abundance.

The targeted SQ-score analysis provided enzyme-level evidence for SQ-related transcriptional signal in the dataset, but SQ-score distributions alone did not provide strong separation between IBS and control samples.

### 3. Species-level taxonomic overview

The species-level MetaPhlAn table contained 326 samples and 626 species-level taxa.

Exploratory visual analysis showed:

- similar global profiles for the most abundant species in IBS and control groups;
- strong overlap of Shannon diversity distributions between IBS and control samples;
- no clear separation between IBS and control samples in PCoA based on Bray-Curtis distance;
- several preliminary species candidates based on mean abundance and prevalence differences.

Candidate taxa included:

- `Bacteroides ovatus`
- `Blautia faecis`
- `Parabacteroides merdae`
- `Phascolarctobacterium faecium`
- `Parabacteroides distasonis`
- `Bacteroides caccae`
- `Alistipes onderdonkii`

These candidates were treated as exploratory before covariate-adjusted testing.

### 4. MaAsLin2 species-level IBS vs control analysis

The general IBS vs control model did not detect FDR-significant species-level associations with IBS status after correction for age, sex, BMI, diet, and ethnicity.

However, significant associations were detected for other covariates, including diet, sex, and ethnicity. This supports the need to adjust for these variables in microbiome association models.

Nominal IBS-associated species included:

- `Bacteroides xylanisolvens`
- `Gemmiger formicilis`
- `Alistipes finegoldii`
- `Bacteroides caccae`
- `Blautia faecis`
- `Eubacterium rectale`
- `Phocaeicola plebeius`

None of these passed FDR correction in the general IBS vs control model.

### 5. IBS subtype species-level analysis

IBS subtype analysis compared Control, IBS-C, IBS-D, and IBS-M groups. IBS-U samples were excluded because the subtype was unspecified.

The subtype model detected two FDR-significant species-level associations with IBS-M:

| Feature | Contrast | Coefficient | p-value | q-value |
|---|---:|---:|---:|---:|
| `Blautia faecis` | IBS-M vs Control | 1.696 | 0.000528 | 0.201 |
| `Vescimonas coprocola` | IBS-M vs Control | 1.445 | 0.001376 | 0.221 |

Both associations were positive, indicating higher abundance in IBS-M compared with controls after adjustment for age, sex, BMI, diet, and ethnicity.

### 6. Metabolomics IBS vs control analysis

The metabolomics MaAsLin2 model tested 601 metabolites across 368 samples.

For the IBS vs control contrast:

- 52 metabolites were nominally associated with IBS at p < 0.05;
- 43 metabolites were FDR-significant at q ≤ 0.25.

Metabolites decreased in IBS included:

- `N-delta-acetylornithine`
- `riboflavin (Vitamin B2)`
- `3-phenylpropionate (hydrocinnamate)`
- `1-methyladenine`
- `phenol sulfate`
- `indolepropionate`
- `indolelactate`
- `phenyllactate`

Metabolites increased in IBS included:

- `palmitate (16:0)`
- `androstenediol (3beta,17beta) disulfate (2)`
- `nervonate (24:1n9)`
- `I-urobilinogen`
- `margarate (17:0)`
- `N-palmitoylglycine`

Sulfur-related metabolites detected among IBS-associated results included:

| Metabolite | Direction in IBS | Coefficient | q-value |
|---|---:|---:|---:|
| `androstenediol (3beta,17beta) disulfate (2)` | Higher | 0.650 | 0.0777 |
| `phenol sulfate` | Lower | -0.640 | 0.229 |

These findings suggest changes in sulfur-associated metabolism, but they do not directly prove altered SQ degradation.

### 7. Cross-omics analysis of sulfur-related metabolites and taxa

Two sulfur-related metabolites were selected for species-metabolite association analysis:

- `phenol sulfate`
- `androstenediol (3beta,17beta) disulfate (2)`

The matched cross-omics dataset contained 234 samples.

For `phenol sulfate`, six species-level taxa passed FDR correction at q ≤ 0.25:

- `Mediterraneibacter faecis`
- `Coprococcus comes`
- `Fusicatenibacter saccharivorans`
- `Blautia massiliensis`
- `Clostridiaceae bacterium`
- `Oscillospiraceae bacterium CLA AA H250`

All significant associations were positive.

For `androstenediol disulfate`, no taxa passed FDR correction, although several nominal positive associations were observed.

### 8. Random Forest multi-block classification

Random Forest classification was used to compare the predictive value of different feature blocks:

- clinical features;
- SQ-score features;
- MetaPhlAn species-level taxa;
- metabolomics features;
- combinations of these blocks.

The matched dataset contained 234 samples.

Final model configuration:

- 50 repeated runs;
- 5-fold cross-validation;
- 1000 trees per model;
- 6 CPU threads.

The strongest classification performance was achieved by metabolomics-containing blocks:

| Feature block | ROC-AUC mean | ROC-AUC SD | Balanced accuracy |
|---|---:|---:|---:|
| metabolites + clinical | 0.898 | 0.046 | 0.825 |
| metabolites + SQ | 0.898 | 0.045 | 0.825 |
| taxa + metabolites + clinical | 0.898 | 0.045 | 0.825 |
| metabolites only | 0.898 | 0.046 | 0.826 |
| all taxa + SQ + metabolites | 0.898 | 0.046 | 0.825 |
| metabolites + SQ + clinical | 0.898 | 0.046 | 0.826 |

SQ-score alone performed poorly:

| Feature block | ROC-AUC mean | Balanced accuracy |
|---|---:|---:|
| SQ only | 0.473 | 0.488 |
| taxa only | 0.560 | 0.548 |
| clinical only | 0.753 | 0.708 |

The results suggest that metabolomics carried the strongest IBS classification signal in this matched dataset, while SQ-score alone was not sufficient for robust IBS classification.

## Key conclusions

1. Direct reconstruction of known SQ degradation pathways was not detected in HUMAnN pathway-level output.
2. General species-level IBS vs control analysis did not identify FDR-significant IBS-associated taxa after covariate adjustment.
3. IBS subtype analysis revealed two FDR-significant positive associations with IBS-M: `Blautia faecis` and `Vescimonas coprocola`.
4. Metabolomics showed the strongest IBS-associated signal, with 43 metabolites passing FDR correction.
5. Two sulfur-related metabolites, `phenol sulfate` and `androstenediol disulfate`, were associated with IBS, suggesting broader sulfur-associated metabolic differences.
6. Cross-omics analysis linked `phenol sulfate` to several bacterial taxa, including `Mediterraneibacter faecis`, `Coprococcus comes`, and `Fusicatenibacter saccharivorans`.
7. Random Forest models confirmed that metabolomics features provided the strongest predictive signal for IBS classification, while SQ-score alone had weak classification performance.

## Limitations

- The analysis is based on public cross-sectional data, so causal conclusions cannot be made.
- HUMAnN did not reconstruct known SQ pathways, which limits direct pathway-level interpretation of SQ degradation.
- Some newer MetaCyc SQ pathway IDs were not present in the local HUMAnN MetaCyc database version.
- Metabolomics and metatranscriptomics were available for only a matched subset of samples.
- Several important clinical covariates from the original study, such as anxiety scores, were not included in all current models.
- Some metabolomics-based Random Forest features may reflect diet, sweetener intake, or host metabolic background rather than IBS-specific biology.


## Main output files

Selected final outputs are stored in:

```text
results/tables/
results/figures/
```

The corresponding full analysis outputs on the HPC cluster were organized under:

```text
results/humann_metatranscriptome/
results/metaphlan_metatranscriptome/
results/metabolomics/
results/crossomics/
results/random_forest/
results/presentation/
results/presentation_figures/
```

Important tables include:

```text
results/humann_metatranscriptome/joined_all/sample_metadata_qc_326.tsv
results/humann_metatranscriptome/joined_all/SQ_pathway_uniref90_hits.tsv
results/humann_metatranscriptome/maaslin2_pathways_311_tss_log/all_results.tsv
results/humann_metatranscriptome/maaslin2_pathways_311_tss_log/significant_results.tsv
results/metaphlan_metatranscriptome/maaslin2_species_326_tss_log_diet_ethnicity_refwhite/significant_results.tsv
results/metaphlan_metatranscriptome/maaslin2_species_305_subtypes_tss_log_diet_ethnicity_refwhite/phenotype_group_significant_results.tsv
results/metabolomics/maaslin2_metabolites_368_group_tss_log_diet_batch/group_IBS_q025_results_with_original_names.tsv
results/crossomics/species_vs_sulfur_metabolites_maaslin2/species_vs_phenol_sulfate/phenol_sulfate_q025_taxa_results.tsv
results/crossomics/species_vs_sulfur_metabolites_maaslin2/summary_plots/crossomics_sulfur_q025_results_combined.tsv
results/random_forest/ibs_multiblock_python/rf_block_summary.tsv
```

Main targeted SQ output files:

```text
results/tables/targeted_sq/filtered_diamond_hits.tsv
results/tables/targeted_sq/diamond_filtering_qc.tsv
results/tables/targeted_sq/sample_gene_abundance.tsv
results/tables/targeted_sq/sq_scores.tsv
results/tables/targeted_sq/sq_scores_wide.tsv
results/tables/targeted_sq/sq_score_group_summary.tsv
results/tables/targeted_sq/sq_score_ibs_hc_stats.tsv
```

## Contributions

This repository contains the analysis performed for the student bioinformatics project. If any scripts, SQ-score files, or intermediate results were produced by collaborators, this should be explicitly indicated in the corresponding script header, table README, or repository documentation.

## References

### Dataset

Jacobs, J. P., Lagishetty, V., Hauer, M. C., Labus, J. S., Dong, T. S., Toma, R., Vuyisich, M., Naliboff, B. D., Lackner, J. M., Gupta, A., Tillisch, K., & Mayer, E. A. (2023).  
**Multi-omics profiles of the intestinal microbiome in irritable bowel syndrome and its bowel habit subtypes.**  
*Microbiome*, 11, 5.  
https://doi.org/10.1186/s40168-022-01450-5

NCBI BioProject / SRA accession used for metatranscriptomic sequencing data: **PRJNA812699**.

### Sulfoquinovose metabolism background

Hanson, B. T., Kits, K. D., Löffler, J., Burrichter, A. G., Fiedler, A., Denger, K., Frommeyer, B., Herbold, C. W., Rattei, T., Karcher, N., Segata, N., Schleheck, D., & Loy, A. (2021).  
**Sulfoquinovose is a select nutrient of prominent bacteria and a source of hydrogen sulfide in the human gut.**  
*The ISME Journal*, 15, 2779–2791.  
https://doi.org/10.1038/s41396-021-00968-0

Wei, Y., Tong, Y., & Zhang, Y. (2022).  
**New mechanisms for bacterial degradation of sulfoquinovose.**  
*Bioscience Reports*, 42(10), BSR20220314.  
https://doi.org/10.1042/BSR20220314

Krasenbrink, J., Hanson, B. T., Weiss, A. S., Borusak, S., Tanabe, T. S., Lang, M., Aichinger, G., Hausmann, B., Berry, D., Richter, A., Marko, D., Mussmann, M., Schleheck, D., Stecher, B., & Loy, A. (2025).  
**Sulfoquinovose is exclusively metabolized by the gut microbiota and degraded differently in mice and humans.**  
*Microbiome*, 13, 184.  
https://doi.org/10.1186/s40168-025-02175-x

### Software and documentation

KneadData documentation:  
https://huttenhower.sph.harvard.edu/kneaddata/

MetaPhlAn documentation:  
https://huttenhower.sph.harvard.edu/metaphlan/

HUMAnN documentation:  
https://huttenhower.sph.harvard.edu/humann/

MaAsLin2 documentation:  
https://github.com/biobakery/Maaslin2


