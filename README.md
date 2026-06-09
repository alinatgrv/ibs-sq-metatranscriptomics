# Multi-omics analysis of sulfoquinovose metabolism in Irritable Bowel Syndrome

## Table of contents

- [Project overview](#project-overview)
- [Aim](#aim)
- [Objectives](#objectives)
- [Dataset](#dataset)
- [Prerequisites and computational requirements](#prerequisites-and-computational-requirements)
- [Repository structure](#repository-structure)
- [Workflow](#workflow)
- [How to reproduce the analysis](#how-to-reproduce-the-analysis)
- [Main results](#main-results)
- [Key conclusions](#key-conclusions)
- [Limitations](#limitations)
- [References](#references)

## Project overview

This repository contains the analysis workflow and selected results for a metatranscriptomic and multi-omics study of sulfoquinovose (SQ) metabolism in the gut microbiome of patients with irritable bowel syndrome (IBS).

Sulfoquinovose is a plant-derived sulfosugar that can be degraded by gut bacteria through several microbial pathways. Depending on the microbial community and downstream cross-feeding interactions, SQ degradation may contribute to different metabolic outputs, including short-chain fatty acid production or sulfur-associated metabolites such as hydrogen sulfide. The central motivation of this project was to investigate whether SQ-related microbial metabolism is detectable in public IBS gut microbiome data and whether it is associated with IBS status or IBS subtypes.

The analysis was based on public metatranscriptomic data from NCBI BioProject **PRJNA812699**, together with available metadata and metabolomics data from the same study.

## Aim

The main aim of the project was to analyze gut microbiome metatranscriptomic data in order to evaluate SQ-related microbial activity and its potential association with IBS.

## Objectives

1. Assess whether SQ-related metabolism is detectable in IBS metatranscriptomes.
2. Characterize SQ-associated enzymes and pathway models using targeted DIAMOND searches.
3. Identify taxonomic and metabolomic features associated with IBS.
4. Investigate links between sulfur-related metabolites and microbial taxa.
5. Evaluate the predictive value of SQ-related and multi-omics features for IBS classification.

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

Large raw sequencing files and full intermediate tables are not stored in this repository. The repository contains scripts, selected summary tables, selected figures. Full raw and intermediate data were stored on the HPC cluster.

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


## Repository structure

This repository is a cleaned reporting version of the analysis project.

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
├── scripts/
│   ├── preprocessing/
│   ├── profiling/
│   ├── differential_abundance/
│   ├── crossomics/
│   ├── machine_learning/
│   └── visualization/
│   └── targeted_sq/
└─── results/
    ├── figures/
    └── tables/
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


