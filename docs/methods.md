## Methods

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
