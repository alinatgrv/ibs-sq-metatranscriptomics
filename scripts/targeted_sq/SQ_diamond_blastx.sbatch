#!/bin/bash

# Run DIAMOND blastx search against the SQ homolog database
# for preprocessed paired-end metatranscriptomic reads.
#
# Input:
#   - paired_rnaseq_path_list.txt
#   - SQ_homolog_prot_db.dmnd
#
# Output:
#   - DIAMOND tabular alignments (*.tsv)

#SBATCH --job-name=SQ_diamond_blastx
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --array=1-326%10
#SBATCH --output=logs/out/diamond_%A_%a.out
#SBATCH --error=logs/err/diamond_%A_%a.err

set -euo pipefail

source "${HOME}/.pyenv/versions/miniconda3-latest/etc/profile.d/conda.sh"

conda activate diamond_env

PROJECT_DIR=${PROJECT_DIR:-$(pwd)}

SAMPLES_LIST=${SAMPLES_LIST:-"${PROJECT_DIR}/data/clean_and_trim/paired_rnaseq_path_list.txt"}
DB=${DB:-"${PROJECT_DIR}/db/SQ_homolog/SQ_homolog_prot_db.dmnd"}
OUTDIR=${OUTDIR:-"${PROJECT_DIR}/diamond_out"}

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_LIST")
read -r r1 r2 <<< "$line"

base=$(basename "$r1" _paired_1.fastq)

cat "$r1" "$r2" | diamond blastx \
    --db "$DB" \
    --query - \
    --out "$OUTDIR/${base}.tsv" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    --evalue 1e-5 \
    --threads "${SLURM_CPUS_PER_TASK}"