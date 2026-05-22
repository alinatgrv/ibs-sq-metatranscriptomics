#!/bin/bash
#SBATCH --job-name=SQ_diamond_blastx
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --array=1-326%10
#SBATCH --output=/beegfs/eoshchekuteva/IBS/logs/out/diamond_%A_%a.out
#SBATCH --error=/beegfs/eoshchekuteva/IBS/logs/err/diamond_%A_%a.err

set -euo pipefail

source /home/eoshchekuteva/.pyenv/versions/miniconda3-latest/etc/profile.d/conda.sh

conda activate diamond_env

cd /home/alina_tgrv/beegfs/IBS_SQ/qc_kneaddata

SAMPLES_LIST="/home/eoshchekuteva/beegfs/IBS/data/clean_and_trim/paired_rnaseq_path_list.txt"
DB="/beegfs/eoshchekuteva/IBS/db/SQ_homolog/SQ_homolog_prot_db.dmnd"
OUTDIR="/beegfs/eoshchekuteva/IBS/diamond_out"

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