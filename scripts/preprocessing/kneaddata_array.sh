#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=compute
#SBATCH --job-name=knead_IBS
#SBATCH --output=/beegfs/alina_tgrv/IBS_SQ/logs/knead_%A_%a.out
#SBATCH --error=/beegfs/alina_tgrv/IBS_SQ/logs/knead_%A_%a.err

set -euo pipefail

source /home/alina_tgrv/.pyenv/versions/miniconda3-3.12-24.7.1-0/etc/profile.d/conda.sh
conda activate ibs_env

cd /beegfs/alina_tgrv/IBS_SQ

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/paired_runs.txt)

R1="reads_folder/${SAMPLE}_1.fastq.gz"
R2="reads_folder/${SAMPLE}_2.fastq.gz"
OUTDIR="qc_kneaddata/${SAMPLE}"
LOGFILE="logs/${SAMPLE}_kneaddata.log"

mkdir -p qc_kneaddata logs

kneaddata \
  --input1 "$R1" \
  --input2 "$R2" \
  --reference-db databases/kneaddata \
  --output "$OUTDIR" \
  --output-prefix "$SAMPLE" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --processes 1 \
  --sequencer-source none \
  --run-fastqc-start \
  --run-fastqc-end \
  --run-trim-repetitive \
  --remove-intermediate-output \
  --log "$LOGFILE"