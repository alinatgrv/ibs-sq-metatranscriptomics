#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=compute
#SBATCH --job-name=HUMAnN_mtx
#SBATCH --output=/home/alina_tgrv/beegfs/IBS_SQ/logs/humann_mtx_%A_%a.out
#SBATCH --error=/home/alina_tgrv/beegfs/IBS_SQ/logs/humann_mtx_%A_%a.err
#SBATCH --array=1-325%3

set -u
set -o pipefail

source /home/alina_tgrv/.pyenv/versions/miniconda3-3.12-24.7.1-0/etc/profile.d/conda.sh
conda activate /home/alina_tgrv/beegfs/conda_envs/humann39_env_fix

BASE=/home/alina_tgrv/beegfs/IBS_SQ
READS=$BASE/qc_kneaddata_metatranscriptome
META=$BASE/metadata
LOGS=$BASE/logs

METAPHLAN_RESULTS=$BASE/results/metaphlan_metatranscriptome
HUMANN_RESULTS=$BASE/results/humann_metatranscriptome

STATUS_DIR=$HUMANN_RESULTS/status
RUNLOG_DIR=$HUMANN_RESULTS/run_logs
TMPBASE=$BASE/tmp/humann_metatranscriptome

SAMPLE_LIST=$META/humann_samples_ready.txt
THREADS=${SLURM_CPUS_PER_TASK:-4}

mkdir -p "$LOGS" "$HUMANN_RESULTS" "$STATUS_DIR" "$RUNLOG_DIR" "$TMPBASE"

if [ ! -s "$SAMPLE_LIST" ]; then
  echo "[ERROR] Sample list not found or empty: $SAMPLE_LIST" >&2
  exit 1
fi

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

if [ -z "${SAMPLE:-}" ]; then
  echo "[ERROR] Could not resolve sample for task $SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

READS_DIR=$READS/$SAMPLE
R1=$READS_DIR/${SAMPLE}_paired_1.fastq
R2=$READS_DIR/${SAMPLE}_paired_2.fastq
TAX_PROFILE=$METAPHLAN_RESULTS/profiles/${SAMPLE}_profile.tsv

OUTDIR=$HUMANN_RESULTS/$SAMPLE
STATUS_FILE=$STATUS_DIR/${SAMPLE}.status
RUNLOG=$RUNLOG_DIR/${SAMPLE}.log

mkdir -p "$OUTDIR"

{
  echo "=================================================="
  echo "[INFO] Date: $(date)"
  echo "[INFO] Host: $(hostname)"
  echo "[INFO] JobID: ${SLURM_JOB_ID:-NA}"
  echo "[INFO] ArrayTaskID: ${SLURM_ARRAY_TASK_ID:-NA}"
  echo "[INFO] Sample: $SAMPLE"
  echo "[INFO] R1: $R1"
  echo "[INFO] R2: $R2"
  echo "[INFO] Tax profile: $TAX_PROFILE"
  echo "[INFO] Output dir: $OUTDIR"
  echo "[INFO] Threads: $THREADS"
  echo "[INFO] Conda env: /home/alina_tgrv/beegfs/conda_envs/humann39_env_fix"
  echo "[INFO] humann path: $(which humann || echo 'not_found')"
  echo "[INFO] humann version:"
  humann --version || true
  echo "=================================================="
} > "$RUNLOG"

if [ -s "$OUTDIR/${SAMPLE}_genefamilies.tsv" ]; then
  echo "[INFO] HUMAnN output already exists, skipping $SAMPLE" | tee -a "$RUNLOG"
  echo "DONE: exists" > "$STATUS_FILE"
  exit 0
fi

if [ ! -f "$R1" ]; then
  echo "[ERROR] Missing R1: $R1" | tee -a "$RUNLOG" >&2
  echo "FAILED: missing_r1" > "$STATUS_FILE"
  exit 1
fi

if [ ! -f "$R2" ]; then
  echo "[ERROR] Missing R2: $R2" | tee -a "$RUNLOG" >&2
  echo "FAILED: missing_r2" > "$STATUS_FILE"
  exit 1
fi

if [ ! -s "$R1" ]; then
  echo "[ERROR] Empty R1: $R1" | tee -a "$RUNLOG" >&2
  echo "FAILED: empty_r1" > "$STATUS_FILE"
  exit 1
fi

if [ ! -s "$R2" ]; then
  echo "[ERROR] Empty R2: $R2" | tee -a "$RUNLOG" >&2
  echo "FAILED: empty_r2" > "$STATUS_FILE"
  exit 1
fi

if [ ! -f "$TAX_PROFILE" ]; then
  echo "[ERROR] Missing tax profile: $TAX_PROFILE" | tee -a "$RUNLOG" >&2
  echo "FAILED: missing_tax_profile" > "$STATUS_FILE"
  exit 1
fi

if [ ! -s "$TAX_PROFILE" ]; then
  echo "[ERROR] Empty tax profile: $TAX_PROFILE" | tee -a "$RUNLOG" >&2
  echo "FAILED: empty_tax_profile" > "$STATUS_FILE"
  exit 1
fi

if [ -n "${SLURM_TMPDIR:-}" ] && [ -d "${SLURM_TMPDIR:-}" ]; then
  WORKDIR=$SLURM_TMPDIR/humann_${SAMPLE}
else
  WORKDIR=$TMPBASE/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SAMPLE}
fi

mkdir -p "$WORKDIR"

cleanup() {
  rm -rf "$WORKDIR"
}
trap cleanup EXIT

LOCAL_R1=$WORKDIR/${SAMPLE}_paired_1.fastq
LOCAL_R2=$WORKDIR/${SAMPLE}_paired_2.fastq
LOCAL_COMBINED=$WORKDIR/${SAMPLE}_combined.fastq
LOCAL_TAX=$WORKDIR/${SAMPLE}_profile.tsv
LOCAL_OUTDIR=$WORKDIR/output

mkdir -p "$LOCAL_OUTDIR"

echo "[INFO] Copying input files to local tmp" | tee -a "$RUNLOG"

if ! cp "$R1" "$LOCAL_R1"; then
  echo "[ERROR] Failed to copy R1 to tmp" | tee -a "$RUNLOG" >&2
  echo "FAILED: copy_r1" > "$STATUS_FILE"
  exit 1
fi

if ! cp "$R2" "$LOCAL_R2"; then
  echo "[ERROR] Failed to copy R2 to tmp" | tee -a "$RUNLOG" >&2
  echo "FAILED: copy_r2" > "$STATUS_FILE"
  exit 1
fi

if ! cp "$TAX_PROFILE" "$LOCAL_TAX"; then
  echo "[ERROR] Failed to copy tax profile to tmp" | tee -a "$RUNLOG" >&2
  echo "FAILED: copy_tax_profile" > "$STATUS_FILE"
  exit 1
fi

echo "[INFO] Combining paired FASTQ files into one input for HUMAnN" | tee -a "$RUNLOG"

if ! cat "$LOCAL_R1" "$LOCAL_R2" > "$LOCAL_COMBINED"; then
  echo "[ERROR] Failed to create combined FASTQ" | tee -a "$RUNLOG" >&2
  echo "FAILED: combine_fastq" > "$STATUS_FILE"
  exit 1
fi

if [ ! -s "$LOCAL_COMBINED" ]; then
  echo "[ERROR] Combined FASTQ is empty: $LOCAL_COMBINED" | tee -a "$RUNLOG" >&2
  echo "FAILED: empty_combined_fastq" > "$STATUS_FILE"
  exit 1
fi

echo "[INFO] Starting HUMAnN for $SAMPLE" | tee -a "$RUNLOG"

if humann \
    --input "$LOCAL_COMBINED" \
    --output "$LOCAL_OUTDIR" \
    --threads "$THREADS" \
    --taxonomic-profile "$LOCAL_TAX" >> "$RUNLOG" 2>&1
then
  echo "[INFO] HUMAnN finished for $SAMPLE" | tee -a "$RUNLOG"
else
  EXIT_CODE=$?
  echo "[ERROR] HUMAnN failed for $SAMPLE with exit code $EXIT_CODE" | tee -a "$RUNLOG" >&2
  echo "FAILED: humann_exit_${EXIT_CODE}" > "$STATUS_FILE"
  exit 1
fi

RAW_GF=$LOCAL_OUTDIR/${SAMPLE}_combined_genefamilies.tsv
RAW_PA=$LOCAL_OUTDIR/${SAMPLE}_combined_pathabundance.tsv
RAW_PC=$LOCAL_OUTDIR/${SAMPLE}_combined_pathcoverage.tsv

if [ ! -s "$RAW_GF" ]; then
  echo "[ERROR] Expected genefamilies output not found: $RAW_GF" | tee -a "$RUNLOG" >&2
  echo "FAILED: missing_genefamilies" > "$STATUS_FILE"
  exit 1
fi

echo "[INFO] Renaming HUMAnN outputs to sample-based names" | tee -a "$RUNLOG"

mv "$RAW_GF" "$LOCAL_OUTDIR/${SAMPLE}_genefamilies.tsv"
mv "$RAW_PA" "$LOCAL_OUTDIR/${SAMPLE}_pathabundance.tsv"
mv "$RAW_PC" "$LOCAL_OUTDIR/${SAMPLE}_pathcoverage.tsv"

echo "[INFO] Copying results back" | tee -a "$RUNLOG"

if ! cp -r "$LOCAL_OUTDIR"/. "$OUTDIR"/; then
  echo "[ERROR] Failed to copy HUMAnN outputs back" | tee -a "$RUNLOG" >&2
  echo "FAILED: copy_outputs_back" > "$STATUS_FILE"
  exit 1
fi

echo "DONE: success" > "$STATUS_FILE"
echo "[INFO] DONE $SAMPLE" | tee -a "$RUNLOG"