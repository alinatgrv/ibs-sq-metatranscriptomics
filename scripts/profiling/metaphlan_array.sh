#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=compute
#SBATCH --job-name=MetaPhlAn_IBS
#SBATCH --output=/home/alina_tgrv/beegfs/IBS_SQ/logs/metaphlan_%A_%a.out
#SBATCH --error=/home/alina_tgrv/beegfs/IBS_SQ/logs/metaphlan_%A_%a.err
#SBATCH --array=1-327%5

set -u
set -o pipefail

# conda inside sbatch
source /home/alina_tgrv/.pyenv/versions/miniconda3-3.12-24.7.1-0/etc/profile.d/conda.sh
conda activate /home/alina_tgrv/beegfs/conda_envs/humann39_env_fix

BASE=/home/alina_tgrv/beegfs/IBS_SQ
READS=$BASE/qc_kneaddata_metatranscriptome
META=$BASE/metadata
LOGS=$BASE/logs
RESULTS=$BASE/results/metaphlan_metatranscriptome
PROFILES=$RESULTS/profiles
STATUS=$RESULTS/status
RUNLOGS=$RESULTS/run_logs
TMPBASE=$BASE/tmp/metaphlan
DB=$BASE/databases/metaphlan
DB_INDEX=mpa_vJun23_CHOCOPhlAnSGB_202307

mkdir -p "$LOGS" "$PROFILES" "$STATUS" "$RUNLOGS" "$TMPBASE"

SAMPLE_LIST=$META/metaphlan_samples_metatranscriptome.txt

if [ ! -s "$SAMPLE_LIST" ]; then
  echo "[ERROR] Sample list not found or empty: $SAMPLE_LIST" >&2
  exit 1
fi

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

if [ -z "${SAMPLE:-}" ]; then
  echo "[ERROR] Could not resolve sample for task $SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

R1=$READS/$SAMPLE/${SAMPLE}_paired_1.fastq
R2=$READS/$SAMPLE/${SAMPLE}_paired_2.fastq

PROFILE_OUT=$PROFILES/${SAMPLE}_profile.tsv
STATUS_FILE=$STATUS/${SAMPLE}.status
RUNLOG=$RUNLOGS/${SAMPLE}.log

THREADS=${SLURM_CPUS_PER_TASK:-4}

{
  echo "=================================================="
  echo "[INFO] Date: $(date)"
  echo "[INFO] Host: $(hostname)"
  echo "[INFO] JobID: ${SLURM_JOB_ID:-NA}"
  echo "[INFO] ArrayTaskID: ${SLURM_ARRAY_TASK_ID:-NA}"
  echo "[INFO] Sample: $SAMPLE"
  echo "[INFO] R1: $R1"
  echo "[INFO] R2: $R2"
  echo "[INFO] DB: $DB"
  echo "[INFO] DB index: $DB_INDEX"
  echo "[INFO] Threads: $THREADS"
  echo "[INFO] Conda env: /home/alina_tgrv/beegfs/conda_envs/humann39_env_fix"
  echo "[INFO] metaphlan path: $(which metaphlan || echo 'not_found')"
  echo "[INFO] MetaPhlAn version:"
  metaphlan --version || true
  echo "=================================================="
} > "$RUNLOG"

if [ -s "$PROFILE_OUT" ]; then
  echo "[INFO] Output already exists, skipping $SAMPLE" | tee -a "$RUNLOG"
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

if [ ! -d "$DB" ]; then
  echo "[ERROR] MetaPhlAn DB directory not found: $DB" | tee -a "$RUNLOG" >&2
  echo "FAILED: missing_db_dir" > "$STATUS_FILE"
  exit 1
fi

if [ ! -f "$DB/${DB_INDEX}.pkl" ]; then
  echo "[ERROR] MetaPhlAn DB pickle not found: $DB/${DB_INDEX}.pkl" | tee -a "$RUNLOG" >&2
  echo "FAILED: missing_db_pkl" > "$STATUS_FILE"
  exit 1
fi

if [ ! -f "$DB/${DB_INDEX}.1.bt2l" ]; then
  echo "[ERROR] MetaPhlAn DB bowtie2 index not found: $DB/${DB_INDEX}.1.bt2l" | tee -a "$RUNLOG" >&2
  echo "FAILED: missing_db_bt2l" > "$STATUS_FILE"
  exit 1
fi

if [ -n "${SLURM_TMPDIR:-}" ] && [ -d "${SLURM_TMPDIR:-}" ]; then
  WORKDIR=$SLURM_TMPDIR/metaphlan_${SAMPLE}
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
LOCAL_PROFILE=$WORKDIR/${SAMPLE}_profile.tsv
LOCAL_BT2=$WORKDIR/${SAMPLE}.bowtie2.bz2

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

echo "[INFO] Starting MetaPhlAn for $SAMPLE" | tee -a "$RUNLOG"

if metaphlan "${LOCAL_R1},${LOCAL_R2}" \
    --input_type fastq \
    --bowtie2db "$DB" \
    --index "$DB_INDEX" \
    --bowtie2out "$LOCAL_BT2" \
    --nproc "$THREADS" \
    -o "$LOCAL_PROFILE" >> "$RUNLOG" 2>&1
then
  echo "[INFO] MetaPhlAn finished for $SAMPLE" | tee -a "$RUNLOG"
else
  EXIT_CODE=$?
  echo "[ERROR] MetaPhlAn failed for $SAMPLE with exit code $EXIT_CODE" | tee -a "$RUNLOG" >&2
  echo "FAILED: metaphlan_exit_${EXIT_CODE}" > "$STATUS_FILE"
  exit 1
fi

if [ ! -s "$LOCAL_PROFILE" ]; then
  echo "[ERROR] Empty output profile: $LOCAL_PROFILE" | tee -a "$RUNLOG" >&2
  echo "FAILED: empty_profile" > "$STATUS_FILE"
  exit 1
fi

echo "[INFO] Copying results back" | tee -a "$RUNLOG"

if ! cp "$LOCAL_PROFILE" "$PROFILE_OUT"; then
  echo "[ERROR] Failed to copy profile back" | tee -a "$RUNLOG" >&2
  echo "FAILED: copy_profile_back" > "$STATUS_FILE"
  exit 1
fi

echo "DONE: success" > "$STATUS_FILE"
echo "[INFO] DONE $SAMPLE" | tee -a "$RUNLOG"