#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=compute
#SBATCH --job-name=PRJNA812699
#SBATCH --output=/home/alina_tgrv/beegfs/IBS_SQ/logs/PRJNA812699_%j.out
#SBATCH --error=/home/alina_tgrv/beegfs/IBS_SQ/logs/PRJNA812699_%j.err

set -euo pipefail

# conda inside sbatch
source /home/alina_tgrv/.pyenv/versions/miniconda3-3.12-24.7.1-0/etc/profile.d/conda.sh
conda activate ibs_env

BASE=/home/alina_tgrv/beegfs/IBS_SQ
READS=$BASE/reads_folder
META=$BASE/metadata
TMP=$BASE/tmp

mkdir -p "$READS" "$META" "$TMP"
cd "$BASE"

# runs list (создаём, если ещё нет)
if [ ! -s "$META/runs.txt" ]; then
  echo "[INFO] Creating runinfo + runs.txt"
  esearch -db sra -query "PRJNA812699[BioProject]" | efetch -format runinfo > "$META/runinfo.csv"

  awk -F',' '
  NR==1 {for(i=1;i<=NF;i++) if($i=="Run") col=i; next}
  {if(col && $col ~ /^SRR/) print $col}
  ' "$META/runinfo.csv" > "$META/runs.txt"
fi

echo "[INFO] Runs found: $(wc -l < "$META/runs.txt")"

# download all runs with retries + resume
while read -r RUN; do
  [ -z "$RUN" ] && continue

  if ls "$READS/${RUN}"_*.fastq.gz >/dev/null 2>&1; then
    echo "[INFO] Skip $RUN (fastq.gz exists)"
    continue
  fi

  echo "[INFO] Processing $RUN"

  ok=0
  for attempt in 1 2 3; do
    echo "[INFO] prefetch attempt $attempt for $RUN"
    if prefetch --output-directory "$READS" "$RUN"; then
      ok=1
      break
    fi
    sleep 30
  done

  if [ "$ok" -ne 1 ]; then
    echo "[ERROR] prefetch failed for $RUN" | tee -a "$META/failed_prefetch.txt" >&2
    rm -rf "$READS/$RUN"
    continue
  fi

  if ! fasterq-dump \
      --split-files \
      --threads 4 \
      --outdir "$READS" \
      --temp "$TMP" \
      "$READS/$RUN/$RUN.sra"; then
    echo "[ERROR] fasterq-dump failed for $RUN" | tee -a "$META/failed_fasterq.txt" >&2
    rm -rf "$READS/$RUN"
    continue
  fi

  pigz -p 4 "$READS/${RUN}"_*.fastq
  rm -rf "$READS/$RUN"

done < "$META/runs.txt"

echo "[INFO] DONE"