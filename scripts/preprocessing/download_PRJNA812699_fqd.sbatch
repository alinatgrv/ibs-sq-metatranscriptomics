#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=compute
#SBATCH --job-name=PRJNA812699_fqd
#SBATCH --output=/home/alina_tgrv/beegfs/IBS_SQ/logs/PRJNA812699_fqd_%j.out
#SBATCH --error=/home/alina_tgrv/beegfs/IBS_SQ/logs/PRJNA812699_fqd_%j.err

set -euo pipefail

source /home/alina_tgrv/.pyenv/versions/miniconda3-3.12-24.7.1-0/etc/profile.d/conda.sh
conda activate ibs_env

BASE=/home/alina_tgrv/beegfs/IBS_SQ
READS=$BASE/reads_folder_fqd
META=$BASE/metadata
TMP=$BASE/tmp_fqd

mkdir -p "$READS" "$META" "$TMP"
cd "$BASE"

if [ ! -s "$META/runs.txt" ]; then
  echo "[INFO] Creating runinfo + runs.txt"
  esearch -db sra -query "PRJNA812699[BioProject]" | efetch -format runinfo > "$META/runinfo.csv"

  awk -F',' '
  NR==1 {for(i=1;i<=NF;i++) if($i=="Run") col=i; next}
  {if(col && $col ~ /^SRR/) print $col}
  ' "$META/runinfo.csv" > "$META/runs.txt"
fi

echo "[INFO] Runs found: $(wc -l < "$META/runs.txt")"

while read -r RUN; do
  [ -z "$RUN" ] && continue

  if ls "$READS/${RUN}"_*.fastq.gz >/dev/null 2>&1; then
    echo "[INFO] Skip $RUN (fastq.gz exists)"
    continue
  fi

  echo "[INFO] Processing $RUN"

  ok=0
  for attempt in 1 2 3; do
    echo "[INFO] fasterq-dump attempt $attempt for $RUN"

    if fasterq-dump \
        --split-files \
        --threads 4 \
        --outdir "$READS" \
        --temp "$TMP" \
        --progress \
        "$RUN"; then
      ok=1
      break
    fi

    rm -f "$READS/${RUN}"_*.fastq "$READS/${RUN}"_*.fastq.gz 2>/dev/null || true
    sleep 30
  done

  if [ "$ok" -ne 1 ]; then
    echo "[ERROR] fasterq-dump failed for $RUN" | tee -a "$META/failed_fqd.txt" >&2
    continue
  fi

  pigz -p 4 "$READS/${RUN}"_*.fastq

done < "$META/runs.txt"

echo "[INFO] DONE"
