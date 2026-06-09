#!/bin/bash
# Generate a sample list for MetaPhlAn processing
# from paired-end FASTQ files.

set -euo pipefail

BASE=${BASE:-$(pwd)}

READS="${BASE}/qc_kneaddata"
META="${BASE}/metadata"

mkdir -p "$META"

OUT="${META}/metaphlan_samples.txt"

if [ ! -d "$READS" ]; then
  echo "[ERROR] Reads directory not found: $READS" >&2
  exit 1
fi

find "$READS" -mindepth 2 -maxdepth 2 -type f -name "*_paired_1.fastq" \
  | sed -E 's#.*/([^/]+)_paired_1\.fastq#\1#' \
  | sort -u > "$OUT"

echo "[INFO] Sample list written to: $OUT"
echo "[INFO] Total samples: $(wc -l < "$OUT")"
echo "[INFO] First samples:"
head "$OUT"