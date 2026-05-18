#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=compute
#SBATCH --output=/beegfs/alina_tgrv/IBS_SQ/logs/multiqc_%A.out
#SBATCH --error=/beegfs/alina_tgrv/IBS_SQ/logs/multiqc_%A.err
#SBATCH --job-name=multiqc-kneaddata-sq

source /home/eoshchekuteva/.pyenv/versions/miniconda3-latest/etc/profile.d/conda.sh

conda activate my_env

cd /beegfs/alina_tgrv/IBS_SQ/qc_kneaddata

multiqc -l paired_list.txt \
  --outdir multiqc_report_paired \
  --title "SQ MultiQC Report" \
  --filename SQ_data_multiqc_report_2.html \
  --force \
  --verbose

echo "MultiQC завершен"
