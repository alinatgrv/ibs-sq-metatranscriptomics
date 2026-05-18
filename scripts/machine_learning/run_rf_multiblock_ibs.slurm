#!/bin/bash
#SBATCH --job-name=RF_IBS
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=logs/rf_ibs_%j.out
#SBATCH --error=logs/rf_ibs_%j.err

set -euo pipefail

cd /home/alina_tgrv/beegfs/IBS_SQ

export RF_N_REPEATS=50
export RF_NUM_TREES=1000
export RF_N_JOBS=${SLURM_CPUS_PER_TASK}

echo "Starting Random Forest IBS multiblock analysis"
echo "Date: $(date)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Repeats: ${RF_N_REPEATS}"
echo "Trees: ${RF_NUM_TREES}"

 /home/alina_tgrv/beegfs/conda_envs/rf_py_env/bin/python -u scripts/run_rf_multiblock_ibs.py

echo "Done"
echo "Date: $(date)"
