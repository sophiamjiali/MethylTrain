#!/bin/bash
#SBATCH --job-name=cohort
#SBATCH --output=/ddn_exa/campbell/sli/methyltrain/logs/cohorts/%x_%j.out
#SBATCH --error=/ddn_exa/campbell/sli/methyltrain/logs/cohorts/%x_%j.err
#SBATCH --time=12:00:00

#SBATCH --partition=gpu
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:0

set -euo pipefail

# Initialize command-line arguments
COHORT_ID="$1"

if [[ -z "${COHORT_ID}" ]]; then
    echo "Usage: sbatch slurm/prepare_cohort.sh <COHORT_ID>"
    exit 1
fi

# Activate Conda Environment
source ~/miniforge3/etc/profile.d/conda.sh
conda activate methyltrain-env

# Navigate to the project directory
cd /ddn_exa/campbell/sli/methyltrain

# Run the sweep
srun python -u scripts/prepare_cohort.py \
    --config "/ddn_exa/campbell/sli/methyltrain/config/${COHORT_ID}_config.yaml" \
    --verbose
