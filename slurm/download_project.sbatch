#!/bin/bash
#SBATCH --job-name=download_project
#SBATCH --output=/ddn_exa/campbell/sli/methyltrain/logs/download_project_%j.out
#SBATCH --error=/ddn_exa/campbell/sli/methyltrain/logs/download_project_%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=gpu
#SBATCH --nodelist=gpu3

set -euo pipefail

# Initialize command-line arguments
PROJECT_ID="$1"

if [[ -z "${PROJECT_ID}" ]]; then
    echo "Usage: sbatch download_project.sbatch <PROJECT_ID>"
    exit 1
fi

echo "Preparing project ${PROJECT_ID}"

# Activate Conda Environment
source ~/miniforge3/etc/profile.d/conda.sh
conda activate methyltrain-env

# Navigate to the project directory
cd /ddn_exa/campbell/sli/methyltrain

# Run the sweep
srun python -u slurm/download_project.py \
    --config "/ddn_exa/campbell/sli/methyltrain/config/${PROJECT_ID}_config.yaml" \
    --verbose \
    --clean-raw-data
