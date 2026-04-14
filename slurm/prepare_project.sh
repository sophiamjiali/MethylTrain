#!/bin/bash
#SBATCH --job-name=project
#SBATCH --output=/ddn_exa/campbell/sli/methyltrain/logs/projects/%x_%j.out
#SBATCH --error=/ddn_exa/campbell/sli/methyltrain/logs/projects/%x_%j.err
#SBATCH --time=12:00:00

#SBATCH --partition=gpu
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:0

set -euo pipefail

# Initialize command-line arguments
PROJECT_ID="$1"

if [[ -z "${PROJECT_ID}" ]]; then
    echo "Usage: sbatch slurm/prepare_project.sh <PROJECT_ID>"
    exit 1
fi

# Activate Conda Environment
source ~/miniforge3/etc/profile.d/conda.sh
conda activate methyltrain-env

# Navigate to the project directory
cd /ddn_exa/campbell/sli/methyltrain

# Run the sweep
srun python -u scripts/prepare_project.py \
    --config "/ddn_exa/campbell/sli/methyltrain/config/${PROJECT_ID}_config.yaml" \
    --verbose \
    --clean-raw-datas
