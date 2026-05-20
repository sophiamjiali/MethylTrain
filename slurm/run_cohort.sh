#!/usr/bin/env bash

set -euo pipefail

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <cohort_name>"
    exit 1
fi

COHORT_NAME="$1"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=methyltrain_${COHORT_NAME}
#SBATCH --output=logs/cohorts/${COHORT_NAME}_%j.out
#SBATCH --error=logs/cohorts/${COHORT_NAME}_%j.err
#SBATCH --time=03:00:00
#SBATCH --nodes=2
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sophiamjia.li@mail.utoronto.ca

source ~/miniconda3/etc/profile.d/conda.sh
conda activate methyltrain-env

python -m methyltrain.cli cohort \
    --name ${PROJECT_NAME}
EOF