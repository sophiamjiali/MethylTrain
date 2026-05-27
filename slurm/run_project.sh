#!/usr/bin/env bash

set -euo pipefail

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <project_name>"
    exit 1
fi

PROJECT_NAME="$1"
LOGS_DIR="/cluster/projects/kumargroup/sophia/logs"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=methyltrain_${PROJECT_NAME}
#SBATCH --output=${LOGS_DIR}/projects/${PROJECT_NAME}_%j.out
#SBATCH --error=${LOGS_DIR}/projects/${PROJECT_NAME}_%j.err
#SBATCH --time=01:55:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=build
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sophiamjia.li@mail.utoronto.ca

module load python3/3.12.11
source /cluster/home/t144807uhn/envs/methyltrain-env/bin/activate

python -m methyltrain.cli project \
    --name ${PROJECT_NAME}
EOF