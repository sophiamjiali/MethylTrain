#!/bin/bash
#SBATCH --job-name=download_project
#SBATCH --output=/ddn_exa/campbell/sli/methyltrain/logs/download_project_%j.>
#SBATCH --error=/ddn_exa/campbell/sli/methyltrain-project/logs/download_project_%j.e>
#SBATCH --time=24:30:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodelist=gpu3

# Activate Conda Environment
source ~/miniforge3/etc/profile.d/conda.sh
conda activate methyltrain-env

# Navigate to the project directory
cd /ddn_exa/campbell/sli/methyltrain

# Run the sweep
srun python slurm/download_project.py \
    --config /ddn_exa/campbell/sli/methyltrain/config/TCGA-UVM_config.yaml \
    --verbose True