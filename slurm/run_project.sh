#!/bin/bash
#SBATCH --output=/cluster/home/t144807uhn/logs/projects/%x/%x_%j.out
#SBATCH --error=/cluster/home/t144807uhn/logs/projects/%x/%x_%j.err
#SBATCH --time=01:55:00
#SBATCH --partition=build
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sophiamjia.li@mail.utoronto.ca

# Make the project-specific logs directory
mkdir -p /cluster/home/t144807uhn/logs/projects/$1

# Activate the virtual environment
module load python3/3.12.11
source /cluster/home/t144807uhn/envs/methyltrain-env/bin/activate

# Resolve the configurations based on the provided job name
python -m methyltrain.cli project --name $1
