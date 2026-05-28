#!/bin/bash
#SBATCH --output=/cluster/home/t144807uhn/logs/%x_%j.out
#SBATCH --error=/cluster/home/t144807uhn/logs//%x_%j.err
#SBATCH --time=01:55:00
#SBATCH --partition=build
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sophiamjia.li@mail.utoronto.ca

module load python3/3.12.11
source /cluster/home/t144807uhn/envs/methyltrain-env/bin/activate

python -m methyltrain.cli project --name $1
