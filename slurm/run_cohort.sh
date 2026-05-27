#!/bin/bash
#SBATCH --time=01:55:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=build
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sophiamjia.li@mail.utoronto.ca


#!/bin/bash
#SBATCH --job-name=methyltrain_${COHORT_NAME}
#SBATCH --output=${LOGS_DIR}/cohorts/${COHORT_NAME}_%j.out
#SBATCH --error=${LOGS_DIR}/cohorts/${COHORT_NAME}_%j.err

module load python3/3.12.11
source /cluster/home/t144807uhn/envs/methyltrain-env/bin/activate

python -m methyltrain.cli cohort --name $1