#!/usr/bin/env bash
PROJECT_NAME=$1
LOGS_DIR="/cluster/projects/kumargroup/sophia/logs/projects"

sbatch --job-name="methylP_${PROJECT_NAME}" \
       --output="${LOGS_DIR}/${PROJECT_NAME}_%j.out" \
       --error="${LOGS_DIR}/${PROJECT_NAME}_%j.err" \
       run_project.sh ${PROJECT_NAME}