#!/usr/bin/env bash
COHORT_NAME=$1
LOGS_DIR="/cluster/projects/kumargroup/sophia/logs/cohorts"

sbatch --job-name="methylC_${COHORT_NAME}" \
       --output="${LOGS_DIR}/${COHORT_NAME}_%j.out" \
       --error="${LOGS_DIR}/${COHORT_NAME}_%j.err" \
       --time=02:00:00 \
       --partition=build \
       run_cohort.sh ${COHORT_NAME}