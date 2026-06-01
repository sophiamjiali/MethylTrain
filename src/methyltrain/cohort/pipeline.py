# ==============================================================================
# Script:           cohort/pipeline.py
# Purpose:          Public orchestration API for the cohort pipeline
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-15
# ==============================================================================

import logging

from methyltrain.cohort.layout import CohortLayout
from methyltrain.cohort.results import CohortResult
from methyltrain.cohort.construction import define_cohort
from methyltrain.cohort.features import define_feature_space
from methyltrain.cohort.alignment import align_distribution
from methyltrain.cohort.normalization import normalize, apply_normalization
from methyltrain.cohort.split import split

logger = logging.getLogger(__name__)


def prepare_cohort(config: dict, layout: CohortLayout) -> CohortResult:
    """
    Upper-level orchestration API for the full cohort pipeline. Downloads and 
    preprocesses DNA Methylation data for the project specified in the 
    configurations.

    Parameters
    ----------
    config : dict
        Configuration dictionary controlling workflow steps.
    layout : CohortLayout
        Object representing a cohort dataset directory layout.

    Returns
    -------
    CohortResult
        The wrapped cohort pipeline results.
    """

    # Fetch aliases for commonly accessed configuration subgroups
    cohort_cfg = config.get('cohort', {})
    paths_cfg = config.get('paths', {})
    toggles_cfg = config.get('toggles', {})
    project_list = config.get('project_list', [])

    logger.info("=====| MethylTrain: Cohort Processing Pipeline |=====\n")
    logger.info("~~~~~| Cohort Details |~~~~~")
    logger.info(f"Cohort ID: {cohort_cfg.get('id', '')}")
    logger.info(f"Seed: {cohort_cfg.get('seed', '')}")
    logger.info(f"Output Directory: {paths_cfg.get('output_dir', '')}")
    logger.info(f"Input Project Directory: {paths_cfg.get('project_dir', '')}")
    logger.info("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    logger.info("~~~~~| Input Project Details |~~~~~")
    logger.info(f"Loaded {len(project_list)} projects:")
    for p in project_list: logger.info(f"- {p}")
    logger.info("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    logger.info("-----| Beginning Pipeline |-----\n")

    # Define the cohort, optionally aggregating to the gene level
    cohort = define_cohort(project_list, config, layout)

    # Align the distributions across projects (batch correction)
    cohort = align_distribution(cohort, config)

    # Define the feature space, returning the final probe set
    cohort, probe_set = define_feature_space(cohort, config)

    if toggles_cfg.get('split', True):
        
        # Split the cohort into train-val-test splits (if toggled)
        train, val, test = split(cohort, config)

        # Apply normalization (winsorization and scaling)
        train, norm_state = normalize(train, config)
        val = apply_normalization(val, norm_state)
        test = apply_normalization(test, norm_state)
    else:
        
        # Else fit scaling on the full cohort
        train = val = test = None
        cohort = normalize(cohort, config)[0]

    logger.info("=====================================================")

    return CohortResult(
        cohort_adata = cohort,
        train_adata = train,
        val_adata = val,
        test_adata = test,
        probe_set = probe_set
    )

# [END]