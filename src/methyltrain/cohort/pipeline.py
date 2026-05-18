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
from methyltrain.io.read import load_annotation
from methyltrain.io.datasets import load_processed_project
from methyltrain.cohort.aggregation import aggregate_cohort, aggregate_genes
from methyltrain.cohort.filtering import mad_probe_filtering
from methyltrain.cohort.correction import batch_correction
from methyltrain.cohort.split import split
from methyltrain.cohort.features import winsorize, scale, extract_probe_set

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
    logger.info("Loaded %d projects:\n- %s",
                 len(project_list), "\n- ".join(project_list))
    logger.info("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    logger.info("-----| Beginning Pipeline |-----\n")

    projects = [load_processed_project(path) for path in project_list]
    cohort = aggregate_cohort(projects, config, layout)

    # Perform MAD probe filtering to reduce probe dimensionality
    if toggles_cfg.get('MAD_probe_filtering', True):
        cohort = mad_probe_filtering(cohort, config)
    
    # Perform batch effect correction across datasets
    if toggles_cfg.get('batch_correction', True):
        cohort = batch_correction(cohort, config)

    # Aggregate to the gene-level based on TSS200/1500 and/or gene bodies
    if toggles_cfg.get('gene_aggregation', True):

        # Load the appropriate annotation object
        annotation = load_annotation(
            platform = cohort.uns['data_source']['platform'], 
            reference_genome = cohort.uns['data_source']['reference_genome'], 
        )
        cohort = aggregate_genes(cohort, annotation, config)

    # Winsorize to increase model training fidelity
    if toggles_cfg.get('winsorize', True):
        cohort = winsorize(cohort, config)

    # Scale the numerical range to increase model training fidelity
    if toggles_cfg.get('scale', True):
        cohort = scale(cohort, config)

    # Split the cohort into train-val-test splits
    if toggles_cfg.get('split', True):
        train, val, test = split(cohort, config)
    else:
        train, val, test = None, None, None

    # Extract the final probe set and its order (as the index)
    probe_set = extract_probe_set(cohort)

    logger.info("=====================================================")

    return CohortResult(
        cohort_adata = cohort,
        train_adata = train,
        val_adata = val,
        test_adata = test,
        probe_set = probe_set
    )

# [END]