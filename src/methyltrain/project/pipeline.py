# ==============================================================================
# Script:           project/pipeline.py
# Purpose:          Public orchestration API for the project pipeline
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-13
# ==============================================================================

import logging

from methyltrain.project.layout import ProjectLayout
from methyltrain.project.results import ProjectResult
from methyltrain.io.datasets import load_raw_project
from methyltrain.project.download import download_methylation
from methyltrain.project.metadata import prepare_metadata
from methyltrain.project.clean import clean_data
from methyltrain.project.qc import quality_control
from methyltrain.project.preprocess import preprocess

logger = logging.getLogger(__name__)


def prepare_project(config: dict, layout: ProjectLayout) -> ProjectResult:
    """
    Upper-level orchestration API for the full project pipeline. Downloads and 
    preprocesses DNA Methylation data for the project specified in the 
    configurations.

    Parameters
    ----------
    config : dict
        Configuration dictionary controlling workflow steps.
    layout : ProjectLayout
        Object representing a project dataset directory layout.

    Returns
    -------
    ProjectResult
        The wrapped project pipeline results.
    """

    # Fetch aliases for commonly accessed configuration subgroups
    project_cfg = config.get('project', {})
    paths_cfg = config.get('paths', {})

    logger.info("=====| MethylTrain: Project Processing Pipeline |=====\n")
    logger.info("~~~~~| Project Details |~~~~~")
    logger.info(f"Project ID: {project_cfg.get('id', '')}")
    logger.info(f"Seed: {project_cfg.get('seed', '')}")
    logger.info(f"GDC API Client Path: {paths_cfg.get('gdc_client', '')}")
    logger.info(f"Output Directory: {paths_cfg.get('output_dir', '')}")
    logger.info("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    logger.info("-----| Beginning Pipeline |-----\n")

    manifest, download_report = download_methylation(config, layout)

    # Extract the IDs of files who's beta values were downloaded
    full_ids = [(d['file_id'], d['file_name']) for d in download_report
                if d.get('download_status') == 1]
    file_ids = [d[0] for d in full_ids]

    metadata, meta_report, bio_report = prepare_metadata(config, file_ids)
    cleaning_report = clean_data(layout, full_ids)

    # Fetch metadata
    adata = load_raw_project(metadata, config, layout)
    adata, qc_report = quality_control(adata, config, layout)
    adata = preprocess(adata, config)

    logger.info("=====================================================\n")

    return ProjectResult(
        adata = adata,
        manifest = manifest,
        metadata = metadata,
        download_report = download_report,
        metadata_report = meta_report,
        biospecimen_report = bio_report,
        cleaning_report = cleaning_report,
        qc_report = qc_report
    )

# [END]