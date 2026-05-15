# ==============================================================================
# Script:           project/pipeline.py
# Purpose:          Public orchestration API for the project pipeline
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-13
# ==============================================================================

import logging

import anndata as ad
from typing import Dict

from methyltrain.project.layout import ProjectLayout
from methyltrain.audit.audit_store import AuditStore
from methyltrain.utils.logging import configure_logger
from methyltrain.io.write import save_manifest, save_metadata
from methyltrain.io.datasets import load_raw_project

from methyltrain.project.download import download_methylation
from methyltrain.project.metadata import prepare_metadata
from methyltrain.project.clean import clean_data
from methyltrain.project.qc import quality_control
from methyltrain.project.preprocess import preprocess


def prepare_project(config: Dict) -> ad.AnnData:
    """
    Upper-level orchestration API for the full project pipeline. Downloads and 
    preprocesses DNA Methylation data for the project specified in the 
    configurations.

    Parameters
    ----------
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    ad.AnnData
        The processed DNA Methylation dataset.
    """

    # Fetch aliases for commonly accessed configuration subgroups
    project_cfg = config.get('project', {})
    paths_cfg = config.get('paths', {})
    
    # Initialize the project's default output layout
    layout = ProjectLayout.from_config(config)
    layout.initialize()
    layout.validate()

    # Initialize the AuditStore for download status logging
    audit = AuditStore(layout.audit_store.with_suffix(".db")) 

    # Initialize a logger to capture verbose output
    logger = configure_logger(level = logging.INFO)

    logger.info("=====| MethylTrain: Project Processing Pipeline |=====\n")
    logger.info("~~~~~| Project Details |~~~~~")
    logger.info(f"Project ID: {project_cfg.get('id', '')}")
    logger.info(f"Seed: {project_cfg.get('seed', '')}")
    logger.info(f"GDC API Client Path: {paths_cfg.get('gdc_client', '')}")
    logger.info(f"Output Directory: {paths_cfg.get('output_dir', '')}")
    logger.info("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    logger.info("-----| Beginning Pipeline |-----\n")

    # ~~~~~| 1. Download and Clean Data |~~~~~

    # Download and clean the raw DNA Methylation data and metadata
    manifest, download_report = download_methylation(config, layout)
    metadata, metadata_report, biospec_report = prepare_metadata(config, audit)
    clean_data(layout, audit)

    # Load the raw data into a CpG matrix
    adata = load_raw_project(config, layout)

    # ~~~~~| 2. Quality Control |~~~~~
    adata, qc_report = quality_control(adata, config, layout)

    # ~~~~~| 3. Preprocessing |~~~~~
    adata = preprocess(adata, config)

    # ~~~~~| 4. I/O |~~~~~

    # Save all metadata to CSV files as specified in the layout
    save_manifest(manifest, layout)
    save_metadata(metadata, layout)

    # ~~~~~| 5. Auditing |~~~~~
    audit.initialize(manifest.index.tolist())
    audit.apply_download_report(download_report)
    audit.apply_metadata_report(metadata_report)
    audit.apply_biospecimen_report(biospec_report)
    audit.apply_qc_report(qc_report)

    # Export the AuditStore object to a CSV file


    # ~~~~~| 5. Cleanup |~~~~~

    # Optionally remove all raw data; processed data and metadata persist
    if config.get('clean_raw_data', False):
        for file in layout.raw_dir.glob("*.parquet"): file.unlink()
        logger.info("Successfully cleaned all raw data.")


    logger.info("=====================================================")

    return adata

# [END]