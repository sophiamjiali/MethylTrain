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
from methyltrain.utils.logging import configure_logger

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

    # Initialize the project's default output layout
    layout = ProjectLayout.from_config(config)
    layout.validate()

    # Initialize a logger to capture verbose output
    logger = configure_logger(level = logging.INFO)


    logger.info(f"Preparing Project: {config.get('project', {}).get('id', '')}")



    # handle all I/O here, including initializing the auditstore, saving 
    # metadata into csv, etc.