# ==============================================================================
# Script:           aggregation.py
# Purpose:          Internal aggregation functions and implementation logic
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-17
# ==============================================================================

import copy

import anndata as ad
import pandas as pd
import numpy as np

from typing import List, Dict

from methyltrain.cohort.layout import CohortLayout
from methyltrain.constants.annotation import PLATFORM_PRIORITY

# =====| Public API |===========================================================

def aggregate_cohort(projects: List[ad.AnnData], 
                     config: Dict,
                     layout: CohortLayout) -> ad.AnnData:
    """
    Aggregates multiple project AnnData objects at the CpG probe x sample 
    matrix level into a single cohort AnnData object. Takes the common set of 
    CpG probes from all projects.

    Resolves dataset-level metadata such that `.uns` is a dictionary with 
    project names as keys, and project-level metadata as their values. 
    Cohort-level metadata is stored in a flat structure.

    Parameters
    ----------
    projects : List[ad.AnnData]
        List of project AnnData objects at the CpG probe x sample level, each representing a single project.
    config : dict
        Configuration dictionary controlling workflow steps.
    layout : CohortLayout
        Object representing a cohort dataset directory layout.

    Returns
    -------
    ad.AnnData
        Aggregated cohort AnnData object at the CpG probe x sample level.
    """
    
    # Concatenate all projects together, keeping the common set of probes
    cohort = ad.concat(
        projects,
        join = "inner",
        label = "project_id",
        keys = [p.uns['provenance']['project_id'] for p in projects]
    )

    # Assert that all projects are the same conversion (beta or M-values)
    conversion = _fetch_conversion(projects)

    # Set the array type for annotation as the highest resolution available
    platform = _fetch_highest_resolution(projects)

    # Fetch the reference genome of all the projects
    reference_genome = _fetch_reference_genome(projects)

    # Initialize the cohort metadata, initializing the projects nested dict.
    cohort.uns = {
        'provenance': {
            'cohort_id': layout.cohort_name,
            'level': 'cohort',
            'data_type': 'cpg_matrix',
            'aggregation_method': 'concatenation',
            'conversion': conversion,
            'pipeline_version': config['version'],
        },
        'data_source': {
            'platform': platform,
            'reference_genome': reference_genome,
        },
        'projects': {
            p.uns['provenance']['project_id']: copy.deepcopy(p.uns)
            for p in projects
        },
        'pipeline': {
            'state': 'processed',
            'steps': ['aggregation']
        }
    }


# =====| Internal Helpers |=====================================================

def _fetch_conversion(projects: List[ad.AnnData]) -> str:
    """
    Returns the conversion type of all projects. Raises an error if all projects do not have the same conversion type.

    Parameters
    ----------
    projects : List[ad.AnnData]
        List of project AnnData objects at the CpG probe x sample level, each representing a single project.

    Returns
    -------
    str
        The conversion type shared by all projects

    Raises
    -------
    ValueError
        If all projects do not have the same conversion type.
    """

    conversion = projects[0].uns['provenance']['conversion']
    for project in projects:
        if project.uns['provenance']['conversion'] != conversion:
            raise ValueError("All projects must be the same conversion type.")
        
    return conversion
        

def _fetch_highest_resolution(projects: List[ad.AnnData]) -> str:
    """
    Returns the highest resolution platform present in the projects such that 
    the annotation later selected will always be the highest resolution.

    Parameters
    ----------
    projects : List[ad.AnnData]
        List of project AnnData objects at the CpG probe x sample level, each representing a single project.

    Returns
    -------
    str
        The highest resolution platform present in the projects provided.
    """

    project_array_types = [p['data_source']['platform'] for p in projects]
    highest_array = PLATFORM_PRIORITY[-1]

    for array_type in PLATFORM_PRIORITY:
        if array_type in project_array_types: 
            highest_array = array_type
            break

    return highest_array


def _fetch_reference_genome(projects: List[ad.AnnData]) -> str:
    """
    Returns the reference genome of all projects. Raises an error if all 
    projects do not share the same reference genome.

    Parameters
    ----------
    projects : List[ad.AnnData]
        List of project AnnData objects at the CpG probe x sample level, each representing a single project.

    Returns
    -------
    str
        The reference genome type shared by all projects

    Raises
    -------
    ValueError
        If all projects do not have the same reference genome.
    """

    reference_genome = projects[0].uns['data_source']['reference_genome']
    for project in projects:
        if project.uns['data_source']['reference_genome'] != reference_genome:
            raise ValueError("All projects must be the same reference genome.")
        
    return reference_genome
    
# [END]