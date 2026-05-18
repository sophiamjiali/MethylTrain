# ==============================================================================
# Script:           aggregation.py
# Purpose:          Internal aggregation functions and implementation logic
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-17
# ==============================================================================

import copy
import logging

import anndata as ad
import numpy as np
import pandas as pd

from typing import List, Dict

from methyltrain.cohort.layout import CohortLayout
from methyltrain.constants.annotation import PLATFORM_PRIORITY

logger = logging.getLogger(__name__)

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

    layout.validate()

    logger.info("=====| Attempting Cohort Aggregation |=====")
    
    # Concatenate all projects together, keeping the common set of probes
    cohort = ad.concat(
        projects,
        join = "inner",
        label = "project_id",
        keys = [p.uns['provenance']['project_id'] for p in projects]
    )

    # Assert that all projects are the same conversion (beta or M-values)
    conversion = _fetch_conversion(projects)
    logger.info("Successfully fetched project conversion types.")

    # Set the array type for annotation as the highest resolution available
    platform = _fetch_highest_resolution(projects)
    logger.info("Successfully fetched the highest array type resolution.")

    # Fetch the reference genome of all the projects
    reference_genome = _fetch_reference_genome(projects)
    logger.info("Successfully fetched all projects' reference genome.")

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

    logger.info("Successfully initialized the cohort metadata.")
    logger.info("=====| Successfully Aggregated the Cohort |=====")

    return cohort


def aggregate_genes(cohort: ad.AnnData, 
                    annotation: pd.DataFrame,
                    config: Dict) -> ad.AnnData:
    """
    Aggregates a project or cohort AnnData object at the CpG probe x sample 
    matrix level to the gene-level. Regions are used to select probes that 
    contribute, but the final matrix has only genes as columns.

    Parameters
    ----------
    adata : ad.AnnData
        Probe-level DNA methylation AnnData (probes x samples).
    annotation : pd.DataFrame
        Simplified annotation table with the following columns:
        - 'probe_id'
        - 'gene_symbol'
        - 'TSS200', 'TSS1500', 'gene_body' (bool)
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    adata : ad.AnnData
        Aggregated AnnData object to the gene-level.
    """

    # Fetch the regions to aggregate (TSS200, TS1500, gene body)
    regions = config.get('preprocessing', {}).get('gene_aggregation', [])

    # Align annotations to the AnnData probes
    annotation = annotation.set_index("probe_id").loc[cohort.var_names]
    
    # Flag probes if they are in the contributing regions list
    annotation['keep'] = annotation[regions].any(axis = 1)
    annotation = annotation[annotation['keep']]

    cpg_matrix = pd.DataFrame(
        np.array(cohort.X), 
        index = cohort.obs_names, 
        columns = cohort.var_names
    )

    # Explode multi-gene probes
    annotation['gene_symbol'] = annotation['gene_symbol'].str.split(';')
    annotation = (annotation.explode('gene_symbol')
                  .dropna(subset = ['gene_symbol']))
    
    cpg_matrix = cpg_matrix[annotation.index]
    cpg_matrix.columns = annotation['gene_symbol'].values
    
    # Aggregate by gene: take the mean across all probes per gene
    gene_matrix = cpg_matrix.T.groupby(level = 0).mean().T
    

    # Update the AnnData object in-place
    cohort.X = gene_matrix.values
    cohort.var = pd.DataFrame(index = gene_matrix.columns)

    cohort.uns['pipeline']['state'] = 'processed'
    cohort.uns['pipeline']['steps'].append("gene_aggregation")
    cohort.uns['provenance']['data_type'] = 'gene_matrix'
    cohort.uns['gene_aggregation'] = {
        'regions': regions
    }

    return cohort


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