# ==============================================================================
# Script:           aggregate.py
# Purpose:          Internal aggregation functions and implementation logic
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-17
# ==============================================================================

import anndata as ad
import pandas as pd
import numpy as np

from typing import List, Dict

from ..constants import PLATFORM_PRIORITY


def cohort_aggregation(adatas: List[ad.AnnData]) -> ad.AnnData:
    """
    Aggregates multiple project AnnData objects at the CpG probe x sample 
    matrix level into a single cohort AnnData object. Takes the common set of 
    CpG probes from all projects.

    Resolves dataset-level metadata such that `.uns` is a dictionary with 
    project names as keys, and project-level metadata as their values. 
    Cohort-level metadata is stored in a flat structure.

    Parameters
    ----------
    adatas : List of ad.AnnData
        List of project AnnData objects at the CpG probe x sample level, each representing a single project. Each object must have:
        - .uns['project_id'] : unique project name
        - .uns['level'] : 'project' or 'cohort'
        - .uns['data_type'] : 'cpg_matrix' or 'gene_matrix'
        - .uns['state'] : 'raw' or 'processed'
        - .uns['preprocessing_steps'] : QC and processing steps performed

    Returns
    -------
    ad.AnnData
        Aggregated cohort AnnData object at the CpG probe x sample level.
    """

    # Collect per-project .uns metadata
    projects_uns = {}
    for adata in adatas:

        project_name = adata.uns.get('project_id')
        project_metadata = {k: v for k, v in adata.uns.items() 
                            if k != 'project_id'}
        
        projects_uns[project_name] = project_metadata

    # Concatenate projects, keeping the common set of CpG probes/genes
    cohort_adata = ad.concat(
        adatas,
        join = "inner",
        label = "project_id",
        keys = [adata.uns['project_id'] for adata in adatas]
    )

    # Set the array type for annotation use as the highest resolution available
    project_arrays = {
        name: meta.get('platform', None) 
        for name, meta in projects_uns.items()
    }

    highest_array = PLATFORM_PRIORITY[-1]
    for array_type in PLATFORM_PRIORITY:
        if array_type in project_arrays.values():
            highest_array = array_type
            break

    cohort_adata.uns['platform'] = highest_array
    cohort_adata.uns['reference_genome'] = adatas[0].uns['reference_genome']
    cohort_adata.uns['projects'] = projects_uns
    cohort_adata.uns['data_type'] = adatas[0].uns['data_type']
    cohort_adata.uns['conversion'] = adatas[0].uns['conversion']

    return cohort_adata


def gene_aggregation(adata: ad.AnnData, 
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
    regions = config.get('gene_aggregation', {}).get('regions', [])

    # Align annotations to the AnnData probes
    annotation = annotation.set_index("probe_id").loc[adata.var_names]
    
    # Flag probes if they are in the contributing regions list
    annotation['keep'] = annotation[regions].any(axis = 1)
    annotation = annotation[annotation['keep']]

    cpg_matrix = pd.DataFrame(
        np.array(adata.X), 
        index = adata.obs_names, 
        columns = adata.var_names
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
    adata.X = gene_matrix.values
    adata.var = pd.DataFrame(index = gene_matrix.columns)
    adata.uns['data_type'] = 'gene_matrix'
    adata.uns['preprocessing_steps'].append("gene_aggregation")

    return adata