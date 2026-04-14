# ==============================================================================
# Script:           prepare.py.py
# Purpose:          Exposes an end-to-end workflow wrapper function
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-07
# ==============================================================================

from typing import Dict, Tuple

import anndata as ad
import pandas as pd

from ..fs.layout import ProjectLayout, CohortLayout
from .steps import (
    download, 
    clean_data, 
    quality_control, 
    preprocess, 
    aggregate_cohort, 
    filter_by_mad,
    cohort_batch_correction,
    aggregate_genes,
    clip_and_scale,
    split,
    load_raw_project,
    load_processed_project,
    save_cohort
)

# =====| Exposed Functions |====================================================

def prepare_dataset(config: Dict, 
                    layout: ProjectLayout,
                    verbose = False) -> Tuple[ad.AnnData, pd.DataFrame]:
    """
    Run the full DNA methylation preprocessing workflow on a given project.

    Parameters
    ----------
    config : dict
        Configuration dictionary controlling workflow steps.
    layout : ProjectLayout
        Object representing a project dataset directory layout.

    Returns
    -------
    ad.AnnData
        The processed dataset.
    """

    # Ensure directories exist
    layout.initialize()

    # Download, clean, and load the project data as an AnnData object
    if verbose: print("Attempting to download data")
    audit_table = download(config, layout, verbose)
    if verbose: print("Successfully downloaded data")

    if verbose: print("Attempting to clean the data")
    audit_table = clean_data(audit_table, layout)
    if verbose: print("Successfully cleaned the data")

    if verbose: print("Attempting to load the raw project data")
    adata = load_raw_project(config, layout)
    if verbose: print("Successfully loaded the raw project data")

    # Perform QC and preprocessing based on user configurations
    if verbose: print("Attempting to perform quality control")
    adata, audit_table = quality_control(adata, audit_table, config, 
                                         layout, verbose)
    if verbose: print("Successfully performed quality control")

    if verbose: print("Attempting to preprocess the data")
    adata = preprocess(adata, config, verbose)
    if verbose: print("Successfully preprocessed the data")

    return (adata, audit_table)
    

def prepare_cohort(config: Dict, 
                   layout: CohortLayout,
                   verbose = False) -> Tuple[ad.AnnData, ad.AnnData, 
                                                  ad.AnnData]:
    """
    Aggregate the full DNA methylation preprocessing workflow outputs
    from multiple project(s) into a single cohort and split the dataset
    into training, validation, and testing splits.

    Parameters
    ----------
    config : dict
        Configuration dictionary controlling workflow steps.
    processed_paths : list of Path
        Paths to the processed output files (*.h5ad) for each project.
    cohort_dir : Path
        Directory for cohort data.
    """

    # Ensure directories exist
    layout.initialize()

    # Load each processed project AnnData object 
    if verbose: print("Attempting to load all project AnnData objects")
    project_adatas = [load_processed_project(path) 
                      for path in layout.project_list]
    if verbose: print("Successfully loaded all project AnnData objects")
    
    # Aggregate the projects into a cohort AnnData object
    if verbose: print("Attempting to aggregate the cohort")
    cohort_adata = aggregate_cohort(project_adatas, layout)
    if verbose: print("Successfully aggregated the cohort")

    # Perform MAD probe filtering to reduce probe dimensionality
    if config.get('toggles', {}).get('MAD_probe_filtering', True):
        if verbose: print("Attempting to perform MAD probe filtering")
        cohort_adata = filter_by_mad(cohort_adata, config)
        if verbose: print("Successfully performed MAD probe filtering")

    # Perform batch effect correction across datasets if toggled
    if config.get('toggles', {}).get('batch_correction', True):
        if verbose: print("Attempting to perform batch correction")
        cohort_adata = cohort_batch_correction(cohort_adata, config)
        if verbose: print("Successfully performed batch correction")

    # Optionally aggregate to the gene-level if toggled by user configurations
    if config.get('toggles', {}).get('gene_aggregation', True):
        if verbose: print("Attempting to perform gene aggregation")
        cohort_adata = aggregate_genes(cohort_adata, config)
        if verbose: print("Successfully aggregated to the gene-level")

    # Winsorize and scale M-values to [-1, 1] promote downstream ML stability
    if config.get('toggles', {}).get('clip_and_scale', True):
        if verbose: print("Attempting to perform winsorization and Min-max scaling")
        cohort_adata = clip_and_scale(cohort_adata, config)
        if verbose: print("Successfully winsorized and scaled the cohort")

    # Split the cohort AnnData object into train-val-test splits
    if config.get('toggles', {}).get('split', True):
        if verbose: print("Attempting to split the cohort into training sets")
        train_adata, val_adata, test_adata = split(cohort_adata, config)
        if verbose: print("Successfully split the cohort into training sets")
    else:
        if verbose: print("Did not split the cohort into training sets (not toggled)")
        train_adata, val_adata, test_adata = None, None, None

    # Save the processed cohort adata prior to splitting
    if verbose: print("Attempting to save the cohort AnnData object")
    save_cohort(cohort_adata, layout)
    if verbose: print("Successfully saved the cohort AnnData object")

    return train_adata, val_adata, test_adata

# [END]