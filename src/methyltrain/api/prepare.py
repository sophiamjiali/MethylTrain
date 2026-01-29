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
    cohort_batch_correction,
    aggregate_genes,
    winsorize,
    split,
    load_raw_project,
    load_processed_project
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
    audit_table = download(config, layout)
    if verbose: print("Successfully downloaded data")

    if verbose: print("Attempting to clean the data")
    audit_table = clean_data(audit_table, layout)
    if verbose: print("Successfully cleaned the data")

    if verbose: print("Attempting to load the raw project data")
    adata = load_raw_project(config, layout)
    if verbose: print("Successfully loaded the raw project data")

    # Perform QC and preprocessing based on user configurations
    if verbose: print("Attempting to perform quality control")
    adata, audit_table = quality_control(adata, audit_table, config, layout)
    if verbose: print("Successfully performed quality control")

    if verbose: print("Attempting to preprocess the data")
    adata = preprocess(adata, config)
    if verbose: print("Successfully preprocessed the data")

    return (adata, audit_table)
    

def prepare_cohort(config: Dict, 
                   layout: CohortLayout) -> Tuple[ad.AnnData, ad.AnnData, 
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
    project_adatas = [load_processed_project(path) 
                      for path in layout.project_list]
    
    # Aggregate the projects into a cohort AnnData object
    cohort_adata = aggregate_cohort(project_adatas, layout)

    # Perform batch effect correction across datasets if toggled
    if config.get('toggles', {}).get('batch_correction', True):
        cohort_adata = cohort_batch_correction(cohort_adata, config)

    # Optionally aggregate to the gene-level if toggled by user configurations
    if config.get('toggles', {}).get('gene_aggregation', True):
        cohort_adata = aggregate_genes(cohort_adata, config)

    # Winsorize values of 0 or 1 to promote downstream ML stability
    if config.get('toggles', {}).get('winsorization', True):
        cohort_adata = winsorize(cohort_adata, config)

    # Split the cohort AnnData object into train-val-test splits
    train_adata, val_adata, test_adata = split(cohort_adata, config)

    return train_adata, val_adata, test_adata
