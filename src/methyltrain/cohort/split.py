# ==============================================================================
# Script:           split.py
# Purpose:          Train-val-test split implementation
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-13
# ==============================================================================

import anndata as ad
import numpy as np

from typing import Dict
from sklearn.model_selection import train_test_split

# =====| Public API |===========================================================

def split(cohort: ad.AnnData, 
          config: Dict) -> tuple[ad.AnnData, ad.AnnData, ad.AnnData]:
    """
    Split a cohort AnnData object into stratified train, validation, and test 
    sets based on project using the ratio provided in the configurations.
    
    Parameters
    ----------
    cohort : ad.AnnData
        AnnData object representing a project's quality-controlled, preprocessed, 
        and batch effect corrected DNA methylation data at the gene matrix level.
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    tuple[ad.AnnData, ad.AnnData, ad.AnnData]
        The train, validation, and test stratified splits as ad.AnnData objects.
    """

    # revise below


    # Ensure obs_names are unique
    if not cohort.obs_names.is_unique:
        cohort.obs_names = cohort.obs['file_id']
        cohort.obs_names_make_unique()

    # Convert obs_names to integer positions for splitting
    all_idx = np.arange(cohort.n_obs)

    # Split into train+validation and test
    train_val_idx, test_idx = train_test_split(
        all_idx,
        test_size = config.get('split', [])[2],
        stratify = cohort.obs['project_id'],
        random_state = config.get('seed', 42),
        shuffle = True
    )

    # Parse train_val_idx to split again
    stratify_array = cohort.obs['project_id'].to_numpy()[train_val_idx]

    # Split into train and validation
    train_idx, val_idx = train_test_split(
        train_val_idx,
        test_size = config.get('split', [])[1],
        stratify = stratify_array,
        random_state = config.get('seed', 42),
        shuffle = True
    )

    # Slice the cohort AnnData object according to the defined splits
    train_adata = cohort[train_idx].copy()
    val_adata = cohort[val_idx].copy()
    test_adata = cohort[test_idx].copy()



    # revise below




    # Update each AnnData object's global metadata
    train_adata.uns['split'] = "training"
    val_adata.uns['split'] = "validation"
    test_adata.uns['split'] = "test"

    train_adata.uns['split_percentage'] = config.get('split', [])[0]
    val_adata.uns['split_percentage'] = config.get('split', [])[1]
    test_adata.uns['split_percentage'] = config.get('split', [])[2]

    return train_adata, val_adata, test_adata