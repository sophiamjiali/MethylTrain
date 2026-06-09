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

    seed = config['cohort']['seed']

    # Convert obs_names to integer positions for splitting
    all_idx = np.arange(cohort.n_obs)

    # Fetch the splitting ratios and their relative ratio
    train_ratio, val_ratio, test_ratio = config['preprocessing']['split']
    if not np.isclose(train_ratio + val_ratio + test_ratio, 1.0):
        raise ValueError("Split ratios must sum to 1.0.")

    # Split into train+validation and test
    train_val_idx, test_idx = train_test_split(
        all_idx,
        test_size = test_ratio,
        stratify = cohort.obs['project_id'],
        random_state = seed,
        shuffle = True
    )

    # Parse train_val_idx to split again using the relative proportion
    val_relative = val_ratio / (train_ratio + val_ratio)
    stratify_array = cohort.obs['project_id'].to_numpy()[train_val_idx]

    # Split into train and validation
    train_idx, val_idx = train_test_split(
        train_val_idx,
        test_size = val_relative,
        stratify = stratify_array,
        random_state = seed,
        shuffle = True
    )

    # Slice the cohort AnnData object according to the defined splits
    train_adata = cohort[train_idx].copy()
    val_adata = cohort[val_idx].copy()
    test_adata = cohort[test_idx].copy()

    # Update each AnnData object's global metadata
    split = config['preprocessing']['split']
    train_adata = _update_metadata(train_adata, 'training', split[0])
    val_adata = _update_metadata(val_adata, 'validation', split[1])
    test_adata = _update_metadata(test_adata, 'testing', split[2])

    return train_adata, val_adata, test_adata


# =====| Internal Helpers |=====================================================

def _update_metadata(adata: ad.AnnData, split: str, percent: str) -> ad.AnnData:
    adata.uns['pipeline']['state'] = 'processed'
    adata.uns['pipeline']['steps'].append('split')
    adata.uns.setdefault('split', {})
    adata.uns['split']['type'] = split
    adata.uns['split']['percentage'] = percent
    return adata