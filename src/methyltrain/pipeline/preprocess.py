# ==============================================================================
# Script:           preprocess.py
# Purpose:          Internal preprocessing functions and implementation logic
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-16
# ==============================================================================

import anndata as ad
import numpy as np

from inmoose.pycombat import pycombat_norm
from typing import Dict

def impute(adata: ad.AnnData):
    """
    Imputes missing beta values in a CpG matrix AnnData object per probe using 
    the mean value across all samples. Returns the imputed object with udpated 
    metadata suitable for further preprocessing.

    Imputation is intended to be performed after sample- and probe-level 
    quality control such that the majority of missing values, identified as 
    Okaytechnical artifacts, are already removed.
    
    Parameters
    ----------
    adata : ad.AnnData
        AnnData object representing a project's DNA methylation data at the 
        CpG matrix level.

    Returns
    -------
    ad.AnnData
        AnnData object representing a project's sample-level imputed 
        DNA methylation data at the CpG matrix level with updated metadata.
    """

    # Compute per-probe mean, ignoring NaNs
    X = np.array(adata.X, copy = True)

    missing_rate = np.isnan(X).mean(axis = 0)
    col_mean = np.nanmean(X, axis = 0)

    # Impute missing values using column means
    mask = np.isnan(X)
    X[mask] = np.take(col_mean, np.where(mask)[1])
    adata.X = X

    # Store metadata
    adata.var['frac_imputed'] = missing_rate
    adata.var['impute_value'] = col_mean

    adata.uns['preprocessing_steps'].append('impute')

    return adata



def batch_correction(adata: ad.AnnData) -> ad.AnnData:
    """
    Performs batch correction on an AnnData object using the InMoose ComBat 
    function (pycombat_norm).

    Correction is applied using the `project_id` column in .obs. The corrected matrix is stored in in-place for simplicity (overwrites .X).

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object representing a project or cohort's optionally quality 
        controlled and/pr preprocessed DNA methylation data at the CpG matrix 
        level. The object must include:
        - .obs['project_id] : project of origin for each sample
        - .X : probe-level beta values (dense, not sparse)
    
    Returns
    -------
    adata : ad.AnnData
        AnnData object representing a project's batch-corrected DNA methylation data at the CpG probe matrix level with updated metadata.

    Raises
    ------
    ValueError
        If NaNs exist in the AnnData object (batch correction requires none).
    """

    X = np.array(adata.X)

    # Ensure no NaNs remain in the AnnData object
    if np.isnan(X).any():
        raise ValueError("Batch correction requires no Nans. Please impute "
                         "or remove missing values first.")
    
    # Apply ComBat from the inmoose package
    batch = adata.obs['project_id'].values
    X = pycombat_norm(X, batch)
    adata.X = X

    adata.uns['preprocessing_steps'].append("batch_correction")

    return adata



def winsorize_data(adata: ad.AnnData, config: Dict) -> ad.AnnData:
    """
    Winsorizes gene-level beta values 
    
    """

    return ad.AnnData()