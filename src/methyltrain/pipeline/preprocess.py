# ==============================================================================
# Script:           preprocess.py
# Purpose:          Internal preprocessing functions and implementation logic
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-16
# ==============================================================================

import anndata as ad
import numpy as np

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



def batch_correction():
    return
# batch effect correction using ComBat (AFTER COHORT AGGREGATION)


def winsorize_data(adata: ad.AnnData, config: Dict) -> ad.AnnData:
    """
    Winsorizes gene-level beta values 
    
    """

    return ad.AnnData()