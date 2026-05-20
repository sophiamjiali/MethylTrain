# ==============================================================================
# Script:           correction.py
# Purpose:          Batch correction implementation
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-13
# ==============================================================================

import logging

import anndata as ad
import numpy as np

from typing import Dict, List
from inmoose.pycombat import pycombat_norm

logger = logging.getLogger(__name__)

# =====| Internal Helpers |=====================================================

def align_distribution(cohort: ad.AnnData, config: Dict) -> ad.AnnData:
    """
    Aligns the distributions of the cohort across projects.
    
    Parameters
    ----------
    adata : ad.AnnData
        AnnData object with DNA methylation M-values in .X.
    config : dict
        Configuration dictionary controlling workflow steps.
    
    Returns
    -------
    adata : ad.AnnData
        AnnData object representing a project's batch-corrected DNA methylation 
        data at the CpG probe matrix level with updated metadata.
    """

    logger.info("=====| Attempting to Align Distributions |=====")

    # Perform batch effect correction across datasets
    if config.get('toggles', {}).get('batch_correction', True):

        # Fetch the batch key and covariates from the configurations
        batch_cfg = config.get('preprocessing', {}).get('batch_correction', {})
        batch_key = batch_cfg.get('batch_key', 'batch_id')
        covariates = batch_cfg.get('covariates', [])

        cohort = _batch_correction(cohort, batch_key, covariates)

    logger.info("=====| Successfully Aligned Distributions |=====")
    return cohort


# =====| Public API |===========================================================

def _batch_correction(adata: ad.AnnData, 
                      batch_key: str, 
                      covariates: List[str]) -> ad.AnnData:
    """
    Performs ComBat batch correction on an AnnData object using InMoose's 
    pycombat_norm function. Automatically checks for singleton batches and 
    collapses to portion- or plate-level if necessary.

    Batch correction should be applied after cohort aggregation.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object with DNA methylation M-values in .X.
    batch_key : str
        The batch key to use for PyCombat batch correction.
    covariates : List[str]
        The covariates to use for PyCombat batch correction.
    
    Returns
    -------
    adata : ad.AnnData
        AnnData object representing a project's batch-corrected DNA methylation 
        data at the CpG probe matrix level with updated metadata.

    Raises
    ------
    ValueError
        If NaNs exist in the AnnData object (batch correction requires none), 
        if data is presented as beta values instead of M-values, or if batch 
        and covariate columns don't exist.
    """

    X = np.array(adata.X, dtype = np.float32)

    # Ensure no NaNs remain in the AnnData object
    if np.isnan(X).any():
        raise ValueError("Batch correction requires no Nans. Please impute "
                         "or remove missing values first.")
    
    # Ensure data type is M-values, not beta values
    if adata.uns['provenance']['conversion'] != "m_value":
        raise ValueError("Batch correction should be performed on M-values. "
                         "Please convert from beta values to M-values first.")

    for col in [batch_key] + covariates:
        if col not in adata.obs.columns: 
            raise ValueError(f"Column {col} not found in adata.obs. Required "
                             "for batch correction.")
        
    # Prepare the batch vector and covariate matrix
    batch = adata.obs[batch_key].astype(str)
    mod = adata.obs[covariates] if covariates else None

    adata.X = pycombat_norm(X.T, batch = batch.values, covar_mod = mod).T

    adata.uns['pipeline']['state'] = 'processed'
    adata.uns['pipeline']['steps'].append("batch_correction")
    adata.uns.setdefault('alignment', {})
    adata.uns['alignment']['batch_correction'] = {
        'method': 'pycombat_norm',
        'batch_key': batch_key,
        'covariates': covariates
    }

    return adata

# [END]