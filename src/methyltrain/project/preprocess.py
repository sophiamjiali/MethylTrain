# ==============================================================================
# Script:           preprocess.py
# Purpose:          Internal preprocessing functions and implementation logic
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-16
# ==============================================================================

import logging

import anndata as ad
import numpy as np

from typing import Dict

logger = logging.getLogger(__name__)

# =====| Public API |===========================================================

def preprocess(adata: ad.AnnData, config: Dict) -> ad.AnnData:
    """
    Preprocess DNA methylation beta values of a project to a gene-level matrix. 
    Returns a samples x genes AnnData object with aligned metadata suitable for 
    downstream analysis.

    Note that further normalization upon level 3 beta values fetched from the 
    GDC is not performed as they are already SeSAMe-harmonized.

    Steps performed are toggled and configured in the user-provided 
    configurations, including the following options in order:
    
    1. Filter low variance / extreme values
    2. Imputation of missing values
    3. M-value conversion 

    Note that batch effect correction is optionally performed after multiple 
    projects have been aggregated into a cohort.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object representing a project's quality-controlled (or raw if 
        QC was not performed) DNA methylation data at the CpG matrix level.
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    ad.AnnData
        AnnData object representing a project's preprocessed DNA methylation 
        data at the CpG probe matrix level with updated metadata.
    """

    logger.info("=====| Attempting to Perform Preprocessing |=====")

    # Fetch all toggles and perform each if indictated
    apply_filter = config.get('toggles', {}).get('filter_variance', True)
    apply_impute = config.get('toggles', {}).get('impute', True)
    apply_mval = config.get('toggles', {}).get('convert_to_mval', True)

    if apply_filter or apply_impute or apply_mval:
        adata.uns['preprocess'] = {}

    # Apply a variance filter if toggled
    if apply_filter: 

        # Extract relevant configurations
        filter_cfg = config.get('preprocessing', {}).get('filter_variance', {})
        min_var = filter_cfg.get('min_variance', 0.0001)
        min_mean = filter_cfg.get('min_mean', 0.01)
        max_mean = filter_cfg.get('max_mean', 0.99)

        adata = _filter_variance(adata, min_var, min_mean, max_mean)
        logger.info("Successfully applied a variance filter.")

    # Impute missing values if toggled
    if apply_impute:
        adata = _impute(adata)
        logger.info("Successfully applied a variance filter.")

    # Convert beta values to M-Values if toggled
    if apply_mval:

        # Extract relevant configurations
        preproc_cfg = config.get('preprocessing', {})
        epsilon = preproc_cfg.get('mval_conversion', {}).get('epsilon', 0.001)

        adata = _convert_to_mval(adata, epsilon)
        logger.info("Successfully converted beta values to M-Values.")

    if apply_filter or apply_impute or apply_mval:
        adata.uns['pipeline']['state'] = 'processed'
        adata.uns['pipeline']['steps'].append('preprocess')

    logger.info("=====| Successfully Performed Preprocessing |=====")

    return adata

# =====| Internal Helpers |=====================================================

def _filter_variance(adata: ad.AnnData,
                     min_variance: float = 0.0001,
                     min_mean: float = 0.01,
                     max_mean: float = 0.99) -> ad.AnnData:
    """
    Filters CpG probes in an AnnData object based on low variance or extreme 
    mean beta values. Probes with near-constant values across samples are 
    removed, as they carry minimal biological information.

    This step should be applied after normalization and before M-value 
    conversion.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object representing a project's DNA methylation data at the 
        CpG matrix level (beta-values).
    min_variance : float
        Minimum variance required to pass the filter.
    min_mean : float
        Minimum mean required to pass the filter
    max_mean : float
        Maximum mean required to pass the filter.

    Returns
    -------
    ad.AnnData
        AnnData object with low-variance or extreme probes removed and updated 
        metadata.
    """
    X = np.array(adata.X, copy = True).astype(np.float32)

    # Compute per-probe variance and mean
    probe_variance = np.var(X, axis = 0)
    probe_mean = np.mean(X, axis = 0)

    # Identify probes to keep and filter
    keep_mask = ((probe_variance >= min_variance) & 
                 (probe_mean >= min_mean) & 
                 (probe_mean <= max_mean))
    
    # Store metadata; mean already stored in probe QC step
    adata.var['variance'] = probe_variance
    adata.var['mean'] = probe_mean

    adata = adata[keep_mask].copy()

    adata.uns['preprocess']['filter_variance'] = {
        'n_total': len(keep_mask),
        'n_pass': int(keep_mask.sum()),
        'n_fail': int(~keep_mask.sum()),
        'params': {
            'min_variance': min_variance,
            'min_mean': min_mean,
            'max_mean': max_mean
        }
    }

    return adata


def _impute(adata: ad.AnnData) -> ad.AnnData:
    """
    Imputes missing beta values in a CpG matrix AnnData object per probe using 
    the mean value across all samples. Returns the imputed object with updated 
    metadata suitable for further preprocessing.

    Imputation is intended to be performed after sample- and probe-level 
    quality control such that the majority of missing values, identified as 
    technical artifacts, are already removed.
    
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
    X = np.asarray(adata.X)

    missing_rate = np.isnan(X).mean(axis = 0)
    col_mean = np.nanmean(X, axis = 0)

    n_imputed = np.isnan(X).sum()
    frac_imputed = n_imputed / X.size

    # Impute missing values using column means
    mask = np.isnan(X)
    X[mask] = col_mean[np.where(mask)[1]]
    adata.X = X

    # Store metadata
    adata.var['frac_imputed'] = missing_rate
    adata.var['impute_value'] = col_mean

    adata.uns['preprocess']['impute'] = {
        'method': 'mean_probe',
        'total_imputed': n_imputed,
        'frac_imputed': frac_imputed
    }
    
    return adata


def _convert_to_mval(adata: ad.AnnData, epsilon: float = 1e-3) -> ad.AnnData:
    """
    Converts beta values in a CpG matrix AnnData object to M-values using the 
    logit-transformation: M = log2(beta / (1 - beta)). This transformation 
    produces unbounded, approximately homoscedastic values suitable for 
    downstream modeling and batch correction.

    Conversion should be performed after normalization and before batch 
    correction.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object representing a project's DNA methylation data at the 
        CpG matrix level.

    Returns
    -------
    ad.AnnData
        AnnData object with M-values in .X and updated preprocessing metadata.
    """

    # Clip beta values to avoid division by zero or log(0)
    beta = np.array(adata.X, copy = True).astype(np.float32)
    beta_min, beta_max = np.nanmin(beta), np.nanmax(beta)

    beta_clipped = np.clip(beta, epsilon, 1 - epsilon)
    mval = np.log2(beta_clipped / (1 - beta_clipped))

    m_min, m_max = np.nanmin(mval), np.nanmax(mval)

    adata.X = mval.copy()

    # Store metadata
    adata.uns['provenance']['conversion'] = 'm_value'
    adata.uns['preprocess']['mval_conversion'] = {
        'method': 'log2it',
        'formula': "log2(beta / (1 - beta))",
        "input_scale": "beta",
        "output_scale": "m_value",
        'epsilon': epsilon,
        'range': {
            'input': {
                "min": float(beta_min),
                "max": float(beta_max),
            },
            "output": {
                "min": float(m_min),
                "max": float(m_max),
            },
        }
    }

    return adata

# [END]