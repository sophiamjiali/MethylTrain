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

def normalize(adata: ad.AnnData) -> ad.AnnData:
    """
    Performs minimal beta-scale normalization and harmonization on DNA 
    methylation values in a CpG matrix AnnData object. This aligns the 
    distributions of beta-values across samples without altering biological 
    signal.

    Normalization should be performed after sample- and probe-level quality 
    control and before filtering low-variance probes.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object representing a project's DNA methylation data at the 
        CpG matrix level.

    Returns
    -------
    ad.AnnData
        AnnData object representing a project's normalized DNA methylation data 
        at the CpG matrix level with updated metadata.
    """

    X = np.array(adata.X, copy = True).astype(np.float32)

    # Compute per-sample median and IQR
    sample_median = np.nanmedian(X, axis = 1)
    sample_iqr = (
        np.nanpercentile(X, 75, axis = 1) -
        np.nanpercentile(X, 25, axis = 1)
    )
    
    # Guard against invalid scaling
    invalid = np.isnan(sample_median) | np.isnan(sample_iqr) | (sample_iqr == 0)
    sample_iqr[invalid] = 1.0
    sample_median[invalid] = 0.0

    # Compute global median and IQR across samples
    global_median = np.nanmedian(X)
    global_iqr = (
        np.nanpercentile(X, 75) -
        np.nanpercentile(X, 25)
    )
    if np.isnan(global_iqr) or global_iqr == 0: global_iqr = 1.0

    # Scale each sample to match the global median and IQR
    X = ((X.T - sample_median) / sample_iqr * global_iqr + global_median).T

    # Clip to [0, 1] to retain valid beta-value bounds
    X = np.clip(X, 0, 1)

    adata.X = X

    # Store metadata
    adata.uns['preprocessing_steps'].append('normalize')

    return adata


def filter_variance(adata: ad.AnnData, config: Dict) -> ad.AnnData:
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
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    ad.AnnData
        AnnData object with low-variance or extreme probes removed and updated 
        metadata.
    """

    # Fetch user-provided configurations with defaults
    filter_cfg = config.get('preprocessing', {}).get('filter_variance', {})

    min_variance = filter_cfg.get('min_variance', 0.0001)
    min_mean = filter_cfg.get('min_mean', 0.01)
    max_mean = filter_cfg.get('max_mean', 0.99)
    
    X = np.array(adata.X, copy = True).astype(np.float32)

    # Compute per-probe variance and mean
    probe_variance = np.var(X, axis = 0)
    probe_mean = np.mean(X, axis = 0)

    # Identify probes to keep and filter
    keep_mask = ((probe_variance >= min_variance) & 
                 (probe_mean >= min_mean) & 
                 (probe_mean <= max_mean))
    
    adata._inplace_subset_var(keep_mask)

    # Store metadata; mean already stored in probe QC step
    adata.var['variance'] = probe_variance[keep_mask]
    adata.uns['preprocessing_steps'].append("filter_variance")

    return adata


def impute(adata: ad.AnnData):
    """
    Imputes missing beta values in a CpG matrix AnnData object per probe using 
    the mean value across all samples. Returns the imputed object with updated 
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


def convert_to_mval(adata: ad.AnnData, epsilon: float = 1e-6) -> ad.AnnData:
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
    X = np.array(adata.X, copy = True).astype(np.float32)
    X = np.clip(X, epsilon, 1 - epsilon)
    adata.X = np.log2(X / (1 - X))

    # Store metadata
    adata.uns['preprocessing_steps'].append('convert_to_mval')
    adata.uns['conversion'] = 'm_value'

    return adata


def batch_correction(adata: ad.AnnData, config: Dict) -> ad.AnnData:
    """
    Performs ComBat batch correction on an AnnData object using InMoose's 
    pycombat_norm function. Automatically checks for singleton batches and 
    collapses to portion- or plate-level if necessary.

    Batch correction should be applied after cohort aggregation.

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

    Raises
    ------
    ValueError
        If NaNs exist in the AnnData object (batch correction requires none), 
        if data is presented as beta values instead of M-values, or if batch 
        and covariate columns don't exist.
    """

    X = np.array(adata.X)

    # Ensure no NaNs remain in the AnnData object
    if np.isnan(X).any():
        raise ValueError("Batch correction requires no Nans. Please impute "
                         "or remove missing values first.")
    
    # Ensure data type is M-values, not beta values
    if adata.uns['conversion'] != "m_value":
        raise ValueError("Batch correction should be performed on M-values. "
                         "Please convert from beta values to M-values first.")
    
    # Ensure batch and covariate columns exist
    batch_cfg = config.get('batch_correction', {})
    batch_key = batch_cfg.get('batch_key', 'batch_id')
    covariates = batch_cfg.get('covariates', ['project_id'])

    for col in [batch_key] + covariates:
        if col not in adata.obs.columns: 
            raise ValueError(f"Column {col} not found in adata.obs. Required "
                             "for batch correction.")
        
    # Prepare the batch vector and covariate matrix
    batch = adata.obs[batch_key].astype(str).copy()
    mod = adata.obs[covariates].copy()

    # Detect singleton batches
    batch_counts = batch.value_counts()
    singletons = batch_counts[batch_counts == 1].index.tolist()

    if singletons:

        # Use the portion as the batch key
        batch_simple = batch.str.split('-').str[0]
        batch_counts_simple = batch_simple.value_counts()
        singletons_simple = batch_counts_simple[batch_counts_simple == 1].index.tolist()

        # Use plate-level as a backup if collapsing to portion fails
        if singletons_simple:
            batch_simple = batch.str.split('-').str[:2].str.join('-')
        
        batch = batch_simple

    # Apply ComBat via the InMoose package
    adata.X = pycombat_norm(counts = X.T, batch = batch.values, covar_mod = mod).T

    adata.uns['preprocessing_steps'].append("batch_correction")

    return adata
