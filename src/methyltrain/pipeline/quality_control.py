# ==============================================================================
# Script:           quality_control.py
# Purpose:          Internal quality-control functions and implementation logic
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-08
# ==============================================================================

import numpy as np
import anndata as ad
import pandas as pd

from typing import Dict

from ..utils.utils import iqr_bounds

def sample_qc(adata: ad.AnnData, config: Dict) -> ad.AnnData:
    """
    Performs sample-level quality control upon DNA methylation values presented as a CpG matrix AnnData object. Returns the quality-controlled CpG matrix AnnData object with updated metadata suitable for further probe-level quality control and/or preprocessing.

    Steps performed are toggled and configured in the user-provided 
    configurations, including the following options in order:

    1. Remove high missingness (above the provided threshold)
    2. Remove outliers from distribution (IQR-based)
    
    Parameters
    ----------
    adata : ad.AnnData
        AnnData object representing a project's raw DNA methylation data at the 
        CpG matrix level.
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    ad.AnnData
        AnnData object representing a project's sample-level quality-controlled 
        DNA methylation data at the CpG matrix level with updated metadata.
    """

    # Fetch relevant thresholds
    qc_cfg = config.get('quality_control', {}).get('sample_qc', {})
    missing_threshold = qc_cfg.get('missing_threshold', 0.20)
    outlier_threshold = qc_cfg.get('outlier_threshold', 1.5)

    # Filter samples for missingness above the provided threshold
    cpg_matrix = np.asarray(adata.X)
    missing_rate = np.isnan(cpg_matrix).mean(axis = 1)
    fail_missing = missing_rate > missing_threshold

    # Filter samples for global signal intensity using a mean beta filter
    sample_mean = np.nanmean(cpg_matrix, axis = 1)
    mean_lo, mean_hi = iqr_bounds(sample_mean, outlier_threshold)
    fail_mean = (sample_mean < mean_lo) | (sample_mean > mean_hi)

    # Generate failure mask using missingness and IQR
    fail_any = fail_missing | fail_mean

    # Set global metadata
    adata.obs['missing_rate'] = missing_rate
    adata.obs['mean_beta'] = sample_mean

    adata._inplace_subset_obs(~fail_any)

    adata.uns['preprocessing_steps'].append("sample_qc")

    return adata


def probe_qc(adata: ad.AnnData, 
             annotation: pd.DataFrame, 
             config: Dict) -> ad.AnnData:
    """
    Performs probe-level quality control upon DNA methylation values presented as a CpG matrix AnnData object. Returns the quality-controlled CpG matrix AnnData object with updated metadata suitable for further preprocessing.

    Steps performed are toggled and configured in the user-provided 
    configurations, including the following options in order:

    1. Remove structural probes
    2. Remove cross-reactive probes
    3. Remove SNP-associated probes
    4. Remove multi-mapped probes
    5. Remove sex chromosome probes
    6. Remove high missingness (above the provided threshold)
    
    Parameters
    ----------
    adata : ad.AnnData
        AnnData object representing a project's sample-level quality-controlled 
        DNA methylation data (or raw if QC was not performed) at the CpG matrix 
        level.
    annotation : pd.DataFrame
        Standardized Illumina annotations provided by the package for the 
        project's array type and genome build.
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    ad.AnnData
        AnnData object representing a project's probe-level quality-controlled 
        DNA methylation data at the CpG matrix level with updated metadata.
    """

    # Fetch relevant thresholds
    probe_cfg = config.get('quality_control', {}).get('probe_qc', {})
    missing_threshold = probe_cfg.get('missing_threshold', 0.05)

    # Remove non-standard probes from the CpG matrix and annotation
    keep_standard = adata.var_names.str.startswith('cg')
    adata._inplace_subset_var(keep_standard)

    # Filter by missingness
    missing_rate = np.isnan(np.asarray(adata.X)).mean(axis = 0)
    keep_missing = missing_rate <= missing_threshold

    # Filter by annotation flags provided in the annotation object
    annotation = annotation.set_index('probe_id').reindex(adata.var_names)
    keep_annotation = pd.Series(True, index = annotation.index)

    if probe_cfg.get('remove_sex_chromosome', True):
        keep_annotation &= ~annotation['is_sex_chr']
    if probe_cfg.get('remove_SNP_associated', True):
        keep_annotation &= ~annotation['has_cpg_snp']
        keep_annotation &= ~annotation['has_sbe_snp']
        keep_annotation &= ~annotation['has_probe_snp']
    if probe_cfg.get('remove_cross_reactive', True):
        keep_annotation &= ~annotation['is_cross_reactive']
    if probe_cfg.get('remove_multi_mapped', True):
        keep_annotation &= ~annotation['is_multi_mapped']

    # Combine filters and subset the matrix
    pass_qc = keep_missing & keep_annotation

    # Set global metadata
    adata.var['missing_rate'] = missing_rate
    adata.uns['preprocessing_steps'].append("probe_qc")

    adata._inplace_subset_var(pass_qc)

    return adata