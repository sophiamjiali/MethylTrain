# ==============================================================================
# Script:           qc.py
# Purpose:          Quality control protocol for DNA methylation data
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-20
# ==============================================================================

import logging

import anndata as ad
import pandas as pd
import numpy as np

from methyltrain.project.layout import ProjectLayout
from methyltrain.project.audit_store import AuditStore
from methyltrain.io.read import load_annotation

logger = logging.getLogger(__name__)

# =====| Public API |===========================================================

def quality_control(adata: ad.AnnData,
                    config: dict,
                    layout: ProjectLayout) -> tuple[ad.AnnData, list[dict]]:
    """
    Performes probe and/or sample quality control upon DNA methylation values 
    presented as a CpG matrix AnnData object. Returns the quality-controlled 
    CpG matrix AnnData object with updated metadata suitable for downstream 
    preprocessing and analysis, and the updated audit table with QC status.

    Note that quality control is intended to be performed upon raw project data 
    (as loaded by `load_raw_project()`) and followed by preprocessing (as per 
    `preprocess()`).

    Steps performed are toggled and configured in the user-provided 
    configurations, including the following options in order:

    1. Sample-level quality control
        a. Remove high missingness (above the provided threshold)
        b. Remove outliers from distribution (above the provided number of SD)
    2. Probe-level quality control 
        a. Remove cross-reactive probes
        b. Remove SNP-associated probes
        c. Remove multi-mapped probes
        d. Remove sex chromosome probes
        e. Remove high missingness (above the provided threshold)

    Generates processed metadata and manifest files for the project after 
    quality control samples may be filtered from the dataset. If no quality 
    control was configured to occur in the user-configurations but this 
    function was still called, the processed files will include identical 
    information to the raw files.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object representing a project's raw DNA methylation data at the 
        CpG matrix level.
    config : dict
        Configuration dictionary controlling workflow steps.
    layout : ProjectLayout
        Object representing a project dataset directory layout.

    Returns
    -------
    tuple[ad.AnnData, list[dict]]
        AnnData object representing a project's quality-controlled DNA 
        methylation data at the CpG matrix level with updated metadata. A 
        quality control report for auditing.
    """

    layout.validate()

    logger.info("~~~~~| Attempting to Perform Quality Control |~~~~~")
    qc_report = [{}]

    # Fetch all toggles and perform each if indicated
    apply_sample_qc = config.get('toggles', {}).get('sample_qc', True)
    apply_probe_qc = config.get('toggles', {}).get('probe_qc', True)

    if apply_sample_qc or apply_probe_qc: adata.uns['qc'] = {}

    # Perform sample-level quality control if toggled
    if apply_sample_qc:

        # Fetch relevant thresholds
        qc_cfg = config.get('quality_control', {}).get('sample_qc', {})
        missing_thresh = qc_cfg.get('missing_threshold', 0.20)
        outlier_thresh = qc_cfg.get('outlier_threshold', 1.5)

        adata, qc_report = _sample_qc(adata, missing_thresh, outlier_thresh)
        logger.info("Successfully performed sample quality control.")

    # Perform probe-level quality control if toggled
    if apply_probe_qc:

        # Load the appropriate annotation provided by the package
        platform = adata.uns['data_source']['platform']
        reference_genome = adata.uns['data_source']['reference_genome']

        annotation = load_annotation(platform, reference_genome)
        logger.info(f"Loaded annotation for platform {platform}"
                    f" and reference genome {reference_genome}.")
        
        # Fetch relevant thresholds
        toggles_cfg = config.get('quality_control', {}).get('probe_qc', {})
        missing_thresh = toggles_cfg.get('missing_threshold', 0.05)

        adata = _probe_qc(adata, annotation, missing_thresh, toggles_cfg)
        logger.info("Successfully performed sample quality control.")


    # Update the AnnData object metadata to reflect processing
    if apply_probe_qc or apply_probe_qc:
        adata.uns['pipeline']['state'] = 'processed'
        adata.uns['pipeline']['steps'].append('qc')

    logger.info("~~~~~| Successfully Performed Quality Control |~~~~~\n")

    return (adata, qc_report)

# =====| Internal Helpers |=====================================================

def _sample_qc(adata: ad.AnnData, 
               missing_threshold: float = 0.20, 
               outlier_threshold: float = 1.5) -> tuple[ad.AnnData, list[dict]]:
    """
    Performs sample-level quality control upon DNA methylation values presented 
    as a CpG matrix AnnData object. Returns the quality-controlled CpG matrix 
    AnnData object with updated metadata suitable for further probe-level 
    quality control and/or preprocessing.

    Steps performed are toggled and configured in the user-provided 
    configurations, including the following options in order:

    1. Remove high missingness (above the provided threshold)
    2. Remove outliers from distribution (IQR-based)
    
    Parameters
    ----------
    adata : ad.AnnData
        AnnData object representing a project's raw DNA methylation data at the 
        CpG matrix level.
    missing_threshold : float
        Threshold of missing values per samples.
    outlier_threhsold : float
        Threshold for the number of standard deviations to mark outliers for 
        removal.

    Returns
    -------
    tuple[ad.AnnData, list[dict]]
        AnnData object representing a project's sample-level quality-controlled 
        DNA methylation data at the CpG matrix level with updated metadata. A 
        quality control status report for auditing.
    """

    # Compute quality control metrics
    cpg_matrix = np.asarray(adata.X)
    missing_rate = np.isnan(cpg_matrix).mean(axis = 1)
    sample_mean = np.nanmean(cpg_matrix, axis = 1)
    
    # Filter samples for missingness and global signal intensity
    fail_missing = missing_rate > missing_threshold
    mean_lo, mean_hi = _iqr_bounds(sample_mean, outlier_threshold)
    fail_mean = (sample_mean < mean_lo) | (sample_mean > mean_hi)

    # Generate failure mask using missingness and IQR
    fail_any = fail_missing | fail_mean
    pass_any = ~fail_any

    # Set global metadata
    adata.obs['missing_rate'] = missing_rate
    adata.obs['mean'] = sample_mean
    adata.obs['qc_pass'] = pass_any

    # Build an auditing report and apply the subsetting
    report = [
        {'file_id': fid, 'qc_pass': bool(pass_any[i])}
        for i, fid in enumerate(adata.obs.index.to_numpy())
    ]
    adata = adata[pass_any].copy()

    adata.uns.setdefault('qc', {})
    adata.uns['qc']['sample_qc'] = {
        'n_total': len(pass_any),
        'n_pass': int(pass_any.sum()),
        'n_fail': int(fail_any.sum()),
        'params': {
            'missing_threshold': missing_threshold,
            'outlier_threshold': outlier_threshold
        }
    }

    return (adata, report)



def _probe_qc(adata: ad.AnnData, 
              annotation: pd.DataFrame,
              missing_threshold: float,
              toggles_cfg: dict
              ) -> ad.AnnData:
    """
    Performs probe-level quality control upon DNA methylation values presented 
    as a CpG matrix AnnData object. Returns the quality-controlled CpG matrix 
    AnnData object with updated metadata suitable for further preprocessing.

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
    missing_threshold : float
        Threshold for the number of missing values to mark a probe for removal.
    toggles_cfg : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    ad.AnnData
        AnnData object representing a project's probe-level quality-controlled 
        DNA methylation data at the CpG matrix level with updated metadata.
    """

    # Remove probes that don't exist in the annotation (alignment)
    annotation = annotation.set_index('probe_id')
    common_probes = adata.var_names.intersection(annotation.index)
    adata._inplace_subset_var(adata.var_names.isin(common_probes))
    annotation = annotation.loc[common_probes]

    # Remove non-standard probes from the CpG matrix and annotation
    keep_standard = adata.var_names.str.startswith('cg')
    adata._inplace_subset_var(keep_standard)
    annotation = annotation.loc[adata.var_names]

    # Filter by missingness
    missing_rate = np.isnan(np.asarray(adata.X)).mean(axis = 0)
    keep_missing = pd.Series(
        missing_rate <= missing_threshold,
        index = adata.var_names
    )

    # Filter by annotation flags provided in the annotation object
    keep_annotation = pd.Series(True, index = annotation.index)

    remove_sex = toggles_cfg.get('remove_sex_chromosome', True)
    remove_SNP = toggles_cfg.get('remove_SNP_associated', True)
    remove_cr = toggles_cfg.get('remove_cross_reactive', True)
    remove_multi = toggles_cfg.get('remove_multi_mapped', True)

    if remove_sex: keep_annotation &= ~annotation['is_sex_chr']
    if remove_SNP:
        keep_annotation &= ~annotation['has_cpg_snp']
        keep_annotation &= ~annotation['has_sbe_snp']
        keep_annotation &= ~annotation['has_probe_snp']
    if remove_cr: keep_annotation &= ~annotation['is_cross_reactive']
    if remove_multi: keep_annotation &= ~annotation['is_multi_mapped']

    # Combine filters and subset the matrix
    pass_qc = keep_missing & keep_annotation

    # Set global metadata
    adata.var['missing_rate'] = missing_rate
    adata.var['qc_pass'] = pass_qc

    adata = adata[:, pass_qc].copy()

    adata.uns.setdefault('qc', {})
    adata.uns['qc']['probe_qc'] = {
        'n_total': len(pass_qc),
        'n_pass': int(pass_qc.sum()),
        'n_fail': int(~pass_qc.sum()),
        'params': {
            'missing_threshold': missing_threshold,
            'filter_toggles': {
                'remove_sex_chromosome': remove_sex,
                'remove_SNP_associated': remove_SNP,
                'remove_cross_reactive': remove_cr,
                'remove_multi_mapped': remove_multi
            }
        }
    }

    return adata


def _iqr_bounds(x, k):
    q1, q3 = np.nanpercentile(x, [25, 75])
    iqr = q3 - q1
    return q1 - k * iqr, q3 + k * iqr

# [END]