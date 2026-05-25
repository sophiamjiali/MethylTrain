# ==============================================================================
# Script:           features.py
# Purpose:          Feature engineering implementation for the cohort workflow
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-13
# ==============================================================================

import logging

import anndata as ad
import pandas as pd
import numpy as np

from typing import Dict

logger = logging.getLogger(__name__)

# =====| Public API |===========================================================

def define_feature_space(cohort: ad.AnnData, 
                         config: Dict) -> tuple[ad.AnnData, pd.DataFrame]:
    """
    Defines the feature space for the cohort, optionally performing stratified 
    MAD probe filtering and returning the final probe set as a pandas DataFrame.

    Parameters
    ----------
    cohort : ad.AnnData
        The Cohort AnnData to perform MAD probe filtering upon.
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    tuple[ad.AnnData, pd.DataFrame]
        The optionally stratified filtered cohort AnnData object and the final 
        probe set.
    """

    logger.info("=====| Attempting to Define the Feature Space |=====")

    # Perform stratified MAD probe filtering to reduce probe dimensionality
    if config.get('toggles', {}).get('MAD_probe_filtering', True):

        # Extract relevant configurations for stratified MAD filtering
        mad_cfg = config.get('preprocessing', {})['stratified_MAD_filter']
        n_p = mad_cfg.get('n_probes', 30000)
        top_n_per_class = mad_cfg.get('top_n_per_class', 15000)

        cohort = _stratified_mad_probe_filter(cohort, n_p, top_n_per_class)
        logger.info("Successfully performed stratified MAD probe filtering.")

    # Extract the final probe set for documentation
    probe_set = _extract_probe_set(cohort)
    logger.info("Successfully extracted the final probe set.")

    cohort.uns['pipeline']['state'] = 'processed'
    logger.info("=====| Successfully Defined the Feature Space |=====")
    
    return (cohort, probe_set)
    

# =====| Internal Helpers |=====================================================

def _stratified_mad_probe_filter(cohort: ad.AnnData, 
                                 n_probes: int,
                                 top_n_per_class: int) -> ad.AnnData:
    """
    Selects the top N most variable probes across the cohort using stratified 
    within-project Median Absolute Deviation (MAD) to prevent tissue-of-origin dominance.

    Parameters
    ----------
    cohort : ad.AnnData
        The Cohort AnnData to perform stratified filtering upon. Must contain 
        'project_id' within cohort.obs.
    n_probes : int
        The final target number of features to retain globally.
    top_n_per_class : int
        The number of highly variable probes to select within each unique project.

    Returns
    -------
    ad.AnnData
        The stratified filtered Cohort AnnData object.
    """

    unique_projects = cohort.obs['project_id'].unique()
    probe_names = cohort.var_names.to_numpy()

    # Track cross-project presence for every probe
    occurrence_counts = np.zeros(cohort.n_vars, dtype = np.int32)
    logger.info(f"Computing stratified MAD profiles across {len(unique_projects)} unique project groups.")
    
    for project in unique_projects:

        # Extract the current project from the cohort
        project_view = cohort[cohort.obs['project_id'] == project]
        X_project = np.asarray(project_view.X, dtype = np.float32)

        # Calculate within-class median absolute deviation (MAD)
        med = np.median(X_project, axis = 0)
        mad = np.median(np.abs(X_project - med), axis = 0)
        mad = np.nan_to_num(mad, nan = 0.0, posinf = 0.0, neginf = 0.0)

        # Identify local highly variable indices
        top_idx = np.argpartition(mad, - top_n_per_class)[-top_n_per_class:]
        occurrence_counts[top_idx] += 1

    # Select globally relevant probes by frequency of high-variance inclusion
    final_idx = np.argpartition(occurrence_counts, -n_probes)[-n_probes:]
    final_idx = final_idx[np.argsort(occurrence_counts[final_idx])[::-1]]

    min_project_overlap = int(np.min(occurrence_counts[final_idx]))
    logger.info(f"Stratified feature selection bounding criteria met. Minimum project occurrence threshold: {min_project_overlap}")
    
    # Slice and document pipeline transformations
    filtered_cohort = cohort[:, final_idx].copy()
    filtered_cohort.var['occurrence_score'] = occurrence_counts[final_idx]

    filtered_cohort.uns['pipeline']['state'] = 'processed'
    filtered_cohort.uns['pipeline']['steps'].append("stratified_MAD_filter")
    filtered_cohort.uns.setdefault('features', {})
    filtered_cohort.uns['features']['stratified_MAD_filter'] = {
        'n_selected': n_probes,
        'n_original': cohort.n_vars,
        'top_n_per_class': top_n_per_class,
        'minimum_project_overlap': min_project_overlap,
        'method': 'stratified_MAD'
    }

    return filtered_cohort


def _extract_probe_set(cohort: ad.AnnData) -> pd.DataFrame:
    """
    Extracts the final probe set from an AnnData object and returns it as a pandas DataFrame for reproducibility and downstream reference.

    Parameters
    ----------
    cohort : ad.AnnData
        The Cohort AnnData to extract the probe set from.

    Returns
    -------
    pd.DataFrame
        The final probe set in order.
    """

    return pd.DataFrame({
        "probe_order": np.arange(cohort.n_vars),
        "probe_id": np.asarray(cohort.var_names).astype(str)
    })

# [END]