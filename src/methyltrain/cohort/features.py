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
    Defines the feature space for the cohort, optionally performing MAD probe 
    filtering and returning the final probe set as a pandas DataFrame.

    Parameters
    ----------
    cohort : ad.AnnData
        The Cohort AnnData to perform MAD probe filtering upon.
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    tuple[ad.AnnData, pd.DataFrame]
        The optionally filtered cohort AnnData object and the final probe set.
    """

    logger.info("=====| Attempting to Define the Feature Space |=====")

    # Perform MAD probe filtering to reduce probe dimensionality
    if config.get('toggles', {}).get('MAD_probe_filtering', True):
        n_p = config.get('preprocessing', {}).get('MAD_probe_filtering', 30000)
        cohort = _mad_probe_filtering(cohort, n_p)
        logger.info("Successfully performed MAD probe filtering.")

    # Extract the final probe set for documentation
    probe_set = _extract_probe_set(cohort)
    logger.info("Successfully extracted the final probe set.")

    cohort.uns['pipeline']['state'] = 'processed'
    logger.info("=====| Successfully Defined the Feature Space |=====")
    
    return (cohort, probe_set)
    

# =====| Internal Helpers |=====================================================

def _mad_probe_filtering(cohort: ad.AnnData, n_probes: int) -> ad.AnnData:
    """
    Selects the top N most variable probes using median absolute deviation 
    (MAD) across samples.

    Parameters
    ----------
    cohort : ad.AnnData
        The Cohort AnnData to perform MAD probe filtering upon.
    n_probes : int
        The number of probes to produce after filtering.

    Returns
    -------
    ad.AnnData
        The Cohort AnnData with MAD probe filtering performed.
    """

    X = np.asarray(cohort.X, dtype = np.float32)
    med = np.median(X, axis = 0)
    mad = np.median(np.abs(X - med), axis = 0)

    # Sanitize the numerical edge cases
    mad = np.nan_to_num(mad, nan = 0.0, posinf = 0.0, neginf = 0.0)
    top_idx = np.argpartition(mad, -n_probes)[-n_probes:]
    top_idx = top_idx[np.argsort(mad[top_idx])[::-1]]

    filtered_cohort = cohort[:, top_idx].copy()
    filtered_cohort.var['mad_score'] = mad[top_idx]

    filtered_cohort.uns['pipeline']['state'] = 'processed'
    filtered_cohort.uns['pipeline']['steps'].append("MAD_probe_filtering")
    filtered_cohort.uns.setdefault('features', {})
    filtered_cohort.uns['features']['mad_probe_filtering'] = {
        'n_selected': n_probes,
        'n_original': cohort.n_vars,
        'method': 'MAD'
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