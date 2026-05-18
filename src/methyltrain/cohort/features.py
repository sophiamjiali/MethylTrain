# ==============================================================================
# Script:           features.py
# Purpose:          Feature engineering implementation for the cohort workflow
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-13
# ==============================================================================

import logging

import anndata as ad
import numpy as np

from typing import Dict

logger = logging.getLogger(__name__)

# =====| Public API |===========================================================

def mad_probe_filtering(cohort: ad.AnnData, config: Dict) -> ad.AnnData:
    """
    Selects the top N most variable probes using median absolute deviation 
    (MAD) across samples.

    Parameters
    ----------
    cohort : ad.AnnData
        The Cohort AnnData to perform MAD probe filtering upon.
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    ad.AnnData
        The Cohort AnnData with MAD probe filtering performed.
    """

    n_probes = config.get('preprocessing', {}).get('MAD_probe_filtering', 30000)

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
    filtered_cohort.uns.setdefault('filtering', {})
    filtered_cohort.uns['filtering']['mad_probe_filtering'] = {
        'n_selected': n_probes,
        'n_original': cohort.n_vars,
        'method': 'MAD'
    }

    return filtered_cohort

# [END]