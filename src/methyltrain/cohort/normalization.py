# ==============================================================================
# Script:           normalization.py
# Purpose:          Internal normalization functions and implementation logic
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-19
# ==============================================================================

import logging

import anndata as ad
import numpy as np

from typing import Dict, Optional, Tuple
from dataclasses import dataclass

logger = logging.getLogger(__name__)

# =====| Internal Class Definition |============================================

@dataclass
class NormalizationState:
    winsorize_enabled: bool
    scale_enabled: bool

    # Winsorization
    lower_bounds: Optional[np.ndarray] = None
    upper_bounds: Optional[np.ndarray] = None

    # Scaling
    scale_min: Optional[np.ndarray] = None
    scale_max: Optional[np.ndarray] = None



# =====| Public API |===========================================================

def normalize(adata: ad.AnnData, 
              config: Dict, 
              norm_state = None) -> Tuple[ad.AnnData, NormalizationState]:
    """
    Applies global winsorization and Min-Max scaling to map M-values to the 
    [-1, 1] range optimized for VAE/Diffusion gradients.

    Winsorization should only be performed at the data level the machine 
    learning model will encounter (i.e. gene-level), as clipping extreme values 
    is only meant to improve model fidelity, rather than a standard 
    preprocessing practice.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object representing a project's quality-controlled, 
        preprocessed, and batch effect corrected DNA methylation data at the 
        CpG probe- or gene-level.
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    ad.AnnData
        The AnnData object with its values clipped and scaled.
    """

    logger.info("=====| Attempting to Normalize the Cohort |=====")

    if norm_state is None:
        norm_state = fit_normalization(adata, config)
        logger.info("Successfully fit for nornalization.")

    
    adata = apply_normalization(adata, norm_state)
    logger.info("Successfully applied normalization.")

    logger.info("=====| Successfully Normalized the Cohort |=====")
    return adata, norm_state


def fit_normalization(adata: ad.AnnData, config: Dict) -> NormalizationState:
    apply_winsorization = config.get('toggles', {}).get('winsorize', True)
    apply_scaling = config.get('toggles', {}).get('scale', True)

    X = np.asarray(adata.X)

    lower = upper = None
    min_ = max_ = None

    # Optionally fit and apply winsorization (for later use by scaling)
    if apply_winsorization:

        # Fetch the user configurations for winsorization
        winsorize_cfg = config['preprocessing']['winsorization']
        q_low = winsorize_cfg.get('lower', 0.01)
        q_high = winsorize_cfg.get('upper', 0.99)

        lower, upper = _fit_winsorize(X, q_low, q_high)
        X_proc = _apply_winsorize(X, lower, upper)
    else: X_proc = X

    # Optionally fit the Min-Max scaling to the winsorized data
    if apply_scaling:
        min_, max_ = _fit_scale(X_proc)

    return NormalizationState(
        winsorize_enabled = apply_winsorization,
        scale_enabled = apply_scaling,
        lower_bounds = lower,
        upper_bounds = upper,
        scale_min = min_,
        scale_max = max_
    )


def apply_normalization(adata: ad.AnnData, 
                        state: NormalizationState) -> ad.AnnData:
    X = np.asarray(adata.X)

    if state.winsorize_enabled or state.scale_enabled:
        adata.uns.setdefault('normalization', {})

    if state.winsorize_enabled:

        if state.lower_bounds is None or state.upper_bounds is None:
            raise ValueError(
                "Winsorization bounds are missing from normalization state."
            )
        
        adata.uns['normalization']['winsorization'] = {
            'lower_bound': state.lower_bounds,
            'upper_bound': state.upper_bounds
        }
        
        X = _apply_winsorize(X, state.lower_bounds, state.upper_bounds)

    if state.scale_enabled:

        if state.scale_min is None or state.scale_max is None:
            raise ValueError(
                "Scaling parameters are missing from normalization state."
            )
        
        adata.uns.setdefault('normalization', {})
        adata.uns['normalization']['scaling'] = {
            'min': state.scale_min,
            'max': state.scale_max
        }

        X = _apply_scale(X, state.scale_min, state.scale_max)

    adata = adata.copy()
    adata.X = X
    return adata

# =====| Internal Helpers |=====================================================

def _fit_winsorize(X: np.ndarray,
                   q_low: float, 
                   q_high: float) -> Tuple[np.ndarray, np.ndarray]:
    lower = np.quantile(X, q_low, axis = 0)
    upper = np.quantile(X, q_high, axis = 0)
    return upper, lower


def _fit_scale(X: np.ndarray,) -> Tuple[np.ndarray, np.ndarray]:
    min_ = np.min(X, axis = 0)
    max_ = np.max(X, axis = 0)
    max_ = np.where(max_ == min_, min_ + 1e-8, max_)
    return min_, max_


def _apply_winsorize(X: np.ndarray,
                     lower: np.ndarray,
                     upper: np.ndarray) -> np.ndarray:
    return np.clip(X, lower, upper)


def _apply_scale(X: np.ndarray,
                 min_: np.ndarray,
                 max_: np.ndarray) -> np.ndarray:
    
    # Masc constant features
    valid = (max_ - min_) > 0
    X_norm = np.zeros_like(X)
    X_norm[:, valid] = (X[:, valid] - min_[valid]) / (max_[valid] - min_[valid])
    return X_norm

# [END]