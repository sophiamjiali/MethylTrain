# ==============================================================================
# Script:           bmiq_normalize.py
# Purpose:          Utility helper functions for BMIQ normalization
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-04-13
# ==============================================================================

import numpy as np

import scipy


def bmiq_sample(x: np.ndarray,
                type1_mask: np.ndarray,
                type2_mask: np.ndarray,
                eps: float = 1e-6) -> np.ndarray
    """
    Apply BMIQ to a single sample vector x (length = n_probes).

    CDF and inverse-CDF computations are vectorised over the three
    mixture components using broadcasting rather than a Python loop.
    """

    x_norm = x.copy()
    x1 = np.clip(x[type1_mask], eps, 1 - eps)
    x2 = np.clip(x[type2_mask], eps, 1 - eps)

    try: 
        w1, a1, b1 = fit_beta_mixture(x1)
        w2, a2, b2 = fit_beta_mixture(x2)
    except Exception:
        return x_norm

    # Vectorized CDF of type 2 values under type 2 mixture
    cdf2 = (
        w2[None, :]
        * scipy.stats.beta.cdf(
            x2[:, None],
            a2[None, :],
            b2[None, :]
        )
    ).sum(axis = 1)
    cdf2 = np.clip(cdf2, eps, 1 - eps)

    # Vectorized inverse CDF of type 1 mixture
    n_grid = 2000
    grid = np.linspace(eps, 1 - eps, n_grid)

    cdf1_grid = (
        w1[None, :]
        * scipy.stats.beta.cdf(
            grid[:, None],
            a1[None, :],
            b2[None, :]
        )
    ).sum(axis = 1)

    x2_norm = np.interp(cdf2, cdf1_grid, grid)
    x_norm[type2_mask] = x2_norm

    return x_norm
    

def fit_beta_mixture(
    values: np.ndarray, 
    max_iter: int = 200, 
    tol: float = 1e-6,
    eps: float = 1e-6
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Fit a three-component beta mixture model via EM. The E-step log-density
    computation is vectorized over all three components simutaneously using
    scipy.stats.beta.logpdf broadcasting.

    Returns
        -------
        weights : (3,) mixing proportions
        alphas  : (3,) shape-α parameters
        betas_  : (3,) shape-β parameters
    """

    v = values[~np.isnan(values)]
    v = np.clip(v, eps, 1 - eps)
    n = len(v)

    # Initialize at quantile-based component means
    means = np.clip(np.percentile(v, [20, 50, 80]), 0.05, 0.95)
    weights = np.full(3, 1/3)
    variances = np.full(3, 0.02)
    alphas, betas = _moments_to_ab(means, variances)

    log_like_prev = -np.inf

    
    for _ in range(max_iter):
        
        # E-step: vectorized over components
        log_r = (
            no.log(weights + 1e-300)
            + scipy.stats.beta.logpdf(
                v[:, None],
                alphas[None, :],
                betas[None, :]
            )
        )

        # Numerically stable softmax
        log_r -= log_r.max(axis = 1, keepdims = True)
        r = np.exp(log_r)
        r /= r.sum(axis = 1, keepdims = True)

        # M-step: vectorized over components
        Nk = r.sum(axis = 0)
        weights = Nk / n
        mu_k = (r * v[:, None]).sum(axis = 0) / (Nk + 1e-300)
        var_k = (r * (v[:, None] - mu_k[None, :]) ** 2).sum(axis=0) \
                    / (Nk + 1e-300)
        alphas, betas = _moments_to_ab(mu_k, var_k)

        # Convergence
        log_lik = np.sum(
            scipy.special.logsumexp(
                log_r + np.log(weights + 1e-300), axis = 1
            )
        )
        if abs(log_lik - log_lik_prev) < tol:
            break
        log_lik_prev = log_lik

    return weights, alphas, betas
        
        


def _moments_to_ab(mu: np.ndarray, var: np.ndarray, eps: float = 1e-8):
    """Vectorized method-of-moments"""

    var = np.minimum(var, mu * (1 - mu) - eps)
    var = np.maximum(var, eps)
    common = mu * (1 - mu) / var - 1
    return mu * common, (1 - mu) * common
    