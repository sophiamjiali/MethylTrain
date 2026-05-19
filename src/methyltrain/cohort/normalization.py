
def normalize()
    
    # fit stage, then apply/transform stage


def _winsorize_fit(cohort: ad.AnnData, config: Dict) -> ad.AnnData:
    """
    Applies global winsorization to improve numerical stability for VAE/
    Diffusion gradients.

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

    clip_values = config.get('clip_and_scale', [0.01, 0.99])
    X = np.asarray(cohort.X, dtype = np.float32)

    # Perform global winsorization across the full matrix
    lower_val = np.percentile(X, clip_values[0] * 100)
    upper_val = np.percentile(X, clip_values[1] * 100)
    X = np.clip(X, lower_val, upper_val)
    
    cohort.X = X.astype(np.float32)

    # Update the object metadata
    cohort.uns['pipeline']['state'] = 'processed'
    cohort.uns.setdefault('features', {})

    adata.uns['features']['winsorization'] = {
        'lower_bound_mval': lower_val,
        'upper_bound_mval': upper_val
    }
    return cohort
# [END]