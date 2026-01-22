# ==============================================================================
# Script:           clean.py
# Purpose:          Internal cleaning functions and implementation logic
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-08
# ==============================================================================

import pandas as pd
import anndata as ad

from ..fs.layout import ProjectLayout

# =====| Clean Metadata |=======================================================

def clean_metadata(adata: ad.AnnData, layout: ProjectLayout) -> None:
    """
    Generates a processed metadata file by filtering the raw metadata file for 
    samples present in the quality-controlled project DNA methylation AnnData 
    object.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object representing a project's raw DNA methylation data at the 
        CpG matrix level.
    layout : ProjectLayout
        Object representing a project dataset directory layout.

    Raises
    ------
    FileNotFoundError
        If the raw metadata file does not exist.
    ValueError
        If the raw or processed metadata file does not have extension `.csv`.
    """

    # Ensure the raw and processed metadata files are valid
    if not layout.metadata.exists():
        raise FileExistsError(f"Metadata file {layout.metadata} does "
                               "not exist.")
    if layout.metadata.suffix != ".csv":
        raise ValueError(f"Metadata field {layout.metadata} must have "
                         "extension `.csv`.")
    
    # Filter out samples that didn't survive quality control
    metadata = pd.read_csv(layout.metadata)
    metadata = metadata[metadata['file_id'].isin(adata.obs['file_id'])]
    metadata.to_csv(layout.metadata)
