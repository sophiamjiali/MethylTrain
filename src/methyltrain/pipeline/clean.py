# ==============================================================================
# Script:           clean.py
# Purpose:          Internal cleaning functions and implementation logic
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-08
# ==============================================================================

import pandas as pd
import anndata as ad

from typing import Dict

from ..fs.layout import ProjectLayout


# for clean_data
def convert_to_parquet(config: Dict, layout: ProjectLayout) -> None:
    return

# =====| Clean Raw Data |=======================================================

def normalize_metadata():

    # for downloading logic


    # clean it to keep the most recent duplicate of file_id (double check actual field)
    # if not already done in query, get only one array type


    # clean the metadata and manifest for samples successfully downloaded
    
    return

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

    raw_metadata = layout.raw_metadata
    processed_metadata = layout.processed_metadata

    # Ensure the raw and processed metadata files are valid
    if not raw_metadata.exists():
        raise FileExistsError(f"Raw metadata file {raw_metadata} does "
                               "not exist.")
    if raw_metadata.suffix != ".csv":
        raise ValueError(f"Raw metadata field {raw_metadata} must have "
                         "extension `.csv`.")
    if processed_metadata.suffix != ".csv":
        raise ValueError(f"Processed metadata field {processed_metadata} must "
                         "have extension `.csv`.")

    # Load the raw metadata file
    raw_metadata = pd.read_csv(layout.raw_metadata)

    # Filter for samples based on the `id` field (standardized in TCGA)
    processed_metadata = raw_metadata[raw_metadata['id'].isin(adata.obs['id'])]

    # Save the filtered metadata
    processed_metadata.to_csv(layout.processed_metadata)
    

def clean_manifest(adata: ad.AnnData, layout: ProjectLayout) -> None:
    """
    Generates a processed manifest file by filtering the raw manifest file for 
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
        If the raw manifest file does not exist.
    ValueError
        If the raw or processed manifest file does not have extension `.csv`.
    """

    raw_manifest = layout.raw_metadata
    processed_manifest = layout.processed_metadata

    # Ensure the raw and processed metadata files are valid
    if not raw_manifest.exists():
        raise FileExistsError(f"Raw manifest file {raw_manifest} does "
                               "not exist.")
    if raw_manifest.suffix != ".csv":
        raise ValueError(f"Raw manifest field {raw_manifest} must have "
                         "extension `.csv`.")
    if processed_manifest.suffix != ".csv":
        raise ValueError(f"Processed manifest field {processed_manifest} must "
                         "have extension `.csv`.")

    # Load the raw metadata file
    raw_manifest = pd.read_csv(layout.raw_manifest)

    # Filter for samples based on the `id` field (standardized in TCGA)
    processed_manifest = raw_manifest[raw_manifest['id'].isin(adata.obs['id'])]

    # Save the filtered metadata
    processed_manifest.to_csv(layout.processed_manifest)
