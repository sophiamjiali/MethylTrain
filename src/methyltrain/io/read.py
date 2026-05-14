# ==============================================================================
# Script:           read.py
# Purpose:          Reading utility functions for the package
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-08
# ==============================================================================

import pandas as pd

from pathlib import Path

from methyltrain.project.layout import ProjectLayout

def _load_metadata(layout: ProjectLayout) -> pd.DataFrame:
    # Loads the metadata table with file_id as index

    layout.validate()
    metadata = pd.read_csv(layout.metadata, sep = '\t', index_col = 0)
    return metadata


def _load_manifest(layout: ProjectLayout) -> pd.DataFrame:
    # Loads the manifest with file_id as index

    layout.validate()
    manifest = pd.read_csv(layout.manifest, sep = '\t', index_col = 0)
    return manifest


def _load_sample(file: Path) -> pd.Series:
    """
    Loads the raw DNA methylation data of a given sample, provided as a 
    `.parquet` with CpG probe ID as the index and `beta_value` as the column 
    name.

    Parameters
    ----------
    file : Path
        Path to a .parquet file containing the beta values of a sample.

    Returns
    -------
    pd.Series
        Beta values of a sample loaded as a Series.

    Raises
    ------
    FileNotFoundError
        If the file path does not exist or is not `.parquet`.
    """

    # Verify the file exists and is a `.parquet` file
    if not file.exists():
        raise FileNotFoundError(f"File was not found: {file}")
    if file.suffix != ".parquet":
        raise FileNotFoundError(f"File must be a `.parquet`: {file}")
    
    sample = pd.read_parquet(file)
    sample = pd.Series(sample['beta_value'].values, 
                       index = sample['probe_id'], 
                       name = str(file.name))
    return sample

# [END]