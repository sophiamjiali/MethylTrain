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

from methyltrain.constants.annotation import (
    PLATFORM_TYPES,
    PLATFORM_PRIORITY,
    REFERENCE_GENOME_TYPES,
    ANNOTATION_hg19_PATHS,
    ANNOTATION_hg38_PATHS
)

def load_metadata(layout: ProjectLayout) -> pd.DataFrame:
    # Loads the metadata table with file_id as index

    layout.validate()
    metadata = pd.read_csv(layout.metadata, sep = '\t', index_col = 0)
    return metadata


def load_manifest(layout: ProjectLayout) -> pd.DataFrame:
    # Loads the manifest with file_id as index

    layout.validate()
    manifest = pd.read_csv(layout.manifest, sep = '\t', index_col = 0)
    return manifest


def load_annotation(platform: str, reference_genome: str) -> pd.DataFrame:
    """
    Load an Illumina annotation based on array type and genome build.

    Parameters
    ----------
    platform : str
        The project's DNA methylation array type (e.g. Illumina 27K, 450K, EPIC)
    reference_genome : str
        The project's genome build (e.g. hg19, hg38)

    Returns
    -------
    pd.DataFrame
        Illumina annotation for the given array type and genome build.

    Raises
    ------
    ValueError
        If the array type or genome build provided in the user-configurations 
        is not valid.
    """

    # Verify the array type and genome build provided are valid
    if platform not in PLATFORM_TYPES:
        raise ValueError(f"Platform {platform} was not recognized from the "
                         f"supported types: {PLATFORM_TYPES}")

    if reference_genome not in REFERENCE_GENOME_TYPES:
        raise ValueError(f"Reference genome {reference_genome} was not "
                         f"recognized from the supported types: "
                         f"{REFERENCE_GENOME_TYPES}")
    
    # Load the appropriate genome build annotation path (provided by package)
    if reference_genome == "GRCh37":
        anno_path = ANNOTATION_hg19_PATHS[platform]
    elif reference_genome == "GRCh38":
        anno_path = ANNOTATION_hg38_PATHS[platform]
    else:
        raise ValueError(f"Reference genome {reference_genome} was not "
                         f"recognized from the supported types: "
                         f"{REFERENCE_GENOME_TYPES}")

    return pd.read_parquet(anno_path)


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