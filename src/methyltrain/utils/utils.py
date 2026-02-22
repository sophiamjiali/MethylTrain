# ==============================================================================
# Script:           utils.py
# Purpose:          General utility functions for the package
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-08
# ==============================================================================

import hashlib
import shutil

import numpy as np
import pandas as pd

from pathlib import Path
from typing import Dict, Any, Tuple

from ..constants import (
    PLATFORM_TYPES, 
    REFERENCE_GENOME_TYPES,
    ANNOTATION_hg19_PATHS,
    ANNOTATION_hg38_PATHS
)

from .. import constants

# ======| File I/O Utilities |==================================================

def load_sample(file: Path) -> pd.Series:
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
    if platform not in constants.PLATFORM_TYPES:
        raise ValueError(f"Platform {platform} was not recognized from the "
                         f"supported types: {constants.PLATFORM_TYPES}")

    if reference_genome not in constants.REFERENCE_GENOME_TYPES:
        raise ValueError(f"Reference genome {reference_genome} was not "
                         f"recognized from the supported types: "
                         f"{constants.REFERENCE_GENOME_TYPES}")
    
    # Load the appropriate genome build annotation path (provided by package)
    if reference_genome == "GRCh37":
        anno_path = constants.ANNOTATION_hg19_PATHS[platform]
    elif reference_genome == "GRCh38":
        anno_path = constants.ANNOTATION_hg38_PATHS[platform]
    else:
        raise ValueError(f"Reference genome {reference_genome} was not "
                         f"recognized from the supported types: "
                         f"{constants.REFERENCE_GENOME_TYPES}")

    # # Verify the array type and genome build provided are valid
    # if platform not in PLATFORM_TYPES:
    #     raise ValueError(f"Platform {platform} was not recognized from the "
    #                      f"supported types: {PLATFORM_TYPES}")

    # if reference_genome not in REFERENCE_GENOME_TYPES:
    #     raise ValueError(f"Reference genome {reference_genome} was not "
    #                      f"recognized from the supported types: "
    #                      f"{REFERENCE_GENOME_TYPES}")
    
    # # Load the appropriate genome build annotation path (provided by package)
    # if reference_genome == "GRCh37":
    #     anno_path = ANNOTATION_hg19_PATHS[platform]
    # elif reference_genome == "GRCh38":
    #     anno_path = ANNOTATION_hg38_PATHS[platform]
    # else:
    #     raise ValueError(f"Reference genome {reference_genome} was not "
    #                      f"recognized from the supported types: "
    #                      f"{REFERENCE_GENOME_TYPES}")

    return pd.read_parquet(anno_path)


# ======| Dictionary Utilities |================================================

def merge_dicts(base: Dict[str, Any], override: Dict[str, Any]):
    """
    Recursively merge two dictionaries. Values from `override` take precedence.

    Parameters
    ----------
    base : dict
        Base dictionary (e.g., defaults).
    override : dict
        Dictionary with overriding values (e.g., user config).

    Returns
    -------
    merged : dict
        Merged dictionary with user values overriding defaults.
    """

    result = base.copy()
    for key, value in override.items():
        if (key in result and isinstance(result[key], dict) 
                          and isinstance(value, dict)):
            result[key] = merge_dicts(result[key], value)
        else:
            result[key] = value
    return result


def check_dict(default: Dict[str, Any], 
               user: Dict[str, Any], 
               path: str = "") -> None:
    """
    Recursively verify that a user-provided configuration dictionary matches
    the structure, types, and constraints of a default configuration dictionary.

    Parameters
    ----------
    default : dict
        The default configuration dictionary serving as the schema. Keys and
        value types define what is valid.
    user : dict
        The user-provided configuration dictionary to validate.
    path : str, optional
        Dot-separated path of keys used internally to indicate the location
        in the nested dictionary for informative error messages. Default is "".

    Raises
    ------
    KeyError
        If a required key from `default` is missing in `user`.
    TypeError
        If the type of a value in `user` does not match the type in `default`.
    ValueError
        If a value in `user` violates a constraint (e.g., allowed enum values).
    """

    for key, val in default.items():
        full_key = f"{path}.{key}" if path else key

        # Verify the key exist in the user configurations
        if key not in user:
            raise KeyError(f"Missing key in configuration: {full_key}")
        
        user_val = user[key]

        # If the value itself is a dictionary in the default, recurse
        if isinstance(val, dict):
            if not isinstance(user_val, dict):
                raise TypeError(f"Key '{full_key}' must be a dictionary")
            
            # Do not check metadata as fields will vary
            if key != "metadata":
                check_dict(val, user_val, path = full_key)
                
        else:

            # Verify the value type
            if not isinstance(user_val, type(val)):
                raise TypeError(f"Key '{full_key}' must be of type "
                                f"{type(val).__name__}")
            
            # Verify value ranges individually
            if key == "split":
                if sum(user_val) != 1:
                    raise ValueError(f"Key '{full_key}' must sum to 1.0")
                
            elif key == "missing_threshold" or key == "outlier_threshold":
                if user_val < 0:
                    raise ValueError(f"Key '{full_key}' must be "
                                     "greater or equal to zero")
            
            elif key == "clip_values":
                if user_val[0] < 0 or user_val[1] > 1:
                    raise ValueError(f"Key '{full_key}' must range between "
                                     "zero and one")
                
            elif key == "array_type":
                if user_val not in PLATFORM_TYPES:
                    raise ValueError(
                        f"Key '{full_key}' must be one of the following "
                        f"supported types: {PLATFORM_TYPES}"
                    )
                
            elif key == "genome_build":
                if user_val not in REFERENCE_GENOME_TYPES:
                    raise ValueError(
                        f"Key '{full_key}' must be one of the following "
                        f"supported types: {REFERENCE_GENOME_TYPES}"
                    )
                
# ======| API Utilities |=======================================================

def verify_gdc_client(gdc_client_path) -> None:
    # Verifies that the gdc-client was downloaded

    if shutil.which(gdc_client_path) is None:
        raise RuntimeError(
            f"gdc-client executable not found at '{gdc_client_path}'.\n"
            "Please install the official GDC Data Transfer Tool from:\n"
            "https://gdc.cancer.gov/access-data/gdc-data-transfer-tool\n"
            "and ensure it is available in your PATH."
        )

def extract_project_id(cases):
    # Helper to extract nested metadata fields
    try:
        return (cases[0]['project']['project_id'] if cases 
                and cases[0].get('project') else pd.NA)
    except Exception:
        return pd.NA
    
def extract_file_name(files):
    # Helper to extract nested metadata fields
    try:
        return (files[0]['file_name'] if files 
                and files[0].get('file_name') else pd.NA)
    except Exception:
        return pd.NA

def extract_sample_type(cases):
    # Helper to extract nested metadata fields
    try:
        return (cases[0]['samples'][0]['sample_type'] if cases 
                and cases[0].get('samples') else pd.NA)
    except Exception:
        return pd.NA

def extract_submitter_id(cases):
    # Helper to extract nested metadata fields
    try:
        return cases[0]['submitter_id'] if cases else pd.NA
    except Exception:
        return pd.NA
    
def extract_aliquot_id(cases):
    # Helper to extract nested metadata fields
    field = (cases[0]['samples'][0]['portions'][0]['analytes'][0]
             ['aliquots'][0]['aliquot_id'])
    try:
        return (field if cases else pd.NA)
    except Exception:
        return pd.NA

def extract_batch_id(barcode: str):
    # Returns the portion and the plate ID together
    
    if pd.isna(barcode): return None
    parts = barcode.split('-')
    if len(parts) >= 6:
        return f"{parts[3]}-{parts[5]}"  # portion + plate
    return None

# ======| Computation Utilities |===============================================

def iqr_bounds(x, k):
    q1, q3 = np.nanpercentile(x, [25, 75])
    iqr = q3 - q1
    return q1 - k * iqr, q3 + k * iqr
