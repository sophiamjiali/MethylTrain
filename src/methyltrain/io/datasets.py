# ==============================================================================
# Script:           qc.py
# Purpose:          I/O functionality for projects and cohorts
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-14
# ==============================================================================

import logging

import pandas as pd
import anndata as ad
import numpy as np

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

from methyltrain.project.layout import ProjectLayout
from methyltrain.cohort.layout import CohortLayout
from methyltrain.io.read import load_metadata, _load_sample

logger = logging.getLogger(__name__)

# =====| Public API |===========================================================

# ~~~~~| Loading |~~~~~
def load_raw_project(config: dict, layout: ProjectLayout):
    """
    Load the raw DNA methylation data of a project as an AnnData object from 
    `.parquet` files in the raw data directory. Column metadata is initialized 
    as the sample ID field specified in the user-provided configurations.

    All raw DNA methylation data files are loaded in parallel using for more 
    efficient loading. Note that `ThreadPoolExecutor.map()` explicitly 
    preserves input order such that the native order of `project_dir` can be 
    used to align sample IDs using the metadata.

    Assumes metadata is perfectly alligned with the data available in the 
    project raw data directory (as per the download() function).

    Default behaviour resolves case-level duplicates (aliquots) by retaining 
    only the first replicate. Performing mean aggregation across aliquots is not
    advised.
    
    Parameters
    ----------
    config: dict
        Configuration dictionary controlling workflow steps.
    layout : ProjectLayout
        Object representing a project dataset directory layout.

    Returns
    -------
    ad.AnnData
        Raw DNA methylation and metadata of the specified project loaded as an 
        AnnData object.

    Raises
    ------
    FileNotFoundError
        If the project directory path does not exist or is empty.
    """

    # Verify all directories exist and are accessible
    layout.validate()
    _validate_directories(layout)

    # Build the raw beta value matrix and align its metadata
    cpg_matrix = _build_beta_matrix(layout)
    metadata = _build_metadata(layout)

    # Initialize the CpG matrix as an AnnData object with aligned metadata
    adata = ad.AnnData(
        X = cpg_matrix.T.values.astype(np.float32),
        obs = metadata,
        var = pd.DataFrame(index = cpg_matrix.index.astype(str))
    )

    # Initialize global metadata for the project
    adata.uns = {
        "provenance": {
            "project_id": layout.project_name,
            'level': 'project',
            "data_type": "cpg_matrix",
            "conversion": "beta_value",
            "pipeline_verson": config['version'],
        },
        "data_source": {
            "platform": config.get("download", {}).get("platform"),
            "reference_genome": config.get("reference_genome"),
        },
        "pipeline": {
            "state": "raw",
            "steps": [],
        },
    }
    
    return adata

def load_processed_project(path: str) -> ad.AnnData:
    """
    Loads a processed project AnnData object.

    Parameters
    ----------
    processed_file : Path or str
        Full path to the processed .h5ad file.

    Returns
    -------
    ad.AnnData
        The loaded processed dataset.

    Raises
    ------
    FileNotFoundError
        If the processed file path does not exist.
    """
    processed_file = Path(path)

    if not processed_file.exists():
        raise FileNotFoundError(f"Processed file not found: {processed_file}")
    return ad.read_h5ad(processed_file)


# ~~~~~| Saving |~~~~~
def save_project(adata: ad.AnnData, layout: ProjectLayout) -> None:
    """
    Saves a project AnnData object.

    The data is not converted to a sparse matrix as beta values will likely see 
    little performance benefit.

    Parameters
    ----------
    adata : ad.AnnData
        Project AnnData object containing DNA methylation data and metadata.
    layout : ProjectLayout
        Object representing a project dataset directory layout.
    """

    layout.validate()
    adata.write_h5ad(layout.adata, compression = "gzip")
    return


def save_cohort(adata: ad.AnnData, layout: CohortLayout):
    """
    Saves a project AnnData object.

    The data is not converted to a sparse matrix as beta values will likely see 
    little performance benefit.

    Parameters
    ----------
    adata : ad.AnnData
        Project AnnData object containing DNA methylation data and metadata.
    layout : ProjectLayout
        Object representing a project dataset directory layout.
    """

    layout.validate()
    adata.write_h5ad(layout.cohort_adata, compression = "gzip")
    return

def save_cohort_train(adata: ad.AnnData | None, layout: CohortLayout):
    layout.validate()
    if adata is not None: 
        adata.write_h5ad(layout.train_adata, compression = "gzip")

def save_cohort_val(adata: ad.AnnData | None, layout: CohortLayout):
    layout.validate()
    if adata is not None: 
        adata.write_h5ad(layout.val_adata, compression = "gzip")

def save_cohort_test(adata: ad.AnnData | None, layout: CohortLayout):
    layout.validate()
    if adata is not None: 
        adata.write_h5ad(layout.test_adata, compression = "gzip")


# =====| Internal Helpers |=====================================================

def _validate_directories(layout: ProjectLayout) -> None:
    """
    Validates the existence and accessibility of the raw data directory used to 
    load the raw beta value matrix.
    """

    raw_dir = layout.raw_dir
    if not raw_dir.is_dir():
        raise FileExistsError(f"Project directory was not found: {raw_dir}")
    if not any(raw_dir.iterdir()):
        raise FileExistsError(f"Project directory is empty: {raw_dir}")
    return


def _build_beta_matrix(layout: ProjectLayout) -> pd.DataFrame:
    """
    Loads all beta values in parallel, concatenating upon the index to build a 
    matrix of raw balues.
    """

    # Load all beta values in parallel as a list of Pandas DataFrames
    files = [f for f in layout.raw_dir.iterdir() if f.suffix == ".parquet"]

    with ThreadPoolExecutor() as ex:
        sample_beta_values = list(ex.map(_load_sample, files))

    # Concatenate on the index to build a matrix
    cpg_matrix = pd.concat(sample_beta_values, axis = 1, join = "outer")
    cpg_matrix = cpg_matrix.sort_index()

    return cpg_matrix


def _build_metadata(layout: ProjectLayout) -> pd.DataFrame:
    """
    Loads the project's metadata, aligning it with the raw beta value matrix.
    """
    metadata = load_metadata(layout)
    metadata = metadata.sort_values(by = 'file_name')
    metadata = metadata.loc[metadata['status'] == 'success']

    return metadata



