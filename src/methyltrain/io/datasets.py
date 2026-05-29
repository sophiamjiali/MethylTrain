# ==============================================================================
# Script:           qc.py
# Purpose:          I/O functionality for projects and cohorts
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-14
# ==============================================================================

import logging
import gc

import pandas as pd
import anndata as ad
import numpy as np

from pathlib import Path

from methyltrain.project.layout import ProjectLayout
from methyltrain.cohort.layout import CohortLayout

logger = logging.getLogger(__name__)

# =====| Public API |===========================================================

# ~~~~~| Loading |~~~~~
def load_raw_project(metadata: pd.DataFrame, 
                     config: dict, 
                     layout: ProjectLayout) -> ad.AnnData:
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

    Memory usage is optimized for a ~8 GB constraint.
    
    Parameters
    ----------
    metadata : pd.DataFrame
        Project metadata queried from the GDC API.
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

    logger.info("~~~~~| Attempting to Load Raw Project |~~~~~")

    # Verify all directories exist and are accessible
    layout.validate()
    _validate_directories(layout)

    # Build the raw beta value matrix and align its metadata
    cpg_matrix, probe_ids = _build_beta_matrix(layout)
    logger.info("Successfully loaded the raw beta value files.")
    metadata = metadata.sort_values(by = 'file_name').reset_index(drop = True)
    
    # Initialize the CpG matrix as an AnnData object with aligned metadata
    X = cpg_matrix.T

    adata = ad.AnnData(
        X = X,
        obs = metadata,
        var = pd.DataFrame(index = probe_ids.astype(str))
    )

    logger.info("Successfully constructed the AnnData object.")

    del cpg_matrix
    gc.collect()

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
            "reference_genome": config['project'].get("reference_genome"),
        },
        "pipeline": {
            "state": "raw",
            "steps": [],
        },
    }

    logger.info("~~~~~| Successfully Loaded Raw Project |~~~~~")

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


def _build_beta_matrix(layout: ProjectLayout):
    """
    Loads all beta values in parallel, concatenating upon the index to build a 
    matrix of raw balues.
    """

    # Load all beta values in parallel as a list of Pandas DataFrames
    files = sorted(layout.raw_dir.glob("*parquet"))

    logger.info(f"Detected {len(files)} raw .parquet files.")

    # Build a global CpG index
    probe_set = set()
    for f in files:
        df = pd.read_parquet(f, columns = ['probe_id'])
        probe_set.update(df['probe_id'].values)
        del df

    global_index = pd.Index(sorted(probe_set))

    # Initialize a look-up table for CpG to row position
    n_rows, n_cols = len(global_index), len(files)
    logger.info(f"Identified {n_rows} probes and {n_cols} samples.")

    probe_to_row = pd.Series(
        np.arange(n_rows, dtype = np.int32),
        index = global_index
    )

    logger.info("Successfully built the CpG-to-row look-up mapping.")
    
    # Preallocate the final matrix
    matrix = np.empty((n_rows, n_cols), dtype = np.float32)
    matrix[:] = np.nan

    # Fill the matrix column-by-column
    for j, f in enumerate(files):
        df = pd.read_parquet(f, columns = ['probe_id', 'beta_value'])

        probes = df['probe_id'].to_numpy(copy = False)
        values = df['beta_value'].to_numpy(dtype = np.float32, copy = False)

        # Vectorized index mapping
        row_idx = probe_to_row.loc[probes].to_numpy()
        matrix[row_idx, j] = values

        del df
        
    return matrix, global_index

# [END]
