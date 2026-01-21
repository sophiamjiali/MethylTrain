# ==============================================================================
# Script:           steps.py
# Purpose:          Workflow steps exposed to the user
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-08
# ==============================================================================

import pandas as pd
import anndata as ad
import numpy as np

from concurrent.futures import ThreadPoolExecutor
from sklearn.model_selection import train_test_split
from typing import Dict, List
from pathlib import Path

from ..pipeline.download import (
    build_manifest, 
    build_metadata,
    download_methylation
)
from ..pipeline.quality_control import sample_qc, probe_qc
from ..pipeline.clean import clean_metadata, clean_manifest
from ..pipeline.preprocess import impute, batch_correction
from ..pipeline.aggregate import cohort_aggregation, gene_aggregation
from ..fs.layout import ProjectLayout, CohortLayout
from ..utils.utils import load_sample, load_annotation


# =====| Workflow |=============================================================

def download(config: Dict, layout: ProjectLayout) -> pd.DataFrame:
    """
    Downloads DNA methylation data as beta values from the TCGA GDC based on the project specified in the configurations object. A manifest of file IDs is created to use the `gdc-client` API, preserving metadata.

    A retry loop will re-run the download command if it fails. The manifest is temporarily split into chunks to break up long API requests (main source of issues when downloading large amounts of data).
    """

    # ==================
    # under construction
    # ==================

    layout.validate()

    # Build a manifest of available files from TCGA GDC
    manifest = build_manifest(config)
    status_log = download_methylation(manifest, config, layout)
    metadata = build_metadata(config)

    # Save manifest, status log, and metadata
    manifest.to_csv(layout.raw_manifest, sep = '\t', header=True, index=False)
    status_log.to_csv(layout.status_log, sep = '\t', header=True, index=False)
    metadata.to_csv(layout.raw_metadata, sep = '\t', header=True, index=False)

    return status_log


def clean_data(layout: ProjectLayout) -> None:

    # set sample ID to config.downloading.sample_id instead of filename

    # ==================
    # under construction
    # ==================

    # - if different aliquots map to the same sample/case, kept a representative
    #   profile by selecting the most recently generated file

    # use convert_to_parquet in clean.py


    return


def quality_control(adata: ad.AnnData, 
                    config: Dict, 
                    layout: ProjectLayout) -> ad.AnnData:
    """
    Performes probe and/or sample quality control upon DNA methylation values 
    presented as a CpG matrix AnnData object. Returns the quality-controlled 
    CpG matrix AnnData object with updated metadata suitable for downstream 
    preprocessing and analysis.

    Note that quality control is intended to be performed upon raw project data (as loaded by `load_raw_project()`) and followed by preprocessing (as per `preprocess()`).

    Steps performed are toggled and configured in the user-provided 
    configurations, including the following options in order:

    1. Sample-level quality control
        a. Remove high missingness (above the provided threshold)
        b. Remove outliers from distribution (above the provided number of SD)
    2. Probe-level quality control 
        a. Remove cross-reactive probes
        b. Remove SNP-associated probes
        c. Remove multi-mapped probes
        d. Remove sex chromosome probes
        e. Remove high missingness (above the provided threshold)

    Generates processed metadata and manifest files for the project after 
    quality control samples may be filtered from the dataset. If no quality 
    control was configured to occur in the user-configurations but this 
    function was still called, the processed files will include identical 
    information to the raw files.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object representing a project's raw DNA methylation data at the 
        CpG matrix level.
    config : dict
        Configuration dictionary controlling workflow steps.
    layout : ProjectLayout
        Object representing a project dataset directory layout.

    Returns
    -------
    ad.AnnData
        AnnData object representing a project's quality-controlled DNA 
        methylation data at the CpG matrix level with updated metadata.

    """

    # Load the appropriate annotation provided by the package, else raises error
    annotation = load_annotation(
        platform = adata.uns['platform'], 
        reference_genome = adata.uns['reference_genome']
    )

    # Perform each quality-control step if toggled by the user-configurations
    toggles = config.get('toggles', {})
    
    if toggles.get('sample_qc', True):
        adata = sample_qc(adata, config)
    
    if toggles.get('probe_qc', True):
        adata = probe_qc(adata, annotation, config)

    # Clean the metadata and manifest; if no QC performed, same as raw
    clean_metadata(adata, layout)
    clean_manifest(adata, layout)

    adata.uns['state'] = 'processed'

    return adata


def preprocess(adata: ad.AnnData, config: Dict) -> ad.AnnData:
    """
    Preprocess DNA methylation beta values of a project to a gene-level matrix. 
    Returns a samples x genes AnnData object with aligned metadata suitable for 
    downstream analysis.

    Steps performed are toggled and configured in the user-provided 
    configurations, including the following options in order:
    
    1. Imputation of missing values

    Note that batch effect correction is optionally performed after multiple 
    projects have been aggregated into a cohort.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object representing a project's quality-controlled (or raw if 
        QC was not performed) DNA methylation data at the CpG matrix level.
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    ad.AnnData
        AnnData object representing a project's preprocessed DNA methylation 
        data at the CpG probe matrix level with updated metadata.
    """

    # Perform each preprocessing step if toggled by the user-configurations
    if config.get('toggles', {}).get('sample_qc', True):
        adata = impute(adata)

    adata.uns['state'] = 'processed'

    return adata


def aggregate_cohort(layout: CohortLayout, 
                     adatas: List[ad.AnnData]) -> ad.AnnData:
    """
    Aggregates multiple project AnnData objects at the CpG probe x sample 
    matrix level into a single cohort AnnData object. Takes the common set of 
    CpG probes from all projects.

    Resolves dataset-level metadata such that `.uns` is a dictionary with 
    project names as keys, and project-level metadata as their values. 
    Cohort-level metadata is stored in a flat structure.

    Parameters
    ----------
    layout : CohortLayout
        Layout object detailing a cohort.
    adatas : List of ad.AnnData
        List of project AnnData objects at the CpG probe x sample level, each representing a single project. Each object must have:
        - .uns['project_id'] : unique project name
        - .uns['level'] : 'project' or 'cohort'
        - .uns['data_type'] : 'cpg_matrix' or 'gene_matrix'
        - .uns['state'] : 'raw' or 'processed'
        - .uns['preprocessing_steps'] : QC and processing steps performed

    Returns
    -------
    ad.AnnData
        Aggregated cohort AnnData object at the CpG probe x sample level with:
        - cohort.obs['project_id'] : indicating original project per sample
        - cohort.uns['projects'] : dict of per-project .uns metadata
        - cohort.uns['cohort_id'] : unique cohort name
        - cohort.uns['level'] : 'cohort'
        - cohort.uns['data_type'] : 'cpg_matrix' or 'gene_matrix'
        - cohort.uns['state'] : 'raw' or 'processed'
        - cohort.uns['preprocessing_steps'] : processing steps performed
    """

    layout.validate()

    cohort_adata = cohort_aggregation(adatas)

    # Populate cohort-level metadata
    cohort_adata.uns['cohort_id'] = layout.cohort_name
    cohort_adata.uns['level'] = 'cohort'
    cohort_adata.uns['state'] = 'raw'
    cohort_adata.uns['preprocessing_steps'] = []

    return cohort_adata


def cohort_batch_correction(adata: ad.AnnData, config: Dict) -> ad.AnnData:
    """
    Corrects batch effects across multiple projects aggregated into a CpG 
    matrix cohort.
    """

    # Optionally perform batch correction on the cohort if toggled
    if config.get('toggles', {}).get('batch_correction', True):
        adata = batch_correction(adata)

    adata.uns['state'] = 'processed'

    return adata


def aggregate_genes(adata: ad.AnnData, config: Dict) -> ad.AnnData:
    """
    Aggregates a project or cohort AnnData object at the CpG probe x sample 
    matrix level to the gene-level. Regions are used to select probes that 
    contribute, but the final matrix has only genes as columns.

    Parameters
    ----------
    adata : ad.AnnData
        Probe-level DNA methylation AnnData (probes x samples).
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    adata : ad.AnnData
        Aggregated AnnData object to the gene-level.
    """

    # Load the appropriate annotation object
    annotation = load_annotation(array_type = adata.uns['array_type'], 
                                 genome_build = adata.uns['genome_build'])
    
    if config.get('toggles', {}).get('gene_aggregation', True):
        adata = gene_aggregation(adata, annotation, config)

    adata.uns['state'] = 'processed'

    return adata



def winsorize(adata: ad.AnnData) -> ad.AnnData:
    """
    Clips data of an AnnData object either at the CpG probe-level or gene-level.
    
    """

    return ad.AnnData()


def split(adata: ad.AnnData, config: Dict) -> tuple[ad.AnnData, ad.AnnData, 
                                                    ad.AnnData]:
    """
    Split a cohort AnnData object into stratified train, validation, and test 
    sets based on project using the ratio provided in the configurations.
    
    Parameters
    ----------
    adata : ad.AnnData
        AnnData object representing a project's quality-controlled, preprocessed, and batch effect corrected DNA methylation data at the gene matrix level.
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    tuple[ad.AnnData, ad.AnnData, ad.AnnData]
        The train, validation, and test stratified splits as ad.AnnData objects.
    """

    # Split into train+validation and test
    train_val_idx, test_idx = train_test_split(
        adata.obs_names, 
        test_size = config.get('split', [])[2],
        stratify = adata.obs['project'],
        random_state = config.get('seed', 42),
        shuffle = True
    )

    # Parse train_val_idx to split again
    stratify_array = np.array(
        adata.obs['project'].reindex(list(train_val_idx)), 
        dtype = str
    )

    # Split into train and validation
    train_idx, val_idx = train_test_split(
        train_val_idx,
        test_size = config.get('split', [])[1],
        stratify = stratify_array,
        random_state = config.get('seed', 42),
        shuffle = True
    )

    # Slice the cohort AnnData object according to the defined splits
    train_adata = adata[train_idx].copy()
    val_adata = adata[val_idx].copy()
    test_adata = adata[test_idx].copy()

    # Update each AnnData object's global metadata
    train_adata.uns['split'] = "training"
    val_adata.uns['split'] = "validation"
    test_adata.uns['split'] = "test"

    train_adata.uns['split_percentage'] = config.get('split', [])[0]
    val_adata.uns['split_percentage'] = config.get('split', [])[1]
    test_adata.uns['split_percentage'] = config.get('split', [])[2]

    return train_adata, val_adata, test_adata


# =====| Project I/O |==========================================================

def load_raw_project(config: Dict, layout: ProjectLayout) -> ad.AnnData:
    """
    Load the raw DNA methylation data of a project as an AnnData object from 
    `.parquet` files in the raw data directory. Column metadata is initialized as the sample ID field specified in the user-provided configurations.

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
    project_dir : Path
        Path to a project directory.
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

    # Verify the project raw data directory exists and is not empty
    raw_dir = layout.raw_dir

    if not raw_dir.is_dir():
        raise FileExistsError(f"Project directory was not found: {raw_dir}")
    if not any(raw_dir.iterdir()):
        raise FileExistsError(f"Project directory is empty: {raw_dir}")

    # Load all beta values in parallel as a list of Pandas DataFrames
    files = [f for f in raw_dir.iterdir() if f.is_file()]

    with ThreadPoolExecutor() as ex:
        sample_beta_values = list(ex.map(load_sample, files))

    # Concatenate on the index to build a matrix
    cpg_matrix = pd.concat(sample_beta_values, axis = 1, join = "outer")
    cpg_matrix = cpg_matrix.sort_index()

    # Load the metadata and align the sample IDs
    metadata = pd.read_csv(layout.raw_metadata)
    metadata = metadata.set_index('file_name').sort_index()

    # Initialize the CpG matrix as an AnnData object with aligned metadata
    adata = ad.AnnData(
        X = cpg_matrix.T.values,
        obs = metadata,
        var = pd.DataFrame(index = cpg_matrix.index)
    )

    # Initialize global metadata for the project
    adata.uns['project_id'] = layout.project_name
    adata.uns['array_type'] = config.get('array_type', '')
    adata.uns['genome_build'] = config.get('genome_build', '')
    adata.uns['level'] = "project"
    adata.uns['data_type'] = "cpg_matrix"
    adata.uns['state'] = "raw"
    adata.uns['preprocessing_steps'] = []

    # Default train-val-test split metadata to NA
    adata.uns['split'] = None
    adata.uns['split_percentage'] = None
    
    return adata


def load_processed_project(processed_file: Path | str) -> ad.AnnData:
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
    processed_file = Path(processed_file)

    if not processed_file.exists():
        raise FileNotFoundError(f"Processed file not found: {processed_file}")
    return ad.read_h5ad(processed_file)


def save_project(path: Path) -> None:



    # ==================
    # under construction
    # ==================







    # Saves any anndata object at a checkpoint

    # convert to sparse if not already
    return