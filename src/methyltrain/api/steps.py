# ==============================================================================
# Script:           steps.py
# Purpose:          Workflow steps exposed to the user
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-08
# ==============================================================================

import shutil

import pandas as pd
import anndata as ad
import numpy as np

from concurrent.futures import ThreadPoolExecutor
from sklearn.model_selection import train_test_split
from typing import Dict, List, Tuple
from pathlib import Path

from ..pipeline.download import (
    build_manifest, 
    build_metadata,
    download_methylation
)
from ..pipeline.preprocess import (
    normalize,
    filter_variance,
    convert_to_mval,
    impute, 
    batch_correction
)
from ..pipeline.audit import initialize_audit_table
from ..pipeline.quality_control import sample_qc, probe_qc
from ..pipeline.clean import clean_metadata
from ..pipeline.aggregate import cohort_aggregation, gene_aggregation
from ..utils.utils import load_sample, load_annotation
from ..fs.layout import ProjectLayout, CohortLayout


# =====| Workflow |=============================================================

def download(config: Dict, 
             layout: ProjectLayout) -> pd.DataFrame:
    """
    Downloads DNA methylation data as beta values from the TCGA GDC based on the project specified in the configurations object. An audit table is built to standardize all attempted files, download status, and metadata status.

    The manifest is created using a GDC API query and resolved at the sample 
    level. Metadata fetching is attempted only for files successfully downloaded.

    Parameters
    ----------
    config : dict
        Configuration dictionary controlling workflow steps.
    layout : ProjectLayout
        Object representing a project dataset directory layout.

    Returns
    -------
    pd.DataFrame
        Audit table containing all atttempted files, download status, and 
        metadata status. Initializes downstream preprocessing and quality 
        control flags. All binary flags include:
            - `downloaded`: 1 if DNA methylation file downloaded, else 0
            - `download_attempts`: number of times a file was attempted 
            - `download_timestamp`: timestamp of successfull download
            - `download_status`: text description for granular logging
            - `metadata_fetched`: 1 if metadata successfully retrieved, else 0
            - `metadata_status`: text description for granular logging
            - `parquet_path`: path to the raw .parquet file (originally .txt)
            - `qc_pass`: 1 if the sample passed QC, else 0
            - `notes`: human-readable notes for verbose documentation
    """

    layout.validate()

    # Build a manifest of available files and attempt to download them
    manifest = build_manifest(config)
    status_log = download_methylation(manifest, config, layout)

    # Build an audit table to hold file status flags to query for metadata
    audit_table = initialize_audit_table(manifest, status_log)
    metadata = build_metadata(audit_table, config)

    # Update the audit table with metadata attempts (rare)
    audit_table = audit_table.merge(metadata[['status']], how = 'left',
                                    left_index = True, right_index = True)
    audit_table['metadata_fetched'] = (audit_table['status']
                                       .eq('success').astype(int))
    audit_table['metadata_status'] = audit_table['status']
    audit_table = audit_table.drop(columns = 'status')

    # Clean the metadata of verbose output to the standard format
    metadata = metadata.loc[metadata['status'] == 'success']

    # Save manifest, status log, metadata, and audit table
    manifest.to_csv(layout.manifest, sep = '\t', header=True, index=False)
    status_log.to_csv(layout.status_log, sep = '\t', header=True, index=False)
    metadata.to_csv(layout.metadata, sep = '\t', header=True, index=False)
    audit_table.to_csv(layout.audit_table, sep = '\t', header=True, index=True)

    return audit_table


def clean_data(audit_table: pd.DataFrame, 
               layout: ProjectLayout) -> pd.DataFrame:
    """
    Cleans raw TCGA DNA methylation beta value .txt files by converting them 
    to .parquet, flattening directory structure and removing accessory files 
    and raw files upon success.

    Updates the audit table with paths to the generated parquets. Renames .
    parquet files to each file's ID (UUID).

    Parameters
    ----------
    audit_table : pd.DataFrame
        Running audit table maintaining the status of all files.
    layout : ProjectLayout
        Object representing a project dataset directory layout.

    Returns
    -------
    pd.DataFrame
        Auit table with the `parquet_path` column filled.
    """

    # Verify the raw data directory exists
    layout.validate()

    # Query for raw files that were successfully downloaded
    downloaded = audit_table.loc[audit_table['downloaded'] == 1].copy()
    for idx, row in downloaded.iterrows():
        file_id = str(idx)
        file_name = row['filename']

        if pd.notna(row['parquet_path']): continue

        # Name the .parquet file with its UUID, not its old file name
        txt_path = layout.raw_dir / file_id / file_name
        parquet_path = layout.raw_dir / f"{file_id}.parquet"

        if not txt_path.exists():
            raise FileNotFoundError(f"Missing raw file: {txt_path}")
        
        # Read the beta values (TCGA standard: probe_id, value)
        txt = pd.read_csv(txt_path, sep = '\t', header = 0, dtype={0: str})
        txt.columns = ['probe_id', 'beta_value']
        txt['beta_value'] = pd.to_numeric(txt['beta_value'], 
                                            errors = "coerce")
        
        txt.to_parquet(parquet_path, index = False)

        # Remove raw files and directory
        txt_path.unlink(missing_ok = True)
        if txt_path.parent.exists(): shutil.rmtree(txt_path.parent)

        audit_table.loc[file_id, "parquet_path"] = str(parquet_path)

    return audit_table


def quality_control(adata: ad.AnnData, 
                    audit_table: pd.DataFrame,
                    config: Dict, 
                    layout: ProjectLayout) -> Tuple[ad.AnnData, pd.DataFrame]:
    """
    Performes probe and/or sample quality control upon DNA methylation values 
    presented as a CpG matrix AnnData object. Returns the quality-controlled 
    CpG matrix AnnData object with updated metadata suitable for downstream 
    preprocessing and analysis, and the updated audit table with QC status.

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
    audit_table : pd.DataFrame
        Running audit table maintaining the status of all files.
    config : dict
        Configuration dictionary controlling workflow steps.
    layout : ProjectLayout
        Object representing a project dataset directory layout.

    Returns
    -------
    Tuple(ad.AnnData, pd.DataFrame)
        AnnData object representing a project's quality-controlled DNA 
        methylation data at the CpG matrix level with updated metadata. Audit 
        table with updated QC flags.

    """

    # Load the appropriate annotation provided by the package, else raises error
    annotation = load_annotation(
        platform = adata.uns['platform'], 
        reference_genome = adata.uns['reference_genome']
    )

    # Perform each quality-control step if toggled by the user-configurations
    toggles = config.get('toggles', {})
    
    if toggles.get('sample_qc', True):

        before_qc = set(adata.obs['file_id'])
        adata = sample_qc(adata, config)
        after_qc = set(adata.obs['file_id'])

        fail_qc = before_qc - after_qc

        # Update the audit table for files present in the adata (not prev. fail)
        audit_table.loc[audit_table.index.isin(fail_qc), 'qc_pass'] = 0
        audit_table.loc[audit_table.index.isin(after_qc), 'qc_pass'] = 1
    
    if toggles.get('probe_qc', True):
        adata = probe_qc(adata, annotation, config)

    # Clean the metadata and manifest; if no QC performed, same as raw
    clean_metadata(adata, layout)

    adata.uns['state'] = 'processed'

    return (adata, audit_table)


def preprocess(adata: ad.AnnData, config: Dict) -> ad.AnnData:
    """
    Preprocess DNA methylation beta values of a project to a gene-level matrix. 
    Returns a samples x genes AnnData object with aligned metadata suitable for 
    downstream analysis.

    Steps performed are toggled and configured in the user-provided 
    configurations, including the following options in order:
    
    1. Normalization / Harmonization
    2. Filter low variance / extreme values
    3. M-value conversion
    4. Imputation of missing values

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
    if config.get('toggles', {}).get('normalize', True):
        adata = normalize(adata)

    if config.get('toggles', {}).get('sample_qc', True):
        adata = impute(adata)

    if config.get('toggles', {}).get('filter_variance', True):
        adata = filter_variance(adata, config)

    if config.get('toggles', {}).get('convert_to_mval', True):
        adata = convert_to_mval(adata)

    if config.get('toggles', {}).get('impute', True):
        adata = impute(adata)

    adata.uns['state'] = 'processed'

    return adata


def aggregate_cohort(adatas: List[ad.AnnData],
                     layout: CohortLayout) -> ad.AnnData:
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
    Performs ComBat batch correction across multiple projects aggregated into a CpG matrix cohort. Automatically checks that the input data is M-values (requirement for parametric batch correction). 

    Parameters
    ----------
    adata : ad.AnnData
        Cohort AnnData object with DNA methylation M-values.
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    adata : ad.AnnData
        AnnData object representing the cohort's batch-corrected DNA 
        methylation data at the CpG probe matrix level with updated metadata.

    Raises
    ------
    ValueError
        If NaNs exist in the AnnData object (batch correction requires none), 
        if data is presented as beta values instead of M-values, or if batch 
        and covariate columns don't exist.
    """
    
    if config.get('toggles', {}).get('batch_correction', True):
        adata = batch_correction(adata, config)

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
    annotation = load_annotation(
        platform = adata.uns['platform'], 
        reference_genome = adata.uns['reference_genome']
    )
    
    adata = gene_aggregation(adata, annotation, config)
    adata.uns['state'] = 'processed'

    return adata


def winsorize(adata: ad.AnnData, config: Dict) -> ad.AnnData:
    """
    Clips data of an AnnData object either at the CpG probe-level or gene-level.

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
        The AnnData object with its beta values clipped per user configurations.
    """

    if config.get('toggles', {}).get('winsorize', True):
        clip_values = config.get('preprocessing', {}).get('clip_values', [0, 0])
        adata.X = np.clip(np.array(adata.X), clip_values[0], clip_values[1])

    adata.uns['state'] = 'processed'
    
    return adata


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
    files = [f for f in raw_dir.iterdir() if f.suffix == ".parquet"]

    with ThreadPoolExecutor() as ex:
        sample_beta_values = list(ex.map(load_sample, files))

    # Concatenate on the index to build a matrix
    cpg_matrix = pd.concat(sample_beta_values, axis = 1, join = "outer")
    cpg_matrix = cpg_matrix.sort_index()

    # Sort by file name as .parquets are loaded in that order
    metadata = pd.read_csv(layout.metadata, sep = '\t')
    metadata = metadata.set_index('file_name', drop = True).sort_index()

    # Initialize the CpG matrix as an AnnData object with aligned metadata
    adata = ad.AnnData(
        X = cpg_matrix.T.values,
        obs = metadata,
        var = pd.DataFrame(index = cpg_matrix.index)
    )

    # Initialize global metadata for the project
    adata.uns['project_id'] = layout.project_name
    adata.uns['platform'] = config.get('download', {}).get('platform', '')
    adata.uns['reference_genome'] = config.get('reference_genome', '')
    adata.uns['level'] = "project"
    adata.uns['data_type'] = "cpg_matrix"
    adata.uns['conversion'] = 'beta_value'
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
    adata.write_h5ad(layout.project_adata, compression = "gzip")


# =====| Metadata I/O |=========================================================

def load_audit_table(layout: ProjectLayout) -> pd.DataFrame:
    # Loads the audit table with file_id as index

    layout.validate()
    audit_table = pd.read_csv(layout.audit_table, sep = '\t', index_col = 0)
    return audit_table


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


def load_status_log(layout: ProjectLayout) -> pd.DataFrame:
    # Loads the status log with file_id as index
    layout.validate()
    status_log = pd.read_csv(layout.status_log, sep = '\t', index_col = 0)
    return status_log

# [END]