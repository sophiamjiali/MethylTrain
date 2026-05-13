# ==============================================================================
# Script:           download.py
# Purpose:          Internal downloading functions and implementation logic
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-20
# ==============================================================================

import sys
import json
import requests
import subprocess
import time
import logging

import pandas as pd

from typing import Dict, List

from methyltrain.audit.audit_store import AuditStore
from methyltrain.project.layout import ProjectLayout
from methyltrain.constants import (
    GDC_QUERY_URL, 
    MAX_RETRIES, 
    GDC_QUERY_BATCH_URL
)
from methyltrain.utils.utils import (
    verify_gdc_client,
    extract_project_id,
    extract_sample_type,
    extract_submitter_id,
    extract_aliquot_id,
    extract_batch_id
)

logger = logging.getLogger(__name__)

# =====| Public API |===========================================================

def download_methylation(config: Dict, 
                         layout: ProjectLayout, 
                         audit: AuditStore) -> pd.DataFrame:
    """
    Downloads DNA methylation data of a TCGA project as beta values from the 
    TCGA GDC based on the project specified in the provided configuration 
    object. An audit store is created to report attempted files, download 
    statIt us, and metadata status.

    A manifest is created using the GDC API and resolved at the sample level. 
    Metadata fetching is attempted only for files who's DNA Methylation data 
    was successfully downloaded.

    Parameters
    ----------
    config : dict
        Configuration dictionary controlling workflow steps.
    audit: AuditStore
        Metadata for downloading and preprocessing fidelity of the project.
    layout : ProjectLayout
        Object representing a project dataset directory layout.

    Returns
    -------
    pd.DataFrame
        The manifest returned by the GDC API used to dowload the project's DNA 
        Methylation Data.
    """

    layout.validate()

    logger.info("=====| Attempting Project Methylation Download |=====")

    # Query the GDC API for the project manifest and initialize the audit store
    manifest = _build_manifest(config)
    logger.info("Successfully queried for the manifest.")

    # Initialize the AuditStore with queried file_id
    audit.initialize(manifest.index.tolist())
    logger.info("Successfully initialized the audit store.")

    # Download the methylation data and update the audit store
    _download_methylation(manifest, audit, config, layout)
    logger.info("Successfully downloaded methylation data.")

    logger.info("~~~~~| Successfully Downloaded Methylation Data |~~~~~\n")
    return manifest


def prepare_metadata(config: Dict, audit: AuditStore) -> pd.DataFrame:
    """
    Downloads biospecimen and metadata of a TCGA project as a CSV file from the 
    TCGA GDC based on the project specified in the provided configuration 
    object. Biospecimen data is merged with the metadata to return a 
    consolidated dataframe.

    Parameters
    ----------
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    pd.DataFrame
        Consolidated metadata and biospecimen data for the project.
    """

    logger.info("~~~~~| Attempting Project Metadata Download |~~~~~")

    # Query the GDC API for the project metadata and biospecimen data
    metadata = _build_metadata(config, audit)
    logger.info("Successfully queried for the metadata.")
    biospecimen = _build_biospecimen(metadata, audit)
    logger.info("Successfully queried for the biospecimen data.")

    # Merge biospecimen data into metadata to consolidate
    metadata = metadata.join(biospecimen[['barcode']], how = "left")
    metadata['batch_id'] = metadata['barcode'].apply(extract_batch_id)
    logger.info("Successfully consolidated metadata and biospecimen data.")

    logger.info("~~~~~| Successfully Downloaded Project Metadata |~~~~~\n")

    return metadata


# =====| Internal Helpers |=====================================================

# ~~~~~| Download Project |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def _build_manifest(config: Dict) -> pd.DataFrame:
    """
    Build a GDC client manifest for a specific project and data type. Queries 
    the GDC API for files matching the user-provided configurations and 
    constructs a manifest compatible with `gdc-client`.

    Filtering is performed at the API level using the following criteria:
    - Case-level fields: Project ID, sample type, and open access
    - File-level fields: Data category, experimental strategy, data type, 
      platform, and reference 

    The user-provided configurations specifies a single platform, reference 
    genome, and sample type, thus the returned manifest is already resolved to 
    the sample-level: each entry corresponds to a unique sample that meets the 
    specified criteria. No additional deduplication by sample is required.

    Assertions are performed on the retrieved metadata to ensure that all files 
    strictly match the requested platform, reference genome, data category, 
    data type, experimental strategy, and sample type.

    Parameters
    ----------
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    pd.DataFrame
        A minimal manifest DataFrame with columns:
        - 'id' : GDC file UUID
        - 'filename' : file name
    """

    dc = config.get('download', {})
    project_id = config.get('project', {}).get('project_id', '')

    # Initialize query filters based on user configurations and defaults
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id",
                                    "value": [project_id]}},
            {"op": "in", "content": {"field": "files.data_category",
                                    "value": [dc['data_category']]}},
            {"op": "in", "content": {"field": "files.experimental_strategy",
                                    "value": [dc['experimental_strategy']]}},
            {"op": "in", "content": {"field": "files.data_type",
                                    "value": [dc['data_type']]}},
            {"op": "in", "content": {"field": "files.platform",
                                    "value": [dc['platform']]}},
            {"op": "in", "content": {"field": "cases.samples.sample_type",
                                    "value": [dc['sample_type']]}},
            {"op": "in", "content": {"field": "files.access",
                                    "value": ["open"]}}
        ]
    }

    # Fetch additional temporary parameters to assert manifest correctness
    params = {
        "filters": json.dumps(filters),
    "fields": ",".join(['file_id', 'file_name', 'cases.project.project_id',
                        'data_category', 'experimental_strategy',
                        'data_type', 'platform', 'cases.samples.sample_type']),
        "size": 10000
    } 

    # Query the GDC API for methylation files with the filters
    response = requests.get(GDC_QUERY_URL, params = params)
    response.raise_for_status()
    hits = response.json()['data']['hits']

    df = pd.DataFrame(hits)

    # Verify the integrity of the manifest and requested data
    assert len(hits) < params["size"], "Query may be truncated."
    assert not df.empty, "No files returned from the GDC query."

    assert df['data_category'].eq(dc['data_category']).all()
    assert df['data_type'].eq(dc['data_type']).all()
    assert df['platform'].eq(dc['platform']).all()
    assert df['experimental_strategy'].eq(dc['experimental_strategy']).all()
    assert df['cases'].apply(lambda cases: any(
        c['sample_type'] == dc['sample_type'] for c in cases[0]['samples'])
    ).all()
    assert df["cases"].apply(lambda cases: any(
        case["project"]["project_id"] == project_id for case in cases)
    ).all()

    # Clean the manifest for unused fields for GDC query
    manifest = df.drop(columns = ['id'])
    manifest = manifest[['file_id', 'file_name']]
    manifest = manifest.set_index('file_id', drop = True)

    return manifest


def _download_methylation(manifest: pd.DataFrame, 
                          audit: AuditStore,
                          config: Dict,
                          layout: ProjectLayout) -> None:
    """
    Downloads DNA methylation files from the GDC using a prevalidated manifest. 
    Includes additional safety and reproducibility checks, ensuring the GDC 
    Data Transfer Tool (`gdc-client`) is installed and available. Uses batch 
    manifest-based downloading with resume support. Asserts file 
    existence after download and logs per-file status.

    Logs per-file download status and timestamp. Failures in download will 
    raise an exception (md5 are typically not available).

    Note that this function does not bundle the `gdc-client` API in the 
    package; users must install it from the official GDC site.

    Parameters
    ----------
    manifest : pd.DataFrame
        Prevalidated manifest that is filtered for the desired platform, 
        reference genome, sample type, and data category. Each row corresponds 
        to a unique sample.
    audit : AuditStore
        SQL Audit store for logging download status and fidelity.
    config : dict
        Configuration dictionary controlling workflow steps.
    layout : ProjectLayout
        Object representing a project dataset directory layout.
    """

    # Verify the `gdc-client` is properly installed on the user's device
    gdc_client_path = config.get('paths', {}).get('gdc_client', '')
    verify_gdc_client(gdc_client_path)

    # Pre-filter the temporary manifest for already downloaded files
    missing_files = []
    for idx, row in manifest.iterrows():
        
        # Check for both raw and cleaned file types
        raw_path = layout.raw_dir / str(idx)
        clean_path = layout.raw_dir / (str(idx) + '.parquet')
        if not raw_path.exists() and not clean_path.exists():
            missing_files.append(row)

    remaining_files = pd.DataFrame(missing_files).reset_index()
    remaining_files = remaining_files.rename(columns = {'index': 'id', 
                                                        'file_name': 'filename'})
    tmp_manifest = layout.raw_dir / "temp_manifest.txt"

    attempt = 0
    while not remaining_files.empty and attempt < MAX_RETRIES:
        attempt += 1

        # Save a temporary manifest
        remaining_files.to_csv(tmp_manifest, sep = '\t', index = False)

        try:
            # Run batch download for remaining files
            cmd = [gdc_client_path,
                    "download", 
                    "-m", str(tmp_manifest), 
                    "-d", str(layout.raw_dir)]

            # Silence standard output, but keep errors for debugging
            subprocess.run(cmd, 
                            check = True, 
                            stdout = subprocess.DEVNULL, 
                            stderr = sys.stderr)
            
        except subprocess.CalledProcessError:
            # Log failure for this batch attempt, will retry failed files
            time.sleep(2 ** attempt)
            continue
        
        # Check existence per file
        still_remaining = [
            row for _, row in remaining_files.iterrows()
            if not (layout.raw_dir / row['id']).exists()
        ]

        # Prepare new manifest with only failed files
        if still_remaining:
            remaining_files = pd.DataFrame(still_remaining)
            remaining_files.to_csv(tmp_manifest, sep = '\t', index = False)
            time.sleep(2 ** attempt)
        else:
            break  # all files downloaded

    # Delete the temporary manifest
    if tmp_manifest.exists(): tmp_manifest.unlink()

    # Log the download status of all files in the audit store
    for idx, row in manifest.iterrows():
        filepath = layout.raw_dir / str(idx) / row['file_name']
        status = 1 if filepath.exists() else 0
        file_name = row.get('file_name', '')

        audit.set_download_status(str(idx), file_name, status)

    return


# ~~~~~| Prepare Metadata |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def _build_metadata(config: Dict, 
                    audit: AuditStore,
                    batch_size: int = 20) -> pd.DataFrame:
    """
    Query the GDC API for metadata corresponding to successfully DNA 
    Mehtlyation data downloaded files in the audit table. Nested fields are 
    flattened to a single level.

    Parameters
    ----------
    audit: AuditStore
        Metadata for downloading and preprocessing fidelity of the project.
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    pd.DataFrame
        Verbose metadata with one row per attempted file. Includes all 
        requested GDC fields as defined in the user-provided configurations, 
        along with a status column indicating query success or failure.
    """

    # Fetch the IDs of all files that were successfully downloaded
    file_ids = audit.get_ids_by_download_status(status = 1)
    if not file_ids: return pd.DataFrame()

    # Prepare the GDC API request
    metadata_fields = config.get('metadata', [])
    results: List[dict] = []

    # Query the metadata in batches to avoid API failures
    for i in range(0, len(file_ids), batch_size):
        batch = file_ids[i:i + batch_size]

        filters = {'op': 'in','content': {'field': 'file_id', 'value': batch}}
        params = {
            'filters': json.dumps(filters),
            'fields': ','.join(metadata_fields),
            'format': 'JSON',
            'size': len(batch)
        }

        # Return an empty dataframe if API request fails
        try:
            response = requests.get(GDC_QUERY_URL, params = params)
            response.raise_for_status()
            hits = response.json()['data']['hits']

            # Log the download status to the audit store
            for hit in hits: audit.set_metadata_status(hit['file_id'], 1)

        except Exception:
            for fid in batch: audit.set_metadata_status(fid, 0)
            continue
    
    # Drop the internal query UUID from the table
    metadata = pd.DataFrame(results).set_index('file_id', drop = True)
    metadata = metadata.drop(columns = ['id'], errors = 'ignore')

    # Extract nested metadata values
    if 'cases' in metadata.columns:
        metadata['project_id'] = (metadata['cases'].apply(extract_project_id))
        metadata['submitter_id'] = metadata['cases'].apply(extract_submitter_id)
        metadata['sample_type'] = metadata['cases'].apply(extract_sample_type)
        metadata['aliquot_id'] = metadata['cases'].apply(extract_aliquot_id)
        metadata = metadata.drop(columns = ['cases'])
        
    return metadata


def _build_biospecimen(metadata: pd.DataFrame,
                       audit: AuditStore,
                       batch_size: int = 20) -> pd.DataFrame:
    """
    Query the GDC API for biospecimen metadata corresponding to the 
    successfully downloaded files (based on metadata, and subsequently DNA 
    Methylation data download status)in the audit table. Nested fields are 
    flattened to a single level.

    Parameters
    ----------
    audit: AuditStore
        Metadata for downloading and preprocessing fidelity of the project.
    config : dict
        Configuration dictionary controlling workflow steps.
    metadata : pd.DataFrame
        Metadata for the project as returned by the GDC API.

    Returns
    -------
    pd.DataFrame
        Verbose biospecimen data with one row per attempted file. Includes all 
        requested GDC fields as defined in the user-provided configurations, 
        along with a status column indicating query success or failure.
    """

    # Map file IDs to case IDs through aliquot IDs
    file_ids = audit.get_ids_by_download_status(status = 1)
    if not file_ids: return pd.DataFrame()

    file_meta = metadata.loc[metadata.index.isin(file_ids),
                             ['aliquot_id', 'submitter_id']].dropna()
    
    # Build a lookup table between levels
    file_to_submitter = file_meta['submitter_id'].to_dict()
    file_to_aliquot = file_meta['aliquot_id'].to_dict()

    submitter_ids = file_meta['submitter_id'].dropna().unique().tolist()
    if not submitter_ids: return pd.DataFrame()

    aliquot_map = {}
    
    # Query the GDC client for the biospecimen data
    for i in range(0, len(submitter_ids), batch_size):
        batch = submitter_ids[i:i + batch_size]

        filters = {
            "op": "in",
            "content": {"field": "submitter_id", "value": batch}
        }
    
        params = {
            "filters": json.dumps(filters),
            "expand": "samples.portions.analytes.aliquots",
            "fields": ",".join([
                'submitter_id', 'samples.portions.submitter_id',
                'samples.portions.analytes.aliquots.aliquot_id',
                'samples.portions.analytes.aliquots.submitter_id'
            ]),
            "format": "JSON",
            "size": len(batch)
        }

        # Query the `cases` endpoint; return an empty dataframe if request fails
        try:
            response = requests.get(GDC_QUERY_BATCH_URL, params = params)
            response.raise_for_status()
            hits = response.json()['data']['hits']

            # Extract the nested barcode values
            for case in hits:
                for sample in case.get('samples', []):
                    for portion in sample.get('portions', []):
                        portion_barcode = portion.get('submitter_id')
                        for analyte in portion.get('analytes', []):
                            for aliquot in analyte.get('aliquots', []):

                                # Use the aliquot-level barcode
                                aliquot_barcode = aliquot.get('submitter_id', 
                                                              portion_barcode)
                                aliquot_map[aliquot.get('aliquot_id')] = {
                                    'submitter_id': case['submitter_id'],
                                    'barcode': aliquot_barcode
                                }
                                            
        except Exception:
            continue

    # Map the barcodes back to file ID
    results = []
    for file_id in file_ids:
        aliquot_id = str(file_to_aliquot[file_id] if file_id in
                         file_to_aliquot else None)
        submitter_id = file_to_submitter.get(file_id, None)
        bio = aliquot_map.get(aliquot_id, None)

        results.append({
            'file_id': file_id,
            'aliquot_id': aliquot_id,
            'submitter_id': submitter_id,
            'barcode': bio['barcode'] if bio else None,
        })

        # Log the download success; initialized as failure, boolean zero
        status = int(bio is not None)
        audit.set_biospecimen_status(file_id, status)

    # Keep only rows corresponding to aliquots in adata
    biospecimen = pd.DataFrame(results).set_index('file_id', drop = True)
    
    return biospecimen

# [END]