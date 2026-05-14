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

from typing import Dict

from methyltrain.audit.audit_store import AuditStore
from methyltrain.project.layout import ProjectLayout
from methyltrain.old_constants import GDC_QUERY_URL, MAX_RETRIES
from methyltrain.utils.utils import verify_gdc_client

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


# =====| Internal Helpers |=====================================================

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

# [END]