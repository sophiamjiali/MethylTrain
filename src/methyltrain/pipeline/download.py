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
import datetime

import pandas as pd

from typing import Dict

from ..fs.layout import ProjectLayout
from ..constants import GDC_QUERY_URL, MAX_RETRIES
from ..utils.utils import (
    verify_gdc_client,
    extract_project_id,
    extract_sample_type,
    extract_submitter_id
)

# =====| Build Manifest |=======================================================

def build_manifest(config: Dict) -> pd.DataFrame:
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

    # Initialize query filters based on user configurations and defaults
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id",
                                    "value": [config['project_id']]}},
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
        "fields": ",".join(['cases.project.project_id', 
                            'file_id', 'file_name', 'platform',
                            'data_category', 'experimental_strategy', 
                            'data_type', 'cases.samples.sample_type']),
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
        case["project"]["project_id"] == config["project_id"] for case in cases)
    ).all()

    # Clean the manifest for unused fields for GDC query
    manifest = (df[['file_id', 'file_name']]
                .rename(columns = {'file_id': 'id', 'file_name': 'filename'}))

    return manifest


# =====| Download Methylation |=================================================

def download_methylation(manifest: pd.DataFrame, 
                         config: Dict,
                         layout: ProjectLayout,
                         verbose = False) -> pd.DataFrame:
    """
    Downloads DNA methylation files from the GDC using a prevalidated manifest. 
    Includes additional safety and reproducibility checks, ensuring the GDC 
    Data Transfer Tool (`gdc-client`) is installed and available. Uses batch manifest-based downloading with resume support. Asserts file 
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
    config : dict
        Configuration dictionary controlling workflow steps.
    layout : ProjectLayout
        Object representing a project dataset directory layout.

    Returns
    -------
    pd.DataFrame
        Status log recording per-file download status ('success' or 'failed').
    """

    # Verify the `gdc-client` is properly installed on the user's device
    gdc_client_path = config.get('gdc_client', '')
    verify_gdc_client(gdc_client_path)

    # Pre-filter the temporary manifest for already downloaded files
    missing_files = []
    for _, row in manifest.iterrows():
        filepath = layout.raw_dir / row['filename']
        if not filepath.exists():
            missing_files.append(row)

    remaining_files = pd.DataFrame(missing_files)
    tmp_manifest = layout.raw_dir/f"{layout.project_name}_temp_manifest.txt"

    if verbose: print(f"=====| Beginning to download {layout.project_name} DNA "
                      f"Methylation Data |=====")

    attempt = 0
    while not remaining_files.empty and attempt < MAX_RETRIES:
        attempt += 1

        # Save a temporary manifest
        remaining_files[['id', 'filename']].to_csv(tmp_manifest, sep = '\t',
                                                   index = False)

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
            print(f"Attempt {attempt} failed, retrying remaining files...")
            time.sleep(2 ** attempt)
        
        # Check existence per file
        still_remaining = []
        for _, row in remaining_files.iterrows():
            filepath = layout.raw_dir / row['filename']

            if not filepath.exists():
                still_remaining.append(row)

        # Prepare new manifest with only failed files
        if still_remaining:
            remaining_files = pd.DataFrame(still_remaining)
            remaining_files[['id', 'filename']].to_csv(tmp_manifest, 
                                                       sep = '\t', 
                                                       index = False)
            time.sleep(2 ** attempt)

        else:
            break  # all files downloaded

    # Delete the temporary manifest
    if tmp_manifest.exists(): tmp_manifest.unlink()

    # Log any files that failed after MAX_RETRIES
    status_log = []
    for _, row in manifest.iterrows():
        filepath = layout.raw_dir / row['id'] / row['filename']

        status_log.append({
            'id': row['id'],
            'filename': row['filename'],
            'status': 'success' if filepath.exists() else 'failed',
            'attempts': attempt,
            'timestamp': datetime.datetime.now(datetime.timezone.utc)
        })

    if verbose: print(f"=====| Finished downloading {layout.project_name} DNA "
                      f"Methylation Data |=====")

    return pd.DataFrame(status_log)


# =====| Build Metadata |=======================================================

def build_metadata(audit_table: pd.DataFrame, config: Dict,) -> pd.DataFrame:
    """
    Query GDC for metadata corresponding to successfully downloaded files in the audit table. Nested fields are flattened to a single level.

    Parameters
    ----------
    audit_table : pd.DataFrame
        Running audit table maintaining the status of all files.
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    pd.DataFrame
        Verbose metadata with one row per attempted file. Includes all 
        requested GDC fields as defined in the user-provided configurations, 
        along with a status column indicating query success or failure.
    """

    # Fetch files that were successfully downloaded
    file_ids = audit_table.loc[audit_table['downloaded'] == 1].index.tolist()
    if not file_ids: return pd.DataFrame()

    # Prepare the GDC API request
    filters = {'op': 'in','content': {'field': 'file_id', 'value': file_ids}}
    metadata_fields = config.get('download', {}).get('metadata', [])

    params = {
        'filters': json.dumps(filters),
        'fields': ','.join(metadata_fields),
        'format': 'JSON',
        'size': len(file_ids)
    }

    # Return an empty dataframe if API request fails
    try:
        response = requests.get(GDC_QUERY_URL, params = params)
        response.raise_for_status()
        hits = response.json()['data']['hits']

    except Exception:
        return pd.DataFrame({
            'file_id': file_ids,
            'status': ['failed'] * len(file_ids)
        })
    
    metadata = pd.DataFrame(hits)
    metadata = metadata.drop(columns = ['id'])

    # Extract nested metadata values
    if 'cases' in metadata.columns:
        metadata['project_id'] = metadata['cases'].apply(extract_project_id)
        metadata['submitter_id'] = metadata['cases'].apply(extract_submitter_id)
        metadata['sample_type'] = metadata['cases'].apply(extract_sample_type)

    metadata = metadata.drop(columns = ['cases'])

    # Update fetched status based on GDC response failures
    metadata['status'] = 'success'

    fetched_ids = set(metadata['file_id'])
    missing_ids = set(file_ids) - fetched_ids
    if missing_ids:
        missing = pd.DataFrame({
            'file_id': list(missing_ids),
            'status': ['failed'] * len(missing_ids)
        })
        metadata = pd.concat([metadata, missing], ignore_index = True)

    return metadata