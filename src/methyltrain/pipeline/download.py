# ==============================================================================
# Script:           download.py
# Purpose:          Internal downloading functions and implementation logic
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-20
# ==============================================================================

import json
import requests
import subprocess
import time
import datetime

import pandas as pd

from typing import Dict, Union

from ..fs.layout import ProjectLayout
from ..constants import GDC_QUERY_URL, MAX_RETRIES
from ..utils.utils import verify_md5, verify_gdc_client

# =====| Build Manifest |=======================================================

def build_manifest(config: Dict) -> pd.DataFrame:
    """
    Build a GDC client manifest for a specific project and data type. Queries the GDC API for files matching the user-provided configurations and constructs a manifest compatible with `gdc-client`.

    Filtering is performed at the API level using the following criteria:
    - Case-level fields: Project ID, sample type, and open access
    - File-level fields: Data category, experimental strategy, data type, 
      platform, and reference 

    The user-provided configurations specifies a single platform, reference genome, and sample type, thus the returned manifest is already resolved to the sample-level: each entry corresponds to a unique sample that meets the specified criteria. No additional deduplication by sample is required.

    Assertions are performed on the retrieved metadata to ensure that all files strictly match the requested platform, reference genome, data category, data type, experimental strategy, and sample type.

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
        - 'md5' : GDC-provided checksum

    """

    pc = config.get('project', {})
    dc = config.get('download', {})

    # Initialize query filters based on user configurations and defaults
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id",
                                     "value": [pc['project_id']]}},
            {"op": "in", "content": {"field": "files.data_category",
                                     "value": [dc['data_category']]}},
            {"op": "in", "content": {"field": "files.experimental_strategy",
                                     "value": [dc['experimental_strategy']]}},
            {"op": "in", "content": {"field": "files.data_type",
                                     "value": [dc['data_type']]}},
            {"op": "in", "content": {"field": "files.platform",
                                     "value": [dc['platform']]}},
            {"op": "in", "content": {"field": "files.reference_genome",
                                     "value": [dc['reference_genome']]}},
            {"op": "in", "content": {"field": "cases.samples.sample_type",
                                     "value": [dc['sample_type']]}},
            {"op": "in", "content": {"field": "files.access",
                                     "value": ["open"]}}
        ]
    }

    # Fetch additional temporary parameters to assert manifest correctness
    params = {
        "filters": json.dumps(filters),
        "fields": ",".join(['file_id', 'file_name', 'md5', 'platform',
                            'reference_genome', 'data_category', 
                            'experimental_strategy', 'data_type', 
                            'cases.samples.sample_type']),
        "size": 10000
    } 

    # Query the GDC API for methylation files with the filters
    response = requests.get(GDC_QUERY_URL, params = params)
    response.raise_for_status()
    files_metadata = response.json()['data']['hits']
    df = pd.DataFrame(files_metadata)

    # Verify the integrity of the manifest and requested data
    assert len(files_metadata) < params["size"], "Query may be truncated."
    assert not df.empty, "No files returned from the GDC query."

    assert df['data_category'].eq(dc['data_category']).all()
    assert df['data_type'].eq(dc['data_type']).all()
    assert df['reference_genome'].eq(dc['reference_genome']).all()
    assert df['platform'].eq(dc['platform']).all()
    assert df['experimental_strategy'].eq(dc['experimental_strategy']).all()
    assert df['cases'].apply(lambda c: any(
        s['sample_type'] == dc['sample_type'] for s in c[0]['samples'])
    ).all()
    assert df['cases'].apply(lambda c: any(
        s['project_id'] == dc['project_id'] for s in c[0]['project'])
    ).all()

    # Clean the manifest for unused fields for GDC query
    manifest = (df[['file_id', 'file_name', 'md5']]
                .rename(columns = {'file_id': 'id', 'file_name': 'filename'}))

    return manifest


# =====| Download Methylation |=================================================

def download_methylation(manifest: pd.DataFrame, 
                         config: Dict,
                         layout: ProjectLayout) -> pd.DataFrame:
    """
    Downloads DNA methylation files from the GDC using a prevalidated manifest. 
    Includes additional safety and reproducibility checks, ensuring the GDC 
    Data Transfer Tool (`gdc-client`) is installed and available. Handles 
    downloading, MD5 verification, and retry logic.

    Verifies MD5 checksums for integrity and retries failed downloads up to 
    `max_retries` with exponential backoff. Logs per-file download status and 
    timestamp. Failures in download or MD5 verification will raise an exception.

    Note that this function does not bundle the `gdc-client` API in the \package; users must install it from the official GDC site.

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
        Status log recording per-file download status.
    """

    # Verify the `gdc-client` is properly installed on the user's device
    gdc_client_path = config.get('gdc_client', '')
    verify_gdc_client(gdc_client_path)

    # Download per file in the manifest using the GDC API
    status_log = []

    for _, row in manifest.iterrows():
        file_id = row['id']
        filename = row['filename']
        md5 = row['md5']
        filepath = layout.raw_dir / filename

        # Attempt at most MAX_RETRIES to download the file
        attempt, status, success, timestamp = 0, 'pending', False, -1

        while attempt < MAX_RETRIES and not success:
            attempt += 1
            timestamp = datetime.datetime.now(datetime.timezone.utc)

            try:
                 # Spawn a subprocess to run the API
                subprocess.run([gdc_client_path,
                                "download",
                                "-d", str(layout.raw_dir),
                                "-f", str(file_id)], 
                                check = True)

                # Verify MD5; if failed, remove the corrupt file and retry
                if verify_md5(filepath, md5):
                    success = True
                    status = 'success'
                else:
                    status = 'failed_md5'
                    filepath.unlink(missing_ok = True)
            
            # Download failed, clean up partial file to avoid persistence
            except subprocess.CalledProcessError:
                status = 'failed_download'
                filepath.unlink(missing_ok = True)

            # If failure, exponential backoff before next attempt
            if not success:
                time.sleep(2 ** attempt)

        status_log.append({
            'id': file_id,
            'filename': filename,
            'md5': md5,
            'status': status,
            'timestamp': timestamp
        })

        if not success:
            raise RuntimeError(f"File {filename} ({file_id}) filed to download "
                               f"successfully after {MAX_RETRIES} attempts.")

    return pd.DataFrame(status_log)




# =====| Build Metadata |=======================================================

def build_metadata(config: Dict) -> pd.DataFrame:
    # Fetches metadata using gdc-client API

    return pd.DataFrame()