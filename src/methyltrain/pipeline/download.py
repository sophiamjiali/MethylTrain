# ==============================================================================
# Script:           download.py
# Purpose:          Internal downloading functions and implementation logic
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-20
# ==============================================================================

import json
import requests

import pandas as pd

from typing import Dict

from ..fs.layout import ProjectLayout
from ..constants import GDC_QUERY_URL

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
    # Returns the status log
    return pd.DataFrame()




# =====| Build Metadata |=======================================================

def build_metadata(config: Dict) -> pd.DataFrame:
    # Fetches metadata using gdc-client API

    return pd.DataFrame()