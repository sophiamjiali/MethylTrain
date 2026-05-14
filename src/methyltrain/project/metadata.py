# ==============================================================================
# Script:           metadata.py
# Purpose:          Queries for project metadata from the GDC API
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-20
# ==============================================================================

import json
import logging
import requests

import pandas as pd

from typing import Dict, List

from methyltrain.audit.audit_store import AuditStore
from methyltrain.constants import GDC_QUERY_URL, GDC_QUERY_BATCH_URL
from methyltrain.utils.utils import (
    extract_batch_id,
    extract_project_id,
    extract_sample_type,
    extract_submitter_id,
    extract_aliquot_id,
    extract_batch_id
)

logger = logging.getLogger(__name__)


# =====| Public API |===========================================================

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