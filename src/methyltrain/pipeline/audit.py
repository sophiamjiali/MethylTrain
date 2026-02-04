# ==============================================================================
# Script:           audit.py
# Purpose:          Internal auditing functions and implementation logic
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-21
# ==============================================================================

import pandas as pd

def initialize_audit_table(manifest: pd.DataFrame,
                           status_log: pd.DataFrame) -> pd.DataFrame:
    """
    Initializes an audit table from the raw manifest and download status log.

    Renames the `id` column from the manifest, as expected by the `gdc-client` 
    API, to `file_id`, the standard internal column name.

    Parameters
    ----------
    manifest : pd.DataFrame
        Raw manifest with columns ['id', 'filename'].
    status_log : pd.DataFrame
        Download status log with columns ['id', 'status', 'attempts', 
        'timestamp']

    Returns
    -------
    pd.DataFrame
        Audit table with one row per file, including initial download flags and 
        placeholders for metadata and quality control stages.
    """

    # Merge the manifest with the status log on the file ID (index)
    audit_table = manifest.merge(
        status_log[['status', 'attempts', 'timestamp']],
        how = 'left', left_index = True, right_index = True
    )

    # Set standard flags
    audit_table['downloaded'] = (audit_table['status'] == 'success').astype(int)
    audit_table['download_status'] = audit_table['status']
    audit_table['download_attempts'] = audit_table['attempts'].astype(int)
    audit_table['download_timestamp'] = audit_table['timestamp']

    audit_table = audit_table.drop(columns = ['status', 'attempts', 
                                              'timestamp'])

    # Initialize future-stage columns
    audit_table['metadata_fetched'] = pd.NA
    audit_table['metadata_status'] = pd.NA
    audit_table['biospecimen_fetched'] = pd.NA
    audit_table['biospecimen_status'] = pd.NA
    audit_table['qc_pass'] = pd.NA
    audit_table['parquet_path'] = pd.NA
    audit_table['parquet_path'] = audit_table['parquet_path'].astype("object")
    audit_table['notes'] = ""

    return audit_table


def update_metadata(audit_table: pd.DataFrame, 
                    metadata: pd.DataFrame) -> pd.DataFrame:
    # Updates the audit_table with the metadata

    audit_table = audit_table.merge(metadata[['status']], how = 'left',
                                    left_index = True, right_index = True)
    audit_table['metadata_fetched'] = (audit_table['status']
                                       .eq('success').astype(int))
    audit_table['metadata_status'] = audit_table['status']
    audit_table = audit_table.drop(columns = 'status')

    return audit_table