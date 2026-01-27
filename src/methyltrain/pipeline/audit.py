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

    # Merge the manifest with the status log on the file ID
    audit_table = manifest.merge(
        status_log[['id', 'status', 'attempts', 'timestamp']],
        on = 'id',
        how = 'left'
    )

    audit_table = audit_table.rename(columns = {'id': 'file_id'})

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
    audit_table['qc_pass'] = pd.NA
    audit_table['parquet_path'] = pd.NA
    audit_table['parquet_path'] = audit_table['parquet_path'].astype("object")
    audit_table['notes'] = ""

    # Initialize table index to the file identifier
    audit_table = audit_table.set_index("file_id", drop = True)

    return audit_table
