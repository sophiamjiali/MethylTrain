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

    Parameters
    ----------
    manifest : pd.DataFrame
        Raw manifest with columns ['id', 'filename', 'md5'].
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

    # Set standard flags
    audit_table['downloaded'] = (audit_table['status'] == 'success').astype(int)
    audit_table['download_status'] = audit_table['status']
    audit_table['download_attempts'] = audit_table['attempts'].astype(int)
    audit_table['download_timestamp'] = audit_table['timestamp']

    audit_table.drop(columns = ['status', 'attempts', 'timestamp'])

    # Initialize future-stage columns
    audit_table['metadata_fetched'] = pd.NA
    audit_table['metadata_status'] = pd.NA
    audit_table['qc_pass'] = pd.NA
    audit_table['notes'] = ""

    return audit_table


def update_audit_table(audit_table: pd.DataFrame) -> pd.DataFrame:

    # Update the audit table with new stage results

    return pd.DataFrame()