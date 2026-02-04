# ==============================================================================
# Script:           load_utils.py
# Purpose:          Loading utility functions for the package
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-08
# ==============================================================================

import pandas as pd

from ..fs.layout import ProjectLayout

# =====| Loading |==============================================================

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

# =====| Saving |===============================================================

def save_audit_table(audit_table: pd.DataFrame, layout: ProjectLayout) -> None:
    audit_table.to_csv(layout.audit_table, sep = '\t', header=True, index=True)

def save_metadata(metadata: pd.DataFrame, layout: ProjectLayout) -> None:
    metadata.to_csv(layout.metadata, sep = '\t', header=True, index=True)

def save_manifest(manifest: pd.DataFrame, layout: ProjectLayout) -> None:
    manifest.to_csv(layout.manifest, sep = '\t', header=True, index=True)

def save_status_log(status_log: pd.DataFrame, layout: ProjectLayout) -> None:
    status_log.to_csv(layout.status_log, sep = '\t', header=True, index=True)

# [END]