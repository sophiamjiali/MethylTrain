# ==============================================================================
# Script:           load_utils.py
# Purpose:          Loading utility functions for the package
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-08
# ==============================================================================

import pandas as pd

from ..fs.layout import ProjectLayout

# =====| Metadata I/O |=========================================================

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

def load_biospecimen(layout: ProjectLayout) -> pd.DataFrame:
    # Loads the status log with file_id as index
    layout.validate()
    biospecimen = pd.read_csv(layout.biospecimen, sep = '\t', index_col = 0)
    return biospecimen