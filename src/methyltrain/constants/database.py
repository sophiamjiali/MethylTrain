# ==============================================================================
# Script:           database.py
# Purpose:          Defines constants for the AuditStore
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-14
# ==============================================================================

AUDIT_NAME = "audit"

AUDIT_FIELDS = {
    "file_id",
    "file_name",
    "download_status",
    "metadata_status",
    "biospecimen_status",
    "qc_status",
    "download_timestamp",
    "metadata_timestamp",
    "biospecimen_timestamp",
    "qc_timestamp",
    "raw_data_path",
    "updated_at"
}