# ==============================================================================
# Script:           project/pipeline.py
# Purpose:          Cleans raw TCGA DNA Methylation files from the GDC API
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-14
# ==============================================================================

import shutil
import logging

import pandas as pd

from pathlib import Path

from methyltrain.project.layout import ProjectLayout
from methyltrain.project.audit_store import AuditStore

logger = logging.getLogger(__name__)

# =====| Public API |===========================================================

def clean_data(layout: ProjectLayout, file_ids: list[str]) -> list[dict]:
    """
    Cleans raw TCGA DNA methylation beta value .txt files by converting them 
    to .parquet, flattening directory structure and removing accessory files 
    and raw files upon success.

    Updates the audit table with paths to the generated parquets. Renames .
    parquet files to each file's ID (UUID).

    Parameters
    ----------
    layout : ProjectLayout
        Object representing a project dataset directory layout.
    file_ids : list[str]
        List of file IDs to clean (based on download status).

    Returns
    -------
    list[dict]
        A cleaning report for auditing.
    """

    # Verify the raw data directory exists
    layout.validate()

    logger.info("=====| Attempting to Clean Raw Data |=====")

    report = []
    for file_id, file_name, raw_data_path in file_ids:
        if raw_data_path is not None: continue

        # Name the .parquet file with its UUID, not its old file name
        txt_path = layout.raw_dir / file_id / file_name
        parquet_path = layout.raw_dir / f"{file_id}.parquet"

        if not txt_path.exists():
            raise FileNotFoundError(f"Missing raw file: {txt_path}")
        
        _convert_txt_to_parquet(txt_path, parquet_path)
        _remove_raw_artifact(txt_path)
        
        # Update the AuditStore with the parquet path
        report.append({'file_id': file_id, 'raw_data_path': parquet_path})

    logger.info("=====| Successfully Cleaned Raw Data |=====")

    return report

# =====| Internal Helpers |=====================================================

def _convert_txt_to_parquet(txt_path: Path, parquet_path: Path) -> None:
    # Read the beta values (TCGA standard: probe_id, value)
    txt = pd.read_csv(txt_path, sep = '\t', header = 0, dtype={0: str})
    txt.columns = ['probe_id', 'beta_value']
    txt['beta_value'] = pd.to_numeric(txt['beta_value'], 
                                        errors = "coerce")
    
    txt.to_parquet(parquet_path, index = False)
    return

def _remove_raw_artifact(txt_path: Path) -> None:
    txt_path.unlink(missing_ok = True)
    if txt_path.parent.exists(): shutil.rmtree(txt_path.parent)
    return

# [END]