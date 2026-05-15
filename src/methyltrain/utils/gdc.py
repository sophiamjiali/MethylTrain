# ==============================================================================
# Script:           gdc.py
# Purpose:          Verifies the GDC Client API
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-08
# ==============================================================================

import shutil

import pandas as pd

# ======| API Utilities |=======================================================

def verify_gdc_client(gdc_client_path) -> None:
    # Verifies that the gdc-client was downloaded

    if shutil.which(gdc_client_path) is None:
        raise RuntimeError(
            f"gdc-client executable not found at '{gdc_client_path}'.\n"
            "Please install the official GDC Data Transfer Tool from:\n"
            "https://gdc.cancer.gov/access-data/gdc-data-transfer-tool\n"
            "and ensure it is available in your PATH."
        )
    
def extract_project_id(cases):
    # Helper to extract nested metadata fields
    try:
        return (cases[0]['project']['project_id'] if cases 
                and cases[0].get('project') else pd.NA)
    except Exception:
        return pd.NA
    
def extract_file_name(files):
    # Helper to extract nested metadata fields
    try:
        return (files[0]['file_name'] if files 
                and files[0].get('file_name') else pd.NA)
    except Exception:
        return pd.NA

def extract_sample_type(cases):
    # Helper to extract nested metadata fields
    try:
        return (cases[0]['samples'][0]['sample_type'] if cases 
                and cases[0].get('samples') else pd.NA)
    except Exception:
        return pd.NA

def extract_submitter_id(cases):
    # Helper to extract nested metadata fields
    try:
        return cases[0]['submitter_id'] if cases else pd.NA
    except Exception:
        return pd.NA
    
def extract_aliquot_id(cases):
    # Helper to extract nested metadata fields
    field = (cases[0]['samples'][0]['portions'][0]['analytes'][0]
             ['aliquots'][0]['aliquot_id'])
    try:
        return (field if cases else pd.NA)
    except Exception:
        return pd.NA

def extract_batch_id(barcode: str):
    # Returns the portion and the plate ID together
    
    if pd.isna(barcode): return None
    parts = barcode.split('-')
    if len(parts) >= 6:
        return f"{parts[3]}-{parts[5]}"  # portion + plate
    return None