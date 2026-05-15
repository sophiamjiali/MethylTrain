# ==============================================================================
# Script:           paths.py
# Purpose:          Global constants for the package provided paths
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-14
# ==============================================================================

from pathlib import Path

# =====| Default Paths |========================================================

# Project root directory
PROJECT_ROOT = Path(__file__).parent.parent.parent
DEFAULT_CONFIG_DIR = PROJECT_ROOT / "configs"

# Query link for the GDC API for fetching DNA methylation beta values
GDC_QUERY_URL = "https://api.gdc.cancer.gov/files"
GDC_QUERY_BATCH_URL = "https://api.gdc.cancer.gov/cases"

ANNOTATION_DIR = PROJECT_ROOT / "resources"