# ==============================================================================
# Script:           constants.py
# Purpose:          Global constants for the package
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-08
# ==============================================================================

from pathlib import Path

# =====| Default Paths |========================================================

# Project root directory
PROJECT_ROOT = Path(__file__).parent.parent.parent

# Query link for the GDC API for fetching DNA methylation beta values
GDC_QUERY_URL = "https://api.gdc.cancer.gov/files"
GDC_QUERY_BATCH_URL = "https://api.gdc.cancer.gov/cases"

# =====| Downloading Logic |====================================================

# Maximum number of attempts to download a file using the `gdc-client`
MAX_RETRIES = 5

# =====| Supported types |======================================================

# Supported platform types; aligns with TCGA metadata case and spelling
PLATFORM_TYPES = [
    "Illumina Human Methylation 27",
    "Illumina Human Methylation 450",
    "Illumina Human Methylation Epic" # later add v2
]

# Defines resolution hierarchy
PLATFORM_PRIORITY = [
    "Illumina Human Methylation Epic",
    "Illumina Human Methylation 450",
    "Illumina Human Methylation 27"
]

# Supported Genome Builds; aligns with TCGA metadata case and spelling
REFERENCE_GENOME_TYPES = ["GRCh37", "GRCh38"]   # hg19 or hg38

# =====| Resource Paths |=======================================================

ANNOTATION_DIR = PROJECT_ROOT / "resources"

ANNOTATION_hg19_PATHS = {
    "Illumina Human Methylation 27": ANNOTATION_DIR / "illumina27k_annotation_hg19.parquet",

    "Illumina Human Methylation 450": ANNOTATION_DIR / "illumina450k_annotation_hg19.parquet",

    "Illumina Human Methylation Epic": ANNOTATION_DIR / "illuminaEPIC_annotation_hg19.parquet"
}

ANNOTATION_hg38_PATHS = {
    "Illumina Human Methylation 27": ANNOTATION_DIR / "illumina27k_annotation_hg38.parquet",

    "Illumina Human Methylation 450": ANNOTATION_DIR / "illumina450k_annotation_hg38.parquet",

    "Illumina Human Methylation Epic": ANNOTATION_DIR / "illuminaEPIC_annotation_hg38.parquet"
}