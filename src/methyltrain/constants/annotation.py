# ==============================================================================
# Script:           annotation.py
# Purpose:          Global constants for the package provided annotations
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-14
# ==============================================================================

from methyltrain.constants.paths import ANNOTATION_DIR
# =====| Supported types |======================================================

# Supported platform types; aligns with TCGA metadata case and spelling
PLATFORM_TYPES = [
    "Illumina Human Methylation 27",
    "Illumina Human Methylation 450",
    "Illumina Methylation Epic" # later add v2
]

# Defines resolution hierarchy
PLATFORM_PRIORITY = [
    "Illumina Methylation Epic",
    "Illumina Human Methylation 450",
    "Illumina Human Methylation 27"
]

# Supported Genome Builds; aligns with TCGA metadata case and spelling
REFERENCE_GENOME_TYPES = ["GRCh37", "GRCh38"]   # hg19 or hg38

# =====| Resource Paths |=======================================================

ANNOTATION_hg19_PATHS = {
    "Illumina Human Methylation 27": ANNOTATION_DIR / "illumina27k_annotation_hg19.parquet",

    "Illumina Human Methylation 450": ANNOTATION_DIR / "illumina450k_annotation_hg19.parquet",

    "Illumina Methylation Epic": ANNOTATION_DIR / "illuminaEPIC_annotation_hg19.parquet"
}

ANNOTATION_hg38_PATHS = {
    "Illumina Human Methylation 27": ANNOTATION_DIR / "illumina27k_annotation_hg38.parquet",

    "Illumina Human Methylation 450": ANNOTATION_DIR / "illumina450k_annotation_hg38.parquet",

    "Illumina Methylation Epic": ANNOTATION_DIR / "illuminaEPIC_annotation_hg38.parquet"
}