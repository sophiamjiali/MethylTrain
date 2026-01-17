# ==============================================================================
# Script:           constants.py
# Purpose:          Global constants for the package
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-08
# ==============================================================================

# =====| Supported types |======================================================

# Supported array types; aligns with TCGA metadata case and spelling
ARRAY_TYPES = [
    "Illumina Human Methylation 27",
    "Illumina Human Methylation 450",
    "Illumina Human Methylation Epic" # later add v2
]

# Supported Genome Builds; aligns with TCGA metadata case and spelling
GENOME_BUILD_TYPES = ["hg19", "hg38"]

# =====| Resource Paths |=======================================================

ANNOTATION_hg19_PATHS = {
    "Illumina Human Methylation 27": "resources/illumina27k_annotation_hg19.parquet",

    "Illumina Human Methylation 450": "resources/illumina450k_annotation_hg19.parquet",

    "Illumina Human Methylation Epic": "resources/illuminaEPIC_annotation_hg19.parquet"
}

ANNOTATION_hg38_PATHS = {
    "Illumina Human Methylation 27": "resources/illumina27k_annotation_hg38.parquet",

    "Illumina Human Methylation 450": "resources/illumina450k_annotation_hg38.parquet",

    "Illumina Human Methylation Epic": "resources/illuminaEPIC_annotation_hg38.parquet"
}