# ==============================================================================
# Script:           defaults.py
# Purpose:          Defines default configurations for the workflow
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-07
# ==============================================================================

# =====| Default Configurations |===============================================

DEFAULT_CONFIG = {

    # Path to the gdc-client executable
    "gdc_client": "/Volumes/FBI_Drive/gdc-client",

    "project_id":   "TCGA-KIRP",      # Project name (as on TCGA)

    "seed": 42,
    "split": [0.60, 0.20, 0.20],

    "download": {

        # Fields to filter the DNA methylation data query
        "data_category":         "DNA Methylation",
        "experimental_strategy": "Methylation Array",
        "data_type":             "Methylation Beta Value",
        "platform":              "Illumina Human Methylation 450" ,
        "reference_genome":      "GRCh38",
        "sample_type":           "Primary Tumor",

        # Fields to request from the query
        "metadata": [
            'file_id',
            'file_name',
            'cases.project.project_id',
            'cases.submitter_id',
            'cases.samples.sample_type',
            'data_category',
            'data_type',
            'experimental_strategy',     # e.g. add mutational/clinical metadata
            'platform',                  #  | to query for here for facilitating
            'reference_genome'           #  | further downstream analysis
        ]
    },

    "toggles": {

        "sample_qc":        True,
        "probe_qc":         True,
        "imputation":       True,    # Perform imputation, else default to 0
        "batch_correction": True,    # Perform batch correction upon cohort
        "gene_aggregation": True,    # Aggregates cohort beta values to genes
        "winsorization":    True     # Clip extreme values for ML stability

    },

    "quality_control": {
        
        "sample_qc": {
            "missing_threshold": 0.05,
            "outlier_threshold": 1.5     
        },

        "probe_qc": {
            "remove_cross_reactive": True,
            "remove_SNP_associated": True,
            "remove_sex_chromosome": True,
            "remove_multi_mapped": True,
            "missing_threshold": 0.05, 
        }
    },

    "preprocessing": {
        "clip_values": [0.001, 0.999]
    },

    "gene_aggregation": {
        "regions": ["TSS200", "TSS1500"]    # TSS200, TSS1500, and/or gene_body
    }
}