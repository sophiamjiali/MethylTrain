#!/usr/bin/env python3
# ==============================================================================
# Script:           prepare_project.py
# Purpose:          Downloads and processes and TCGA Project dataset
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-26
# ==============================================================================

import argparse

from methyltrain.api.prepare import prepare_dataset
from methyltrain.config.loader import load_config
from methyltrain.fs.layout import ProjectLayout
from methyltrain.api.steps import save_project
from methyltrain.utils.load_utils import save_audit_table

def main():

    # Parse the arguments provided to the entry-point script
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type = str, required = True)
    parser.add_argument("--verbose", action = "store_true", 
                        help = "Enable verbose logging")
    parser.add_argument("--clean-raw-data", action = "store_true", 
                        help = "Delete raw data")
    args = parser.parse_args()

    # Load the user-provided configurations
    config = load_config(args.config)
    project = config.get('project_id', '')

    if args.verbose: 
        print(f"=====| Beginning to process project {project} |=====")

    # Initialize the default project layout
    layout = ProjectLayout(
        project_name = project,
        root_dir = "./data",
        raw_dir = f"./data/raw/{project}",
        audit_table = f"./data/metadata/{project}/{project}_audit_table.csv",
        metadata = f"./data/metadata/{project}/{project}_metadata.csv",
        manifest = f"./data/metadata/{project}/{project}_manifest.csv",
        status_log = f"./data/metadata/{project}/{project}_status_log.csv",
        project_adata = f"./data/processed/{project}_adata.h5ad"
    )
    layout.initialize()
    layout.validate()

    adata, audit_table = prepare_dataset(config, layout, args.verbose)
    save_project(adata, layout)

    # Save the updated audit table
    save_audit_table(audit_table, layout)

    # Optional clean-up: remove raw data, manifest and audit table persists
    if args.clean_raw_data: 
        for parquet_file in layout.raw_dir.glob("*.parquet"):
            parquet_file.unlink()

    if args.verbose: 
        print(f"=====| Finished processing project {project} |=====")
        

if __name__ == "__main__":
    main()