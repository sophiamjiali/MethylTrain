#!/usr/bin/env python3
# ==============================================================================
# Script:           prepare_cohort.py
# Purpose:          Downloads and processes and TCGA Project dataset
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-03-17
# ==============================================================================

import argparse
import os

import anndata as ad

from methyltrain.api.prepare import prepare_cohort
from methyltrain.config.loader import load_config
from methyltrain.fs.layout import CohortLayout
from methyltrain.api.steps import save_cohort

def main():

    # Parse the arguments provided to the entry-point script
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type = str, required = True)
    parser.add_argument("--verbose", action = "store_true", 
                        help = "Enable verbose logging")
    args = parser.parse_args()

    # Load the user-provided configurations
    config = load_config(args.config)
    cohort = config.get('project_id', '')

    if args.verbose: 
        print(f"=====| Beginning to process cohort {cohort} |=====")

    # Initialize the projects provided in the configurations
    project_dir = config.get('project_dir', '')
    all_projects = config.get('projects', [])
    project_list = []
    for project in all_projects:
        project_list.append(os.path.join(project_dir, f"{project}_adata.h5ad"))
    
    if args.verbose: print("Cohort contains projects: ", ", ".join(all_projects))

    # Initialize the default cohort layout
    layout = CohortLayout(
        cohort_name = cohort,
        root_dir = "./data",
        project_list = project_list,
        cohort_adata = f"../data/cohort/{cohort}_cohort_adata.h5ad",
        train_adata = f"../data/cohort/{cohort}_train_adata.h5ad",
        val_adata = f"../data/cohort/{cohort}_val_adata.h5ad",
        test_adata = f"../data/cohort/{cohort}_test_adata.h5ad"
    )
    layout.initialize()
    layout.validate()

    train_adata, val_adata, test_adata = prepare_cohort(config, layout, args.verbose)

    # Save the training splits if toggled
    if config.get('toggles', {}).get('split', True):
        train_adata.write_h5ad(layout.train_adata)
        val_adata.write_h5ad(layout.val_adata)
        test_adata.write_h5ad(layout.test_adata)
        if args.verbose: print("Saved")

    if args.verbose: 
        print(f"=====| Finished processing cohort {cohort} |=====")
        

if __name__ == "__main__":
    main()
