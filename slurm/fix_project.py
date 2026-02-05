#!/usr/bin/env python3

import argparse
from sys import audit

from methyltrain.api.prepare import prepare_dataset
from methyltrain.config.loader import load_config
from methyltrain.fs.layout import ProjectLayout
from methyltrain.api.steps import save_project
from methyltrain.pipeline.download import build_metadata
from methyltrain.utils.load_utils import (
    save_audit_table, load_audit_table, save_metadata
)

def main():

    # Parse the arguments provided to the entry-point script
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type = str, required = True)
    args = parser.parse_args()

    # Load the user-provided configurations
    config = load_config(args.config)
    project = config.get('project_id', '')

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

    audit_table = load_audit_table(layout)
    metadata = build_metadata(audit_table, config)

    save_metadata(metadata, layout)
        

if __name__ == "__main__":
    main()