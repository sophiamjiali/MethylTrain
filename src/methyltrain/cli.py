# ==============================================================================
# Script:           cli.py
# Purpose:          Main stable entry for both project and cohort workflows
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-15
# ==============================================================================

import logging
import argparse

from pathlib import Path

from methyltrain.project.layout import ProjectLayout
from methyltrain.audit.audit_store import AuditStore
from methyltrain.utils.logging import configure_logger
from methyltrain.config.loader import load_config
from methyltrain.project.pipeline import prepare_project
from methyltrain.cohort.pipeline import prepare_cohort
from methyltrain.io.datasets import save_project, save_cohort
from methyltrain.io.write import save_manifest, save_metadata

from methyltrain.constants.paths import DEFAULT_CONFIG_DIR


# =====| CLI Entry-Point |======================================================

def main():
    args = _build_parser().parse_args()

    config_path = _resolve_config_path(args.command, args.config, args.name)
    config = load_config(config_path)

    if args.command == "project": run_project_command(config)
    elif args.command == "cohort": run_cohort_command(config)


def run_project_command(config: dict) -> None:
    """
    Wraps calling and saving the outputs of the project workflow.
    """

    # Initialize the project's default output layout
    layout = ProjectLayout.from_config(config)
    layout.initialize()
    layout.validate()

    # Initialize a logger and AuditStore to capture verbose output
    logger = configure_logger(level = logging.INFO)
    results = prepare_project(config, layout)

    # Capture all audit reports and export it as a CSV
    audit = AuditStore(layout.audit_store.with_suffix(".db"))
    audit.initialize(results.manifest.index.tolist())
    audit.apply_download_report(results.download_report)
    audit.apply_metadata_report(results.metadata_report)
    audit.apply_biospecimen_report(results.biospecimen_report)
    audit.apply_cleaning_report(results.cleaning_report)
    audit.apply_qc_report(results.qc_report)
    audit.export_csv(layout.audit_store)

    # Save all metadata to CSV files as specified in the layout
    save_manifest(results.manifest, layout)
    save_metadata(results.metadata, layout)
    save_project(results.adata, layout)

def run_cohort_command(config: dict) -> None:
    """
    Wraps calling and saving the outputs of the cohort workflow.
    """

    prepare_cohort(config)

# =====| Internal Helpers |=====================================================

def _build_parser():
    parser = argparse.ArgumentParser("MethylTrain")

    subparsers = parser.add_subparsers(dest = "command")

    project = subparsers.add_parser("project")
    project.add_argument("--config", type = str, default = None)
    project.add_argument("--name", type = str, default = None)

    cohort = subparsers.add_parser("cohort")
    cohort.add_argument("--config", type = str, default = None)
    cohort.add_argument("--name", type = str, default = None)

    return parser


def _resolve_config_path(command: str,
                         config: str | None, 
                         name: str | None) -> Path:
    """
    Resolve configuration from either an explicitly provided path or a name.
    """

    if config and name:
        raise ValueError("Provide only one of --config or --name")

    if config:
        path = Path(config)
        if not path.exists():
            raise FileNotFoundError(f"Config file not found: {config}")
        return path
    
    if name:
        candidate = DEFAULT_CONFIG_DIR / command / f"{name}.yaml"
        if not candidate.exists():
            raise FileNotFoundError(f"Named config not found: {name}")
        return candidate

    raise ValueError("Must provide either --config or --name")

if __name__ == "__main__":
    main()