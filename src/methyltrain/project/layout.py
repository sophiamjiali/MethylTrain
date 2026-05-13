# ==============================================================================
# Script:           layout.py
# Purpose:          Defines and manages the filesystem layout for a dataset
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-13
# ==============================================================================

from pathlib import Path
from typing import List


class ProjectLayout:
    """
    Encapsulates the directory structure for a DNA methylation dataset 
    associated with a single TCGA project.

    This object centralizes all raw, metadata, manifest, processed, and 
    training directories. It is initialized from the user-specified 
    configurations object, which defines a root directory that within will hold 
    the projects (and cohorts) folders.

    Individual TCGA projects are defined in a self-titled (ID) folder within 
    the output_dir/projects/ folder.

    Parameters
    ----------
    config : dict
        Configuration dictionary controlling workflow steps.

    Attributes
    ----------
    project_name : str
    root_dir : Path
    project_dir: Path
    raw_dir : Path
    metadata_dir : Path
    manifest: Path
    metadata : Path
    audit_store : Path
    adata : Path
    dir_paths: List(Path)
    file_paths: List(Path)
    """

    def __init__(self, project_name: str, root_dir: str):
        self.project_name = project_name
        
        # Initialize directory structure
        self.project_dir = Path(root_dir)
        self.raw_dir = self.project_dir / "raw"
        self.metadata_dir = self.project_dir / "metadata"

        # Define output object destinations
        self.manifest = self.metadata_dir / "manifest.csv"
        self.metadata = self.metadata_dir / "metadata.csv"
        self.audit_store = self.metadata_dir / "audit_store.csv"
        self.adata = self.project_dir / f"{self.project_name}_adata.h5ad"

        # Define shortcut attributes for directories and files
        self.dir_paths = [self.project_dir, self.raw_dir, self.metadata_dir]
        self.file_paths = [self.manifest, self.metadata, self.audit_store]


    @classmethod
    def from_config(cls, config):
        project_id = config.get('project_id', '')
        root = (
            Path(config.get('paths', {}).get('output_dir', ''))
            / "projects"
            / project_id
        )
        return cls(project_name = project_id, root_dir = root)


    def initialize(self):
        for d in self.dir_paths: d.mkdir(parents = True, exist_ok = True)
        for p in self.file_paths: p.parent.mkdir(parents = True, exist_ok=True)


    def validate(self):

        # Validate required directories
        md: List[Path] = [d for d in self.dir_paths if not d.exists()]
        if md: raise FileNotFoundError(f"Missing directories: {md}")

        # Validate required parent directories
        mf: List[Path] = [p for p in self.file_paths if not p.parent.exists()]
        if mf: raise FileNotFoundError(f"Missing parent directories: {mf}")

# [END]