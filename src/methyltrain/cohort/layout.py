# ==============================================================================
# Script:           cohort/layout.py
# Purpose:          Defines and manages the filesystem layout for a cohort
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-14
# ==============================================================================

from pathlib import Path
from typing import List

from methyltrain.constants.cohort_layout import (
    COHORT_ADATA,
    TRAIN_ADATA,
    VAL_ADATA,
    TEST_ADATA
)

class CohortLayout:
    """
    Encapsulates the directory structure for a DNA methylation dataset 
    associated with a multi-project cohort.

    This object centralizes project and training directories. Users can provide either:

    1. A single `root_dir` with default filenames, or
    2. Full paths for each individual directory.

    Parameters
    ----------
    cohort_name : str
        Name of the cohort.
    project_list : List of str or Path
        List of paths to project .h5ad AnnData objects.
    cohort_adata : str or Path, optional
        Path for the full cohort AnnData object.
    train_adata : str or Path, optional
        Path for the training-split cohort AnnData object.
    val_adata : str or Path, optional
        Path for the validation-split cohort AnnData object.
    test_adata : str or Path, optional
        Path for the test-split cohort AnnData object.

    Attributes
    ----------
    cohort_name : str
    project_list : List[Path]
    cohort_dir : Path
    cohort_adata : Path
    train_adata : Path
    val_adata : Path
    test_adata : Path
    dir_paths: List(Path)
    file_paths: List(Path)
    """

    def __init__(self, 
                 cohort_name: str,
                 project_dir: str,
                 cohort_dir: str,
                 project_list: List[str]):
        
        self.cohort_name = cohort_name
        self.project_list = [Path(p) for p in project_list]

        # Initialize directory structure
        self.project_list = Path(project_dir)
        self.cohort_dir = Path(cohort_dir)
        
        # Define output object destinations
        self.cohort_adata = self.cohort_dir / COHORT_ADATA
        self.train_adata = self.cohort_dir / TRAIN_ADATA
        self.val_adata = self.cohort_dir / VAL_ADATA
        self.test_adata = self.cohort_adata / TEST_ADATA

        # Define shortcut attributes for directories and files
        self.dir_paths = [self.cohort_dir]
        self.file_paths = [self.cohort_adata, self.train_adata, 
                           self.val_adata, self.test_adata]


    @classmethod
    def from_config(cls, config):
        cohort_id = config.get('cohort', {}).get('id', '')
        project_dir = config.get('paths', {}).get('project_dir', '')
        root = (
            Path(config.get('paths', {}).get('output_dir', ''))
            / "cohorts"
            / cohort_id
        )
        project_list = config.get('project_list', [])
        return cls(cohort_id, project_dir, root, project_list)
    

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