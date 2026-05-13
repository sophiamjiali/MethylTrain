# ==============================================================================
# Script:           layout.py
# Purpose:          Defines and manages the filesystem layout for a dataset
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-07
# ==============================================================================

from pathlib import Path
from typing import Optional, Union, List

# Typing alias for Path-like objects
StrPath = Union[str, Path]


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
    project_list : List of Path
    cohort_adata : Path
    train_adata : Path
    val_adata : Path
    test_adata : Path
    """

    def __init__(self, 
                 cohort_name: str = "",
                 root_dir: Optional[StrPath] = None,
                 project_list: List = [StrPath],
                 cohort_adata: Optional[StrPath] = None,
                 train_adata: Optional[StrPath] = None,
                 val_adata: Optional[StrPath] = None,
                 test_adata: Optional[StrPath] = None):
        
        self.cohort_name = cohort_name
        self.project_list = project_list
        
        # Use the current working directory if `root_dir` is not provided
        root_path: Path = Path(root_dir) if root_dir is not None else Path.cwd()
        
        # Files are user-specified paths if provided, else defaults
        self.cohort_adata: Path = (
            Path(cohort_adata) if cohort_adata is not None else 
            root_path / f"{cohort_name}_cohort_adata.h5ad"
        )

        self.train_adata: Path = (
            Path(train_adata) if train_adata is not None else 
            root_path / f"{cohort_name}_train_adata.h5ad"
        )

        self.val_adata: Path = (
            Path(val_adata) if val_adata is not None else 
            root_path / f"{cohort_name}_val_adata.h5ad"
        )

        self.test_adata: Path = (
            Path(test_adata) if test_adata is not None else 
            root_path / f"{cohort_name}_test_adata.h5ad"
        )

        self.files = [self.cohort_adata, self.train_adata, self.val_adata, 
                      self.test_adata]


    def initialize(self):
        """
        Create all required directories if they do not exist.
        """

        for p in self.files:
            p.parent.mkdir(parents = True, exist_ok = True)
    
    def validate(self):
        """
        Ensure all required directories exist.

        Raises
        ------
        FileNotFoundError
            If any required directory is missing.
        ValueError
            If any of the files do not have extension `.csv`.
        """

        missing_paths: List[Path] = [p for p in self.files if 
                                     not p.parent.exists()]
        if missing_paths:
            raise FileNotFoundError(f"The following required paths are "
                                    f"not initialized: {missing_paths}")
        
        incorrect_paths: List[Path] = [p for p in self.files if 
                                       p.suffix != ".h5ad"]
        if incorrect_paths:
            raise ValueError(f"The following required paths must contain the `."
                             f"h5ad` extension: {incorrect_paths}")


