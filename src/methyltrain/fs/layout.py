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
    Encapsulates the directory strucdture for a DNA methylation dataset 
    associated with a multi-project cohort.

    This object centralizes project and training directories.
    """

    def __init__(self, 
                 cohort_name: str = "",
                 project_list: List = [],
                 training_dir: StrPath = ""):
        
        self.cohort_name = cohort_name
        self.project_list = project_list
        self.training_dir = training_dir



        # ==================
        # under construction
        # ==================


    def initialize(self):
        return






class ProjectLayout:
    """
    Encapsulates the directory structure for a DNA methylation dataset 
    associated with a single project.

    This object centralizes all raw, metadata, manifest, processed, and 
    training directories. Users can provide either:

    1. A single `root_dir` (subdirectories will default to `raw`, `metadata`, 
       `manifest`, `processed`, and `training`) or
    2. Full paths for each individual directory.

    Parameters
    ----------
    project_name : str
        Name of the project.
    root_dir : str or Path, optional
        Root directory of the dataset. Defaults to the current working 
        directory.
    raw_dir : str or Path, optional
        Directory for raw methylation data. Overrides `root_dir` default if 
        provided.
    raw_metadata : str or Path, optional
        Path for the raw metadata file. Overrides `root_dir` default if 
        provided.
    raw_manifest : str or Path, optional
        Path for the raw manifest file. Overrides `root_dir` default if 
        provided.
    processed_metadata : str or Path, optional
        Path for the processed metadata file. Overrides `root_dir` default if 
        provided.
    processed_manifest : str or Path, optional
        Path for the processed manifest file. Overrides `root_dir` default if 
        provided.
    processed_dir : str or Path, optional
        Directory for final processed data. Overrides `root_dir` default if 
        provided.

    Attributes
    ----------
    project_name : str
    raw_dir : Path
    raw_metadata : Path
    raw_manifest: Path
    processed_metadata : Path
    processed_manifest : Path
    processed_dir : Path
    """

    def __init__(self,
                 project_name: str = "",
                 root_dir: Optional[StrPath] = None,
                 raw_dir: Optional[StrPath] = None,
                 raw_metadata: Optional[StrPath] = None,
                 raw_manifest: Optional[StrPath] = None,
                 processed_metadata: Optional[StrPath] = None,
                 processed_manifest: Optional[StrPath] = None,
                 processed_dir: Optional[StrPath] = None):
        
        self.project_name = project_name
        
        # Use the current working directory if `root_dir` is not provided
        root_path: Path = Path(root_dir) if root_dir is not None else Path.cwd()
        
        # Directories are user-specified paths if provided, else defaults
        self.raw_dir: Path = (Path(raw_dir) if raw_dir is not None else 
                              root_path / "raw")
        
        self.raw_metadata: Path = (
            Path(raw_metadata) if raw_metadata is not None 
            else root_path / f"{[project_name]}_raw_metadata.csv"
        )

        self.raw_manifest: Path = (
            Path(raw_manifest) if raw_manifest is not None 
            else root_path / f"{project_name}_raw_manifest.csv"
        )

        self.processed_metadata: Path = (
            Path(processed_metadata) if processed_metadata is not None 
            else root_path / f"{[project_name]}_processed_metadata.csv"
        )

        self.processed_manifest: Path = (
            Path(processed_manifest) if processed_manifest is not None 
            else root_path / f"{project_name}_processed_manifest.csv"
        )

        self.processed_dir: Path = (Path(processed_dir) if processed_dir is not 
                                   None else root_path / "processed")

        self.paths = [self.raw_dir, self.raw_metadata, self.raw_manifest, 
                      self.processed_metadata, self.processed_manifest, 
                      self.processed_dir]
        
        self.files = [self.raw_metadata, self.raw_manifest, self.
                      processed_metadata, self.processed_manifest]


    def initialize(self):
        """
        Create all required directories if they do not exist.

        Parameters
        ----------
        overwrite : bool, default False
            if True, remove existing directories before creating them.
        """

        for p in self.paths:
            if p.is_dir():
                p.mkdir(parents = True, exist_ok = True)
            elif p.is_file():
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

        missing_paths: List[Path] = [p for p in self.paths 
                                    if (p.is_dir() and not p.exists())
                                    or (p.is_file() and not p.parent.exists())]
        if missing_paths:
            raise FileNotFoundError(f"The following required paths are "
                                    f"not initialized: {missing_paths}")
        
        incorrect_paths: List[Path] = [p for p in self.files if 
                                       p.suffix != ".csv"]
        if incorrect_paths:
            raise ValueError(f"The following required paths must contain the `."
                             f"csv` extension: {incorrect_paths}")