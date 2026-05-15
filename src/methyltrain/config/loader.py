# ==============================================================================
# Script:           loader.py
# Purpose:          Loads user-provided YAML configurations
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-07
# ==============================================================================

import yaml
from pathlib import Path

def load_config(config_path: str):
    """
    Load a configuration file.

    Parameters
    ----------
    config_path : str or Path, optional
        Path to a YAML configuration file.

    Returns
    -------
    config : dict
        Configuration dictionary with defaults overridden by user-provided 
        values.

    Raises
    ------
    FileNotFoundError
        If a configuration path is provided and doesn't exist.
    """
    path = Path(config_path)
    if not path.is_file():
            raise FileNotFoundError(f"Configuration file not found: {path}")
    
    with path.open('r') as f:
        config = yaml.safe_load(f) or {}
    return config

# [END]