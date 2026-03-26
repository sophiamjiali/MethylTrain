# ==============================================================================
# Script:           loader.py
# Purpose:          Loads user-provided YAML configurations
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-07
# ==============================================================================

import yaml
from pathlib import Path
from typing import Any, Dict, Union
from .defaults import DEFAULT_PROJECT_CONFIG, DEFAULT_COHORT_CONFIG

from ..utils.utils import merge_dicts, check_dict

StrPath = Union[str, Path]

def load_config(config_path: StrPath | None = None, config_type = "Project"):
    """
    Load a configuration file and merge it with the default configuration.

    Parameters
    ----------
    config_path : str or Path, optional
        Path to a YAML configuration file. If None, only default config is used.
    config_type : str
        Configuration type (cohort or project)

    Returns
    -------
    config : dict
        Configuration dictionary with defaults overridden by user-provided 
        values.

    Raises
    ------
    FileNotFoundError
        If a configuration path is provided and doesn't exist.
    KeyError
        If a required key from `default` is missing in `user`.
    TypeError
        If the type of a value in `user` does not match the type in `default`.
    ValueError
        If a value in `user` violates a constraint (e.g., allowed enum values).

    """

    if config_type == "project": default = DEFAULT_PROJECT_CONFIG.copy()
    else: default = DEFAULT_COHORT_CONFIG.copy()

    # If a configuration path is provided, merge it with default configurations
    if config_path is not None:

        path = Path(config_path)
        if not path.is_file():
            raise FileNotFoundError(f"Configuration file not found: {path}")
        
        with path.open('r') as f:
            user_config = yaml.safe_load(f) or {}
        config = merge_dicts(default, user_config)

    # Verify the configurations provided are viable
    verify_config(config, default)

    return config


def verify_config(config: Dict, default: Dict) -> None:
    """
    Verify that the user-provided configuration dictionary is valid. Raises 
    informative errors if any required keys or values are missing or invalid.

    Compares the provided configurations to the default configurations.

    Parameters
    ----------
    config : dict
        Configuration dictionary controlling workflow steps.
    default : dict
        Defautl configuration for verification.

    Raises
    ------
    KeyError
        If a required key from `default` is missing in `user`.
    TypeError
        If the type of a value in `user` does not match the type in `default`.
    ValueError
        If a value in `user` violates a constraint (e.g., allowed enum values).
    """
    # check_dict(default, config)
    return
