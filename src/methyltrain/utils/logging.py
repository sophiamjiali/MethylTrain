# ==============================================================================
# Script:           logging.py
# Purpose:          Defines a logger for verbose pipeline output
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-13
# ==============================================================================

import logging
import sys

def configure_logger(level: int = logging.INFO) -> None:
    """
    Configure and return a package-level logger to capture and record pipeline 
    verbose output.
    """

    logging.basicConfig(
        level=level,
        stream=sys.stdout,
        format="[%(asctime)s] %(levelname)s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )
    return

# [END]