# ==============================================================================
# Script:           logging.py
# Purpose:          Defines a logger for verbose pipeline output
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-05-13
# ==============================================================================

import logging
import sys

def configure_logger(name: str = "methyltrain", 
                     level: int = logging.INFO) -> logging.Logger:
    """
    Configure and return a package-level logger to capture and record pipeline 
    verbose output.
    """

    logger = logging.getLogger(name)
    if logger.handlers: return logger

    logger.setLevel(level)

    handler = logging.StreamHandler(sys.stdout)
    logger.addHandler(handler)

    formatter = logging.Formatter(
        fmt = "[%(asctime)s] %(levelname)s - %(message)s",
        datefmt = "%Y-%m-%d %H:%M:%S"
    )
    handler.setFormatter(formatter)

    logger.propagate = False
    return logger
    


