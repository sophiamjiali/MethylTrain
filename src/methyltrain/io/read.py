# ==============================================================================
# Script:           read.py
# Purpose:          Reading utility functions for the package
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-08
# ==============================================================================

import pandas as pd

from methyltrain.project.layout import ProjectLayout

def load_metadata(layout: ProjectLayout) -> pd.DataFrame:
    # Loads the metadata table with file_id as index

    layout.validate()
    metadata = pd.read_csv(layout.metadata, sep = '\t', index_col = 0)
    return metadata


def load_manifest(layout: ProjectLayout) -> pd.DataFrame:
    # Loads the manifest with file_id as index

    layout.validate()
    manifest = pd.read_csv(layout.manifest, sep = '\t', index_col = 0)
    return manifest

# [END]