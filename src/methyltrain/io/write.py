# ==============================================================================
# Script:           write.py
# Purpose:          Writing utility functions for the package
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-08
# ==============================================================================

import pandas as pd

from methyltrain.project.layout import ProjectLayout
from methyltrain.cohort.layout import CohortLayout


def save_metadata(metadata: pd.DataFrame, layout: ProjectLayout) -> None:
    # Saves the metadata table with file_id as the first column
    metadata.to_csv(layout.metadata, sep = '\t', header=True, index=True)

def save_manifest(manifest: pd.DataFrame, layout: ProjectLayout) -> None:
    # Saves the manifest table with file_id as the first column
    manifest.to_csv(layout.manifest, sep = '\t', header=True, index=True)

def save_probe_set(probe_set: pd.DataFrame, layout: CohortLayout) -> None:
    # Note that NO index is included
    probe_set.to_csv(layout.probe_set, sep = '\t', header=True, index=False)

# [END]