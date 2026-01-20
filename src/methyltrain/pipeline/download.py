# ==============================================================================
# Script:           download.py
# Purpose:          Internal downloading functions and implementation logic
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-20
# ==============================================================================

import pandas as pd

from typing import Dict

from ..fs.layout import ProjectLayout
from ..constants import GDC_QUERY_URL

def build_manifest(config: Dict) -> pd.DataFrame:
    # Fetches manifest using gdc-client API

    pc = config.get('project', {})
    dc = config.get('download', {})

    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id",
                                     "value": [pc('project_id')]}},
            {"op": "in", "content": {"field": "data_category",
                                     "value": [dc('data_category')]}},
            {"op": "in", "content": {"field": "experimental_strategy",
                                     "value": [dc('experimental_strategy')]}},
            {"op": "in", "content": {"field": "data_type",
                                     "value": [dc('data_type')]}},
            {"op": "in", "content": {"field": "platform",
                                     "value": [dc('platform')]}},
            {"op": "in", "content": {"field": "reference_genome",
                                     "value": [dc('reference_genome')]}},
            {"op": "in", "content": {"field": "cases.samples.sample_type",
                                     "value": [dc('sample_type')]}},
            {"op": "in", "content": {"field": "files.access",
                                     "value": ["open"]}}
        ]
    }


    return pd.DataFrame()



def download_methylation(manifest: pd.DataFrame, 
                         config: Dict,
                         layout: ProjectLayout) -> pd.DataFrame:
    # Returns the status log
    return pd.DataFrame()






def build_metadata(config: Dict) -> pd.DataFrame:
    # Fetches metadata using gdc-client API

    return pd.DataFrame()