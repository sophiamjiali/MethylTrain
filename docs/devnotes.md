

- use CPTAC-3 as example/testing dataset


- working on downloading and cleaning logic
- refactor cleaning logic to take sample identifiers, so I can clean the download with the same function
    - look at implementation: sample identifier in configurations, map that using the metadata?
    - not implemented for now


- need to implement:
    - download.py: build_manifest()
    - download.py: download_methylation()
    - download.py: download_metadata()
    - clean.py: normalize_metadata() for downloading logic
    - layout.py: CohortLayout class definition