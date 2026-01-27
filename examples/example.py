
# change this to import directlry from methyltrain (add to methyltrain/__init__.py)

from methyltrain.api.steps import (
    download,
    clean_data,
    quality_control,
    preprocess,
    aggregate_cohort,
    cohort_batch_correction,
    aggregate_genes,
    winsorize,
    split,
    load_raw_project,
    load_processed_project,
    save_project
)
from methyltrain.fs.layout import ProjectLayout
from methyltrain.config.loader import load_config

def main():

    # User provides top-level paths
    layout = ProjectLayout(
        raw_dir = "/data/raw",
        metadata_dir = "/data/metadata",
        manifest_dir = "/data/manifest",
        processed_dir = "/data/processed",
        training_dir = "/data/training"
    )
    layout.initialize()

    # User provides a configuration file
    config = load_config("methylation_preproc.yaml")




    # Save the dataset as an AnnData with the default processed object name
    project_name = config.get('project', '')
    project_path = layout.processed_dir / f"{project_name}_anndata.h5ad"
    save_project(adata, project_path)




if __name__ == "__main__":
    main()