# MethylTrain

A Python package for high-throughput downloading and data engineering of DNA methylation data from the [NCI Genomic Data Commons (GDC)](https://portal.gdc.cancer.gov/) for machine learning workflows.

MethylTrain handles the full data preparation pipeline — from raw TCGA methylation array downloads through preprocessing and cohort assembly — producing clean, ML-ready datasets (AnnData format) for use in downstream modelling.

---

## Features

- Automated bulk download of DNA methylation data from NCI GDC
- Per-project and multi-project cohort preparation via SLURM
- Preprocessing of raw methylation arrays (filtering, normalization, quality control)
- Outputs ML-ready AnnData (`.h5ad`) objects compatible with PyTorch and scikit-learn workflows
- Configuration-driven pipeline via YAML

---

## Repository Structure

```
MethylTrain/
├── config/             # Pipeline configuration files
├── docs/               # Documentation
├── notebooks/          # Exploratory Jupyter notebooks
├── resources/          # Reference materials
├── scripts/            # Pipeline entry-points
├── slurm/              # SLURM batch scripts for HPC execution
├── src/
│   └── methyltrain/    # Core package source
├── environment.yaml    # Conda environment
└── pyproject.toml      # Package metadata
```

---

## Installation

```bash
conda env create -f environment.yaml
conda activate methyltrain-env
pip install -e .
```

Key dependencies: Python 3.11, NumPy, Pandas, AnnData, h5py, scikit-learn, scipy, statsmodels, inmoose, zarr.

---

## Usage

MethylTrain is designed to run on HPC clusters via SLURM, though the underlying scripts can be run directly as well.

### Prepare a Single Project

Downloads and preprocesses methylation data for a single TCGA project:

```bash
sbatch slurm/prepare_project.sh <TCGA-PROJECT-CODE>
# e.g.
sbatch slurm/prepare_project.sh TCGA-GBM
```

### Prepare a Cohort

Assembles and preprocesses data across multiple projects into a unified cohort:

```bash
sbatch slurm/prepare_cohort.sh <cohort>
```

---

## Supported Projects

MethylTrain has been used with the following cancer projects:

| Project | Cancer Type |
|---------|-------------|
| TCGA-BLCA | Bladder urothelial carcinoma |
| TCGA-BRCA | Breast invasive carcinoma |
| TCGA-CESC | Cervical squamous cell carcinoma & endocervical adenocarcinoma |
| TCGA-CHOL | Cholangiocarcinoma |
| TCGA-COAD | Colon adenocarcinoma |
| TCGA-ESCA | Oesophageal carcinoma |
| TCGA-GBM | Glioblastoma multiforme |
| TCGA-HNSC | Head and neck squamous cell carcinoma |
| TCGA-KIRC | Kidney renal clear cell carcinoma |
| TCGA-KIRP | Kidney renal papillary cell carcinoma |
| TCGA-LUAD | Lung adenocarcinoma |
| TCGA-OV | Ovarian serous cystadenocarcinoma |
| TCGA-PAAD | Pancreatic adenocarcinoma |
| TCGA-UVM | Uveal melanoma |
| CPTAC-3 | NCI Clinical Proteomic Tumor Analysis Consortium |

---

## Configuration

Pipeline behaviour is controlled via YAML files in `config/`. Parameters include GDC access settings, preprocessing options (filtering thresholds, normalization strategy), and output paths.

---

## License

Distributed under the BSD 3-Clause License. See [`LICENSE`](LICENSE) for details.
