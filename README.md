# MethylTrain
A Python package for high-throughput downloading and data engineering of DNA methylation data from NCI GDC for machine learning workflows

Next Steps
- analysis layer: CpGtools / DMseg / other quick analyses tools
- output some quick and basic visualizations
- basic GUI?


To Prepare a Project:

```
sbatch slurm/prepare_project.sh <TCGA-...>
```

To Prepare a Cohort:

```
sbatch slurm/prepare_cohort.sh <cohort>
```

Current training projects:
- TCGA-BLCA
- TCGA-BRCA
- TCGA-GBM
- TCGA-HNSC
- TCGA-KIRC
