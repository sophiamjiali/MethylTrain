Final shape: (2143 x 18990)

(below fetched from RNA-CDM paper)

Project Code | Cancer Type                                | # Samples | # Download | Preproc. | AnnData Shape | Time to Download
---------------------------------------------------------------------------------------------------------------------
TCGA-LUAD    | Lung adenocarcinoma                        | 520       |         |       | 
TCGA-KIRP    | Kidney renal papillary cell carcinoma      | 298       |         |       | 
TCGA-COAD    | Colon adenocarcinoma                       | 289       |         |       | 
TCGA-CESC    | C. squamous cell carcinoma & Endoc. adeno. | 277       |         |       | 
TCGA-GBM     | Glioblastoma multiforme                    | 212       |         |       | 
TCGA-PAAD    | Pancreatic adenocarcinoma                  | 202       |         |       | 
TCGA-ESCA    | Oesophageal carcinoma                      | 156       |         |       | 
TCGA-OV      | Ovarian serous cystadenocarcinoma          | 83        |         |        | 
TCGA-UVM     | Uveal melanoma                             | 80        | 80         |        |  |
TCGA-CHOL    | Cholangiocarcinoma                         | 36        | 36         | 36       | (36 X 309993) | 15 minutes


Downloaded Status
[] TCGA-LUAD
[] TCGA-KIRP
[] TCGA-COAD
[] TCGA-CESC
[] TCGA-GBM
[] TCGA-PAAD
[] TCGA-ESCA
[X] TCGA-OV 
[X] TCGA-UVM
[X] TCGA-CHOL

[] CPTAC-3



To-Do:
- adding normalization, batch correction, variance filtering to workflow
- M-value transformation as well




- batch_correction (ComBat-seq, M-values) -> AFTER AGGREGATION
    - checks if M-value toggled, use one approach per

- add batch variables from metadata
- steps.py/preprocess(): add everything
- update configs to match default