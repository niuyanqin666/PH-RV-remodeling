# PH-RV Remodeling Analysis

## Overview

This repository contains analysis code for a cross-model, cross-species study of ventricular remodeling in pulmonary hypertension.

The study identifies:

* A conserved ventricular remodeling program across left and right ventricles
* A robust right ventricular (RV) consensus remodeling module using multi-dataset integration
* External validation in an independent CTEPH cohort

---

## Analysis workflow

The analysis is organized by figures:

* **Fig1–3**: Definition of Shared / LV-only / RV-only remodeling programs
* **Fig4**: Cross-model and cross-species replication analysis
* **Fig5**: Anchor-based definition of RV consensus remodeling module
* **Fig6**: External validation in CTEPH dataset

---

## Repository structure

* `scripts/` – R scripts for each figure (Fig1–Fig6)
* `data/` – input data (large files partially excluded, see below)
* `outputs/` – generated figures

---

## Data availability

### 1. Public transcriptomic datasets (GEO)

The following datasets are available from:

* Gene Expression Omnibus

* GSE198618 (human RV)

* GSE240921 (human RV)

* GSE242014 (rat PAB)

* GSE240923 (rat MCT)

* GSE186989 (rat SuHx)

* GSE133402 (rat hypoxia)

* GSE266139 (rat MCT discovery dataset)

---

### 2. CTEPH dataset (external validation, Figure 6)

Data were obtained from the **source data** of:

* Transcriptional changes of the extracellular matrix in chronic thromboembolic pulmonary hypertension govern right ventricle remodeling and recovery

The following supplementary files are used:

* 44161_2025_672_MOESM4_ESM
* 44161_2025_672_MOESM5_ESM
* 44161_2025_672_MOESM6_ESM
* 44161_2025_672_MOESM7_ESM

---

### 3. Proteomics and Olink data

Data were obtained from the **source data** of:

* Identification of LTBP-2 as a plasma biomarker for right ventricular dysfunction in human pulmonary arterial hypertension

The following supplementary files are used:

* 44161_2022_113_MOESM3_ESM
* 44161_2022_113_MOESM4_ESM

---

### Notes

* Some large files (e.g., GSE186989) are not included due to GitHub size limitations
* Please refer to `data/*.txt` files for download instructions where applicable
* All scripts assume data files are placed in the `data/` directory

---

## How to run

1. Download required datasets (see `data/` folder and instructions)
2. Place all files into the `data/` directory
3. Run scripts in R:

```r
source("scripts/Fig4.R")
```

Each script reproduces one figure.

---

## Contact

Yanqin Niu
niuyq@szu.edu.cn
