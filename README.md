# Cardinal Analysis – Differential Gene Expression (DEG)

## Overview

This repository contains code for performing differential gene expression (DGE) analysis on pseudobulk single-cell RNA-seq data generated within the Cardinal project.

The analysis is performed per cell type using donor-level pseudobulk counts and donor-level metadata. The pipeline includes:

- Filtering donors by minimum cell count, ethnicity, duplicates
- Harmonisation of donor identifiers
- PCA and covariate association diagnostics
- Linear modeling using limma/edgeR
- Per-celltype per disease trait result export

This project is designed to run inside the UK Biobank RAP environment.

---

## Project Structure

cardinal_analysis_DGE/code/
- scDEA_LM_CT2_UKB.R # Main DGE script for UKB cohort to analyse pseudobulks at celltype level 2 granularity
- PCA_CT2.R # PCA script for data at celltype level 2 granularity
- scDEA_LM_CT2_G&H.R # Main DGE script for G&H cohort to analyse pseudobulks at celltype level 2 granularity
- scDEA_LM_CT3_UKB.R # Main DGE script for UKB cohort to analyse pseudobulks at celltype level 3 granularity
- PCA_CT3.R # PCA script for data at celltype level 3 granularity
- scDEA_LM_CT3_G&H.R # Main DGE script for G&H cohort to analyse pseudobulks at celltype level 3 granularity
- scDEA_LMM_CT3_G&H.R # Old DGE script to analyse pseudobulks at celltype level 2 granularity using a Linear Mixed Model

- renv.lock # Reproducible R environment snapshot
- README.md

---

## Input Data

### 1. Pseudobulk Counts

Counts files per cell type:

<celltype>.count.agg_sum.tsv


Structure:
- Rows: genes
- Columns: donor IDs
- First column: gene_id

---

### 2. Ncells Files 

<celltype>.ncells_x_donor.tsv

Used to filter donors with low cell numbers.

---

### 3. Donor Metadata - master_table_F3.csv

Donor-level metadata derived from cell-level master table:

- eid
- tranche_id
- pool_id
- donor_uid_tpd_norm
- sex
- age
- bmi
- smoking status
- ancestry

---

## PCA Diagnostics

For each cell type:

- PCA is performed on pseudobulk counts
- R² association with covariates is calculated
- main technical variation sources are regressed in linear model (PC1-4, corresponding to pool_id and tranche_id-associated variation)

---

## Differential Expression

DGE is performed per cell type using:

- edgeR
- limma

Models include relevant covariates (e.g., age, sex, bmi, ancestry/smoking_status).

---

## Environment & Reproducibility

This project uses `renv` for reproducibility.

To restore the environment in RAP:

```r
install.packages("renv")
renv::restore()

---
