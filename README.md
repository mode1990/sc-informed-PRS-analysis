# sc-informed-PRS-analysis
single cell informed polygenic risk scoring 

# scRNAseq Informed Polygenic Risk Scoring

This repository contains the R script analysis pipeline for scRNAseq informed polygenic risk scoring. This analysis is designed to integrate single-cell RNA sequencing (scRNAseq) data with polygenic risk scores (PRS) for predicting clinical outcomes.

## Table of Contents
1. [Installation](#installation)
2. [Pipeline Overview](#pipeline-overview)
3. [Usage](#usage)
4. [Input Data](#input-data)
5. [Output Data](#output-data)
6. [Scripts](#scripts)
7. [License](#license)
8. [Contact](#contact)

## Installation

Ensure you have R and the necessary libraries installed. The required libraries include:

```r
install.packages(c("data.table", "dplyr", "readxl", "lm.beta", "caret", "e1071", "pROC", "glmnet"))
```

Additionally, make sure to set the library path correctly if you are using specific directories for package installations:

```r
.libPaths(c("/data/Common_Folder/R/Single_cell_packages/", .libPaths()))
```

## Pipeline Overview

The pipeline follows these major steps:
1. Load annotated GWAS genes.
2. Load and filter differentially expressed genes (DEGs).
3. Merge DEGs with GWAS genes.
4. Define SNP column format.
5. Load Parkinson's disease GWAS summary statistics.
6. Merge SNPs and prepare PRSice2 input.
7. Calculate PRS using PRSice2.
8. Predict clinical outcomes using linear and logistic regression.
9. Perform k-fold cross-validation.
10. Apply regularization techniques to address overfitting.

## Usage

To run the analysis, follow these steps:

1. Ensure all required input files are available in the specified directories.
2. Set the correct paths for input and output files within the script.
3. Execute the R script using R or RStudio.

```r
source("path/to/your/script.R")
```

## Input Data

- Annotated GWAS genes: `/data/dehestani/scPRS_analysis/Annotated_PD_GWAS/gwasgenes`
- DEG list: `path/to/DEGs`
- PD GWAS summary statistics: `/data/dehestani/scPRS_analysis/Sum_stats/Chang2017_GWAS.tab`
- PRSice2 output: `/data/dehestani/scPRS_analysis/PRsice2_output/bestPRS`
- Clinical outcomes: `/data/dehestani/scPRS_analysis/Clinical_outcomes/clinical_outcomes.csv`
- Covariates file: `/data/dehestani/scPRS_analysis/Clinical_outcomes/Covariates`

## Output Data

- PRSice2 input summary statistics: `/data/dehestani/scPRS_analysis/PRsice2_input_sumstats/ODC_test.txt`
- Various regression models' outputs and cross-validation results.
- ROC curves and AUC values for model performance evaluation.

## Scripts

The main script contains several sections, each performing specific tasks:

1. **Loading and Preprocessing Data**:
    - Load GWAS genes and DEG list.
    - Filter and merge data.
    - Define SNP column format.

2. **PRS Calculation**:
    - Merge GWAS summary statistics with DEGs.
    - Prepare input for PRSice2.

3. **Clinical Outcome Prediction**:
    - Linear regression models for UPDRS-III, MoCA, and BDI-II.
    - Logistic regression for case/control prediction.

4. **Cross-validation and Regularization**:
    - K-fold cross-validation for regression models.
    - Lasso and Ridge regularization.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For any questions or issues, please contact Mo Dehestani at smdehestani@gmail.com


