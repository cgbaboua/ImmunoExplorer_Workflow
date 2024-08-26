# ImmunoExplorer Workflow

## Table of Contents
- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation](#installation)
  - [Clone the Repository](#clone-the-repository)
- [Pipeline Description](#pipeline-description)
  - [Differential Expression Analysis (DEG)](#differential-expression-analysis-deg)
  - [WGCNA Analysis](#wgcna-analysis)
  - [Correlation Analysis](#correlation-analysis)
  - [Lasso Regression](#lasso-regression)
  - [Corr_genes Analysis](#corr_genes-analysis)
  - [Pathways](#pathways)
  - [Protein-Protein Interaction (PPI)](#protein-protein-interaction-ppi)
  - [Random Walk with Restart (RWR)](#random-walk-with-restart-rwr)
- [Usage Instructions](#usage-instructions)
- [Dependencies](#dependencies)


## Overview

The ImmunoExplorer Workflow is designed to investigate autoimmune diseases with a particular focus on Chronic Spontaneous Urticaria (CSU). The pipeline includes transcriptomic analyses, co-expression module identification using WGCNA, Lasso regression for identifying key genes, protein-protein interaction mapping, and exploration of functional networks through Random Walk with Restart (RWR) on multiplex heterogeneous networks.

The pipeline enables the detection of crucial genes and biological pathways implicated in disease pathophysiology, providing insights for further investigation and potential therapeutic targeting.

## System Requirements

- **R version:** 4.3.2 or higher
- **Adequate storage:** Handling of large transcriptomic datasets.

## Installation

### Clone the Repository

Clone the workflow repository using the following command:

```bash
git clone https://github.com/cgbaboua/ImmunoExplorer_Workflow.git
```
## Pipeline Description

### Differential Expression Analysis (DEG)
**Script:** `DEG.Rmd`

The pipeline begins by identifying differentially expressed genes (DEGs) between patient samples and healthy controls using the Limma package. This step highlights the genes that are significantly up- or downregulated in disease states.

- **Inputs:** GEO datasets (GSE72540, GSE57178, GSE72541, GSE185516)
- **Outputs:** CSV files containing DEG lists.

### WGCNA Analysis
**Script:** `WGCNA.Rmd`

WGCNA (Weighted Gene Co-expression Network Analysis) is used to identify modules of co-expressed genes. These modules may represent functionally related groups of genes that are involved in similar biological processes.

- **Inputs:** DEG list from the previous step.
- **Outputs:** Co-expression modules and heatmaps correlating modules with clinical traits.

### Correlation Analysis
**Script:** `Correlation.Rmd`

This step identifies correlations between WGCNA co-expression modules and clinical traits, offering insights into how certain moodules relate to disease severity and other phenotypic factors.

- **Inputs:** Co-expression modules and clinical trait data.
- **Outputs:** Corrplots, corr matrix and lists of correlated modules with clinical data.

### Lasso Regression
**Script:** `Lasso.Rmd`

The Lasso regression method identifies key genes within the modules identified by WGCNA. This regularization method helps narrow down genes to the most influential ones within a co-expression module.

- **Inputs:** Co-expression modules from WGCNA.
- **Outputs:** List of key genes that contribute most significantly to disease traits.

### Corr_genes Analysis
**Script:** `Corr_genes.Rmd`

This script focuses on extracting and filtering genes that show strong correlation with clinical traits, prioritizing genes based on their correlation values.

- **Inputs:** DEG list and clinical data.
- **Outputs:** Filtered list of genes with high correlation scores.

### Pathways Enrichment Analysis
**Script:** `Pathways.Rmd`

The pathways enrichment step involves mapping the key genes to known biological pathways. This allows for the identification of overrepresented pathways, providing insights into the biological processes most affected in the disease state.

- **Inputs:** Key genes from lasso and DEG list.
- **Outputs:** Enriched pathways and their associated genes.


### Protein-Protein Interaction (PPI)
**Script:** `PPI.Rmd`

This step involves identifying interactions between proteins, leveraging databases such as STRING to construct a protein-protein interaction (PPI) network.

- **Inputs:** DEG list.
- **Outputs:** PPI network visualizations and interaction lists.

### Random Walk with Restart (RWR)
**Script:** `RWR.Rmd`

RWR is applied on a multiplex heterogeneous network that combines protein-protein interactions, gene co-expression data, enriched pathways, and clinical correlations. This step helps to prioritize genes and pathways that are central in the disease network.

- **Inputs:** PPI networks,enriched biological pathways,correlations.
- **Outputs:** Network of prioritized genes and pathways with their potential involvement in disease progression.

## Usage Instructions

1. **Set up the environment**: Ensure the required version of R is installed, along with all necessary libraries.
2. **Review outputs**: Check the CSV files, heatmaps, and network visualizations generated during each stage.

## Dependencies

### R Packages:
- `Limma`
- `WGCNA`
- `STRINGdb`
- `caret`
- `DOSE`
- `enrichplot`
- `corrplot`
- `igraph`
- `dnet`
- `RandomWalkRestartMH`

### Additional Libraries:
- `Dataframes`
- `tidyverse`
- Refer to each specific script for more details on additional required packages.
