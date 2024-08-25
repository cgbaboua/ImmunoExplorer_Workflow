# ImmunoExplorer Workflow

## Table of Contents
- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation](#installation)
  - [Clone the Repository](#clone-the-repository)
- [Pipeline Description](#pipeline-description)
  - [Differential Expression Analysis (DEG)](#differential-expression-analysis-deg)
  - [WGCNA Analysis](#wgcna-analysis)
  - [Lasso Regression](#lasso-regression)
  - [Protein-Protein Interaction (PPI)](#protein-protein-interaction-ppi)
  - [Random Walk with Restart (RWR)](#random-walk-with-restart-rwr)
- [Usage Instructions](#usage-instructions)
- [Dependencies](#dependencies)
- [Additional Resources](#additional-resources)

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
