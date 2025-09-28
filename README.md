# WDP (WONDER Data Pipeline) Project

[![Python](https://img.shields.io/badge/Python-3.11+-blue.svg)](https://www.python.org/)
[![PyMC](https://img.shields.io/badge/PyMC-5.8+-green.svg)](https://www.pymc.io/)
[![R](https://img.shields.io/badge/R-4.x-blue.svg)](https://www.r-project.org/)
[![INLA](https://img.shields.io/badge/R--INLA-Latest-orange.svg)](https://www.r-inla.org/)

**WONDER Data Pipeline (WDP)** is a comprehensive data processing and Bayesian modeling system for epidemiological spatiotemporal analysis, focusing on the relationship between pesticide exposure and cancer mortality. The project supports two main modeling approaches: R-INLA based fast approximate inference and PyMC based full Bayesian MCMC sampling.

## Table of Contents

1. [Project Overview](#1-project-overview)
2. [Code Structure](#2-code-structure)
3. [Data Sources](#3-data-sources)
4. [Data Cleaning](#4-data-cleaning)
5. [Principal Component Analysis (PCA)](#5-principal-component-analysis-pca)
6. [PyMC Analysis Pipeline](#6-pymc-analysis-pipeline)
7. [BYM2 Model Specification](#7-bym2-model-specification)
8. [Quick Start](#8-quick-start)
9. [Configuration Management](#9-configuration-management)
10. [Quality Control and Diagnostics](#10-quality-control-and-diagnostics)
11. [Common Issues and Solutions](#11-common-issues-and-solutions)
12. [References and Documentation](#12-references-and-documentation)
13. [Contributors and Contact](#13-contributors-and-contact)

## 1. Project Overview

### 1.1. Core Features
- **Data Integration**: CDC mortality data, USGS pesticide usage, socioeconomic indicators, environmental factors.
- **Spatiotemporal Modeling**: BYM2 Bayesian spatiotemporal models with county-level adjacency effects and temporal trends.
- **Dual Modeling Framework**: R-INLA (fast) and PyMC (flexible) implementations.
- **Covariate Engineering**: PCA-derived Social Vulnerability Index and climate factors.
- **Cluster Computing**: SLURM job scheduling for large-scale batch analysis.

### 1.2. Technology Stack
- **Python**: Data processing, PyMC Bayesian modeling.
- **R**: INLA modeling, spatial statistical analysis.
- **Configuration Management**: Unified `config.yaml` path configuration.
- **Cluster Computing**: SLURM job submission and management.

## 2. Code Structure

```
WDP/
├── Code/
│   ├── Clean/                   # Data cleaning and processing scripts
│   │   ├── CDC_*.py             # CDC data cleaning (AAMR, Death, Location, Urbanization)
│   │   ├── County_Adjacency.py  # Spatial adjacency matrix generation
│   │   ├── ENV_*.py             # Environmental data processing (GEE, LUR, NLDAS)
│   │   ├── PNSP_*.py            # Pesticide data processing (Weight and Density)
│   │   └── SE_*.py              # Socioeconomic data processing
│   ├── Download/                # Data download scripts
│   ├── INLA/                    # R-INLA based analysis scripts
│   │   ├── INLA_Main.R          # Main production script for INLA models
│   │   ├── INLA_Config/         # Configuration files
│   │   ├── INLA_Scripts/        # Core INLA modeling scripts
│   │   └── INLA_Utils/          # Utility functions
│   ├── PYMC/                    # PyMC based analysis system
│   │   ├── main.py              # Command-line interface and workflow
│   │   ├── Utils_Data.py        # Data loading and preprocessing
│   │   ├── Utils_Model.py       # BYM2 model implementation
│   │   ├── Utils_Result.py      # Results extraction and formatting
│   │   └── Utils_Others.py      # Utility functions
│   └── Analysis/                # Analysis and PCA scripts
│       └── PCA.py               # Systematic PCA and covariate generation
├── Data/
│   ├── Original/                # Raw data from various sources
│   │   ├── CDC WONDER/          # CDC mortality data
│   │   ├── BEA/                 # Bureau of Economic Analysis
│   │   ├── CACES LUR/           # Air pollution data
│   │   ├── County Shapeline/    # Geographic boundary files
│   │   ├── GEE/                 # Google Earth Engine exports
│   │   ├── LAUS/                # Labor force statistics
│   │   ├── NLDAS/               # Meteorological data
│   │   ├── SAIPE/               # Poverty and income estimates
│   │   ├── SEER Population/     # Population demographics
│   │   ├── USDA ERS/            # Education and rural classifications
│   │   └── USGS PNSP/           # Pesticide usage estimates
│   └── Processed/               # Cleaned, analysis-ready datasets
│       ├── CDC/                 # Mortality and geographic data
│       ├── Environmental/       # Climate and land use data
│       ├── PCA/                 # Principal component analysis results
│       ├── Pesticide/           # Pesticide usage and density data
│       └── Socioeconomic/       # Economic and demographic indicators
├── Result/
│   ├── Filter/                  # Filtered analysis results
│   ├── PyMC_Results/            # PyMC modeling outputs
│   ├── Figures/                 # Visualization outputs
│   └── Tables/                  # Summary tables and diagnostics
├── config.yaml                  # Central configuration file
└── pymc_task.sh                 # SLURM job submission script
```

## 3. Data Sources

This project integrates data from multiple public sources to construct a comprehensive analytical dataset.

### 3.1. CDC (Outcomes and Geographic Classification)
- **Source**: CDC WONDER (Wide-ranging Online Data for Epidemiologic Research)
- **Content**: County-level mortality data (1999-2020), geographic classifications, and urbanization codes.
- **Key Files**: `Data/Original/ CDC WONDER/`, `Data/Original/CDC WONDER AAMR/`

### 3.2. Pesticide (USGS PNSP)
- **Source**: USGS Pesticide National Synthesis Project (PNSP)
- **Content**: County-level pesticide usage estimates (1992-2019).
- **Key Files**: `Data/Original/USGS PNSP/`

### 3.3. Socioeconomic
- **Sources**:
    - **SAIPE**: Poverty and income data.
    - **LAUS**: Labor force and unemployment data.
    - **BEA**: Economic indicators (GDP, per capita income).
    - **USDA-ERS**: Education levels and rural-urban classifications.
    - **SEER**: Population structure (age, race, sex).
    - **TIGER/Line**: County shapefiles for spatial analysis.
- **Content**: A wide range of county-level socioeconomic and demographic variables.
- **Key Files**: `Data/Original/SAIPE/`, `Data/Original/LAUS/`, `Data/Original/BEA/`, `Data/Original/USDA ERS/`, `Data/Original/SEER Population/`, `Data/Original/County Shapeline/`

### 3.4. Environmental
- **Sources**:
    - **CACES LUR**: Air pollution data (O3, CO, SO2, NO2, PM10, PM2.5).
    - **GEE (Google Earth Engine)**: Land cover (NLCD) and surface water (JRC) data.
    - **NLDAS**: Meteorological data (temperature, precipitation, wind, etc.).
- **Content**: Environmental and climate-related variables at the county level.
- **Key Files**: `Data/Original/CACES LUR/`, `Data/Original/GEE/`, `Data/Original/NLDAS/`

## 4. Data Cleaning

The `Code/Clean/` directory contains scripts for processing raw data into analysis-ready formats. All paths are managed through `config.yaml`.

### 4.1. CDC Data Processing
- **Scripts**: `CDC_AAMR_Merge.py`, `CDC_Death_Merge.py`, `CDC_Location_Urbanization.py`
- **Output**: `Data/Processed/CDC/`
- **Description**: Cleans and merges mortality, location, and urbanization data. `CDC_Data_Integrity_Checker.py` is used for validation.
- **Key Variables**:
    - **AAMR**: `COUNTY_FIPS`, `Year`, `Deaths`, `Population`, `AAMR` (Age-Adjusted Mortality Rate).
    - **Location**: `COUNTY_FIPS`, `HHS_Region`, `Census_Region`, `Census_Division`.
    - **Urbanization**: `COUNTY_FIPS`, `Year`, `Urbanization_Code`, `Urbanization_Type`.

### 4.2. Pesticide Data Processing
- **Scripts**: `PNSP_Merge.py`, `PNSP_Density.py`
- **Output**: `Data/Processed/Pesticide/`
- **Description**: Merges yearly pesticide usage files, maps compounds to categories (`mapping.csv`), and calculates usage density. `PNSP_Unique_Name.py` helps identify all unique compounds.
- **Key Variables**:
    - **PNSP.csv**: `COUNTY_FIPS`, `Year`, `cat{id}_min/avg/max` (kg/county/year), `chem{id}_min/avg/max` (kg/county/year).
    - **PNSP_Density.csv**: `COUNTY_FIPS`, `Year`, `cat{id}_min/avg/max` (kg/km²/year), `chem{id}_min/avg/max` (kg/km²/year).

### 4.3. Socioeconomic Data Processing
- **Scripts**: `SE_SAIPE_Merge.py`, `SE_LAUS_Merge.py`, `SE_BEA_Merge.py`, `SE_USDA_ERS_Education_Merge.py`, `SE_SEER_Population.py`, `County_Adjacency.py`
- **Output**: `Data/Processed/Socioeconomic/`
- **Description**: Processes various socioeconomic indicators and generates the county adjacency matrix required for spatial models.
- **Key Variables**:
    - **Poverty & Income**: `Poverty_Percent_All_Ages`, `Median_Household_Income`.
    - **Unemployment**: `Unemployment_Rate`.
    - **Economy**: `Per_Capita_Income`.
    - **Education**: `Less_Than_High_School_Percent`, `College_Plus_Percent`.
    - **Adjacency**: `county_from`, `county_to` for spatial models.

### 4.4. Environmental Data Processing
- **Scripts**: `ENV_LUR_Merge.py`, `ENV_GEE_Merge.py`, `ENV_NLDAS_Meterology.py`
- **Output**: `Data/Processed/Environmental/`
- **Description**: Processes air pollution, land use, and meteorological data. GEE data is collected using scripts like `ENV_GEE_County.py`, `ENV_GEE_JRC.py`, and `ENV_GEE_NLCD.py`.
- **Key Variables**:
    - **Air Pollution**: `O3`, `CO`, `SO2`, `NO2`, `PM10`, `PM25`.
    - **Land Use**: `nlcd_forest_km2`, `nlcd_urban_km2`, `nlcd_agriculture_km2`.
    - **Meteorology**: `tas_mean_annual` (temperature), `prcp_sum_annual` (precipitation).

## 5. Principal Component Analysis (PCA)

- **Script**: `Code/Analysis/PCA.py` (uses `PCA_Data_Loading.py`)
- **Outputs**:
  - `Data/Processed/PCA/PCA_Master_Covariables.csv` (master covariates with SVI_PC1 and ENV_PC1–PC3 scores)
  - `Result/Tables/PCA_Diagnose.csv` (unified diagnostics: eigenvalues, loadings, VIF)
  - `Result/Figure_Original_Data/*` (scree, loadings, sample scores for plots)
- **Description**: Performs VIF-based selection then PCA for SVI and ENV domains. `PCA_Plot.py` can visualize figures from the saved plot-data CSVs.

### 5.1. Socioeconomic Vulnerability Index (SVI)

- **Input Variables**:
    - `Poverty_Percent_All_Ages` (SAIPE)
    - `Median_Household_Income` (SAIPE)
    - `Unemployment_Rate` (LAUS)
    - `Per_Capita_Income` (BEA)
    - `Less_Than_High_School_Percent` (USDA-ERS)
    - `College_Plus_Percent` (USDA-ERS)
- **Processing**:
    1. VIF-based variable selection (threshold: 10.0).
    2. PCA with Kaiser criterion for component retention.
    3. SVI_PC1 is inverted to ensure higher values represent greater vulnerability.

### 5.2. Environmental Factors

- **Input Variables**:
    - Annual and seasonal temperature, precipitation, wind speed, humidity, and radiation from NLDAS.
    - Air pollutants (O₃, CO, SO₂, NO₂, PM₁₀, PM₂.₅) from CACES-LUR.
- **Output Components**: The analysis generates orthogonal environmental factors. Modeling uses ENV_PC1–PC3 when available (Kaiser criterion: eigenvalue > 1). Exact explained variance and loadings are reported in `Result/Tables/PCA_Diagnose.csv`.

  Based on the current diagnostics (`Result/Tables/PCA_Diagnose.csv`):
  - ENV_PC1 explains 33.1% of the variance
  - ENV_PC2 explains 28.3% of the variance
  - ENV_PC3 explains 13.3% of the variance
  - Cumulative (PC1–PC3): 74.7%

## 6. PyMC Analysis Pipeline

- **Main Script**: `Code/PYMC/main.py`
- **Description**: This is the command-line interface for running the Bayesian spatiotemporal analysis using PyMC. It orchestrates data loading, model fitting, and results extraction.

### 6.1. Key Modules
- `Utils_Data.py`: Handles loading and preprocessing of all data required for the model.
- `Utils_Model.py`: Implements the BYM2 model in PyMC.
- `Utils_Result.py`: Extracts and formats results from the model trace.
- `Utils_Others.py`: Contains helper functions for configuration, validation, and parsing.

### 6.2. Command-Line Interface
The `main.py` script provides a flexible CLI for running analyses.

| Parameter | Description | Default | Examples |
|-----------|-------------|---------|----------|
| `--disease` | ICD disease codes | (required) | `C81-C96`, `C50,C34` |
| `--compound` | Pesticide compounds/categories | (required) | `2`, `cat21`, `Atrazine` |
| `--model` | Model variants | `M0` | `M0`, `M1,M2,M3` |
| `--lag` | Exposure lag years | `5` | `5`, `10` |
| `--measure` | Exposure measures | `Weight` | `Weight`, `Density` |
| `--estimate` | Exposure estimates | `avg` | `min,avg,max` |
| `--sampling-mode` | MCMC configuration | `test` | `test`, `production` |
| `--dry-run` | Preview parameters without running | - | (flag) |

Advanced overrides: you can also pass `--cores`, `--chains`, `--draws`, `--tune`, `--target-accept`, `--output-dir`, and `--config-path` to control sampling and I/O.

### 6.3. SLURM Cluster Submission
- **Script**: `pymc_task.sh`
- **Description**: A shell script for submitting batch jobs to a SLURM cluster. It sets environment variables and calls `Code/PYMC/main.py` with appropriate arguments.

## 7. BYM2 Model Specification

For each county $i$ and time $t$, observed deaths $y_{i,t}$ follow:

$$y_{i,t} \sim \text{Poisson}(\lambda_{i,t})$$

$$\log(\lambda_{i,t}) = \log(E_{i,t}) + \eta_{i,t}$$

$$\eta_{i,t} = \alpha + \beta \log(X_{i,t-\text{lag}} + c) + \sum_j \gamma_j Z_{j,i,t} + \phi_i + v_i + \theta_t$$

Where:
- $\alpha$: Intercept
- $\beta$: Lagged pesticide exposure effect
- $\mathbf{Z}_{i,t}$: PCA-derived covariates (SVI, climate factors)
- $\phi_i$: Structured spatial effects (CAR prior)
- $v_i$: Unstructured spatial effects (IID prior)
- $\theta_t$: Temporal effects (RW1 prior)

### 7.1. Model Variants

| Model | Description | Covariates |
|-------|-------------|------------|
| **M0** | Base model | Pesticide + spatial + temporal effects only |
| **M1** | Social vulnerability | M0 + SVI_PC1 |
| **M2** | Environmental factors | M0 + ENV_PC1, ENV_PC2, ENV_PC3 |
| **M3** | Full model | M0 + SVI_PC1 + ENV_PC1, ENV_PC2, ENV_PC3 |

## 8. Quick Start

### 8.1. Environment Setup

**Python Environment (PyMC)**
```bash
# Create PyMC environment
conda create -n pymc python=3.11 -y
conda activate pymc

# Install dependencies
conda install -c conda-forge \
  pymc=5.8.0 arviz=0.16.0 pytensor=2.17.0 \
  numpy pandas scipy networkx geopandas \
  matplotlib seaborn tqdm pyyaml statsmodels
```

**R Environment (INLA)**
```r
# Install R-INLA
install.packages("INLA", repos=c(getOption("repos"), 
  INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

# Install other dependencies
install.packages(c("dplyr", "data.table", "sf", "spdep"))
```

### 8.2. Basic Usage

**PyMC Analysis**
```bash
# Single analysis
python Code/PYMC/main.py --disease C81-C96 --compound 2 --model M0 --lag 5 --estimate avg

# Batch analysis with multiple parameters
python Code/PYMC/main.py \
  --disease C81-C96,C50 \
  --compound 2,cat21 \
  --model M0,M3 \
  --lag 5,10 \
  --estimate min,avg,max \
  --measure Weight,Density

# Dry run to preview parameter resolution
python Code/PYMC/main.py --disease C81-C96 --compound 2,cat21 --estimate avg,max --dry-run
```

**SLURM Cluster Submission**
```bash
# Basic cluster job
sbatch pymc_task.sh

# Custom parameters
sbatch pymc_task.sh --compound "2,9,cat21" --estimate "avg,max" --model "M0,M3"
```

## 9. Configuration Management

All file paths and parameters are centrally managed in `config.yaml`:

```yaml
# Core data directories
data_directories:
  original: "Data/Original"
  processed: "Data/Processed"

# PyMC analysis configuration  
pymc_analysis:
  data_files:
    pca_covariates: "Data/Processed/PCA/PCA_Master_Covariables.csv"
    cdc_data_template: "Data/Processed/CDC/AAMR_{disease_code}.csv"
    pesticide_data: "Data/Processed/Pesticide/PNSP.csv"
    pesticide_density_data: "Data/Processed/Pesticide/PNSP_Density.csv"
    adjacency_data: "Data/Processed/Socioeconomic/County_Adjacency_List.csv"

  # Model configurations
  models:
    M0: {name: "base", covariates: []}
    M1: {name: "socioeconomic", covariates: ["SVI_std"]}
  M2: {name: "environmental", covariates: ["Climate1_std", "Climate2_std", "Climate3_std"]}  
  M3: {name: "full", covariates: ["SVI_std", "Climate1_std", "Climate2_std", "Climate3_std"]}

  # Sampling configurations
  sampling:
    test: {draws: 1000, tune: 500, chains: 2, cores: 2, target_accept: 0.8}
    production: {draws: 4000, tune: 2000, chains: 4, cores: 4, target_accept: 0.95}
```

## 10. Quality Control and Diagnostics

### 10.1. Data Validation
- **Integrity Checker**: `Code/Clean/CDC_Data_Integrity_Checker.py`
- **Coverage Validation**: Year ranges, geographic completeness.
- **Variable Consistency**: Column presence, data type validation.

### 10.2. Model Diagnostics
- **Convergence**: R-hat statistics, trace plots.
- **Effective Sample Size**: Bulk and tail ESS metrics.
- **Model Comparison**: WAIC, LOO cross-validation.
- **Posterior Predictive**: Goodness-of-fit assessment.

### 10.3. Output Validation
- **Parameter Bounds**: Biologically plausible relative risks.
- **Statistical Significance**: Bayesian p-values with multiple testing awareness.
- **Spatial Patterns**: Residual spatial autocorrelation checks.

## 11. Common Issues and Solutions

### 11.1. Environment Setup
```bash
# PyMC installation issues
conda install -c conda-forge pymc::pymc=5.8.0

# R-INLA installation  
install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable")
```

### 11.2. Memory Management
```bash
# Reduce batch size for large analyses
python Code/PYMC/main.py --measure Weight  # instead of Weight,Density

# Increase SLURM memory allocation
#SBATCH --mem-per-cpu=4G  # in pymc_task.sh
```

### 11.3. Convergence Issues
```python
# Increase sampling parameters
python Code/PYMC/main.py --draws 6000 --tune 3000 --target-accept 0.99

# Check spatial matrix connectivity
# (This is handled internally by the data loading scripts)
```

## 12. References and Documentation

### 12.1. Methodological References
1. Riebler A, et al. (2016). An intuitive Bayesian spatial model for disease mapping. *Statistical Methods in Medical Research* 25(4):1145-1165.
2. Salvatier J, et al. (2016). Probabilistic programming in Python using PyMC3. *PeerJ Computer Science* 2:e55.
3. Rue H, et al. (2009). Approximate Bayesian inference for latent Gaussian models. *Journal of the Royal Statistical Society B* 71(2):319-392.

### 12.2. Data Source Documentation
- **CDC WONDER**: https://wonder.cdc.gov/
- **USGS PNSP**: https://water.usgs.gov/nawqa/pnsp/
- **NLDAS**: https://ldas.gsfc.nasa.gov/nldas/
- **CACES**: https://www.caces.us/data

## 13. Contributors and Contact

**Development Team**: WDP Analysis Team  
**Institution**: Environmental Health Research  
**Last Updated**: September 28, 2025  
**Version**: 2.2

### 13.1. Getting Help
- Check `config.yaml` for path configurations.
- Use `--dry-run` to preview parameter resolution.
- Review model diagnostics in result CSVs.
- Consult convergence statistics for MCMC issues.

### 13.2. Contributing
1. Follow existing code structure and naming conventions.
2. Update `config.yaml` for new data sources or paths.
3. Document changes in commit messages.
4. Test both PyMC and INLA pipelines for compatibility.