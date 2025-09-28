# WDP (WONDER Data Pipeline) Project

[![Python](https://img.shields.io/badge/Python-3.11+-blue.svg)](https://www.python.org/)
[![PyMC](https://img.shields.io/badge/PyMC-5.8+-green.svg)](https://www.pymc.io/)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![INLA](https://img.shields.io/badge/R--INLA-Latest-orange.svg)](https://www.r-inla.org/)

**WONDER Data Pipeline (WDP)** is a comprehensive data processing and Bayesian modeling system for epidemiological spatiotemporal analysis, focusing on the relationship between pesticide exposure and cancer mortality. The project supports two main modeling approaches: R-INLA based fast approximate inference and PyMC based full Bayesian MCMC sampling.

## 🎯 **Project Overview**

### **Core Features**
- **Data Integration**: CDC mortality data, USGS pesticide usage, socioeconomic indicators, environmental factors
- **Spatiotemporal Modeling**: BYM2 Bayesian spatiotemporal models with county-level adjacency effects and temporal trends
- **Dual Modeling Framework**: R-INLA (fast) and PyMC (flexible) implementations
- **Covariate Engineering**: PCA-derived Social Vulnerability Index and climate factors
- **Cluster Computing**: SLURM job scheduling for large-scale batch analysis

### **Technology Stack**
- **Python**: Data processing, PyMC Bayesian modeling
- **R**: INLA modeling, spatial statistical analysis
- **Configuration Management**: Unified config.yaml path configuration
- **Cluster Computing**: SLURM job submission and management

## 📁 **Project Structure**Sources and Processing Guide

This document records the full pipeline from Original Data → Cleaning/Processing → Analysis-ready datasets in the WDP project. As we migrate/standardize scripts, please continue to update this file with each data source’s inputs, outputs, variable definitions, and assumptions to ensure reproducibility and maintainability.

## Contents

- CDC (Outcomes and Geographic Classification)
- Pesticide (USGS PNSP)
- Socioeconomic (SAIPE / LAUS / BEA / USDA-ERS / SEER / County Adjacency)
- Environmental (CACES LUR / GEE / NLDAS)
- PCA-Derived Covariates (SVI and Climate Factors)

---

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
│       └── PCA_Systematic.py    # Systematic PCA and covariate generation
├── Data/
│   ├── Original/                # Raw data from various sources
│   │   ├──  CDC WONDER/         # CDC mortality data (note the leading space)
│   │   ├── BEA/                 # Bureau of Economic Analysis
│   │   ├── CACES LUR/           # Air pollution data
│   │   ├── County Shapeline/    # Geographic boundary files
│   │   ├── GEE/                 # Google Earth Engine exports
│   │   ├── LAUS/                # Labor force statistics
│   │   ├── NLDAS/               # Meteorological data
│   │   ├── SAIPE/               # Poverty and income estimates
│   │   ├── SEER Population/     # Population demographics
│   │   ├── USDA ERS/           # Education and rural classifications
│   │   └── USGS PNSP/          # Pesticide usage estimates
│   └── Processed/               # Cleaned, analysis-ready datasets
│       ├── CDC/                 # Mortality and geographic data
│       ├── Environmental/       # Climate and land use data
│       ├── PCA/                 # Principal component analysis results
│       ├── Pesticide/           # Pesticide usage and density data
│       └── Socioeconomic/       # Economic and demographic indicators
├── Result/
│   ├── Filter/                  # Filtered analysis results
│   ├── PyMC_Results/           # PyMC modeling outputs
│   ├── Figures/                # Visualization outputs
│   └── Tables/                 # Summary tables and diagnostics
├── config.yaml                 # Central configuration file
└── pymc_task.sh               # SLURM job submission script
```

## 🚀 **Quick Start**

### **Environment Setup**

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

### **Basic Usage**

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

## 🔬 **Bayesian Spatiotemporal Modeling**

### **BYM2 Model Specification**

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

### **Model Variants**

| Model | Description | Covariates |
|-------|-------------|------------|
| **M0** | Base model | Pesticide + spatial + temporal effects only |
| **M1** | Social vulnerability | M0 + SVI_PC1 |
| **M2** | Environmental factors | M0 + ENV_PC1, ENV_PC2 |
| **M3** | Full model | M0 + SVI_PC1 + ENV_PC1 + ENV_PC2 |

### **Exposure Parameterization**

The system supports flexible exposure estimation:

- **Compounds**: Individual chemicals (e.g., `2` for 2,4-D) or categories (e.g., `cat21` for herbicides)
- **Estimates**: `min`, `avg`, `max` estimates from USGS PNSP
- **Measures**: `Weight` (kg/county/year) or `Density` (kg/km²/year)
- **Lag periods**: 5-year or 10-year exposure windows

## 📊 **Data Sources and Processing**

### **CDC Mortality Data**
- **Source**: CDC WONDER (Wide-ranging Online Data for Epidemiologic Research)
- **Coverage**: County-level, 1999-2020
- **Processing**: Age-adjusted mortality rates, censored data handling
- **Key Files**:
  - `Data/Processed/CDC/AAMR_{ICD}.csv`: Age-adjusted mortality rates
  - `Data/Processed/CDC/Location.csv`: County geographic classifications
  - `Data/Processed/CDC/Urbanization.csv`: Rural-urban classifications

### **Pesticide Exposure Data**
- **Source**: USGS PNSP (Pesticide National Synthesis Project)
- **Coverage**: County-level, 1992-2019
- **Processing**: Compound mapping, category aggregation, density calculations
- **Key Files**:
  - `Data/Processed/Pesticide/PNSP.csv`: Usage weights (kg/county/year)
  - `Data/Processed/Pesticide/PNSP_Density.csv`: Usage densities (kg/km²/year)
  - `Data/Processed/Pesticide/mapping.csv`: Compound-category relationships

### **Socioeconomic Indicators**
- **Sources**: SAIPE (poverty), LAUS (unemployment), BEA (income), USDA-ERS (education), SEER (demographics)
- **Coverage**: County-level panels, 1999-2020
- **Processing**: VIF-based variable selection, PCA dimension reduction
- **Key Variables**: Poverty rates, unemployment, income, education levels

### **Environmental Factors**
- **Sources**: NLDAS (meteorology), CACES-LUR (air pollution), GEE (land cover)
- **Coverage**: County-level, 1999-2020
- **Processing**: Seasonal aggregation, spatial interpolation, PCA factors
- **Key Variables**: Temperature, precipitation, air pollutants, land use

### **Spatial Data**
- **Source**: TIGER/Line county shapefiles
- **Processing**: Adjacency matrix generation for spatial modeling
- **Output**: `Data/Processed/Socioeconomic/County_Adjacency_List.csv`

## 🧮 **Principal Component Analysis (PCA)**

### **Social Vulnerability Index (SVI)**

**Input Variables**:
- Poverty_Percent_All_Ages (SAIPE)
- Median_Household_Income (SAIPE)  
- Unemployment_Rate (LAUS)
- Per_Capita_Income (BEA)
- Less_Than_High_School_Percent (USDA-ERS)
- College_Plus_Percent (USDA-ERS)

**Processing**:
1. VIF-based variable selection (threshold: 10.0)
2. PCA with Kaiser criterion for component retention
3. SVI_PC1 inversion (higher values = greater vulnerability)

### **Environmental Factors**

**Input Variables**:
- Annual and seasonal temperature (NLDAS)
- Annual and seasonal precipitation (NLDAS)
- Wind speed, humidity, radiation (NLDAS)
- Air pollutants: O₃, CO, SO₂, NO₂, PM₁₀, PM₂.₅ (CACES-LUR)

**Output Components**:
- **ENV_PC1**: Climate/precipitation component (33.1% variance)
- **ENV_PC2**: Temperature/radiation component (28.3% variance)  
- **ENV_PC3**: Wind/seasonal patterns (13.3% variance)

**Master Output**: `Data/Processed/PCA/Master_Covariates.csv`

## 📈 **Results and Output Format**

### **CSV Results Structure**

```csv
Timestamp,Disease,Exposure,Category,Measure,Estimate,Lag,Model,RR_per_SD,RR_per_IQR,RR_Q1_vs_Q1,RR_Q2_vs_Q1,RR_Q3_vs_Q1,RR_Q4_vs_Q1,P_Value,R_hat,ESS_bulk,WAIC,N_Counties,N_Records,Status_Message
2025-09-28 14:30:25,C81-C96,2,4-D,Herbicides,Weight,avg,5,M3,1.0234 (0.9876, 1.0612),1.1234 (1.0123, 1.2456),1.0 (Reference),1.0123 (0.9987, 1.0256),1.0287 (1.0012, 1.0587),1.0456 (1.0089, 1.0845),0.034*,1.0012,1250,2342.12,3108,98765,SUCCESS
```

### **Key Result Metrics**

- **RR_per_SD**: Relative risk per standard deviation increase in log-exposure
- **RR_per_IQR**: Relative risk for interquartile range increase
- **RR_Qx_vs_Q1**: Quartile contrast relative risks
- **P_Value**: Bayesian two-sided p-value with significance markers
- **R_hat**: Convergence diagnostic (< 1.01 preferred)
- **ESS_bulk**: Effective sample size (≥ 400 preferred)
- **WAIC**: Watanabe-Akaike Information Criterion

## ⚙️ **Configuration Management**

All file paths and parameters are centrally managed in `config.yaml`:

```yaml
# Core data directories
data_directories:
  original: "Data/Original"
  processed: "Data/Processed"

# PyMC analysis configuration  
pymc_analysis:
  data_files:
    pca_covariates: "Data/Processed/PCA/Master_Covariates.csv"
    cdc_data_template: "Data/Processed/CDC/AAMR_{disease_code}.csv"
    pesticide_data: "Data/Processed/Pesticide/PNSP.csv"
    pesticide_density_data: "Data/Processed/Pesticide/PNSP_Density.csv"
    adjacency_data: "Data/Processed/Socioeconomic/County_Adjacency_List.csv"

  # Model configurations
  models:
    M0: {name: "base", covariates: []}
    M1: {name: "socioeconomic", covariates: ["SVI_std"]}
    M2: {name: "environmental", covariates: ["Climate1_std", "Climate2_std"]}  
    M3: {name: "full", covariates: ["SVI_std", "Climate1_std", "Climate2_std"]}

  # Sampling configurations
  sampling:
    test: {draws: 1000, tune: 500, chains: 2, cores: 2, target_accept: 0.8}
    production: {draws: 4000, tune: 2000, chains: 4, cores: 4, target_accept: 0.95}
```

## 🔧 **Command Line Interface**

### **PyMC Parameters**

| Parameter | Description | Default | Examples |
|-----------|-------------|---------|----------|
| `--disease` | ICD disease codes | C81-C96 | C81-C96, C50, C34 |
| `--compound` | Pesticide compounds | 2 | 2, cat21, Atrazine |
| `--estimate` | Exposure estimates | avg | min, avg, max |
| `--measure` | Exposure measures | Weight | Weight, Density |
| `--model` | Model variants | M0 | M0, M1, M2, M3 |
| `--lag` | Exposure lag years | 5 | 5, 10 |
| `--sampling-mode` | MCMC configuration | test | test, production |
| `--cores` | Number of cores | auto | 4, 8, auto |
| `--chains` | Number of chains | auto | 4, 8, auto |
| `--draws` | Draws per chain | - | 1000, 2000 |
| `--tune` | Tuning steps | - | 1000, 2000 |
| `--target-accept` | Target acceptance | - | 0.8, 0.9 |
| `--output-dir` | Output directory | Result/PyMC_Results | custom/path |
| `--config-path` | Config file path | config.yaml | custom/config.yaml |
| `--verbose` | Verbose output | - | Flag |
| `--dry-run` | Preview parameters | - | Flag |

### **SLURM Script Parameters**

```bash
# Environment variables (edit pymc_task.sh)
DISEASE_CODE="C81-C96"
COMPOUND="2,9,cat21,cat33"  
ESTIMATE_TYPE="avg,max"
MODEL_TYPE="M0,M3"
LAG_YEARS="5,10"
MEASURE_TYPE="Weight"

# Sampling configuration  
DRAWS="4000"
TUNE="2000"
CHAINS="4" 
CORES="4"
TARGET_ACCEPT="0.95"

# Command line override
sbatch pymc_task.sh --compound "2,cat21" --estimate "min,max"
```

## 🔍 **Quality Control and Diagnostics**

### **Data Validation**
- **Integrity Checker**: `Code/Clean/CDC_Data_Integrity_Checker.py`
- **Coverage Validation**: Year ranges, geographic completeness
- **Variable Consistency**: Column presence, data type validation

### **Model Diagnostics**
- **Convergence**: R-hat statistics, trace plots
- **Effective Sample Size**: Bulk and tail ESS metrics  
- **Model Comparison**: WAIC, LOO cross-validation
- **Posterior Predictive**: Goodness-of-fit assessment

### **Output Validation**
- **Parameter Bounds**: Biologically plausible relative risks
- **Statistical Significance**: Bayesian p-values with multiple testing awareness
- **Spatial Patterns**: Residual spatial autocorrelation checks

## 🚨 **Common Issues and Solutions**

### **Environment Setup**
```bash
# PyMC installation issues
conda install -c conda-forge pymc::pymc=5.8.0

# R-INLA installation  
install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable")
```

### **Memory Management**
```bash
# Reduce batch size for large analyses
python main.py --measure Weight  # instead of Weight,Density

# Increase SLURM memory allocation
#SBATCH --mem-per-cpu=4G  # in pymc_task.sh
```

### **Convergence Issues**
```python
# Increase sampling parameters
--draws 6000 --tune 3000 --target-accept 0.99

# Check spatial matrix connectivity
assert nx.is_connected(adjacency_graph)
```

## 📚 **References and Documentation**

### **Methodological References**
1. Riebler A, et al. (2016). An intuitive Bayesian spatial model for disease mapping. *Statistical Methods in Medical Research* 25(4):1145-1165.
2. Salvatier J, et al. (2016). Probabilistic programming in Python using PyMC3. *PeerJ Computer Science* 2:e55.
3. Rue H, et al. (2009). Approximate Bayesian inference for latent Gaussian models. *Journal of the Royal Statistical Society B* 71(2):319-392.

### **Data Source Documentation**
- **CDC WONDER**: https://wonder.cdc.gov/
- **USGS PNSP**: https://water.usgs.gov/nawqa/pnsp/
- **NLDAS**: https://ldas.gsfc.nasa.gov/nldas/
- **CACES**: https://www.caces.us/data

## 👥 **Contributors and Contact**

**Development Team**: WDP Analysis Team  
**Institution**: Environmental Health Research  
**Last Updated**: September 28, 2025  
**Version**: 2.2

### **Getting Help**
- Check `config.yaml` for path configurations
- Use `--dry-run` to preview parameter resolution
- Review model diagnostics in result CSVs
- Consult convergence statistics for MCMC issues

### **Contributing**
1. Follow existing code structure and naming conventions
2. Update `config.yaml` for new data sources or paths
3. Document changes in commit messages
4. Test both PyMC and INLA pipelines for compatibility

### 1) Location (Static County-Level Geography)

- **Cleaning script**: `Code/Clean/CDC_Location_Urbanization.py`
- **Input directory**: `Data/Original/ CDC WONDER/Location and Urbanization`
  - **Files used**:
    - `Location_HHS_State.csv` (HHS regions)
    - `Location_Region_Division_State.csv` (Census regions/divisions)
- **Output file**: `Data/Processed/CDC/Location.csv`
- **Granularity**: County-level (one row per county; static across years)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `County`: County name with state abbreviation
  - `HHS_Region`: HHS region (e.g., “HHS Region #1 …”)
  - `Census_Region`: Census region (e.g., Northeast/West)
  - `Census_Division`: Census division (e.g., New England)

**Use**: Static join table for any county-level panel; join via `COUNTY_FIPS`.

---

### 2) Urbanization (County × Year Panel)

- **Cleaning script**: `Code/Clean/CDC_Location_Urbanization.py`
- **Input directory**: `Data/Original/ CDC WONDER/Location and Urbanization`
  - **Files used**: `Location_County_Urbanization*.csv` (spanning 1999–2020)
- **Output file**: `Data/Processed/CDC/Urbanization.csv`
- **Granularity**: County × Year (panel)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (1999–2020)
  - `County`: County name with state abbreviation
  - `Urbanization_Code`: Code from “2013 Urbanization Code”
  - `Urbanization_Type`: Category label (e.g., Large Central Metro, NonCore)

**Processing rules**:
- For duplicate county–year records, keep the first occurrence.
- Year is coerced to integer; records with missing year are dropped.

**Use**: Yearly urbanization classification for county-level analyses; join via `COUNTY_FIPS` and `Year`.

---

### 3) Cancer Mortality (ICD Group Merge)

- **Cleaning script**: `Code/Clean/CDC_Death_Merge.py`
- **Input base directory**: `Data/Original/ CDC WONDER/`
  - **Example**: `Data/Original/ CDC WONDER/C00-C97`
- **Output file**: `Data/Processed/CDC/<ICD>.csv` (e.g., `C00-C97.csv`)
- **Granularity**: County × Year (panel; 1999–2019)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int, filtered to 1999–2019)
  - `County`: County name with state abbreviation
  - `Sex`: Sex code
  - `Race`: Race category
  - `Age`: Ten-year age group
  - `Deaths`: Count of deaths (int)
  - `Population`: Population (int)
  - `SD`: Crude rate standard error

**Processing rules and notes**:
- Paths are resolved from `config.yaml` (`data_sources.cdc_wonder.integrity_base_dir`).
- Set `MANUAL_ICD_GROUP` at the top of the script to select the disease group directory.
- Merges all CSVs in the specified ICD folder.

---

### 4) AAMR (Age-Adjusted Mortality Rate, County × Year, Aggregated)

- **Cleaning script**: `Code/Clean/CDC_AAMR_Merge.py`
- **Input directory**: `Data/Original/CDC WONDER AAMR/{ICD}`（例如 `C81-C96`）
- **Output file**: `Data/Processed/CDC/AAMR_{ICD}.csv`（例如 `AAMR_C81-C96.csv`）
- **Granularity**: County × Year（1999–2020）
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int)
  - `County`: County name with state abbreviation
  - `Deaths`: Death count (Int64，可空)
  - `Population`: Population (Int64，可空)
  - `CMR`: Crude Mortality Rate
  - `CMR_SE`: Standard error of CMR
  - `AAMR`: Age-adjusted mortality rate
  - `AAMR_SE`: Standard error of AAMR

**Processing rules**:
- Drop non-county totals and invalid FIPS（移除诸如全国合计 `00nan` 等记录）。
- Standardize FIPS to 5-digit string；`Deaths`/`Population` 保存为可空整型（缺失保留为 NaN）。
- Year coerced to integer; filter to 1999–2020.
- Deduplicate county-year by preferring non-missing SE; then non-missing deaths.

**Use**: Aggregated outcome for spatio-temporal Bayesian smoothing（策略三），在 `Code/Test/BSTM_Run.py` 中可通过 `--aggregated` 使用。

---

### 5) CDC Data Integrity Checker

- **Utility script**: `Code/Clean/CDC_Data_Integrity_Checker.py`
- **Input directory**: Configurable via `MANUAL_ICD_GROUPS` (e.g., "C00-C14")
- **Purpose**: Validates data completeness and structure across CDC WONDER files
- **Function**: Checks for required columns, data ranges, and file consistency

**Key features**:
- Validates presence of key columns (Year, County, County Code, Sex, Race, Age)
- Checks statistical columns (Deaths, Population, Crude Rate Standard Error)
- Reports year ranges and missing data patterns
- Compact summary mode for efficient reporting
- Configurable to check single or multiple ICD disease groups

**Use**: Run before `CDC_Death_Merge.py` to identify data quality issues

---

## Pesticide (USGS PNSP)

### 1) PNSP Compound Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/PNSP_Merge.py`
- **Input directory**: `Data/Original/USGS PNSP`
  - **Files used**: Individual year files (1999-2012, 2018-2019) + combined file (2013-2017)
  - **Mapping file**: `Data/Processed/Pesticide/mapping.csv`
- **Output file**: `Data/Processed/Pesticide/PNSP.csv`
- **Granularity**: County × Year (panel; 1999–2019)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int, 1999–2019)
  - `cat{id}_min/avg/max`: Category-level estimates in kg/county/year (35 categories)
  - `chem{id}_min/avg/max`: Compound-specific estimates in kg/county/year (509 compounds)

**Processing rules and notes**:
- Paths are read from `config.yaml`.
- Merges all PNSP files from 1999-2019.
- Uses `mapping.csv` to map compound names to IDs and categories.

### 2) PNSP Pesticide Density Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/PNSP_Density.py`
- **Input files**:
  - `Data/Processed/Pesticide/PNSP.csv` (Pesticide weights)
  - `Data/Processed/Environmental/NLCD_JRC.csv` (Agricultural land area)
- **Output file**: `Data/Processed/Pesticide/PNSP_Density.csv`
- **Granularity**: County × Year (panel; 1999–2019)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int, 1999–2019)
  - `cat{id}_min/avg/max`: Category-level estimates in kg/km²/year
  - `chem{id}_min/avg/max`: Compound-specific estimates in kg/km²/year

**Processing rules and notes**:
- Calculates pesticide density by dividing the application weight from `PNSP.csv` by the agricultural land area (`nlcd_agriculture_km2`) from `NLCD_JRC.csv`.
- This provides a measure of exposure intensity that accounts for the size of a county's agricultural footprint.
- If agricultural area is zero or missing, the resulting density is zero.

### 3) Compound Name Mapping and Utility Scripts

- **Utility script**: `Code/Clean/PNSP_Unique_Name.py`
- **Input directory**: `Data/Original/USGS PNSP` (all .txt and .csv files)
- **Output file**: `Data/Processed/Pesticide/Unique_Name.txt`
- **Purpose**: Extracts unique pesticide compound names from all PNSP files.

**File**: `Data/Processed/Pesticide/mapping.csv`
- **Contains**: 509 compounds with classification into 35 categories.
- **Used by**: `PNSP_Merge.py` for data reshaping and categorization.

---

## Socioeconomic

### 1) SAIPE Poverty and Income Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/SE_SAIPE_Merge.py`
- **Input directory**: `Data/Original/SAIPE`
- **Output file**: `Data/Processed/Socioeconomic/Poverty_Income.csv`
- **Granularity**: County × Year (panel; 1999–2019)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int)
  - `Poverty_Percent_All_Ages`: Poverty rate for all ages (%)
  - `Median_Household_Income`: Median household income (dollars)

### 2) LAUS Labor Force and Unemployment Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/SE_LAUS_Merge.py`
- **Input directory**: `Data/Original/LAUS`
- **Output file**: `Data/Processed/Socioeconomic/Unemployment.csv`
- **Granularity**: County × Year (panel; 1999–2019)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int)
  - `Labor_Force`: Total labor force (integer)
  - `Employed`: Number of employed persons (integer)
  - `Unemployed`: Number of unemployed persons (integer)
  - `Unemployment_Rate`: Unemployment rate (%)

### 3) BEA Economic Indicators Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/SE_BEA_Merge.py`
- **Input directory**: `Data/Original/BEA`
- **Output file**: `Data/Processed/Socioeconomic/GDP.csv`
- **Granularity**: County × Year (panel; 1999–2019)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int)
  - `Population`: Population count (integer)
  - `Total_GDP_10K_USD`: Total GDP in 10,000 USD (float)
  - `Per_Capita_Income`: Per capita income (integer)

### 4) USDA ERS Education Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/SE_USDA_ERS_Education_Merge.py`
- **Input directory**: `Data/Original/USDA ERS`
- **Output file**: `Data/Processed/Socioeconomic/Education.csv`
- **Granularity**: County × Year (panel; 1999–2020, interpolated)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int)
  - `Less_Than_High_School_Percent`: Less than high school education (%)
  - `High_School_Only_Percent`: High school diploma only (%)
  - `Some_College_Percent`: Some college or associate degree (%)
  - `College_Plus_Percent`: Four-year college or higher (%)
  - `Rural_Urban_Continuum_Code`: Rural-urban continuum classification code
  - `Urban_Influence_Code`: Urban influence classification code

### 5) SEER Population Structure Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/SE_SEER_Population.py`
- **Input directory**: `Data/Original/SEER Population`
  - **File used**: `us.1990_2023.singleages.through89.90plus.adjusted.txt`
- **Output file**: `Data/Processed/Socioeconomic/Population_Structure.csv`
- **Granularity**: County × Year (panel; 1999–2020)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int, 1999–2020)
  - `Race`: Numeric code (1: White, 2: Black, 3: American Indian/Alaska Native, 4: Asian or Pacific Islander)
  - `Origin`: Numeric code (0: Non-Hispanic, 1: Hispanic)
  - `Sex`: Numeric code (1: Male, 2: Female)
  - `Age`: Numeric code (0-89 for single years, 90 for 90+)
  - `Population`: Population count (integer)

### 6) County Adjacency Data (Spatial Relationships)

- **Cleaning script**: `Code/Clean/County_Adjacency.py`
- **Input directory**: `Data/Original/County Shapeline`
  - **File used**: `tl_2015_us_county.shp` (county shapefile)
- **Output files**: 
  - `Data/Processed/Socioeconomic/County_Adjacency_Matrix.csv`
  - `Data/Processed/Socioeconomic/County_Adjacency_List.csv`
- **Granularity**: County spatial relationships (static)
- **Variables**:
  - **Matrix format** (`County_Adjacency_Matrix.csv`):
    - Rows and columns: County GEOIDs (5-digit FIPS codes)
    - Values: Boolean (True = counties are adjacent, False = not adjacent)
  - **Edge list format** (`County_Adjacency_List.csv`):
    - `county_from`: Origin county GEOID (string)
    - `county_to`: Destination county GEOID (string) 
    - `adjacency_weight`: Boolean adjacency indicator (True for adjacent counties)

**Processing rules and notes**:
- Uses GeoPandas to determine spatial adjacency based on touching boundaries
- Matrix is symmetric (if county A is adjacent to county B, then B is adjacent to A)
- Self-adjacency is set to False (counties are not adjacent to themselves)
- Edge list contains only True adjacency relationships (18,962 relationships total)
- Used for spatial modeling in Bayesian analysis (ICAR/BYM models)

---

## Environmental

### 1) CACES LUR Air Pollution Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/ENV_LUR_Merge.py`
- **Input directory**: `Data/Original/CACES LUR`
  - **Files used**: 
    - `1999-2019 O3.csv`
    - `1999-2020 CO SO2 NO2 PM10 PM25.csv`
- **Output file**: `Data/Processed/Environmental/Air_Pollution.csv`
- **Granularity**: County × Year (panel; 1999–2019)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int)
  - `O3`: Ozone concentration (pred_wght)
  - `CO`: Carbon Monoxide concentration (pred_wght)
  - `SO2`: Sulfur Dioxide concentration (pred_wght)
  - `NO2`: Nitrogen Dioxide concentration (pred_wght)
  - `PM10`: Particulate Matter < 10µm concentration (pred_wght)
  - `PM25`: Particulate Matter < 2.5µm concentration (pred_wght)

**Processing rules**:
- Removes geographical coordinates (lat, lon) and state abbreviations
- Converts from long format (pollutant column) to wide format (one column per pollutant)
- Filters years to 1999-2019 for consistency with other datasets

---

### 2) GEE Land Cover and Surface Water Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/ENV_GEE_Merge.py`
- **Input directory**: `Data/Original/GEE`
  - **Files used**: 
    - `county_base*.csv` (county area data)
    - `jrc_water_*.csv` (JRC surface water by year)
    - `nlcd_landuse_*.csv` (NLCD land cover by year)
- **Output file**: `Data/Processed/Environmental/NLCD_JRC.csv`
- **Granularity**: County × Year (panel; 1999–2020, interpolated)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int)
  - `total_area_km2`: Total county area in square kilometers
  - `jrc_permanent_water_km2`: Area of permanent surface water (JRC)
  - `jrc_seasonal_water_km2`: Area of seasonal surface water (JRC)
  - `nlcd_forest_km2`: Forested area (NLCD)
  - `nlcd_water_km2`: Water area (NLCD)
  - `nlcd_urban_km2`: Urban/developed area (NLCD)
  - `nlcd_agriculture_km2`: Agricultural area (NLCD)
  - `nlcd_cropland_km2`: Cropland area (NLCD)
  - `nlcd_pasture_km2`: Pasture/hay area (NLCD)
  - `nlcd_wetland_km2`: Wetland area (NLCD)
  - `nlcd_shrub_km2`: Shrub/scrub area (NLCD)
  - `nlcd_grassland_km2`: Grassland/herbaceous area (NLCD)
  - `nlcd_barren_km2`: Barren land area (NLCD)

**Processing rules**:
- Merges county base data with JRC water and NLCD land cover data
- Creates complete county-year panel for 1999-2020
- Applies linear interpolation to fill missing yearly data
- Values are rounded to 4 decimal places and clipped to non-negative

**GEE Data Collection Scripts** (Google Earth Engine):
- `Code/Clean/ENV_GEE_County.py`: Exports county base data with total area
- `Code/Clean/ENV_GEE_JRC.py`: Exports JRC water data by year (1999-2020)
- `Code/Clean/ENV_GEE_NLCD.py`: Exports NLCD land cover data by available years
- `Code/Clean/ENV_GEE_utils.py`: Utility functions for GEE processing

---

### 3) NLDAS Meteorology Data (County × Year Panel)

- **Cleaning script**: `Code/Clean/ENV_NLDAS_Meterology.py`
- **Input directories**:
  - `Data/Original/NLDAS` (NLDAS NetCDF files)
  - `Data/Original/County Shapeline/` (county boundaries)
- **Output file**: `Data/Processed/Environmental/NLDAS_{START_YEAR}_{END_YEAR}.csv`
- **Granularity**: County × Year (panel; 1999–2019)
- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `year`: Calendar year (int)
  - `tas_mean_annual/DJF/MAM/JJA/SON`: Mean temperature (°C)
  - `wind_mean_annual/DJF/MAM/JJA/SON`: Mean wind speed (m/s)
  - `prcp_sum_annual/DJF/MAM/JJA/SON`: Total precipitation (mm/month)
  - `rh_mean_annual/DJF/MAM/JJA/SON`: Mean relative humidity (%)
  - `swrad_mean_annual/DJF/MAM/JJA/SON`: Mean shortwave radiation (W/m²)
  - `lwrad_mean_annual/DJF/MAM/JJA/SON`: Mean longwave radiation (W/m²)
  - `psurf_mean_annual/DJF/MAM/JJA/SON`: Mean surface pressure (kPa)
  - `cape_mean_annual/DJF/MAM/JJA/SON`: Mean CAPE (J/kg)
  - `potevap_sum_annual/DJF/MAM/JJA/SON`: Total potential evaporation (mm/month)

**Processing rules**:
- Processes 0.125° NLDAS gridded data to county-level aggregates
- Uses spatial weight matrix to map grid cells to counties
- Calculates annual and seasonal (DJF/MAM/JJA/SON) statistics
- Temperature converted from Kelvin to Celsius
- Wind speed calculated from U and V components
- Relative humidity calculated from specific humidity, temperature, and pressure

---

## Conventions

- **Keys**:
  - Use `COUNTY_FIPS` (5-digit string) for county joins.
  - Panel datasets include `Year` (int).
- **Paths**:
  - Original: `Data/Original/...`
  - Processed: `Data/Processed/...`
- **Cleaning assumptions**:
  - Drop all-empty rows at read time.
  - Standardize FIPS to 5-digit string.

---

## Maintenance Notes

After each cleaning script is stabilized, please update this README with:
- Data source description, coverage window, key variables and units.
- Cleaning rules and assumptions.
- Output structure (granularity, row counts if helpful) and join keys to other tables.

## Recent Updates

**Additional Utility Scripts Available**:
- `Code/Clean/CDC_Data_Integrity_Checker.py`: Validates CDC data completeness
- `Code/Clean/PNSP_Unique_Name.py`: Extracts unique pesticide compound names
- `Code/Clean/ENV_GEE_*.py`: Google Earth Engine data collection scripts
- `Code/Test/BSTM_Run.py`: Bayesian spatio-temporal model runner
- `Code/Test/Data_Loading.py`: Model data loading utilities

**File Organization**:
- All cleaning scripts follow consistent path resolution via `config.yaml`
- Environmental data processing includes multiple GEE utility scripts
- Test/modeling scripts are organized in `Code/Test/` directory
- Model results are systematically stored in `Result/Test/` with timestamped summaries

---

## PCA-Derived Covariates (SVI and Climate Factors)

### 1) Master Covariate File

- **Generation script**: `Code/PCA/PCA_Systematic.py`
- **Output file**: `Data/Processed/PCA/Master_Covariates.csv`
- **Granularity**: County × Year (panel; 1999–2020)
- **Description**: This file contains the final, model-ready covariates derived from the systematic PCA workflow. It is the recommended source of covariates for all subsequent analyses, including the Bayesian models.

- **Variables**:
  - `COUNTY_FIPS`: 5-digit county FIPS (string)
  - `Year`: Calendar year (int)
  - `SVI_PCA`: **Socioeconomic Vulnerability Index**. A composite index where higher values indicate higher socioeconomic vulnerability (e.g., higher poverty, higher unemployment, lower education). Derived from a VIF-screened set of socioeconomic variables.
  - `Climate_Factor_1`, `Climate_Factor_2`, etc.: **Composite Climate Factors**. Orthogonal (uncorrelated) climate indices derived from a VIF-screened set of annual and seasonal meteorological variables. The interpretation of each factor (e.g., "Warm & Dry" vs. "Cool & Wet") can be found in the PCA diagnostic tables in `Result/Tables/PCA_Systematic/`.
  - Other demographic and land use variables included for direct use in modeling.

**Use**: This is the primary source file for covariates in all modeling stages. Join with outcome data via `COUNTY_FIPS` and `Year`.