# WDP Project Structure

This document provides an overview of the directory structure for the WONDER Data Pipeline (WDP) project.

## Project Layout

The project is organized into several key directories to separate code, data, and results.

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
│   │   └── Download_*.py        # Automated data download utilities
│   ├── INLA/                    # R-INLA based analysis scripts
│   │   ├── BYM_INLA_Production.R # Main production script for INLA models
│   │   ├── BYM.md               # Technical design document for the INLA analysis
│   │   ├── utils/               # Helper R scripts for data processing, modeling, etc.
│   │   └── ...
│   ├── PCA/
│   │   └── PCA_Systematic.py    # Systematic PCA and covariate generation
│   └── Test/                    # Analysis and modeling scripts (Legacy or Python-based)
│       ├── BSTM_Run.py          # Bayesian Spatio-Temporal Model runner
│       ├── Bayes_Spatio_Temporal_Model.py  # Core Bayesian modeling functions
│       └── Data_Loading.py      # Data loading utilities for modeling
├── Data/
│   ├── Original/                # Raw, unmodified data from various sources
│   │   ├──  CDC WONDER/         # Note the leading space in the directory name
│   │   ├── BEA/
│   │   ├── CACES LUR/
│   │   ├── County Shapeline/
│   │   ├── GEE/
│   │   ├── LAUS/
│   │   ├── NLDAS/
│   │   ├── SAIPE/
│   │   ├── SEER Population/
│   │   ├── USDA ERS/
│   │   └── USGS PNSP/
│   └── Processed/               # Cleaned, analysis-ready datasets
│       ├── CDC/
│       │   ├── AAMR_{ICD}.csv   # Aggregated AAMR by county-year (e.g., AAMR_C81-C96.csv)
│       │   ├── Location.csv     # Static county geography and classifications
│       │   ├── Urbanization.csv # County-year urbanization classifications
│       │   └── {ICD}.csv        # Raw mortality data by ICD group (e.g., C00-C97.csv)
│       ├── Environmental/
│       │   ├── Air_Pollution.csv    # CACES LUR air pollution data
│       │   ├── NLCD_JRC.csv         # GEE land cover and surface water data
│       │   └── NLDAS_*.csv          # NLDAS meteorological data
│       ├── Pesticide/
│       │   ├── PNSP.csv             # Merged USGS PNSP pesticide data (weight)
│       │   ├── PNSP_Density.csv     # Calculated pesticide density data
│       │   ├── mapping.csv          # Compound-to-category mapping
│       │   └── Unique_Name.txt      # List of all unique compound names
│       ├── PCA/
│       │   └── Master_Covariates.csv # PCA-derived SVI and Climate factors
│       └── Socioeconomic/
│           ├── County_Adjacency_Matrix.csv   # Spatial adjacency matrix
│           ├── County_Adjacency_List.csv     # Spatial adjacency edge list
│           ├── Education.csv                 # USDA ERS education data
│           ├── GDP.csv                       # BEA economic indicators
│           ├── Population_Structure.csv      # SEER population demographics
│           ├── Poverty_Income.csv            # SAIPE poverty and income data
│           └── Unemployment.csv              # LAUS unemployment data
├── Result/
│   ├── Figures/
│   │   └── PCA_Systematic/      # Figures from the systematic PCA workflow
│   ├── Filter/                  # Filtered data subsets
│   ├── INLA_Top10/              # Results from INLA model runs
│   ├── Logs/                    # Log files from script executions
│   ├── Tables/
│   │   └── PCA_Systematic/      # Tables from the systematic PCA workflow
│   └── Test/                    # Model results and diagnostics (Legacy)
│       ├── BSTM_Summary_*.yaml  # Bayesian model run summaries
│       └── {Disease}_{Pesticide}_{Model}_lag{X}/  # Model-specific results
├── config.yaml                  # Central configuration file for all paths
├── README.md                    # Main project README
├── README_Data.md               # Detailed guide to data sources and variables
└── README_Structure.md          # This file: project structure overview
```

## Core Modules and Directory Highlights

**Core Modules**:
- **Code/Clean/**: Contains all data cleaning scripts for each data source
- **Code/Download/**: Scripts to download raw data from external sources  
- **Code/Test/**: Bayesian modeling, analysis, and testing scripts
- **Data/Original/**: Raw data from external sources (organized by data source)
- **Data/Processed/**: Cleaned, structured datasets ready for analysis
- **Result/**: Output files including logs, figures, model results, and diagnostics
- **config.yaml**: Central configuration file for all file paths and settings

**Key Scripts by Category**:
- `Code/Clean/CDC_*.py`: CDC mortality, geography, and AAMR data processing
- `Code/Clean/ENV_*.py`: Environmental data processing (LUR, GEE, NLDAS) 
- `Code/Clean/SE_*.py`: Socioeconomic data (SAIPE, LAUS, BEA, USDA, SEER)
- `Code/Clean/PNSP_*.py`: Pesticide data processing and compound mapping
- `Code/Clean/County_Adjacency.py`: Generates adjacency matrix for spatial modeling
- `Code/PCA/PCA_Systematic.py`: Runs the full systematic PCA workflow to generate master covariates
- `Code/Test/BSTM_Run.py`: Entry point for Bayesian spatio-temporal modeling
- `Code/Test/Bayes_Spatio_Temporal_Model.py`: Core Bayesian model implementation
- `Code/Test/Data_Loading.py`: Data loading and preprocessing for modeling

**Analysis Output Directories**:
- `Result/Test/`: Bayesian model outputs, diagnostics, and summary files
- `Result/Test/BSTM_Summary_*.yaml`: Model run summaries with timestamps
- `Result/Test/{Disease}_{Pesticide}_{Model}_lag{X}/`: Model-specific results and diagnostics

## Configuration Management

All file paths and directory settings are managed centrally in the `config.yaml` file located in the project root (`WDP/`).

### Key Principles:

1.  **Centralized Paths**: All scripts should read from `config.yaml` to get the correct paths for data input and output. This avoids hard-coded paths and makes the project portable.
2.  **Relative Paths**: The configuration uses relative paths from the project root, making it easy to run the project from any machine without modification.
3.  **Clear Naming**: The configuration keys in `config.yaml` are named to clearly correspond to the data source or directory they represent (e.g., `data_directories`, `data_sources`). For AAMR processing, see `data_sources.cdc_wonder.aamr_original` and `data_sources.cdc_wonder.aamr_output_pattern`.

## How to Run Scripts

All scripts should be run from the project's root directory (`WDP/`).

### 1. Downloading Data

Download scripts are located in `Code/Download/`.

```bash
# Example: Download USGS PNSP data
python Code/Download/Download_PNSP.py
```

### 2. Cleaning Data

Data cleaning scripts are located in `Code/Clean/`.

```bash
# Example: Run the SEER Population cleaning script
python Code/Clean/SE_SEER_Population.py

# Example: Process CDC AAMR data (set ICD group in script)
python Code/Clean/CDC_AAMR_Merge.py

# Example: Check data integrity before processing
python Code/Clean/CDC_Data_Integrity_Checker.py

# Example: Generate spatial adjacency matrix
python Code/Clean/County_Adjacency.py
```

### 3. Covariate Generation (PCA)

Run the systematic PCA script to generate the master covariate file.

```bash
# Run the full systematic PCA workflow
python Code/PCA/PCA_Systematic.py
```

### 4. Analysis and Modeling

Analysis and modeling scripts are located in `Code/Test/`.

```bash
# Example: Run Bayesian Spatio-Temporal Model
python Code/Test/BSTM_Run.py --aggregated

# Example: Run with specific disease and pesticide categories
python Code/Test/BSTM_Run.py --disease C81-C96 --pesticide herbicide
```

The scripts will automatically find the `config.yaml` file, read the necessary paths, and save the output to the appropriate `Data/Processed/` subdirectory.