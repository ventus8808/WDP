# WDP INLA Analysis System - API Reference

## Overview

The WDP INLA Analysis System is a production-ready framework for Bayesian spatiotemporal analysis of pesticide exposure and health outcomes. This system has been completely restructured for better maintainability, robustness, and reproducibility.

## File Structure

```
Code/INLA/
├── INLA_Main.R                    # Main analysis script
├── INLA_Config/                   # Configuration files
│   └── analysis_config.yaml       # Primary configuration
├── INLA_Utils/                    # Utility functions
│   ├── INLA_Utils_Data.R          # Data processing utilities
│   ├── INLA_Utils_Model.R         # Model fitting utilities  
│   ├── INLA_Utils_Results.R       # Result extraction utilities
│   ├── INLA_Utils_Dashboard.R     # Output formatting utilities
│   ├── INLA_Utils_Validation.R    # Data validation utilities
│   └── INLA_Utils_Logger.R        # Logging and error handling
├── INLA_Scripts/                  # Execution scripts
│   └── run_single_compound.sh     # SLURM submission script
├── INLA_Temp/                     # Temporary files
│   ├── graphs/                    # Spatial graph files
│   ├── models/                    # Model cache files
│   └── logs/                      # Runtime logs
├── INLA_Dependencies/             # Dependency management
│   ├── install_packages.R         # Automatic package installer
│   ├── check_dependencies.R       # Dependency verification
│   └── R_packages/                # Local package cache
└── INLA_Docs/                     # Documentation
    ├── API_Reference.md           # This file
    ├── BYM_Model.md              # Model methodology
    ├── INLA_Installation_Guide.md # Complete INLA setup guide
    └── Troubleshooting.md        # Problem-solving guide
```

## Command Line Interface

### Basic Usage

```bash
Rscript INLA_Main.R [OPTIONS]
```

### Options

| Option | Description | Default | Examples |
|--------|-------------|---------|----------|
| `--config` | Path to configuration file | `INLA_Config/analysis_config.yaml` | Custom config file |
| `--disease-code` | Disease code to analyze | `C81-C96` | `C50`, `C34` |
| `--measure-type` | Exposure measure types | `Weight` | `Weight,Density` |
| `--pesticide-category` | Pesticide selection | `TEST` | `compound:1`, `compound:1,2,3` |
| `--estimate-types` | Exposure estimates | `avg` | `min,avg,max` |
| `--lag-years` | Lag periods in years | `5` | `5,10` |
| `--model-types` | Model complexity levels | `M0,M1,M2,M3` | `M0,M1` |
| `--verbose` | Enable detailed logging | `FALSE` | Flag, no value |
| `--dry-run` | Validate without running | `FALSE` | Flag, no value |

### Examples

#### Single compound analysis
```bash
Rscript INLA_Main.R \
  --disease-code "C81-C96" \
  --pesticide-category "compound:1" \
  --measure-type "Weight,Density" \
  --estimate-types "avg" \
  --lag-years "5,10" \
  --model-types "M0,M1,M2,M3" \
  --verbose
```

#### Quick validation check
```bash
Rscript INLA_Main.R \
  --disease-code "C50" \
  --pesticide-category "compound:5" \
  --dry-run
```

## Core Functions

### Data Processing (`INLA_Utils_Data.R`)

#### `load_all_data(config, disease_code, measure_type)`
Loads and validates all required data files using robust path management.

**Parameters:**
- `config`: Configuration list from YAML
- `disease_code`: Disease code (e.g., "C81-C96")
- `measure_type`: "Weight" or "Density"

**Returns:** List containing all loaded datasets

#### `calculate_lagged_exposure(pesticide_data, pesticide_col_name, lag_years, config)`
Calculates rolling mean exposure with specified lag period.

**Parameters:**
- `pesticide_data`: Pesticide exposure data frame
- `pesticide_col_name`: Column name for specific pesticide
- `lag_years`: Number of years for lag calculation
- `config`: Configuration list

**Returns:** Data frame with lagged exposure values

### Model Fitting (`INLA_Utils_Model.R`)

#### `create_spatial_structure(adjacency_data, counties_in_data, category_id, config)`
Creates spatial adjacency structure for INLA modeling.

**Parameters:**
- `adjacency_data`: County adjacency data frame
- `counties_in_data`: Vector of county FIPS codes
- `category_id`: Unique identifier for file naming
- `config`: Configuration list

**Returns:** List with INLA graph object and file path

#### `fit_inla_model(formula, model_data, config)`
Fits Bayesian spatial model using INLA.

**Parameters:**
- `formula`: Model formula object
- `model_data`: Prepared data frame
- `config`: Configuration list

**Returns:** Fitted INLA model object

### Validation (`INLA_Utils_Validation.R`)

#### `validate_input_data(data_list, config)`
Comprehensive validation of input data structure and completeness.

#### `check_data_completeness(model_data, min_counties, min_records)`
Checks if data meets minimum requirements for analysis.

#### `validate_model_requirements(model_data, formula, config)`
Validates data requirements before model fitting.

### Logging (`INLA_Utils_Logger.R`)

#### `setup_logger(log_level, log_file, console_output)`
Initializes logging system with specified configuration.

#### `log_error_with_trace(error, context, additional_info)`
Enhanced error logging with full stack trace for debugging.

## Configuration File

The system uses YAML configuration files for all settings. Key sections:

### Data Paths
All paths are relative to project root using `here` package:
```yaml
data_paths:
  pca_covariates: "Data/Processed/PCA/Master_Covariates.csv"
  cdc_data_template: "Data/Processed/CDC/{disease_code}.csv"
  pesticide_data: "Data/Processed/Pesticide/PNSP.csv"
  # ... other paths
```

### Output Configuration
```yaml
output:
  base_dir: "Result/INLA_Results"
  filename_template: "Results_{disease}_{exposure}_{timestamp}.csv"
```

### Model Settings
```yaml
model_fitting:
  spatial:
    model_type: "bym2"
    graph_dir: "Code/INLA/INLA_Temp/graphs"
```

## Error Handling

The system implements comprehensive error handling:

1. **Input Validation**: All inputs are validated before processing
2. **Data Completeness**: Minimum thresholds ensure reliable results
3. **Model Convergence**: Automatic checks for model fitting issues
4. **Resource Management**: Proper cleanup of temporary files
5. **Detailed Logging**: Complete error traces for debugging

## Environment Setup

### INLA Installation
For detailed INLA installation instructions, especially on HPC servers, see the complete guide:
```
INLA_Docs/INLA_Installation_Guide.md
```

### Dependencies
Run the automatic installer:
```bash
Rscript INLA_Dependencies/install_packages.R
```

### Verification
Check all dependencies:
```bash
Rscript INLA_Dependencies/check_dependencies.R
```

## Performance Considerations

- **Memory**: Minimum 16GB RAM recommended for large datasets
- **CPU**: Multi-threading support via INLA configuration
- **Storage**: Ensure sufficient space in temp directories
- **Network**: INLA package requires internet for initial installation

## Integration with SLURM

Use the optimized SLURM script:
```bash
sbatch INLA_Scripts/run_single_compound.sh C81-C96 1 Weight,Density avg 5,10 M0,M1,M2,M3
```

## Output Format

Results are saved as CSV files with the following columns:
- `Timestamp`: Analysis timestamp
- `Disease`, `Exposure`, `Category`: Analysis identifiers
- `Measure`, `Estimate`, `Lag`, `Model`: Parameter specifications
- `Q1`, `Q2`, `Q3`, `Q4`, `Q5`: Quintile relative risks
- `P_Value`: Statistical significance
- `N_Counties`, `N_Records`: Sample size information

## Version History

- **v2.0**: Complete restructure with robust path management
- **v1.x**: Original implementation (deprecated)

For detailed troubleshooting, see `Troubleshooting.md`.