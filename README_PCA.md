# WDP Principal Component Analysis (PCA) Documentation

This document provides detailed information about the Principal Component Analysis (PCA) workflow in the WDP project, including the variables used, generated components, and geographic coverage.

## Overview

The PCA analysis in WDP generates two sets of principal components from socio-economic and climate data:
1. **Social Vulnerability Index (SVI)** - Derived from socio-economic variables
2. **Climate Factors** - Derived from meteorological and environmental variables

These components are used as covariates in the Bayesian spatio-temporal modeling pipeline.

## Input Variables

### Social Vulnerability Index (SVI) Analysis

The following socio-economic variables are used in the SVI PCA:

| Variable Name | Description | Source |
|---------------|-------------|--------|
| Poverty_Percent_All_Ages | Percentage of population below poverty level | SAIPE |
| Median_Household_Income | Median household income in dollars | SAIPE |
| Unemployment_Rate | Unemployment rate | LAUS |
| Per_Capita_Income | Per capita income | BEA |
| Less_Than_High_School_Percent | Percentage of population age 25+ with less than high school education | USDA ERS |
| College_Plus_Percent | Percentage of population age 25+ with bachelor's degree or higher | USDA ERS |

### Environment Factors Analysis

The following climate and environmental variables are used in the Environment PCA:

| Variable Name | Description | Source |
|---------------|-------------|--------|
| tas_mean_annual | Annual mean air temperature | NLDAS |
| prcp_sum_annual | Annual total precipitation | NLDAS |
| wind_mean_annual | Annual mean wind speed | NLDAS |
| rh_mean_annual | Annual mean relative humidity | NLDAS |
| swrad_mean_annual | Annual mean shortwave radiation | NLDAS |
| potevap_sum_annual | Annual total potential evaporation | NLDAS |
| tas_mean_DJF | Mean air temperature for December-January-February | NLDAS |
| tas_mean_MAM | Mean air temperature for March-April-May | NLDAS |
| tas_mean_JJA | Mean air temperature for June-July-August | NLDAS |
| tas_mean_SON | Mean air temperature for September-October-November | NLDAS |
| prcp_sum_DJF | Total precipitation for December-January-February | NLDAS |
| prcp_sum_MAM | Total precipitation for March-April-May | NLDAS |
| prcp_sum_JJA | Total precipitation for June-July-August | NLDAS |
| prcp_sum_SON | Total precipitation for September-October-November | NLDAS |

## Generated Components

### Social Vulnerability Index (SVI)

The SVI PCA typically generates 1-2 principal components that capture the main dimensions of social vulnerability. The first component (SVI_PC1) is inverted so that higher values represent higher vulnerability.

### Environmental Factors

The Environmental PCA typically generates 2-3 principal components that capture the main dimensions of climate variability. These components represent combinations of temperature, precipitation, and other meteorological factors.

## Geographic Coverage

The PCA analysis covers all U.S. counties with valid data across the following time period:

- **Years**: 1999-2020
- **Geographic Units**: All U.S. counties (FIPS codes)
- **Spatial Resolution**: County-level

The analysis excludes:
- Counties with missing FIPS codes
- Data outside the 1999-2020 time range
- Counties with insufficient data for PCA computation

## Output Files

The PCA workflow generates the following key output files:

1. **Data/Processed/PCA/PCA_Master_Covariables.csv** - Master covariate file with PCA scores
2. **Result/Tables/PCA_Diagnose.csv** - Comprehensive diagnostics table with VIF values
3. **Result/Figure_Original_Data/** - Raw data for plotting
4. **Result/Figures/PCA_Analysis/** - Visualization outputs

## Analysis Pipeline

The PCA analysis follows these steps:

1. **Data Loading** - Load and merge all processed data files
2. **Variable Selection** - Apply VIF-based variable selection to remove multicollinearity
   - Iteratively remove variables with VIF > 10.0
   - Record VIF values for all variables (kept and removed)
3. **PCA Execution** - Perform PCA using the Kaiser criterion to determine components to keep
4. **Component Interpretation** - Analyze loadings to understand component meaning
5. **SVI Adjustment** - Invert SVI PC1 to align with vulnerability (higher = more vulnerable)
6. **Output Generation** - Create master covariate file and comprehensive diagnostic tables
7. **Plotting Data** - Generate raw data files for visualization scripts

## Diagnostics and Quality Control

The analysis includes several quality control measures and comprehensive diagnostics:

### Variable Selection
- VIF threshold of 10.0 to remove highly collinear variables
- Iterative variable removal based on VIF values
- Final VIF values reported in diagnostics table

### Component Selection
- Kaiser criterion to determine the number of components to retain
- Eigenvalue analysis for component interpretation
- Explained variance ratios for each component

### Data Quality
- Data filtering for the 1999-2020 time period
- Missing data handling through listwise deletion
- Geographic coverage validation

### Diagnostics Output
The **PCA_Diagnose.csv** file contains comprehensive diagnostics including:
- PCA Type (SVI or Climate)
- Component information (PC1, PC2, etc.)
- Eigenvalues and explained variance percentages
- Variable loadings on each component
- **VIF values** for all variables used in the analysis
- Cumulative variance explained