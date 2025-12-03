# Visualization for Typology and LandUse Stratifications

## Overview

This directory contains visualization scripts for creating forest plots from County Economic Typology and Land Use Cluster stratified analysis results.

## Scripts

- **`Code/Visualization/Visualization_Typology.py`**: Creates forest plots for County Economic Typology stratifications
- **`Code/Visualization/Visualization_LandUse.py`**: Creates forest plots for Land Use Cluster stratifications

## Usage

### Typology Visualization

```bash
# Generate plots for a specific cancer type
python Code/Visualization/Visualization_Typology.py --icd C00_C97

# Generate plots for all available cancer types
python Code/Visualization/Visualization_Typology.py --all
```

### LandUse Visualization

```bash
# Generate plots for a specific cancer type
python Code/Visualization/Visualization_LandUse.py --icd C00_C97

# Generate plots for all available cancer types
python Code/Visualization/Visualization_LandUse.py --all
```

## Input Files

The scripts expect result files in `Result/brms_Typology_LandUse/`:
- Typology: `{ICD}_Typology.csv`
- LandUse: `{ICD}_LandUse.csv`

## Output

### Typology Plots
- **Location**: `Result/brms_Typology_Visualization/`
- **Format**: `{ICD}_Typology_{EQI_Period}_{AAMR_Period}_Lag{lag}.png`
- **Panels**: 6 panels (Farming, Mining, Manufacturing, Government, Services, Nonspecialized)

### LandUse Plots
- **Location**: `Result/brms_LandUse_Visualization/`
- **Format**: `{ICD}_LandUse_{EQI_Period}_{AAMR_Period}_Lag{lag}.png`
- **Panels**: 4 panels (Natural, Water-Sensitive, Agricultural, Urban)

## Plot Features

- Forest plot style with grayscale coloring
- Error bars with 95% credible intervals
- 5 scenarios per cancer type (combinations of EQI/AAMR periods and lags)
- 6 EQI domains per panel: Overall, Air, Water, Land, Built, Social
- 4 quintiles displayed per domain: Q2, Q3, Q4, Q5 (Q1 is reference)

## Example

```bash
# Generate all Typology plots
python Code/Visualization/Visualization_Typology.py --all

# Generate all LandUse plots
python Code/Visualization/Visualization_LandUse.py --all
```

This will create plots for all available cancer types across all 5 scenarios.

---

**Last Updated**: 2024-12-02
