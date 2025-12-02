# Typology and LandUse Stratification Analysis

## Overview

This directory contains scripts for running interval-censored Bayesian mixed models stratified by **County Economic Typology** (USDA ERS 2004) and **Land Use Clusters** (NLCD/JRC analysis).

## Files

- **`cmdstan_Typology_LandUse.R`**: R script that performs the Bayesian analysis for both stratifications
- **`submit_cmdstan_Typology_LandUse.sh`**: SLURM array job submission script

## Stratification Variables

### 1. County Economic Typology (`econdep`)
Based on USDA ERS 2004 County Typology classification:

| Code | Type | Description |
|------|------|-------------|
| 1 | Farming | Farming-dependent counties |
| 2 | Mining | Mining-dependent counties |
| 3 | Manufacturing | Manufacturing-dependent counties |
| 4 | Government | Federal/State government-dependent counties |
| 5 | Services | Services-dependent counties |
| 6 | Nonspecialized | Nonspecialized counties |

### 2. Land Use Clusters (`landuse_cluster`)
Based on NLCD/JRC land use analysis (k=4 clustering):

| Code | Type | Description |
|------|------|-------------|
| 0 | Natural | Natural/forest-dominated counties |
| 1 | Water-Sensitive | Counties with significant water/wetland features |
| 2 | Agricultural | Agriculture-dominated counties |
| 3 | Urban | Urban-dominated counties |

## Input Data

**Required file**: `Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Cluster_Climate_Typology_LandUse.csv`

This file contains:
- County-level AAMR data (interval-censored)
- EQI quintiles (overall and domain-specific)
- Smoking rates
- County economic typology codes (`econdep`)
- Land use cluster assignments (`landuse_cluster`)

## Usage

### Local Testing

```bash
# Test on a single cancer type with reduced iterations
Rscript Code/brms/cmdstan_Typology_LandUse.R \
  --cancer-types "C00_C97" \
  --chains 2 \
  --iter 800 \
  --warmup 300 \
  --test
```

### SLURM Cluster Submission

The submission script automatically discovers all cancer types and submits one array job per cancer type (maximum 18 tasks). Each task runs both Typology and LandUse stratifications.

```bash
# Submit all cancer types
bash Code/brms/submit_cmdstan_Typology_LandUse.sh

# Check job status
squeue -u $USER

# Check specific job output
tail -f cmdstan_typology_landuse_<JOB_ID>_<TASK_ID>.out
```

## Analysis Workflow

For each cancer type, the script:

1. **Typology Stratification**:
   - Analyzes 6 strata (Farming, Mining, Manufacturing, Government, Services, Nonspecialized)
   - Runs overall EQI model + multi-domain model for each stratum
   - Outputs results to: `Result/brms_Typology_LandUse/<CANCER>_Typology.csv`

2. **LandUse Stratification**:
   - Analyzes 4 strata (Natural, Water-Sensitive, Agricultural, Urban)
   - Runs overall EQI model + multi-domain model for each stratum
   - Outputs results to: `Result/brms_Typology_LandUse/<CANCER>_LandUse.csv`

## Models Fitted

For each stratum × scenario combination:

1. **Overall EQI Model**: `AAMR ~ Smoking_Rate + EQI_Q2 + EQI_Q3 + EQI_Q4 + EQI_Q5 + (1|State)`
2. **Multi-Domain Models** (5 separate models):
   - `AAMR ~ Smoking_Rate + EQI_Air_Q2-Q5 + (1|State)`
   - `AAMR ~ Smoking_Rate + EQI_Water_Q2-Q5 + (1|State)`
   - `AAMR ~ Smoking_Rate + EQI_Land_Q2-Q5 + (1|State)`
   - `AAMR ~ Smoking_Rate + EQI_Built_Q2-Q5 + (1|State)`
   - `AAMR ~ Smoking_Rate + EQI_Social_Q2-Q5 + (1|State)`

## Output Format

Results are saved as CSV files with the following structure:

| Column | Description |
|--------|-------------|
| `ICD_Code` | Cancer type (ICD-10 code) |
| `EQI_Period` | EQI exposure period (2000_2005 or 2006_2010) |
| `AAMR_Period` | AAMR outcome period (2006_2010, 2011_2015, or 2016_2020) |
| `Lag` | Lag years between exposure and outcome (5, 10, or 15) |
| `Model` | Model identifier (e.g., "Typology_Farming_EQI", "LandUse_Urban_EQI_Air") |
| `Q1` | Reference level (always 0.00) |
| `Q2` | Effect estimate for quintile 2: `mean(95%CI)*` |
| `Q3` | Effect estimate for quintile 3 |
| `Q4` | Effect estimate for quintile 4 |
| `Q5` | Effect estimate for quintile 5 |
| `Q2_p` - `Q5_p` | Two-sided posterior p-values (4 decimals) |
| `Q2_rhat` - `Q5_rhat` | R-hat convergence diagnostics (4 decimals) |
| `Q2_ess_bulk` - `Q5_ess_bulk` | Effective sample size (bulk) |
| `Q2_ess_tail` - `Q5_ess_tail` | Effective sample size (tail) |

Significance markers: `***` p<0.001, `**` p<0.01, `*` p<0.05

## Scenarios

Each stratified analysis runs 5 EQI-AAMR lag scenarios:

1. **EQI 2000-2005 → AAMR 2006-2010** (5-year lag)
2. **EQI 2000-2005 → AAMR 2011-2015** (10-year lag)
3. **EQI 2000-2005 → AAMR 2016-2020** (15-year lag)
4. **EQI 2006-2010 → AAMR 2011-2015** (5-year lag)
5. **EQI 2006-2010 → AAMR 2016-2020** (10-year lag)

## Computational Requirements

- **CPU**: 16 cores per task (recommended)
- **Memory**: 48 GB per task
- **Time**: ~24 hours per cancer type (for all strata and scenarios)
- **Total**: 18 tasks × 24 hours = 432 core-hours (with parallelization: ~24-48 hours wall time)

## Resource Management

The SLURM script is designed to respect cluster limits:
- Maximum 18 concurrent array tasks (one per cancer type)
- Each task is self-contained and runs both stratifications
- Uses 80% of available cores for parallel chains
- Includes proper error handling and logging

## Output Directory

All results are saved to: **`Result/brms_Typology_LandUse/`**

Example output files:
```
Result/brms_Typology_LandUse/
├── C00_C97_Typology.csv       # All-site cancer by economic typology
├── C00_C97_LandUse.csv        # All-site cancer by land use
├── C34_Typology.csv           # Lung cancer by economic typology
├── C34_LandUse.csv            # Lung cancer by land use
└── ...
```

## Troubleshooting

### Common Issues

1. **"Data not found"**: Ensure `EQI_Triangulation_AAMR_df.py` has been run to generate the input CSV

2. **"Missing cols"**: Verify the input CSV contains `econdep` and `landuse_cluster` columns

3. **Convergence warnings**: Check R-hat values in output; values > 1.01 indicate convergence issues

4. **Memory errors**: Increase `--mem` in SLURM script or reduce `--chains`

### Debugging

```bash
# Run in test mode with verbose output
Rscript Code/brms/cmdstan_Typology_LandUse.R \
  --cancer-types "C00_C97" \
  --test \
  2>&1 | tee debug_output.log
```

## Contact

For questions or issues, refer to the main WDP project documentation or consult the project PI.

---

**Last Updated**: 2024-12-02
**Authors**: WDP Project Team