# Delta Model Ridgeline Extraction

## Overview

This workflow extracts posterior MCMC draws from delta models (EQI change vs. AAMR change) for ridgeline plot generation. It supports cluster stratification and multi-domain analysis.

## Files

- `cmdstan_delta_ridgeline.R` - R script for fitting delta models and extracting posterior draws
- `submit_delta_ridgeline.sh` - SLURM array submission script for batch processing

## Key Features

- **Delta modeling**: Analyzes relationship between EQI changes (improved/worsened) and cancer mortality changes
- **Cluster stratification**: Supports k=3 and k=4 county cluster analyses
- **Multi-lag support**: 5, 10, and 15-year lag periods
- **Two model types**:
  - **Overall**: Single overall EQI change effect (improved/worsened categories)
  - **Multi-domain**: Separate effects for Air, Water, Land, Built, and Social domains
- **Complete posterior extraction**: Saves all MCMC draws for ridgeline visualization

## Data Requirements

### Input Files

1. **Delta data** (default: `Data/Processed/df_EQI_AAMR/EQI_AAMR_Delta.csv`)
   - Required columns: `COUNTY_FIPS`, `Cancer_Type`, `Lag_Years`, `Delta_AAMR_Lower`, `Delta_AAMR_Upper`, `Delta_EQI`, `Delta_EQI_Air`, `Delta_EQI_Water`, `Delta_EQI_Land`, `Delta_EQI_Built`, `Delta_EQI_Social`, `Delta_Smoking_Rate`, `RUCC`, `EQI_Change_Category`

2. **Cluster assignments** (default: `Result/Cluster_Visualization/EQI_Clusters_All_K.csv`)
   - Required columns: `COUNTY_FIPS`, `cluster_3`, `cluster_4`

## Usage

### Single Run (Interactive)

Run a single combination manually:

```bash
Rscript Code/brms/cmdstan_delta_ridgeline.R \
  --cancer C00_C97 \
  --k 3 \
  --lag 5 \
  --model overall \
  --chains 4 \
  --iter 2000 \
  --warmup 1000 \
  --seed 1234
```

### Batch Processing (SLURM)

The submission script automatically generates and submits all task combinations:

```bash
# Auto-submit mode (discovers tasks and submits array)
bash Code/brms/submit_delta_ridgeline.sh
```

This generates 12 tasks for C00_C97:
- 2 k values (3, 4)
- 3 lag values (5, 10, 15)
- 2 model types (overall, multi)
- Total: 2 × 3 × 2 = 12 tasks

### Command-Line Options

```bash
--cancer          Cancer type [default: C00_C97]
--k               Number of clusters: 3 or 4 [default: 3]
--lag             Lag years: 5, 10, or 15 [default: 5]
--model           Model type: overall or multi [default: overall]
--data            Input data path [default: Data/Processed/df_EQI_AAMR/EQI_AAMR_Delta.csv]
--cluster-data    Cluster assignment path [default: Result/Cluster_Visualization/EQI_Clusters_All_K.csv]
--output-dir      Output directory [default: Result/Ridgeline_Delta]
--chains          Number of MCMC chains [default: 4]
--iter            Total iterations per chain [default: 2000]
--warmup          Warmup iterations [default: 1000]
--adapt-delta     Adapt delta for Stan [default: 0.95]
--max-treedepth   Max treedepth for Stan [default: 12]
--seed            Random seed [default: 1234]
--test            Test mode: reduced iterations
```

## Output Structure

### Output Files

Results are saved to `Result/Ridgeline_Delta/` as RDS files:

```
Result/Ridgeline_Delta/
├── C00_C97_k3_Lag5_OVERALL.rds
├── C00_C97_k3_Lag5_MULTI.rds
├── C00_C97_k3_Lag10_OVERALL.rds
├── C00_C97_k3_Lag10_MULTI.rds
├── C00_C97_k3_Lag15_OVERALL.rds
├── C00_C97_k3_Lag15_MULTI.rds
├── C00_C97_k4_Lag5_OVERALL.rds
├── C00_C97_k4_Lag5_MULTI.rds
├── C00_C97_k4_Lag10_OVERALL.rds
├── C00_C97_k4_Lag10_MULTI.rds
├── C00_C97_k4_Lag15_OVERALL.rds
└── C00_C97_k4_Lag15_MULTI.rds
```

### RDS Content

Each RDS file contains:

```r
list(
  draws = data.frame(
    iteration,      # MCMC iteration number
    draw,           # Posterior draw value
    parameter,      # Stan parameter name (e.g., "beta[1]")
    covariate,      # Covariate name (e.g., "Improved", "Worsened")
    effect_type,    # "Improved" or "Worsened" (or domain-specific)
    cluster,        # Cluster ID (e.g., "Cluster1")
    model_type,     # "overall" or "multi"
    cancer,         # Cancer type
    k,              # Number of clusters
    lag             # Lag years
  ),
  diagnostics = data.frame(
    variable,       # Parameter name
    mean, median, sd, mad, q5, q95,
    rhat, ess_bulk, ess_tail,
    cluster, model_type
  ),
  beta_mappings = list(...),    # Parameter-to-covariate mappings
  model_summaries = list(...),  # Per-cluster diagnostics
  metadata = list(
    cancer, k, lag, model_type,
    n_obs, n_clusters, clusters,
    chains, iter, warmup, adapt_delta, max_treedepth, seed,
    max_rhat, min_ess, timestamp
  )
)
```

## Loading and Using Results

```r
# Load results
result <- readRDS("Result/Ridgeline_Delta/C00_C97_k3_Lag5_OVERALL.rds")

# Access posterior draws
draws <- result$draws

# Filter for specific effect
improved_draws <- draws[draws$effect_type == "Improved", ]
worsened_draws <- draws[draws$effect_type == "Worsened", ]

# Check diagnostics
max(result$metadata$max_rhat)  # Should be < 1.01
min(result$metadata$min_ess)   # Should be > 100

# Create ridgeline plot (example with ggridges)
library(ggridges)
library(ggplot2)

ggplot(draws, aes(x = draw, y = cluster, fill = effect_type)) +
  geom_density_ridges(alpha = 0.7) +
  facet_wrap(~ covariate) +
  theme_minimal() +
  labs(x = "Posterior Draw", y = "Cluster",
       title = "Delta EQI Effects on Cancer Mortality Change")
```

## Model Specifications

### Overall EQI Model

Outcome: `Delta_AAMR` (interval-censored)

Predictors:
- `Improved`: Binary indicator (EQI_Change_Category == "Improved")
- `Worsened`: Binary indicator (EQI_Change_Category == "Worsened")
- `Delta_Smoking_Rate`: Change in smoking rate
- `RUCC`: Rural-Urban Continuum Code
- Random intercepts by state

### Multi-Domain Model

Outcome: `Delta_AAMR` (interval-censored)

Predictors:
- `Air_Improved`, `Air_Worsened`: Air domain changes
- `Water_Improved`, `Water_Worsened`: Water domain changes
- `Land_Improved`, `Land_Worsened`: Land domain changes
- `Built_Improved`, `Built_Worsened`: Built domain changes
- `Social_Improved`, `Social_Worsened`: Social domain changes
- `Delta_Smoking_Rate`: Change in smoking rate
- `RUCC`: Rural-Urban Continuum Code
- Random intercepts by state

## Diagnostics

The script monitors convergence:
- **Rhat**: Should be < 1.01 for all parameters
- **ESS (Effective Sample Size)**: Should be > 100 (preferably > 400)

Warnings are issued if diagnostics fail thresholds.

## Computational Requirements

- **Memory**: 48 GB recommended (32 GB minimum)
- **CPU**: 16 cores recommended for parallel chains
- **Time**: ~30-60 minutes per task (varies by cluster size and model complexity)
- **Storage**: ~50-200 MB per RDS file

## Workflow Integration

This workflow follows the project conventions:
- Reads paths from `config.yaml` (via relative resolution)
- Outputs to `Result/Ridgeline_Delta/` directory
- Uses cluster assignments from `Result/Cluster_Visualization/EQI_Clusters_All_K.csv`
- Compatible with existing delta analysis pipeline (`Delta_bayesian_Cluster.R`)

## Comparison with Standard Ridgeline Workflow

| Feature | Standard Ridgeline | Delta Ridgeline |
|---------|-------------------|-----------------|
| Input data | `EQI_AAMR_Cluster_Climate.csv` | `EQI_AAMR_Delta.csv` |
| Outcome | AAMR quintiles | Delta AAMR (change) |
| Predictors | EQI quintiles | EQI change (improved/worsened) |
| Stratification | None | Cluster-based (k=3, k=4) |
| Scenarios | 3 lag periods | 3 lag periods × 2 k values |

## Troubleshooting

### Common Issues

1. **"Cluster column not found"**
   - Ensure cluster data file contains `cluster_3` and `cluster_4` columns
   - Check that `--k` matches available cluster columns
   - Default location: `Result/Cluster_Visualization/EQI_Clusters_All_K.csv`

2. **"No data after filtering"**
   - Verify cancer type exists in delta data
   - Check that lag value is present (5, 10, or 15)

3. **High Rhat / Low ESS**
   - Increase iterations: `--iter 4000 --warmup 2000`
   - Increase adapt_delta: `--adapt-delta 0.99`
   - Check for data issues (extreme values, missing clusters)

4. **CmdStan compilation errors**
   - Ensure `TBB_CXX_TYPE=gcc` is set
   - Load devtoolset-8: `module load devtoolset-8`

## Related Scripts

- `Delta_bayesian_Cluster.R` - Main delta analysis (produces coefficient tables)
- `cmdstan_main_ridgeline.R` - Standard ridgeline extraction (non-delta models)
- `submit_Delta_bayesian_Cluster_array.sh` - Batch delta analysis submission

## Citation

If using this workflow, cite the parent WDP project and note the delta modeling approach for environmental health change analysis.