# BRMS Analysis Pipeline - Refactored Architecture

This directory contains a **refactored** Bayesian multilevel regression pipeline using the `brms` package in R. The architecture has been redesigned following the "responsibility inversion" principle for improved maintainability and robustness.

## Architecture Overview

### Before: Complex Orchestrator + Passive Engine
- Python reads data files to generate parameter combinations
- R script accepts 7 positional command-line arguments  
- Tight coupling and fragile interfaces

### After: Simple Orchestrator + Intelligent Engine
- Python only reads configuration and schedules scenarios
- R script receives scenario name and self-configures
- Clean separation of concerns and robust interfaces

## Files

- `01_prepare_data.R`: Data preparation script (unchanged)
- `02_run_brms_model.R`: **Refactored** intelligent modeling script  
- `03_main_runner.py`: **Refactored** simple orchestration script
- `04_process_results.R`: Post-processing script (unchanged)

## Key Improvements

### 1. Configuration-Driven Design
All scenario definitions are now centralized in `config.yaml`:

```yaml
brms_analysis:
  scenarios:
    - name: "LungCancer_TotalEQI_Lag5_AllRUCC"
      active: true
      cancer_type: "C34"
      eqi_period: "0005" 
      time_period: "2006-2010"
      lag_years: 5
      rucc_filter: null
      domain: "total"
      formula_key: "total_eqi_quintile"
      family: "gaussian"
```

### 2. Simple Command Interface
The R script now accepts a single parameter:

```bash
# Old (fragile): 7 positional parameters
Rscript 02_run_brms_model.R scenario_name prefix eqi_period time_period lag rucc domain

# New (robust): scenario name only
Rscript 02_run_brms_model.R --scenario LungCancer_TotalEQI_Lag5_AllRUCC
```

### 3. Intelligent Self-Configuration
The R script now:
- Reads `config.yaml` to find scenario configuration
- Automatically applies all necessary data filters
- Dynamically selects appropriate EQI columns based on domain/RUCC
- Handles fallback strategies for model fitting

## Usage

### Run All Active Scenarios
```bash
python 03_main_runner.py
```

### Run Individual Scenario  
```bash
Rscript 02_run_brms_model.R --scenario LungCancer_TotalEQI_Lag5_AllRUCC
```

### Add New Scenarios
Simply add to `config.yaml`:

```yaml
- name: "MyNewScenario"
  active: true
  cancer_type: "C50"
  eqi_period: "0610"
  time_period: "2011-2015" 
  lag_years: 3
  rucc_filter: [1, 2]  # Urban only
  domain: "air"
  formula_key: "total_eqi_quintile"
  family: "gaussian"
```

## Outputs

Results are saved in multiple formats:
- **Model objects**: `Result/brms/model_fits/`
- **CSV summaries**: `Result/brms/csv_outputs/`  
- **LMM-compatible**: `Result/brms/brms_{ICD}.csv`
- **Diagnostics**: Convergence and ESS metrics

## Benefits of Refactoring

1. **Maintainability**: Configuration changes don't require code modifications
2. **Robustness**: Single parameter interface eliminates argument order errors
3. **Transparency**: All scenario parameters explicitly defined in config
4. **Scalability**: Easy to add new scenarios without touching Python code
5. **Reproducibility**: Complete scenario specification in version control
