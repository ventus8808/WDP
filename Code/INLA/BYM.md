# Bayesian Spatiotemporal Model (BYM2) Implementation with R-INLA
# Production System Documentation and Status

## Executive Summary

This document describes the production implementation of Bayesian spatiotemporal models using R-INLA (Integrated Nested Laplace Approximation) for analyzing pesticide exposure and health outcomes. The system has been fully implemented and tested, supporting both local Docker execution and remote HPC server deployment.

**Current Status**: ✅ **Production Ready** (Last Updated: 2024-09-15)

---

## 1. System Architecture Overview

### 1.1 Core Components

- **Main Script**: `BYM_INLA_Production.R` - Production-level analysis orchestrator
- **Utility Modules**:
  - `utils/data_processing.R` - Data loading, transformation, and quality control
  - `utils/model_fitting.R` - INLA model specification and fitting
  - `utils/result_extraction.R` - Result extraction and formatting
  - `utils/dashboard_printing.R` - Real-time progress monitoring
- **Configuration**: `config/analysis_config.yaml` - Centralized parameter management
- **Containerization**: `Dockerfile.inla_arm64` - Docker image for reproducible environments

### 1.2 Model Specifications

#### Core Models
- **M0**: Base model (pesticide + spatial + temporal effects)
- **M1**: M0 + Socioeconomic covariates (SVI_PCA)
- **M2**: M0 + Environmental covariates (Climate_Factor_1, Climate_Factor_2)
- **M3**: M0 + Full covariates (all of the above)

#### Dose-Response Modeling
- **Linear**: Continuous log-transformed exposure with standardization
- **Non-linear**: Random walk (RW2) on binned exposure values
- **Flexible**: Configuration-driven selection per model type

### 1.3 Key Features

1. **Continuous Dose-Response Analysis**
   - Log transformation with robust outlier detection
   - Standardized exposure metrics (RR per SD)
   - P90 vs P10 comparisons for interpretability

2. **Performance Optimizations**
   - Pre-loaded data for multi-measure analyses
   - Sparse matrix operations for spatial components
   - Efficient memory management

3. **Production Features**
   - Real-time progress dashboard
   - Comprehensive error handling
   - Immediate result writing (crash-resistant)
   - Validation and quality checks

---

## 2. Mathematical Framework

### 2.1 Base Model Structure

For county $i$, year $t$:

$$O_{i,t} \sim \text{Poisson}(\lambda_{i,t})$$

$$\log(\lambda_{i,t}) = \log(E_{i,t}) + \eta_{i,t}$$

**Linear dose-response**:
$$\eta_{i,t} = \alpha + \beta \cdot \log(X_{i,t-lag} + c) + \sum_{j} \gamma_j Z_{j,i,t} + \phi_i + \theta_i + \tau_t$$

**Non-linear dose-response**:
$$\eta_{i,t} = \alpha + f(X_{binned,i,t-lag}) + \sum_{j} \gamma_j Z_{j,i,t} + \phi_i + \theta_i + \tau_t$$

Where:
- $O_{i,t}$ = Observed death counts
- $E_{i,t}$ = Expected death counts (age-standardized offset)
- $X_{i,t-lag}$ = Pesticide exposure with lag
- $c$ = Small constant for log transformation
- $f(\cdot)$ = RW2 smooth function for non-linear models
- $\phi_i$ = Structured spatial effect (CAR)
- $\theta_i$ = Unstructured spatial effect (IID)
- $\tau_t$ = Temporal effect (RW1)

### 2.2 Prior Specifications

- **Fixed effects**: Default INLA priors (proper but vague)
- **Spatial precision**: PC prior with P(σ > 1) = 0.01
- **Temporal precision**: PC prior with P(σ > 1) = 0.01
- **BYM2 mixing parameter**: PC prior centered at 0.5

---

## 3. Data Pipeline

### 3.1 Input Data Structure

**Required Files**:
- `Data/Processed/CDC/{disease_code}.csv` - Health outcome data
- `Data/Processed/Pesticide/PNSP.csv` - Pesticide weight data
- `Data/Processed/Pesticide/PNSP_Density.csv` - Pesticide density data
- `Data/Processed/PCA/Master_Covariates.csv` - Standardized covariates
- `Data/Processed/Socioeconomic/County_Adjacency_List.csv` - Spatial structure

### 3.2 Data Processing Steps

1. **Exposure Transformation**
   ```r
   # Log transformation with small constant
   c_constant <- min(exposure[exposure > 0]) / 2
   exposure_log <- log(exposure + c_constant)
   exposure_std <- scale(exposure_log)
   ```

2. **Outlier Detection**
   - IQR-based method on log scale
   - Flags outliers beyond Q1 - 3×IQR or Q3 + 3×IQR
   - Reports percentage of outliers

3. **Temporal Aggregation**
   - 5-year or 10-year rolling averages
   - Lag application before analysis period

---

## 4. Execution Environments

### 4.1 Local Docker Execution

**Quick Test**:
```bash
./Code/INLA/run_in_docker.sh
```

**Full Analysis**:
```bash
./Code/INLA/run_in_docker.sh "Rscript Code/INLA/BYM_INLA_Production.R \
  --pesticide-category compound:2 \
  --measure-type Weight,Density \
  --estimate-types min,avg,max \
  --lag-years 5 \
  --model-types M0,M1,M2,M3 \
  --disease-code C81-C96 \
  --verbose"
```

### 4.2 Remote HPC Server Execution

**Server Details**:
- Host: cancon.hpccube.com
- Port: 65023
- Username: acf4pijnzl
- Path: `/public/home/acf4pijnzl/WDP_Analysis`

**Deployment**:
```bash
# Deploy and start analysis
./Code/INLA/deploy_to_server_fixed.sh

# Monitor progress
./monitor_analysis.sh

# Download results
./download_results.sh
```

### 4.3 Resource Requirements

- **Memory**: 8-16 GB recommended
- **CPU**: 4+ cores for parallel spatial computations
- **Storage**: ~5 GB for data, ~2 GB for results
- **Time**: 2-5 minutes per model, 30-60 minutes for full run

---

## 5. Output Specifications

### 5.1 Result File Format

**Columns** (CSV format):
- `Timestamp`: Analysis timestamp
- `Disease`: ICD code
- `Exposure`: Compound/category name
- `Category`: Pesticide category
- `Measure`: Weight or Density
- `Estimate`: min/avg/max
- `Lag`: Lag years
- `Model`: M0/M1/M2/M3
- `Dose_Response_Type`: Linear/Non-linear_RW2
- `RR_Per_SD`: Relative risk per SD increase
- `RR_Per_SD_Lower/Upper`: 95% CI bounds
- `RR_P90_vs_P10`: 90th vs 10th percentile comparison
- `RR_P90_vs_P10_Lower/Upper`: 95% CI bounds
- `P_Value`: Significance with stars (***/**/* for p<0.001/0.01/0.05)
- `DIC`: Deviance Information Criterion
- `WAIC`: Watanabe-Akaike Information Criterion
- `N_Counties`: Number of counties in analysis
- `N_Records`: Total observations
- `Status_Message`: SUCCESS or error description

### 5.2 Output Locations

- **Local**: `Result/Filter/Results_{disease}_{exposure}_{timestamp}.csv`
- **Server**: `/public/home/acf4pijnzl/WDP_Results/{analysis_name}/`

---

## 6. Recent Updates (2024-09)

### 6.1 Continuous Dose-Response Implementation

- **Removed**: Arbitrary quintile-based RR calculations
- **Added**: Scientifically rigorous continuous exposure modeling
- **Impact**: Publication-ready statistical methodology

### 6.2 Performance Enhancements

- **Data Pre-loading**: Reduced I/O by ~80% for multi-measure analyses
- **Memory Optimization**: Efficient handling of large spatial matrices
- **Progress Monitoring**: Real-time dashboard with ETA

### 6.3 Server Deployment

- **Automated Scripts**: Full deployment pipeline for HPC execution
- **Monitoring Tools**: Remote progress tracking
- **Result Retrieval**: Automated download scripts

---

## 7. Usage Examples

### 7.1 Single Compound Analysis

```bash
# Analyze 2,4-D (compound ID: 2)
Rscript Code/INLA/BYM_INLA_Production.R \
  --pesticide-category compound:2 \
  --measure-type Density \
  --estimate-types avg \
  --lag-years 5 \
  --model-types M0,M1,M2,M3 \
  --disease-code C81-C96
```

### 7.2 Category Analysis

```bash
# Analyze all herbicides
Rscript Code/INLA/BYM_INLA_Production.R \
  --pesticide-category "Herbicide" \
  --measure-type Weight,Density \
  --estimate-types min,avg,max \
  --lag-years 5,10 \
  --model-types M0,M3
```

### 7.3 Comprehensive Analysis

```bash
# Full analysis with linear and non-linear models
./Code/INLA/run_full_24D_analysis.sh
```

---

## 8. Troubleshooting

### 8.1 Common Issues

1. **Memory Errors**
   - Reduce batch size in config
   - Use single measure type
   - Increase Docker memory allocation

2. **Convergence Warnings**
   - Check data completeness
   - Verify spatial connectivity
   - Review outlier percentages

3. **Missing Covariates**
   - Model automatically adjusts
   - Check warnings in output
   - Verify PCA file completeness

### 8.2 Validation Checks

- **Data Quality**: Check outlier reports
- **Model Fit**: Compare DIC/WAIC across models
- **Results**: Verify RR bounds are reasonable (0.1-10.0)
- **Significance**: Review p-value distribution

---

## 9. Future Enhancements

### 9.1 Planned Features

- [ ] Distributed computing support
- [ ] Interactive result visualization
- [ ] Automated report generation
- [ ] Multi-disease batch processing

### 9.2 Research Extensions

- [ ] Mixture distributions for zero-inflation
- [ ] Space-time interaction terms
- [ ] Multi-pollutant models
- [ ] Exposure uncertainty propagation

---

## 10. References

1. Rue H, Martino S, Chopin N (2009). Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. *J R Stat Soc Series B*, 71:319-392.

2. Simpson D, Rue H, Riebler A, et al (2017). Penalising model component complexity: A principled, practical approach to constructing priors. *Stat Sci*, 32:1-28.

3. Greenland S (1995). Dose-response and trend analysis in epidemiology: alternatives to categorical analysis. *Epidemiology*, 6:356-365.

---

## Contact and Support

For technical questions or issues:
- Review: `docs/WONDER_Code_Modifications_Summary.md`
- Check: `docs/Testing_Instructions.md`
- Scripts: All production scripts include `--help` option

**Last Production Run**: 2024-09-15 (2,4-D comprehensive analysis)  
**Validation Status**: ✅ All tests passing  
**Documentation Version**: 3.0 (2024-09-15)