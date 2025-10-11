#!/usr/bin/env Rscript
# ==============================================================================
# EQI Change → Cancer Mortality Change: Three-Phase Bayesian Analysis
# ==============================================================================
# 
# This script implements a comprehensive Bayesian analysis examining the 
# relationship between changes in Environmental Quality Index (EQI) and 
# changes in cancer mortality rates using interval-censored mixed models.
#
# THREE-PHASE ANALYSIS STRATEGY:
# ------------------------------------------------------------------------------
# Phase 1: NATIONAL LEVEL (全国层面)
#   - Model 1.1: Overall EQI change effect
#   - Model 1.2: Single domain effects (Air, Water, Land, Built, Social separately)
#   - Model 1.3: Multi-domain model (all domains simultaneously)
#
# Phase 2: STRATIFIED BY RUCC (城乡分层)
#   - Repeat all Phase 1 models separately for RUCC=1,2,3,4
#   - Tests whether effects differ by urbanization level
#
# Phase 3: WITHIN-RUCC RELATIVE CHANGE (同类层面相对排名变化)
#   - Model 3.1: Overall EQI relative ranking change within same RUCC category
#   - Model 3.2: Domain-specific relative ranking changes within same RUCC category
#   
#   说明：Phase 3分析的是县级单位在其**同一城乡化层级（RUCC）内部**的相对排名变化。
#   例如：某县在RUCC=1内从环境质量排名靠后变为靠前，这种"相对改善"对癌症死亡率的影响。
#   与Phase 1/2不同，这里使用的是RUCC_EQI_Change_Category（同类内相对变化），而非EQI_Change_Category（绝对变化）。
#
# USAGE EXAMPLES:
# ------------------------------------------------------------------------------
# # All cancer types (full analysis)
# Rscript Code/brms/Delta_bayesian.R
#
# # Specific cancer types
# Rscript Code/brms/Delta_bayesian.R --cancer-types="C00_C97,C18_C21"
#
# # Test mode (reduced iterations)
# Rscript Code/brms/Delta_bayesian.R --test --cancer-types="C00_C97"
#
# # Custom output directory
# Rscript Code/brms/Delta_bayesian.R --output-dir="Result/brms_delta_test"
#
# OUTPUT FORMAT:
# ------------------------------------------------------------------------------
# CSV with 9 columns (no extra domain columns):
#   - Model: Overall_EQI, Air, Water_Multi, RUCC1_Built, Within_RUCC_Social, etc.
#   - ICD_Code: C00_C97, C18_C21, etc.
#   - Lag: 5 or 10
#   - MRD_Q_Improved: Effect of Improved category (Mean(95%CI)*)
#   - MRD_Q_Worsened: Effect of Worsened category
#   - Intercept_Baseline_Change: Baseline delta AAMR
#   - Control_delta_Smoking_Rate: Smoking control effect
#   - N_Counties: Sample size
#   - Rhat_max: Maximum Rhat for convergence check
#
# Note: Multi-Domain models output 5 rows (Air_Multi, Water_Multi, Land_Multi,
#       Built_Multi, Social_Multi) instead of 1 row with extra columns.
#
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(cmdstanr)
  library(posterior)
})

# Suppress global variable warnings
utils::globalVariables(c(
  'COUNTY_FIPS', 'State', 'RUCC', 'Cancer_Type', 'Lag',
  'delta_AAMR_lower', 'delta_AAMR_upper', 'delta_Smoking_Rate',
  'EQI_Change_Category', 'Air_Change_Category', 'Water_Change_Category',
  'Land_Change_Category', 'Built_Change_Category', 'Social_Change_Category',
  'RUCC_EQI_Change_Category', 'RUCC_EQI_Air_Change_Category', 
  'RUCC_EQI_Water_Change_Category', 'RUCC_EQI_Land_Change_Category',
  'RUCC_EQI_Built_Change_Category', 'RUCC_EQI_Social_Change_Category'
))

# ==============================================================================
# COMMAND LINE ARGUMENTS
# ==============================================================================

option_list <- list(
  make_option(c("--data"), type="character", 
              default="Data/Processed/df_EQI_AAMR/EQI_AAMR_Delta.csv", 
              help="Input delta data path"),
  make_option(c("--output-dir"), type="character", 
              default="Result/brms_delta", 
              help="Output directory for results"),
  make_option(c("--cancer-types"), type="character", 
              default=NA, 
              help="Comma separated cancer codes (e.g., 'C00_C97,C18_C21'). If not provided, all will be analyzed."),
  make_option(c("--lag"), type="integer", 
              default=NA, 
              help="Lag years to analyze: 5, 10, or NA for both (default: both)"),
  make_option(c("--chains"), type="integer", 
              default=4, 
              help="Number of MCMC chains"),
  make_option(c("--iter"), type="integer", 
              default=2000, 
              help="Total iterations per chain"),
  make_option(c("--warmup"), type="integer", 
              default=1000, 
              help="Warmup iterations per chain"),
  make_option(c("--adapt-delta"), type="double", 
              default=0.95, 
              help="Target acceptance rate"),
  make_option(c("--max-treedepth"), type="integer", 
              default=12, 
              help="Maximum tree depth"),
  make_option(c("--min-n"), type="integer", 
              default=50, 
              help="Minimum sample size for overall analysis"),
  make_option(c("--min-n-rucc"), type="integer", 
              default=30, 
              help="Minimum sample size for RUCC-stratified analysis"),
  make_option(c("--seed"), type="integer", 
              default=1234, 
              help="Random seed for reproducibility"),
  make_option(c("--test"), action="store_true", 
              default=FALSE, 
              help="Test mode with reduced iterations")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Test mode: reduce iterations
if (opt$test) {
  opt$iter <- min(opt$iter, 800)
  opt$warmup <- min(opt$warmup, 300)
  message("[TEST MODE] Reduced iterations: iter=", opt$iter, " warmup=", opt$warmup)
}

set.seed(opt$seed)

# Configure parallel processing (auto-detect and use 80%)
cores_avail <- parallel::detectCores(logical=TRUE)
cores_used <- max(1, floor(cores_avail * 0.8))
options(mc.cores = cores_used)
message("🖥️  Detected cores: ", cores_avail, " | Auto-configured to use: ", cores_used, " (", 
        round(cores_used/cores_avail*100, 0), "%)")

# ==============================================================================
# STAN MODEL DEFINITION
# ==============================================================================
# Interval-censored mixed model with state-level random effects
# Identical to cmdstan_interval_runner.R but adapted for delta outcomes

stan_code <- "
data {
  int<lower=1> N;              // Number of observations
  int<lower=1> S;              // Number of states
  array[N] int<lower=1,upper=S> state;  // State indicator
  vector[N] y_lower;           // Lower bound of delta AAMR
  vector[N] y_upper;           // Upper bound of delta AAMR
  array[N] int<lower=0,upper=2> cens;   // Censoring indicator (0=exact, 2=interval)
  int<lower=1> K;              // Number of predictors
  matrix[N,K] X;               // Design matrix
}
parameters {
  vector[K] beta;              // Fixed effects coefficients
  vector[S] z_u;               // Non-centered state random effects
  real<lower=0> sigma;         // Residual standard deviation
  real<lower=0> sigma_u;       // State random effect standard deviation
}
transformed parameters {
  vector[S] u = sigma_u * z_u; // Centered state random effects
}
model {
  // Priors
  beta ~ normal(0, 5);
  z_u ~ normal(0, 1);
  sigma ~ exponential(1);
  sigma_u ~ exponential(1);
  
  // Likelihood
  for (i in 1:N) {
    real mu = X[i] * beta + u[state[i]];
    if (cens[i] == 0) {
      // Exact observation
      target += normal_lpdf(y_lower[i] | mu, sigma);
    } else {
      // Interval censored
      real p_upper = normal_cdf(y_upper[i] | mu, sigma);
      real p_lower = normal_cdf(y_lower[i] | mu, sigma);
      real diff = fmax(p_upper - p_lower, 1e-12);
      target += log(diff);
    }
  }
}
"

# Compile Stan model
stan_file <- file.path(tempdir(), "delta_interval_mixed_model.stan")
writeLines(stan_code, stan_file)
message("Compiling Stan model...")
mod <- cmdstan_model(stan_file)
message("✓ Stan model compiled successfully")

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

# Create concise model name from analysis components
create_model_name <- function(analysis_phase, rucc_stratum, model_type, domain_analyzed) {
  # Base name for different model types
  base_name <- switch(model_type,
    "Overall_EQI" = "EQI",
    "Single_Domain" = domain_analyzed,
    "Multi_Domain" = "Multi",
    "Unknown"
  )
  
  # Add prefix based on analysis phase
  if (analysis_phase == "Stratified") {
    return(paste0("RUCC", rucc_stratum, "_", base_name))
  } else if (analysis_phase == "Within_RUCC") {
    return(paste0("Within_RUCC_", base_name))
  } else {
    # National: no prefix
    return(base_name)
  }
}

# Significance marker based on p-value
sig_mark <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("")
}

# Format posterior draws as "Mean(95%CI)***"
format_cell <- function(draws) {
  if (length(draws) == 0 || all(is.na(draws))) return("")
  ci <- quantile(draws, c(0.025, 0.975), na.rm=TRUE)
  # Two-sided p-value
  p <- 2 * min(mean(draws > 0, na.rm=TRUE), mean(draws < 0, na.rm=TRUE))
  sprintf("%.2f(%.2f,%.2f)%s", mean(draws, na.rm=TRUE), ci[1], ci[2], sig_mark(p))
}

# Append results to CSV (create if doesn't exist)
# Note: Always quote fields to handle commas in confidence intervals
append_result <- function(result_df, output_path) {
  if (!file.exists(output_path)) {
    suppressWarnings(write.table(
      result_df, output_path, 
      sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE
    ))
  } else {
    suppressWarnings(write.table(
      result_df, output_path, 
      sep=",", col.names=FALSE, row.names=FALSE, append=TRUE, quote=TRUE
    ))
  }
}

# Wrapper function to fit model and save results with error handling
fit_and_save <- function(data, design_fn, model_type, analysis_phase, rucc_stratum,
                         domain, outcome_var, lag_years, output_file, ...) {
  result <- fit_delta_model(data, design_fn, ...)
  effects <- extract_change_effects(result$fit, result$design_info, model_type)
  row <- format_output_row(effects, analysis_phase, rucc_stratum, model_type,
                          domain, outcome_var, lag_years, result$design_info$n_counties)
  append_result(row, output_file)
  return(effects$Rhat_max)
}

# ==============================================================================
# DATA PREPARATION FUNCTIONS
# ==============================================================================

# Build design matrix for Model 1.1: Overall EQI change
build_design_overall_eqi <- function(d) {
  # Ensure Change_Category is a factor with Stable as reference
  d$EQI_Change_Category <- factor(
    d$EQI_Change_Category, 
    levels = c("Stable", "Improved", "Worsened")
  )
  
  # Create model matrix: Intercept + delta_Smoking_Rate + EQI_Change (2 dummies)
  X <- model.matrix(~ delta_Smoking_Rate + EQI_Change_Category, data=d)
  
  list(
    X = X,
    coef_names = colnames(X),
    n_counties = nrow(d)
  )
}

# Build design matrix for Model 1.2: Single domain
build_design_single_domain <- function(d, domain) {
  # domain: "Air", "Water", "Land", "Built", "Social"
  col_name <- paste0(domain, "_Change_Category")
  
  if (!col_name %in% names(d)) {
    stop("Column not found: ", col_name)
  }
  
  # Ensure factor with Stable as reference
  d[[col_name]] <- factor(
    d[[col_name]], 
    levels = c("Stable", "Improved", "Worsened")
  )
  
  # Create formula dynamically
  formula_str <- paste0("~ delta_Smoking_Rate + ", col_name)
  X <- model.matrix(as.formula(formula_str), data=d)
  
  list(
    X = X,
    coef_names = colnames(X),
    n_counties = nrow(d)
  )
}

# Build design matrix for Model 1.3: Multi-domain (all 5 domains)
build_design_multi_domain <- function(d) {
  # Ensure all domain Change_Category variables are factors
  domain_cols <- c(
    "Air_Change_Category", "Water_Change_Category", "Land_Change_Category",
    "Built_Change_Category", "Social_Change_Category"
  )
  
  for (col in domain_cols) {
    if (!col %in% names(d)) {
      stop("Column not found: ", col)
    }
    d[[col]] <- factor(d[[col]], levels = c("Stable", "Improved", "Worsened"))
  }
  
  # Create model matrix with all domains
  X <- model.matrix(
    ~ delta_Smoking_Rate + Air_Change_Category + Water_Change_Category + 
      Land_Change_Category + Built_Change_Category + Social_Change_Category,
    data = d
  )
  
  list(
    X = X,
    coef_names = colnames(X),
    n_counties = nrow(d)
  )
}

# Build design matrix for Model 3.1: Within-RUCC overall
build_design_rucc_overall <- function(d) {
  d$RUCC_EQI_Change_Category <- factor(
    d$RUCC_EQI_Change_Category,
    levels = c("Stable", "Improved", "Worsened")
  )
  
  X <- model.matrix(~ delta_Smoking_Rate + RUCC_EQI_Change_Category, data=d)
  
  list(
    X = X,
    coef_names = colnames(X),
    n_counties = nrow(d)
  )
}

# Build design matrix for Model 3.2: Within-RUCC single domain
build_design_rucc_single_domain <- function(d, domain) {
  # domain: "Air", "Water", "Land", "Built", "Social"
  col_name <- paste0("RUCC_EQI_", domain, "_Change_Category")
  
  if (!col_name %in% names(d)) {
    stop("Column not found: ", col_name)
  }
  
  d[[col_name]] <- factor(
    d[[col_name]],
    levels = c("Stable", "Improved", "Worsened")
  )
  
  formula_str <- paste0("~ delta_Smoking_Rate + ", col_name)
  X <- model.matrix(as.formula(formula_str), data=d)
  
  list(
    X = X,
    coef_names = colnames(X),
    n_counties = nrow(d)
  )
}

# ==============================================================================
# RESULT EXTRACTION FUNCTIONS
# ==============================================================================

# Extract change effects from posterior draws
extract_change_effects <- function(fit, design_info, model_type) {
  # Get posterior draws as data frame
  draws_df <- as_draws_df(fit$draws())
  
  # Extract coefficients based on model type
  coef_names <- design_info$coef_names
  
  # Initialize result list
  result <- list(
    Intercept = numeric(0),
    MRD_Improved = numeric(0),
    MRD_Worsened = numeric(0),
    Control_Smoking = numeric(0)
  )
  
  # Extract Intercept
  if ("(Intercept)" %in% coef_names) {
    idx <- which(coef_names == "(Intercept)")
    result$Intercept <- draws_df[[paste0("beta[", idx, "]")]]
  }
  
  # Extract Improved effect (look for pattern "Improved")
  improved_col <- coef_names[grepl("Improved", coef_names)]
  if (length(improved_col) > 0) {
    idx <- which(coef_names == improved_col[1])
    result$MRD_Improved <- draws_df[[paste0("beta[", idx, "]")]]
  }
  
  # Extract Worsened effect (look for pattern "Worsened")
  worsened_col <- coef_names[grepl("Worsened", coef_names)]
  if (length(worsened_col) > 0) {
    idx <- which(coef_names == worsened_col[1])
    result$MRD_Worsened <- draws_df[[paste0("beta[", idx, "]")]]
  }
  
  # Extract smoking control
  smoking_col <- coef_names[grepl("delta_Smoking_Rate", coef_names)]
  if (length(smoking_col) > 0) {
    idx <- which(coef_names == smoking_col[1])
    result$Control_Smoking <- draws_df[[paste0("beta[", idx, "]")]]
  }
  
  # For multi-domain models, we need to extract multiple domain effects
  if (model_type == "Multi_Domain") {
    domains <- c("Air", "Water", "Land", "Built", "Social")
    for (dom in domains) {
      # Improved
      improved_col <- coef_names[grepl(paste0(dom, ".*Improved"), coef_names)]
      if (length(improved_col) > 0) {
        idx <- which(coef_names == improved_col[1])
        result[[paste0("MRD_", dom, "_Improved")]] <- draws_df[[paste0("beta[", idx, "]")]]
      }
      # Worsened
      worsened_col <- coef_names[grepl(paste0(dom, ".*Worsened"), coef_names)]
      if (length(worsened_col) > 0) {
        idx <- which(coef_names == worsened_col[1])
        result[[paste0("MRD_", dom, "_Worsened")]] <- draws_df[[paste0("beta[", idx, "]")]]
      }
    }
  }
  
  # Extract Rhat
  rhat_vals <- fit$summary(variables="beta")$rhat
  result$Rhat_max <- max(rhat_vals, na.rm=TRUE)
  
  return(result)
}

# Format results into output row(s)
# For Multi-Domain models, returns 5 rows (one per domain)
# For other models, returns a single-row data frame
format_output_row <- function(effects, analysis_phase, rucc_stratum, model_type,
                               domain_analyzed, outcome_var, lag_years, n_counties) {
  
  icd_code <- gsub("^delta_AAMR_", "", outcome_var)
  
  # Helper function to create a single result row
  create_row <- function(model_name, improved_draws, worsened_draws) {
    data.frame(
      ICD_Code = icd_code,                          # Column 1: ICD first
      Model = model_name,                            # Column 2: Model name
      Lag = lag_years,
      MRD_Q_Improved = format_cell(improved_draws),
      MRD_Q_Worsened = format_cell(worsened_draws),
      Intercept_Baseline_Change = format_cell(effects$Intercept),
      Control_delta_Smoking_Rate = format_cell(effects$Control_Smoking),
      N_Counties = n_counties,
      Rhat_max = round(effects$Rhat_max, 4),
      stringsAsFactors = FALSE
    )
  }
  
  # For Multi-Domain: create 5 rows (one per domain)
  if (model_type == "Multi_Domain") {
    domains <- c("Air", "Water", "Land", "Built", "Social")
    prefix <- if (analysis_phase == "Stratified") paste0("RUCC", rucc_stratum, "_")
              else if (analysis_phase == "Within_RUCC") "Within_RUCC_"
              else ""
    
    rows <- lapply(domains, function(dom) {
      model_name <- paste0(prefix, dom, "_Multi")
      improved <- effects[[paste0("MRD_", dom, "_Improved")]]
      worsened <- effects[[paste0("MRD_", dom, "_Worsened")]]
      create_row(model_name, improved, worsened)
    })
    
    return(do.call(rbind, rows))
  } else {
    # For non-multi-domain: single row
    model_name <- create_model_name(analysis_phase, rucc_stratum, model_type, domain_analyzed)
    return(create_row(model_name, effects$MRD_Improved, effects$MRD_Worsened))
  }
}

# ==============================================================================
# MODEL FITTING WRAPPER
# ==============================================================================

fit_delta_model <- function(data_subset, design_function, ...) {
  # Prepare design matrix
  design_info <- design_function(data_subset, ...)
  
  # Prepare Stan data
  data_subset$State_numeric <- as.integer(factor(data_subset$State))
  n_states <- length(unique(data_subset$State))
  
  # Determine censoring indicator
  data_subset$cens <- ifelse(
    data_subset$delta_AAMR_lower == data_subset$delta_AAMR_upper, 0, 2
  )
  
  stan_data <- list(
    N = nrow(data_subset),
    S = n_states,
    state = data_subset$State_numeric,
    y_lower = data_subset$delta_AAMR_lower,
    y_upper = data_subset$delta_AAMR_upper,
    cens = data_subset$cens,
    K = ncol(design_info$X),
    X = design_info$X
  )
  
  # Fit model
  fit <- mod$sample(
    data = stan_data,
    chains = opt$chains,
    iter_sampling = opt$iter - opt$warmup,
    iter_warmup = opt$warmup,
    adapt_delta = opt$`adapt-delta`,
    max_treedepth = opt$`max-treedepth`,
    seed = opt$seed,
    refresh = 0,
    show_messages = FALSE
  )
  
  return(list(fit = fit, design_info = design_info))
}

# ==============================================================================
# MAIN ANALYSIS PIPELINE
# ==============================================================================

run_analysis <- function() {
  # Load data
  project_root <- normalizePath(".")
  data_path <- file.path(project_root, opt$data)
  
  if (!file.exists(data_path)) {
    stop("❌ Data file not found: ", data_path)
  }
  
  message("📊 Loading data from: ", data_path)
  dt <- fread(data_path)
  
  # Validate required columns
  required_cols <- c(
    "COUNTY_FIPS", "State", "RUCC", "Cancer_Type", "Lag",
    "delta_AAMR_lower", "delta_AAMR_upper", "delta_Smoking_Rate",
    "EQI_Change_Category", "Air_Change_Category", "Water_Change_Category",
    "Land_Change_Category", "Built_Change_Category", "Social_Change_Category",
    "RUCC_EQI_Change_Category", "RUCC_EQI_Air_Change_Category",
    "RUCC_EQI_Water_Change_Category", "RUCC_EQI_Land_Change_Category",
    "RUCC_EQI_Built_Change_Category", "RUCC_EQI_Social_Change_Category"
  )
  
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("❌ Missing required columns: ", paste(missing_cols, collapse=", "))
  }
  
  # Filter complete cases
  dt <- dt[complete.cases(dt[, ..required_cols])]
  message("✓ Loaded ", nrow(dt), " complete observations")
  
  # Filter by lag if specified
  if (!is.na(opt$lag)) {
    if (!opt$lag %in% c(5, 10)) {
      stop("❌ --lag must be 5 or 10")
    }
    dt <- dt[Lag == opt$lag]
    message("✓ Filtered to Lag=", opt$lag, ": ", nrow(dt), " observations")
  } else {
    message("✓ Analyzing both Lag=5 and Lag=10")
  }
  
  # Determine cancer types to analyze
  all_cancers <- sort(unique(dt$Cancer_Type))
  if (is.na(opt$`cancer-types`)) {
    selected_cancers <- all_cancers
  } else {
    selected_cancers <- str_split(opt$`cancer-types`, ",", simplify=TRUE) %>%
      as.vector() %>%
      str_trim()
    invalid <- setdiff(selected_cancers, all_cancers)
    if (length(invalid) > 0) {
      stop("❌ Invalid cancer types: ", paste(invalid, collapse=", "))
    }
  }
  
  message("🎯 Cancer types to analyze: ", paste(selected_cancers, collapse=", "))
  
  # Setup output directory
  out_dir <- file.path(project_root, opt$`output-dir`)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive=TRUE)
    message("📁 Created output directory: ", out_dir)
  }
  
  # Main loop over cancer types
  for (cancer in selected_cancers) {
    message("\n", paste(rep("=", 70), collapse=""))
    message("🔬 Analyzing: ", cancer)
    message(paste(rep("=", 70), collapse=""))
    
    cancer_data <- dt[Cancer_Type == cancer]
    
    if (nrow(cancer_data) < opt$`min-n`) {
      message("⚠️  Skipping ", cancer, " (n=", nrow(cancer_data), " < min_n=", opt$`min-n`, ")")
      next
    }
    
    # Output file: {ICD_Code}_delta.csv
    icd_code <- gsub("^delta_AAMR_", "", cancer)
    output_file <- file.path(out_dir, paste0(icd_code, "_delta.csv"))
    outcome_var <- cancer  # Keep full Cancer_Type for internal use
    
    # Get unique lag values in this cancer's data
    lag_values <- sort(unique(cancer_data$Lag))
    message("📊 Lag values in data: ", paste(lag_values, collapse=", "))
    
    # Loop over each lag value
    for (lag_val in lag_values) {
      message("\n🕐 Analyzing Lag=", lag_val)
      
      # Filter to specific lag
      cancer_lag_data <- cancer_data[Lag == lag_val]
      
      if (nrow(cancer_lag_data) < opt$`min-n`) {
        message("  ⚠️  Skipping Lag=", lag_val, " (n=", nrow(cancer_lag_data), " < min_n=", opt$`min-n`, ")")
        next
      }
      
      lag_years <- lag_val
    
      # ========================================================================
      # PHASE 1: NATIONAL LEVEL ANALYSIS
      # ========================================================================
      message("  📍 Phase 1: National Level Analysis")
      domains <- c("Air", "Water", "Land", "Built", "Social")
      
      # Model 1.1: Overall EQI
      message("    → Model 1.1: Overall EQI")
      tryCatch({
        rhat <- fit_and_save(cancer_lag_data, build_design_overall_eqi, "Overall_EQI",
                            "National", "All", "Overall", outcome_var, lag_years, output_file)
        message("      ✓ Completed (Rhat_max=", round(rhat, 3), ")")
      }, error = function(e) message("      ✗ Failed: ", e$message))
      
      # Model 1.2: Single domains
      for (i in seq_along(domains)) {
        message("    → Model 1.2.", i, ": ", domains[i])
        tryCatch({
          rhat <- fit_and_save(cancer_lag_data, build_design_single_domain, "Single_Domain",
                              "National", "All", domains[i], outcome_var, lag_years, output_file, domains[i])
          message("      ✓ Completed (Rhat_max=", round(rhat, 3), ")")
        }, error = function(e) message("      ✗ Failed: ", e$message))
      }
      
      # Model 1.3: Multi-domain
      message("    → Model 1.3: Multi-domain")
      tryCatch({
        rhat <- fit_and_save(cancer_lag_data, build_design_multi_domain, "Multi_Domain",
                            "National", "All", "All", outcome_var, lag_years, output_file)
        message("      ✓ Completed (Rhat_max=", round(rhat, 3), ")")
      }, error = function(e) message("      ✗ Failed: ", e$message))
    
      # ========================================================================
      # PHASE 2: STRATIFIED BY RUCC
      # ========================================================================
      message("\n  📍 Phase 2: RUCC-Stratified Analysis")
      
      for (rucc_level in 1:4) {
        rucc_data <- cancer_lag_data[RUCC == rucc_level]
        
        if (nrow(rucc_data) < opt$`min-n-rucc`) {
          message("    ⚠️  Skipping RUCC=", rucc_level, " (n=", nrow(rucc_data), 
                  " < min_n=", opt$`min-n-rucc`, ")")
          next
        }
        
        message("    → RUCC=", rucc_level, " (n=", nrow(rucc_data), ")")
        rucc_str <- as.character(rucc_level)
        
        # Overall EQI
        tryCatch({
          fit_and_save(rucc_data, build_design_overall_eqi, "Overall_EQI",
                      "Stratified", rucc_str, "Overall", outcome_var, lag_years, output_file)
          message("      ✓ EQI")
        }, error = function(e) message("      ✗ EQI: ", e$message))
        
        # Single domains
        for (dom in domains) {
          tryCatch({
            fit_and_save(rucc_data, build_design_single_domain, "Single_Domain",
                        "Stratified", rucc_str, dom, outcome_var, lag_years, output_file, dom)
          }, error = function(e) message("      ✗ ", dom, ": ", e$message))
        }
        message("      ✓ Single domains")
        
        # Multi-domain
        tryCatch({
          fit_and_save(rucc_data, build_design_multi_domain, "Multi_Domain",
                      "Stratified", rucc_str, "All", outcome_var, lag_years, output_file)
          message("      ✓ Multi")
        }, error = function(e) message("      ✗ Multi: ", e$message))
      }
    
      # ========================================================================
      # PHASE 3: WITHIN-RUCC RELATIVE CHANGE
      # ========================================================================
      message("\n  📍 Phase 3: Within-RUCC Relative Change")
      
      # Overall EQI ranking change
      message("    → Within_RUCC_EQI")
      tryCatch({
        rhat <- fit_and_save(cancer_lag_data, build_design_rucc_overall, "Overall_EQI",
                            "Within_RUCC", "All", "Overall", outcome_var, lag_years, output_file)
        message("      ✓ Completed (Rhat_max=", round(rhat, 3), ")")
      }, error = function(e) message("      ✗ Failed: ", e$message))
      
      # Single domain ranking changes
      for (i in seq_along(domains)) {
        message("    → Within_RUCC_", domains[i])
        tryCatch({
          rhat <- fit_and_save(cancer_lag_data, build_design_rucc_single_domain, "Single_Domain",
                              "Within_RUCC", "All", domains[i], outcome_var, lag_years, output_file, domains[i])
          message("      ✓ Completed (Rhat_max=", round(rhat, 3), ")")
        }, error = function(e) message("      ✗ Failed: ", e$message))
      }
    } # End lag loop
    
    message("📂 Results saved to: ", output_file)
  } # End cancer loop
  
  message("\n", paste(rep("=", 70), collapse=""))
  message("✅ All analyses complete!")
  message("� Output directory: ", out_dir)
  message(paste(rep("=", 70), collapse=""))
}

# ==============================================================================
# EXECUTE
# ==============================================================================

if (sys.nframe() == 0) {
  run_analysis()
}
