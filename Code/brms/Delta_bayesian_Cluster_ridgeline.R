#!/usr/bin/env Rscript
# ==============================================================================
# EQI Change → Cancer Mortality Change: Ridgeline Plot Data Generator
# ==============================================================================
#
# This script runs the Bayesian analysis for Delta EQI vs Delta AAMR and
# extracts the full posterior draws for "Improved" and "Worsened" effects.
# The output is saved as .rds files (one per cancer type) containing lists
# of posterior draws for National and Cluster-stratified models.
#
# These .rds files are intended to be used for generating Ridgeline plots (Joyplots).
#
# USAGE:
#   Rscript Code/brms/Delta_bayesian_Cluster_ridgeline.R --cancer-types="C00_C97" --k="3"
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

# ==============================================================================
# COMMAND LINE ARGUMENTS
# ==============================================================================

option_list <- list(
  make_option(c("--data"), type="character",
              default="Data/Processed/df_EQI_AAMR/EQI_AAMR_Delta.csv",
              help="Input delta data path"),
  make_option(c("--output-dir"), type="character",
              default="Result/brms_delta_cluster_ridgeline",
              help="Output directory for .rds files"),
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
              help="Minimum sample size for analysis"),
  make_option(c("--k"), type="character",
              default="3,4",
              help="Comma separated k values to analyze (e.g., '3,4'). Default: 3,4"),
  make_option(c("--seed"), type="integer",
              default=1234,
              help="Random seed for reproducibility"),
  make_option(c("--test"), action="store_true",
              default=FALSE,
              help="Test mode with reduced iterations")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Test mode settings
if (opt$test) {
  opt$iter <- min(opt$iter, 800)
  opt$warmup <- min(opt$warmup, 300)
  message("[TEST MODE] Reduced iterations: iter=", opt$iter, " warmup=", opt$warmup)
}

set.seed(opt$seed)

# Parallel processing
cores_avail <- parallel::detectCores(logical=TRUE)
cores_used <- max(1, floor(cores_avail * 0.8))
options(mc.cores = cores_used)
message("🖥️  Detected cores: ", cores_avail, " | Using: ", cores_used)

# ==============================================================================
# STAN MODEL
# ==============================================================================

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
      target += normal_lpdf(y_lower[i] | mu, sigma);
    } else {
      real p_upper = normal_cdf(y_upper[i] | mu, sigma);
      real p_lower = normal_cdf(y_lower[i] | mu, sigma);
      real diff = fmax(p_upper - p_lower, 1e-12);
      target += log(diff);
    }
  }
}
"

stan_file <- file.path(tempdir(), "delta_interval_ridgeline.stan")
writeLines(stan_code, stan_file)
mod <- cmdstan_model(stan_file)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Build design matrix for Overall EQI
build_design_overall_eqi <- function(d) {
  d$EQI_Change_Category <- factor(d$EQI_Change_Category, levels = c("Stable", "Improved", "Worsened"))
  X <- model.matrix(~ delta_Smoking_Rate + EQI_Change_Category, data=d)
  list(X = X, coef_names = colnames(X), n_counties = nrow(d))
}

# Build design matrix for Single Domain
build_design_single_domain <- function(d, domain) {
  col_name <- paste0(domain, "_Change_Category")
  d[[col_name]] <- factor(d[[col_name]], levels = c("Stable", "Improved", "Worsened"))
  formula_str <- paste0("~ delta_Smoking_Rate + ", col_name)
  X <- model.matrix(as.formula(formula_str), data=d)
  list(X = X, coef_names = colnames(X), n_counties = nrow(d))
}

# Extract draws and summary statistics for Improved/Worsened
extract_draws_and_summary <- function(fit, design_info, model_label) {
  draws_df <- as_draws_df(fit$draws("beta"))
  coef_names <- design_info$coef_names

  # Identify indices for Improved and Worsened
  # Pattern: "Improved" and "Worsened" are appended to the variable name by model.matrix
  idx_imp <- grep("Improved", coef_names)
  idx_wor <- grep("Worsened", coef_names)

  if (length(idx_imp) == 0 || length(idx_wor) == 0) {
    warning("Could not find Improved/Worsened coefficients for ", model_label)
    return(NULL)
  }

  # Extract draws
  imp_draws <- draws_df[[paste0("beta[", idx_imp, "]")]]
  wor_draws <- draws_df[[paste0("beta[", idx_wor, "]")]]

  n_draws <- length(imp_draws)

  # Create long format draws for plotting
  draws_long <- data.frame(
    draw_id = rep(1:n_draws, 2),
    category = rep(c("Improved", "Worsened"), each = n_draws),
    effect = c(imp_draws, wor_draws)
  )

  # Calculate summary stats
  summ <- fit$summary("beta")

  # Helper to get stats for a specific index
  get_stats <- function(idx, label) {
    row <- summ[idx, ]
    data.frame(
      category = label,
      mean = row$mean,
      q2.5 = row$q5,   # posterior package default quantiles
      q97.5 = row$q95,
      rhat = row$rhat,
      ess_bulk = row$ess_bulk,
      ess_tail = row$ess_tail
    )
  }

  summary_df <- rbind(
    get_stats(idx_imp, "Improved"),
    get_stats(idx_wor, "Worsened")
  )

  list(
    draws = draws_long,
    summary = summary_df,
    n_obs = design_info$n_counties
  )
}

# Fit model wrapper
fit_model <- function(data_subset, design_fn, ...) {
  design_info <- design_fn(data_subset, ...)

  data_subset$State_numeric <- as.integer(factor(data_subset$State))
  n_states <- length(unique(data_subset$State))
  data_subset$cens <- ifelse(data_subset$delta_AAMR_lower == data_subset$delta_AAMR_upper, 0, 2)

  stan_data <- list(
    N = nrow(data_subset), S = n_states, state = data_subset$State_numeric,
    y_lower = data_subset$delta_AAMR_lower, y_upper = data_subset$delta_AAMR_upper,
    cens = data_subset$cens, K = ncol(design_info$X), X = design_info$X
  )

  fit <- mod$sample(
    data = stan_data, chains = opt$chains,
    iter_sampling = opt$iter - opt$warmup, iter_warmup = opt$warmup,
    adapt_delta = opt$`adapt-delta`, max_treedepth = opt$`max-treedepth`,
    seed = opt$seed, refresh = 0, show_messages = FALSE
  )

  list(fit = fit, design_info = design_info)
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

# Load Data
project_root <- normalizePath(".")
data_path <- file.path(project_root, opt$data)
if (!file.exists(data_path)) stop("Data file not found: ", data_path)

message("📊 Loading data...")
dt <- fread(data_path)

# Load Clusters
cluster_path <- file.path(project_root, "Result/Cluster_Visualization/EQI_Clusters_All_K.csv")
if (!file.exists(cluster_path)) stop("Cluster data not found: ", cluster_path)
cluster_dt <- fread(cluster_path)

# Parse K values
k_values <- as.integer(str_trim(unlist(str_split(opt$k, ","))))
cluster_cols <- paste0("cluster_", k_values)
dt <- merge(dt, cluster_dt[, c("COUNTY_FIPS", cluster_cols), with=FALSE], by = "COUNTY_FIPS", all.x = TRUE)

# Filter complete cases
req_cols <- c("COUNTY_FIPS", "State", cluster_cols, "Cancer_Type", "Lag",
              "delta_AAMR_lower", "delta_AAMR_upper", "delta_Smoking_Rate",
              "EQI_Change_Category", "Air_Change_Category", "Water_Change_Category",
              "Land_Change_Category", "Built_Change_Category", "Social_Change_Category")
dt <- dt[complete.cases(dt[, ..req_cols])]

# Filter Lag
if (!is.na(opt$lag)) dt <- dt[Lag == opt$lag]

# Select Cancers
all_cancers <- sort(unique(dt$Cancer_Type))
if (is.na(opt$`cancer-types`)) {
  selected_cancers <- all_cancers
} else {
  selected_cancers <- str_trim(unlist(str_split(opt$`cancer-types`, ",")))
}

# Output Directory
if (!dir.exists(opt$`output-dir`)) dir.create(opt$`output-dir`, recursive = TRUE)

# ------------------------------------------------------------------------------
# Analysis Loop
# ------------------------------------------------------------------------------

for (cancer in selected_cancers) {
  message("\n", paste(rep("=", 60), collapse=""))
  message("🔬 Processing: ", cancer)
  message(paste(rep("=", 60), collapse=""))

  cancer_data <- dt[Cancer_Type == cancer]
  if (nrow(cancer_data) < opt$`min-n`) {
    message("⚠️ Skipping (insufficient data)")
    next
  }

  # Container for all results for this cancer type
  # Structure: list(Lag5_National_Overall = ..., Lag5_K3_C0_Air = ..., etc.)
  cancer_results <- list()

  lags <- sort(unique(cancer_data$Lag))

  for (lag_val in lags) {
    message("  🕐 Lag: ", lag_val)
    lag_data <- cancer_data[Lag == lag_val]
    if (nrow(lag_data) < opt$`min-n`) next

    # --- Phase 1: National ---
    message("    📍 National Analysis")

    # Overall EQI
    tryCatch({
      res <- fit_model(lag_data, build_design_overall_eqi)
      extracted <- extract_draws_and_summary(res$fit, res$design_info, "National_Overall")
      key <- paste0("Lag", lag_val, "_National_Overall")
      cancer_results[[key]] <- extracted
    }, error = function(e) message("      ✗ Failed National Overall: ", e$message))

    # Single Domains
    domains <- c("Air", "Water", "Land", "Built", "Social")
    for (dom in domains) {
      tryCatch({
        res <- fit_model(lag_data, build_design_single_domain, dom)
        extracted <- extract_draws_and_summary(res$fit, res$design_info, paste0("National_", dom))
        key <- paste0("Lag", lag_val, "_National_", dom)
        cancer_results[[key]] <- extracted
      }, error = function(e) message("      ✗ Failed National ", dom, ": ", e$message))
    }

    # --- Phase 2: Cluster Stratified ---
    message("    📍 Cluster Analysis")

    for (k_val in k_values) {
      cluster_col <- paste0("cluster_", k_val)
      cluster_ids <- 0:(k_val - 1)

      for (cid in cluster_ids) {
        c_data <- lag_data[get(cluster_col) == cid]
        if (nrow(c_data) < opt$`min-n`) next

        prefix <- paste0("Lag", lag_val, "_K", k_val, "_C", cid)

        # Overall EQI
        tryCatch({
          res <- fit_model(c_data, build_design_overall_eqi)
          extracted <- extract_draws_and_summary(res$fit, res$design_info, paste0(prefix, "_Overall"))
          cancer_results[[paste0(prefix, "_Overall")]] <- extracted
        }, error = function(e) message("      ✗ Failed ", prefix, " Overall"))

        # Single Domains
        for (dom in domains) {
          tryCatch({
            res <- fit_model(c_data, build_design_single_domain, dom)
            extracted <- extract_draws_and_summary(res$fit, res$design_info, paste0(prefix, "_", dom))
            cancer_results[[paste0(prefix, "_", dom)]] <- extracted
          }, error = function(e) message("      ✗ Failed ", prefix, " ", dom))
        }
      }
    }
  }

  # Save Results
  if (length(cancer_results) > 0) {
    out_file <- file.path(opt$`output-dir`, paste0(cancer, "_ridgeline.rds"))
    saveRDS(cancer_results, out_file)
    message("✅ Saved results to: ", out_file)
  } else {
    message("⚠️ No results generated for ", cancer)
  }
}

message("\nDone.")
