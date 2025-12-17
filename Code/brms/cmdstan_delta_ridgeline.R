#!/usr/bin/env Rscript
# Delta model ridgeline posterior extraction
# Generates complete MCMC draws for improved/worsened EQI changes
# Supports: cluster stratification (k=3, k=4) × multiple lags × Overall/Multi-domain models
# Usage: Rscript cmdstan_delta_ridgeline.R --cancer <type> --k <3|4> --lag <5|10|15> --model <overall|multi>
# Output: Result/Ridgeline_Delta/{Cancer}_k{k}_Lag{lag}_{Model}.rds

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(cmdstanr)
  library(posterior)
})

utils::globalVariables(c(
  "Delta_EQI", "Delta_EQI_Air", "Delta_EQI_Water", "Delta_EQI_Land",
  "Delta_EQI_Built", "Delta_EQI_Social", "Delta_Smoking_Rate",
  "State_FIPS", "RUCC", "EQI_Change_Category", "Cluster_ID"
))

# ============================================================================
# Parse command-line arguments
# ============================================================================
option_list <- list(
  make_option(c("--cancer"),
    type = "character", default = "C00_C97",
    help = "Cancer type [default: %default]"
  ),
  make_option(c("--k"),
    type = "integer", default = 3,
    help = "Number of clusters (3 or 4) [default: %default]"
  ),
  make_option(c("--lag"),
    type = "integer", default = 5,
    help = "Lag years: 5, 10, or 15 [default: %default]"
  ),
  make_option(c("--model"),
    type = "character", default = "overall",
    help = "Model type: overall or multi [default: %default]"
  ),
  make_option(c("--data"),
    type = "character",
    default = "Data/Processed/df_EQI_AAMR/EQI_AAMR_Delta.csv",
    help = "Input data path [default: %default]"
  ),
  make_option(c("--cluster-data"),
    type = "character",
    default = "Data/Processed/County_Cluster/county_cluster_assignments.csv",
    help = "Cluster assignment data [default: %default]"
  ),
  make_option(c("--output-dir"),
    type = "character", default = "Result/Ridgeline_Delta",
    help = "Output directory [default: %default]"
  ),
  make_option(c("--chains"),
    type = "integer", default = 4,
    help = "Number of MCMC chains [default: %default]"
  ),
  make_option(c("--iter"),
    type = "integer", default = 2000,
    help = "Total iterations per chain [default: %default]"
  ),
  make_option(c("--warmup"),
    type = "integer", default = 1000,
    help = "Warmup iterations [default: %default]"
  ),
  make_option(c("--adapt-delta"),
    type = "double", default = 0.95,
    help = "Adapt delta for Stan [default: %default]"
  ),
  make_option(c("--max-treedepth"),
    type = "integer", default = 12,
    help = "Max treedepth for Stan [default: %default]"
  ),
  make_option(c("--seed"),
    type = "integer", default = 1234,
    help = "Random seed [default: %default]"
  ),
  make_option(c("--test"),
    action = "store_true", default = FALSE,
    help = "Test mode: reduced iterations"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Test mode: reduce iterations
if (opt$test) {
  opt$iter <- 800
  opt$warmup <- 400
  message("[TEST MODE] iter=", opt$iter, " warmup=", opt$warmup)
}

set.seed(opt$seed)

# Validate inputs
if (!opt$k %in% c(3, 4)) {
  stop("Invalid k: must be 3 or 4")
}
if (!opt$lag %in% c(5, 10, 15)) {
  stop("Invalid lag: must be 5, 10, or 15")
}
if (!opt$model %in% c("overall", "multi")) {
  stop("Invalid model: must be 'overall' or 'multi'")
}

CANCER_TYPE <- opt$cancer
K_VALUE <- opt$k
LAG <- opt$lag
MODEL_TYPE <- opt$model

# Core config
cores_avail <- parallel::detectCores(logical = TRUE)
cores_used <- max(1, floor(cores_avail * 0.8))
options(mc.cores = cores_used)

message("========================================")
message("Delta Ridgeline Production Run")
message("========================================")
message("Cancer Type:  ", CANCER_TYPE)
message("K Clusters:   ", K_VALUE)
message("Lag:          ", LAG, " years")
message("Model:        ", toupper(MODEL_TYPE))
message("Iterations:   ", opt$iter, " (warmup: ", opt$warmup, ")")
message("Chains:       ", opt$chains)
message("Cores used:   ", cores_used, " / ", cores_avail)
message("Seed:         ", opt$seed)
message("========================================")

# ============================================================================
# Stan model (interval-censored mixed model for delta)
# ============================================================================
stan_code <- "data {
  int<lower=1> N;
  int<lower=1> S;
  array[N] int<lower=1,upper=S> state;
  vector[N] y_lower;
  vector[N] y_upper;
  array[N] int<lower=0,upper=2> cens;
  int<lower=1> K;
  matrix[N,K] X;
}
parameters {
  vector[K] beta;
  vector[S] z_u;
  real<lower=0> sigma;
  real<lower=0> sigma_u;
}
transformed parameters {
  vector[S] u = sigma_u * z_u;
}
model {
  beta ~ normal(0,5);
  z_u ~ normal(0,1);
  sigma ~ exponential(1);
  sigma_u ~ exponential(1);
  for (i in 1:N) {
    real mu = X[i] * beta + u[state[i]];
    if (cens[i]==0) {
      target += normal_lpdf(y_lower[i] | mu, sigma);
    } else {
      real p_up = normal_cdf(y_upper[i] | mu, sigma);
      real p_lo = normal_cdf(y_lower[i] | mu, sigma);
      real diff = fmax(p_up - p_lo, 1e-12);
      target += log(diff);
    }
  }
}"

stan_file <- file.path(tempdir(), "interval_mixed_delta_ridgeline.stan")
writeLines(stan_code, stan_file)
message("Compiling Stan model...")
mod <- cmdstan_model(stan_file, quiet = TRUE)
message("✓ Model compiled")

# ============================================================================
# Load and filter data
# ============================================================================
project_root <- normalizePath(".")
data_path <- file.path(project_root, opt$data)
cluster_path <- file.path(project_root, opt$cluster_data)

if (!file.exists(data_path)) {
  stop("Data file not found: ", data_path)
}
if (!file.exists(cluster_path)) {
  stop("Cluster data file not found: ", cluster_path)
}

message("Loading data: ", basename(data_path))
dt <- fread(data_path)

message("Loading cluster assignments: ", basename(cluster_path))
cluster_dt <- fread(cluster_path)

# Merge cluster assignments
cluster_col <- paste0("Cluster_k", K_VALUE)
if (!cluster_col %in% names(cluster_dt)) {
  stop("Cluster column not found: ", cluster_col)
}
cluster_dt <- cluster_dt[, .(COUNTY_FIPS = County_FIPS, Cluster_ID = get(cluster_col))]
dt <- merge(dt, cluster_dt, by = "COUNTY_FIPS", all.x = TRUE)

# Check required columns
required_cols <- c(
  "COUNTY_FIPS", "Cancer_Type", "Lag_Years", "Delta_AAMR_Lower", "Delta_AAMR_Upper",
  "Delta_Smoking_Rate", "RUCC", "Delta_EQI", "Delta_EQI_Air", "Delta_EQI_Water",
  "Delta_EQI_Land", "Delta_EQI_Built", "Delta_EQI_Social", "EQI_Change_Category", "Cluster_ID"
)
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols)) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Create State_FIPS if not present
if (!"State_FIPS" %in% names(dt)) {
  dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]
}

# Filter data
message("Filtering: Cancer=", CANCER_TYPE, ", Lag=", LAG)
dt <- dt[Cancer_Type == CANCER_TYPE & Lag_Years == LAG]

if (nrow(dt) == 0) {
  stop("No data after filtering")
}

message("✓ Loaded ", nrow(dt), " observations")

# ============================================================================
# Helper functions for design matrix construction
# ============================================================================

build_design_overall <- function(data) {
  # Overall EQI model: Delta_EQI (improved/worsened) + smoking + RUCC
  data <- data %>%
    mutate(
      Improved = as.integer(EQI_Change_Category == "Improved"),
      Worsened = as.integer(EQI_Change_Category == "Worsened")
    )

  X <- model.matrix(~ Improved + Worsened + Delta_Smoking_Rate + RUCC - 1, data = data)

  list(X = X, design_names = colnames(X))
}

build_design_multi <- function(data) {
  # Multi-domain model: each domain separately
  data <- data %>%
    mutate(
      Air_Improved = as.integer(Delta_EQI_Air < 0),
      Air_Worsened = as.integer(Delta_EQI_Air > 0),
      Water_Improved = as.integer(Delta_EQI_Water < 0),
      Water_Worsened = as.integer(Delta_EQI_Water > 0),
      Land_Improved = as.integer(Delta_EQI_Land < 0),
      Land_Worsened = as.integer(Delta_EQI_Land > 0),
      Built_Improved = as.integer(Delta_EQI_Built < 0),
      Built_Worsened = as.integer(Delta_EQI_Built > 0),
      Social_Improved = as.integer(Delta_EQI_Social < 0),
      Social_Worsened = as.integer(Delta_EQI_Social > 0)
    )

  formula_str <- paste0(
    "~ Air_Improved + Air_Worsened + Water_Improved + Water_Worsened + ",
    "Land_Improved + Land_Worsened + Built_Improved + Built_Worsened + ",
    "Social_Improved + Social_Worsened + Delta_Smoking_Rate + RUCC - 1"
  )

  X <- model.matrix(as.formula(formula_str), data = data)

  list(X = X, design_names = colnames(X))
}

# ============================================================================
# Fit model and extract draws
# ============================================================================

fit_and_extract_draws <- function(data, model_type, cluster_id = NULL) {

  cluster_str <- if (!is.null(cluster_id)) paste0("Cluster", cluster_id) else "AllClusters"

  message("\n--- Fitting ", toupper(model_type), " model: ", cluster_str, " ---")

  # Build design matrix
  if (model_type == "overall") {
    design_info <- build_design_overall(data)
  } else {
    design_info <- build_design_multi(data)
  }

  X <- design_info$X
  design_names <- design_info$design_names

  message("Design matrix: ", nrow(X), " × ", ncol(X))
  message("Predictors: ", paste(design_names, collapse = ", "))

  # State indexing
  states <- sort(unique(data$State_FIPS))
  state_index <- match(data$State_FIPS, states)
  n_states <- length(states)

  message("States: ", n_states)

  # Stan data
  stan_data <- list(
    N = nrow(data),
    S = n_states,
    state = state_index,
    y_lower = data$Delta_AAMR_Lower,
    y_upper = data$Delta_AAMR_Upper,
    cens = rep(1L, nrow(data)),  # All interval-censored for delta
    K = ncol(X),
    X = X
  )

  # Fit model
  message("Running MCMC...")
  fit <- mod$sample(
    data = stan_data,
    chains = opt$chains,
    parallel_chains = min(opt$chains, cores_used),
    iter_warmup = opt$warmup,
    iter_sampling = opt$iter - opt$warmup,
    adapt_delta = opt$adapt_delta,
    max_treedepth = opt$max_treedepth,
    seed = opt$seed,
    refresh = 500,
    show_messages = FALSE
  )

  message("✓ Sampling complete")

  # Extract draws
  draws_df <- as_draws_df(fit$draws())
  beta_cols <- grep("^beta\\[", names(draws_df), value = TRUE)
  n_beta <- length(beta_cols)

  if (n_beta != length(design_names)) {
    warning("Beta count mismatch: ", n_beta, " vs ", length(design_names))
  }

  # Create beta mapping
  beta_mapping <- data.frame(
    parameter = beta_cols,
    covariate = design_names,
    stringsAsFactors = FALSE
  )

  # Diagnostics
  summ <- fit$summary(variables = beta_cols)
  max_rhat <- max(summ$rhat, na.rm = TRUE)
  min_ess <- min(summ$ess_bulk, na.rm = TRUE)

  message("Diagnostics: max Rhat = ", round(max_rhat, 4),
          ", min ESS = ", round(min_ess, 0))

  # Prepare output
  list(
    draws = draws_df,
    beta_cols = beta_cols,
    beta_mapping = beta_mapping,
    design_names = design_names,
    diagnostics = summ,
    cluster = cluster_str,
    model_type = model_type,
    max_rhat = max_rhat,
    min_ess = min_ess
  )
}

# ============================================================================
# Extract draws for improved/worsened effects
# ============================================================================

extract_effect_draws <- function(result, effect_type) {
  # effect_type: "Improved" or "Worsened" (for overall model)
  # For multi-domain: "Air_Improved", "Water_Worsened", etc.

  beta_mapping <- result$beta_mapping
  draws_df <- result$draws

  # Find matching parameter
  idx <- which(grepl(effect_type, beta_mapping$covariate, fixed = TRUE))

  if (length(idx) == 0) {
    return(NULL)
  }

  effect_draws <- list()

  for (i in idx) {
    param_name <- beta_mapping$parameter[i]
    cov_name <- beta_mapping$covariate[i]

    draws <- draws_df[[param_name]]

    effect_draws[[cov_name]] <- data.frame(
      iteration = 1:length(draws),
      draw = draws,
      parameter = param_name,
      covariate = cov_name,
      effect_type = effect_type,
      cluster = result$cluster,
      model_type = result$model_type,
      stringsAsFactors = FALSE
    )
  }

  if (length(effect_draws) > 0) {
    do.call(rbind, effect_draws)
  } else {
    NULL
  }
}

# ============================================================================
# Main analysis loop
# ============================================================================

message("\n" , rep("=", 60), sep = "")
message("Starting analysis for k=", K_VALUE)
message(rep("=", 60), sep = "")

all_results <- list()
all_draws_list <- list()

# Get unique clusters for this k
clusters <- sort(unique(dt$Cluster_ID))
message("Clusters: ", paste(clusters, collapse = ", "))

# Fit models for each cluster
for (cluster_id in clusters) {

  message("\n", rep("-", 60), sep = "")
  message("Processing Cluster ", cluster_id, " / ", max(clusters))
  message(rep("-", 60), sep = "")

  cluster_data <- dt[Cluster_ID == cluster_id]

  message("N = ", nrow(cluster_data), " counties")

  # Fit model
  result <- fit_and_extract_draws(cluster_data, MODEL_TYPE, cluster_id)
  all_results[[paste0("Cluster", cluster_id)]] <- result

  # Extract improved/worsened draws
  if (MODEL_TYPE == "overall") {
    improved_draws <- extract_effect_draws(result, "Improved")
    worsened_draws <- extract_effect_draws(result, "Worsened")

    if (!is.null(improved_draws)) all_draws_list[[length(all_draws_list) + 1]] <- improved_draws
    if (!is.null(worsened_draws)) all_draws_list[[length(all_draws_list) + 1]] <- worsened_draws

  } else {
    # Multi-domain: extract all domain-specific effects
    domains <- c("Air", "Water", "Land", "Built", "Social")
    for (domain in domains) {
      improved_draws <- extract_effect_draws(result, paste0(domain, "_Improved"))
      worsened_draws <- extract_effect_draws(result, paste0(domain, "_Worsened"))

      if (!is.null(improved_draws)) all_draws_list[[length(all_draws_list) + 1]] <- improved_draws
      if (!is.null(worsened_draws)) all_draws_list[[length(all_draws_list) + 1]] <- worsened_draws
    }
  }
}

# Combine all draws
message("\n", rep("=", 60), sep = "")
message("Combining results...")

if (length(all_draws_list) > 0) {
  combined_draws <- do.call(rbind, all_draws_list)
  combined_draws$cancer <- CANCER_TYPE
  combined_draws$k <- K_VALUE
  combined_draws$lag <- LAG
} else {
  combined_draws <- data.frame()
  warning("No draws extracted!")
}

# Combine diagnostics
all_diagnostics <- lapply(all_results, function(x) {
  diag <- x$diagnostics
  diag$cluster <- x$cluster
  diag$model_type <- x$model_type
  diag
})
combined_diagnostics <- do.call(rbind, all_diagnostics)

# Overall diagnostics summary
max_rhat_overall <- max(combined_diagnostics$rhat, na.rm = TRUE)
min_ess_overall <- min(combined_diagnostics$ess_bulk, na.rm = TRUE)

message("Overall diagnostics:")
message("  Max Rhat: ", round(max_rhat_overall, 4))
message("  Min ESS:  ", round(min_ess_overall, 0))

if (max_rhat_overall > 1.01) {
  warning("Rhat > 1.01 detected! Model may not have converged.")
}
if (min_ess_overall < 100) {
  warning("ESS < 100 detected! Increase iterations.")
}

# ============================================================================
# Save output
# ============================================================================

out_dir <- file.path(project_root, opt$output_dir)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  message("Created output directory: ", out_dir)
}

# Clean cancer name for filename
cancer_clean <- gsub("^delta_AAMR_", "", CANCER_TYPE)
model_str <- toupper(MODEL_TYPE)

out_file <- file.path(
  out_dir,
  sprintf("%s_k%d_Lag%d_%s.rds", cancer_clean, K_VALUE, LAG, model_str)
)

# Save results
output <- list(
  draws = combined_draws,
  diagnostics = combined_diagnostics,
  beta_mappings = lapply(all_results, function(x) x$beta_mapping),
  model_summaries = lapply(all_results, function(x) x$diagnostics),
  metadata = list(
    cancer = CANCER_TYPE,
    k = K_VALUE,
    lag = LAG,
    model_type = MODEL_TYPE,
    n_obs = nrow(dt),
    n_clusters = length(clusters),
    clusters = clusters,
    chains = opt$chains,
    iter = opt$iter,
    warmup = opt$warmup,
    adapt_delta = opt$adapt_delta,
    max_treedepth = opt$max_treedepth,
    seed = opt$seed,
    max_rhat = max_rhat_overall,
    min_ess = min_ess_overall,
    timestamp = Sys.time()
  )
)

saveRDS(output, out_file)
message("\n✓ Saved ridgeline data: ", out_file)

# Summary table
message("\n", rep("=", 60), sep = "")
message("Summary of extracted draws:")
message(rep("=", 60), sep = "")

if (nrow(combined_draws) > 0) {
  summary_table <- combined_draws %>%
    group_by(covariate, cluster, effect_type) %>%
    summarise(
      n_draws = n(),
      mean = mean(draw),
      sd = sd(draw),
      q025 = quantile(draw, 0.025),
      q975 = quantile(draw, 0.975),
      .groups = "drop"
    )

  print(summary_table, n = Inf)

  message("\nTotal draws extracted: ", nrow(combined_draws))
} else {
  message("No draws to summarize")
}

message("\n", rep("=", 60), sep = "")
message("✅ Delta ridgeline extraction complete!")
message(rep("=", 60), sep = "")
