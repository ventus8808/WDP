#!/usr/bin/env Rscript
# Production ridgeline posterior extraction for Delta Models (Change in EQI vs Change in AAMR)
# Generates complete MCMC draws for "Improved" and "Worsened" categories (ref: Stable)
# Usage: Rscript Delta_bayesian_ridgeline.R --lag 5 --model overall --k 1
# Output: Result/Ridgeline/Delta_{Cancer}_Lag{lag}_{Cluster}_{Model}.rds

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(cmdstanr)
  library(posterior)
})

utils::globalVariables(c("EQI_Change_Category", "delta_AAMR_Lower", "delta_AAMR_Upper", "delta_Smoking_Rate", "State_FIPS", "Cluster_K", "Lag_Years", "Cancer_Type"))

# ============================================================================
# Parse command-line arguments
# ============================================================================
option_list <- list(
  make_option(c("--lag"), type="integer", default=5, help="Lag years (5, 10, 15) [default: %default]"),
  make_option(c("--model"), type="character", default="overall", help="Model type: overall or domain [default: %default]"),
  make_option(c("--cancer"), type="character", default="C00_C97", help="Cancer type [default: %default]"),
  make_option(c("--k"), type="integer", default=1, help="Cluster K (1=National) [default: %default]"),
  make_option(c("--cluster-id"), type="integer", default=1, help="Cluster ID to analyze (if k>1) [default: %default]"),
  make_option(c("--data"), type="character", default="Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Cluster_Climate.csv", help="Input data path"),
  make_option(c("--cluster-data"), type="character", default="Data/Processed/Cluster/Cluster_Assignment.csv", help="Cluster assignment data path"),
  make_option(c("--output-dir"), type="character", default="Result/Ridgeline", help="Output directory"),
  make_option(c("--chains"), type="integer", default=4),
  make_option(c("--iter"), type="integer", default=2000),
  make_option(c("--warmup"), type="integer", default=1000),
  make_option(c("--adapt-delta"), type="double", default=0.95),
  make_option(c("--max-treedepth"), type="integer", default=12),
  make_option(c("--seed"), type="integer", default=1234),
  make_option(c("--test"), action="store_true", default=FALSE)
)

opt <- parse_args(OptionParser(option_list=option_list))

if (opt$test) {
  opt$iter <- 800
  opt$warmup <- 400
  message("[TEST MODE] iter=", opt$iter, " warmup=", opt$warmup)
}

set.seed(opt$seed)

# Validate inputs
if (!opt$lag %in% c(5, 10, 15)) stop("Invalid lag: must be 5, 10, or 15")
if (!opt$model %in% c("overall", "domain")) stop("Invalid model: must be 'overall' or 'domain'")

CANCER_TYPE <- opt$cancer
LAG <- opt$lag
MODEL_TYPE <- opt$model
K_VAL <- opt$k
CLUSTER_ID <- opt$`cluster-id`

# Core config
cores_avail <- parallel::detectCores(logical = TRUE)
cores_used <- max(1, floor(cores_avail * 0.8))
options(mc.cores = cores_used)

message("========================================")
message("Delta Ridgeline Production Run")
message("========================================")
message("Cancer Type:  ", CANCER_TYPE)
message("Lag:          ", LAG)
message("Model:        ", toupper(MODEL_TYPE))
message("Cluster K:    ", K_VAL, if(K_VAL > 1) paste0(" (ID: ", CLUSTER_ID, ")") else " (National)")
message("Iterations:   ", opt$iter, " (warmup: ", opt$warmup, ")")

# ============================================================================
# Stan Model (Interval Censored)
# ============================================================================
stan_code <- "
data {
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
  beta ~ normal(0, 5);
  z_u ~ normal(0, 1);
  sigma ~ exponential(1);
  sigma_u ~ exponential(1);
  for (i in 1:N) {
    real mu = X[i] * beta + u[state[i]];
    if (cens[i] == 0) {
      target += normal_lpdf(y_lower[i] | mu, sigma);
    } else {
      real p_up = normal_cdf(y_upper[i] | mu, sigma);
      real p_lo = normal_cdf(y_lower[i] | mu, sigma);
      real diff = fmax(p_up - p_lo, 1e-12);
      target += log(diff);
    }
  }
}
"
stan_file <- file.path(tempdir(), "delta_interval_mixed_model.stan")
writeLines(stan_code, stan_file)
mod <- cmdstan_model(stan_file)

# ============================================================================
# Data Loading & Prep
# ============================================================================
project_root <- normalizePath(".")
data_path <- file.path(project_root, opt$data)
if (!file.exists(data_path)) stop("Data not found: ", data_path)

dt <- fread(data_path)

# Filter Cancer and Lag
dt <- dt[Cancer_Type == CANCER_TYPE & Lag_Years == LAG]
if (nrow(dt) == 0) stop("No data found for Cancer=", CANCER_TYPE, " Lag=", LAG)

# Cluster filtering
if (K_VAL > 1) {
  cluster_path <- file.path(project_root, opt$`cluster-data`)
  if (!file.exists(cluster_path)) stop("Cluster data not found: ", cluster_path)
  cluster_dt <- fread(cluster_path)

  # Expected cluster col: Cluster_K{k}
  cluster_col <- paste0("Cluster_K", K_VAL)
  if (!cluster_col %in% names(cluster_dt)) stop("Cluster column ", cluster_col, " not found")

  # Merge
  dt <- merge(dt, cluster_dt[, c("COUNTY_FIPS", cluster_col), with=FALSE], by="COUNTY_FIPS")

  # Filter by specific cluster ID
  dt <- dt[get(cluster_col) == CLUSTER_ID]
  if (nrow(dt) < 50) stop("Insufficient data for Cluster K=", K_VAL, " ID=", CLUSTER_ID, " (n=", nrow(dt), ")")
}

# Ensure required columns exist
req_cols <- c("delta_AAMR_Lower", "delta_AAMR_Upper", "delta_Smoking_Rate", "cens", "State_FIPS")
miss <- setdiff(req_cols, names(dt))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse=", "))

# Helper to build design matrix
build_design <- function(d, cat_col) {
  # Ensure factor with Stable reference
  if (!cat_col %in% names(d)) stop("Column ", cat_col, " not found")

  # Filter complete cases
  d_sub <- d[!is.na(get(cat_col)) & !is.na(delta_Smoking_Rate) & !is.na(delta_AAMR_Lower) & !is.na(delta_AAMR_Upper)]

  d_sub[[cat_col]] <- factor(d_sub[[cat_col]], levels = c("Stable", "Improved", "Worsened"))

  form <- as.formula(paste0("~ delta_Smoking_Rate + ", cat_col))
  X <- model.matrix(form, data=d_sub)

  list(X=X, df=d_sub, names=colnames(X))
}

# Helper to extract draws for Improved/Worsened
extract_draws <- function(draws_df, design_names, cat_col) {
  # Look for columns like "{cat_col}Improved" and "{cat_col}Worsened"
  imp_name <- paste0(cat_col, "Improved")
  wor_name <- paste0(cat_col, "Worsened")

  idx_imp <- match(imp_name, design_names)
  idx_wor <- match(wor_name, design_names)

  if (is.na(idx_imp) || is.na(idx_wor)) return(NULL)

  imp_draws <- draws_df[[paste0("beta[", idx_imp, "]")]]
  wor_draws <- draws_df[[paste0("beta[", idx_wor, "]")]]

  list(Improved = imp_draws, Worsened = wor_draws)
}

extract_diag <- function(draws_df, design_names, cat_col, suffix, label, summ_df) {
  col_name <- paste0(cat_col, suffix)
  idx <- match(col_name, design_names)
  if (is.na(idx)) return(NULL)

  beta_name <- paste0("beta[", idx, "]")

  # summ_df has variable names like "beta[1]", "beta[2]"
  row <- summ_df[summ_df$variable == beta_name, ]

  list(
    category = label,
    mean = mean(draws_df[[beta_name]]),
    q025 = quantile(draws_df[[beta_name]], 0.025),
    q975 = quantile(draws_df[[beta_name]], 0.975),
    rhat = row$rhat,
    ess_bulk = row$ess_bulk
  )
}

# ============================================================================
# Analysis
# ============================================================================

results_list <- list()

if (MODEL_TYPE == "overall") {
  domains <- c("EQI")
  cat_cols <- c("EQI_Change_Category")
} else {
  domains <- c("Air", "Water", "Land", "Built", "Social")
  cat_cols <- paste0(domains, "_Change_Category")
}

all_draws_wide <- list()
all_draws_long <- list()
all_summaries <- list()

for (i in seq_along(domains)) {
  dom <- domains[i]
  cat_col <- cat_cols[i]

  message("Processing domain: ", dom)

  des <- build_design(dt, cat_col)

  # Stan data
  states <- sort(unique(des$df$State_FIPS))
  state_idx <- match(des$df$State_FIPS, states)

  stan_data <- list(
    N = nrow(des$df),
    S = length(states),
    state = state_idx,
    y_lower = des$df$delta_AAMR_Lower,
    y_upper = des$df$delta_AAMR_Upper,
    cens = des$df$cens,
    K = ncol(des$X),
    X = des$X
  )

  # Fit
  init_fun <- function() list(beta=rep(0, stan_data$K), z_u=rep(0, stan_data$S), sigma=1, sigma_u=1)

  fit <- try(mod$sample(
    data = stan_data,
    chains = opt$chains,
    iter_sampling = opt$iter - opt$warmup,
    iter_warmup = opt$warmup,
    adapt_delta = opt$`adapt-delta`,
    max_treedepth = opt$`max-treedepth`,
    parallel_chains = min(opt$chains, cores_used),
    refresh = 0,
    seed = opt$seed,
    init = rep(list(init_fun()), opt$chains)
  ), silent = TRUE)

  if (inherits(fit, "try-error")) {
    message("Model failed for ", dom)
    next
  }

  # Compute summary before converting to simple DF to preserve chain info for Rhat
  summ_df <- posterior::summarize_draws(fit$draws("beta"))

  draws <- as_draws_df(fit$draws("beta"))
  # colnames(draws) are already beta[1], beta[2]... from cmdstanr

  # Extract
  effs <- extract_draws(draws, des$names, cat_col)
  if (is.null(effs)) {
    message("Could not find coefficients for ", dom)
    next
  }

  n_draws <- length(effs$Improved)

  # Wide
  dw <- data.frame(
    draw_id = 1:n_draws,
    Improved = effs$Improved,
    Worsened = effs$Worsened
  )
  all_draws_wide[[dom]] <- dw

  # Long
  dl <- data.frame(
    draw_id = rep(1:n_draws, 2),
    domain = dom,
    category = factor(rep(c("Improved", "Worsened"), each=n_draws), levels=c("Improved", "Worsened")),
    effect = c(effs$Improved, effs$Worsened)
  )
  all_draws_long[[dom]] <- dl

  # Summary
  s1 <- extract_diag(draws, des$names, cat_col, "Improved", "Improved", summ_df)
  s2 <- extract_diag(draws, des$names, cat_col, "Worsened", "Worsened", summ_df)

  sdf <- rbind(as.data.frame(s1), as.data.frame(s2))
  sdf$domain <- dom
  all_summaries[[dom]] <- sdf

  message("✓ ", dom, " complete")
}

# Combine
combined_draws_long <- do.call(rbind, all_draws_long)
combined_summary <- do.call(rbind, all_summaries)

if (is.null(combined_summary)) stop("No models converged successfully")

# ============================================================================
# Save Output
# ============================================================================
out_dir <- file.path(project_root, opt$`output-dir`)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Filename
# Delta_{Cancer}_Lag{lag}_K{k}_C{cluster}_{Model}.rds
cluster_tag <- if (K_VAL == 1) "National" else paste0("K", K_VAL, "_C", CLUSTER_ID)
model_tag <- ifelse(MODEL_TYPE == "overall", "Overall", "MultiDomain")

out_file <- file.path(out_dir, sprintf(
  "Delta_%s_Lag%d_%s_%s.rds",
  CANCER_TYPE, LAG, cluster_tag, model_tag
))

ridge_data <- list(
  metadata = list(
    cancer_type = CANCER_TYPE,
    lag = LAG,
    model = MODEL_TYPE,
    k = K_VAL,
    cluster_id = CLUSTER_ID,
    n_obs = nrow(dt), # Note: this is n_obs for the last processed domain if multi, but they share the same dt subset
    timestamp = Sys.time()
  ),
  draws_wide = all_draws_wide,
  draws_long = combined_draws_long,
  summary = combined_summary
)

saveRDS(ridge_data, out_file)
message("Saved output to: ", out_file)
