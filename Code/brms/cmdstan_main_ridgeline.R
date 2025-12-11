#!/usr/bin/env Rscript
# Production ridgeline posterior extraction for C00_C97
# Generates complete MCMC draws for Q2-Q5 quintiles across scenarios and models
# Supports: 3 lag scenarios (5,10,15) × 2 models (Overall EQI, Multi-domain)
# Usage: Rscript cmdstan_main_ridgeline.R --scenario <1-3> --model <overall|multi>
# Output: Result/Ridgeline/{Cancer}_{EQI}_{AAMR}_Lag{lag}_{Model}.rds

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(cmdstanr)
  library(posterior)
})

utils::globalVariables(c("EQI", "EQI_Air", "EQI_Water", "EQI_Land", "EQI_Built", "EQI_Social", "Smoking_Rate", "State_FIPS", "RUCC"))

# ============================================================================
# Parse command-line arguments
# ============================================================================
option_list <- list(
  make_option(c("--scenario"),
    type = "integer", default = 1,
    help = "Scenario ID: 1=Lag5, 2=Lag10, 3=Lag15 [default: %default]"
  ),
  make_option(c("--model"),
    type = "character", default = "overall",
    help = "Model type: overall or multi [default: %default]"
  ),
  make_option(c("--cancer"),
    type = "character", default = "C00_C97",
    help = "Cancer type [default: %default]"
  ),
  make_option(c("--data"),
    type = "character",
    default = "Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Cluster_Climate.csv",
    help = "Input data path [default: %default]"
  ),
  make_option(c("--output-dir"),
    type = "character", default = "Result/Ridgeline",
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
if (!opt$scenario %in% 1:3) {
  stop("Invalid scenario: must be 1, 2, or 3")
}
if (!opt$model %in% c("overall", "multi")) {
  stop("Invalid model: must be 'overall' or 'multi'")
}

# Map scenario ID to EQI/AAMR periods
scenario_map <- list(
  list(id = 1, eqi = "2000-2005", aamr = "2006-2010", lag = 5),
  list(id = 2, eqi = "2000-2005", aamr = "2011-2015", lag = 10),
  list(id = 3, eqi = "2000-2005", aamr = "2016-2020", lag = 15)
)

scen <- scenario_map[[opt$scenario]]
EQI_PERIOD <- scen$eqi
AAMR_PERIOD <- scen$aamr
LAG <- scen$lag
CANCER_TYPE <- opt$cancer
MODEL_TYPE <- opt$model

# Core config
cores_avail <- parallel::detectCores(logical = TRUE)
cores_used <- max(1, floor(cores_avail * 0.8))
options(mc.cores = cores_used)

message("========================================")
message("Ridgeline Production Run")
message("========================================")
message("Cancer Type:  ", CANCER_TYPE)
message("EQI Period:   ", EQI_PERIOD)
message("AAMR Period:  ", AAMR_PERIOD)
message("Lag:          ", LAG)
message("Model:        ", toupper(MODEL_TYPE))
message("Iterations:   ", opt$iter, " (warmup: ", opt$warmup, ")")
message("Chains:       ", opt$chains)
message("Cores used:   ", cores_used, " / ", cores_avail)
message("Seed:         ", opt$seed)
message("========================================")

# ============================================================================
# Stan model (interval-censored mixed model)
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

stan_file <- file.path(tempdir(), "interval_mixed_ridgeline.stan")
writeLines(stan_code, stan_file)
message("Compiling Stan model...")
mod <- cmdstan_model(stan_file, quiet = TRUE)
message("✓ Model compiled")

# ============================================================================
# Load and filter data
# ============================================================================
project_root <- normalizePath(".")
data_path <- file.path(project_root, opt$data)

if (!file.exists(data_path)) {
  stop("Data file not found: ", data_path)
}

message("Loading data: ", basename(data_path))
dt <- fread(data_path)

# Check required columns
req_cols <- c(
  "COUNTY_FIPS", "EQI_Period", "Time_Period", "Lag_Years",
  "Cancer_Type", "AAMR_Lower", "AAMR_Upper", "Smoking_Rate",
  "RUCC", "EQI", "EQI_Air", "EQI_Water", "EQI_Land", "EQI_Built", "EQI_Social"
)
miss_cols <- setdiff(req_cols, names(dt))
if (length(miss_cols)) {
  stop("Missing required columns: ", paste(miss_cols, collapse = ", "))
}

# Create State_FIPS if not present
if (!"State_FIPS" %in% names(dt)) {
  dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]
}

# Filter data
dt <- dt[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
dt[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]

# RUCC restriction (1-4 only, or NA)
dt <- dt[RUCC %in% 1:4 | is.na(RUCC)]

# Filter to specific scenario
scen_dt <- dt[Cancer_Type == CANCER_TYPE &
  EQI_Period == EQI_PERIOD &
  Time_Period == AAMR_PERIOD]

message("Filtered data:")
message("  Total rows: ", nrow(scen_dt))
message("  Counties:   ", length(unique(scen_dt$COUNTY_FIPS)))
message("  States:     ", length(unique(scen_dt$State_FIPS)))

if (nrow(scen_dt) < 50) {
  stop("Insufficient data (n=", nrow(scen_dt), "). Need at least 50 observations.")
}

# ============================================================================
# Build design matrix based on model type
# ============================================================================
message("Building design matrix for ", toupper(MODEL_TYPE), " model...")

if (MODEL_TYPE == "overall") {
  # Overall EQI model: Intercept + Smoking_Rate + EQI Q2-Q5
  scen_dt <- scen_dt %>%
    mutate(EQI_factor = factor(EQI, levels = 1:5))

  # Keep only complete cases
  scen_dt <- scen_dt[complete.cases(scen_dt[, c(
    "Smoking_Rate", "EQI_factor",
    "AAMR_Lower", "AAMR_Upper",
    "cens", "State_FIPS"
  )]), ]

  # Build model matrix
  mm <- model.matrix(~ Smoking_Rate + EQI_factor, scen_dt,
    contrasts.arg = list(EQI_factor = contr.treatment(5))
  )
  colnames(mm) <- make.names(colnames(mm))
  design_names <- colnames(mm)

  message("Design matrix (Overall):")
  message("  Rows:    ", nrow(mm))
  message("  Columns: ", ncol(mm))
  message("  Names:   ", paste(design_names, collapse = ", "))
} else if (MODEL_TYPE == "multi") {
  # Multi-domain model: Intercept + Smoking_Rate + all 5 domains (Q2-Q5 each)
  scen_dt <- scen_dt %>%
    mutate(
      EQI_Air_factor = factor(EQI_Air, levels = 1:5),
      EQI_Water_factor = factor(EQI_Water, levels = 1:5),
      EQI_Land_factor = factor(EQI_Land, levels = 1:5),
      EQI_Built_factor = factor(EQI_Built, levels = 1:5),
      EQI_Social_factor = factor(EQI_Social, levels = 1:5)
    )

  # Keep only complete cases
  scen_dt <- scen_dt[complete.cases(scen_dt[, c(
    "Smoking_Rate",
    "EQI_Air_factor", "EQI_Water_factor",
    "EQI_Land_factor", "EQI_Built_factor",
    "EQI_Social_factor",
    "AAMR_Lower", "AAMR_Upper",
    "cens", "State_FIPS"
  )]), ]

  # Build model matrix
  form <- as.formula("~ Smoking_Rate + EQI_Air_factor + EQI_Water_factor + EQI_Land_factor + EQI_Built_factor + EQI_Social_factor")
  mm <- model.matrix(form, scen_dt,
    contrasts.arg = list(
      EQI_Air_factor = contr.treatment(5),
      EQI_Water_factor = contr.treatment(5),
      EQI_Land_factor = contr.treatment(5),
      EQI_Built_factor = contr.treatment(5),
      EQI_Social_factor = contr.treatment(5)
    )
  )
  colnames(mm) <- make.names(colnames(mm))
  design_names <- colnames(mm)

  message("Design matrix (Multi-domain):")
  message("  Rows:    ", nrow(mm))
  message("  Columns: ", ncol(mm))
  message("  Domains: Air, Water, Land, Built, Social (Q2-Q5 each)")
}

# Encode state IDs
states <- sort(unique(scen_dt$State_FIPS))
state_index <- match(scen_dt$State_FIPS, states)

# ============================================================================
# Prepare Stan data
# ============================================================================
stan_data <- list(
  N = nrow(scen_dt),
  S = length(states),
  state = state_index,
  y_lower = scen_dt$AAMR_Lower,
  y_upper = scen_dt$AAMR_Upper,
  cens = scen_dt$cens,
  K = ncol(mm),
  X = mm
)

# Custom initial values
init_fun <- function() {
  list(
    beta = rep(0, stan_data$K),
    z_u = rep(0, stan_data$S),
    sigma = 50,
    sigma_u = 10
  )
}

# ============================================================================
# Run MCMC sampling
# ============================================================================
message("Starting MCMC sampling...")
message("  This may take 10-30 minutes depending on data size...")

fit <- mod$sample(
  data = stan_data,
  chains = opt$chains,
  iter_sampling = opt$iter - opt$warmup,
  iter_warmup = opt$warmup,
  adapt_delta = opt$`adapt-delta`,
  max_treedepth = opt$`max-treedepth`,
  parallel_chains = min(opt$chains, cores_used),
  refresh = 100,
  seed = opt$seed,
  init = rep(list(init_fun()), opt$chains)
)

message("✓ Sampling complete")

# ============================================================================
# Extract posterior draws
# ============================================================================
message("Extracting posterior draws...")

draws_df <- as_draws_df(fit$draws("beta"))
beta_cols <- grep("^beta\\[", colnames(draws_df), value = TRUE)
n_beta <- length(beta_cols)

if (n_beta != length(design_names)) {
  warning(
    "Beta columns (", n_beta, ") don't match design matrix columns (",
    length(design_names), ")"
  )
}

# Map beta columns to design names
beta_mapping <- data.frame(
  beta_col = beta_cols,
  param_name = design_names,
  stringsAsFactors = FALSE
)

message("Beta parameter mapping:")
print(beta_mapping)

# ============================================================================
# Extract quintile draws and compute summaries
# ============================================================================
message("Computing summary statistics...")

summ <- posterior::summarize_draws(fit$draws("beta"))

# Helper function to extract draws for a quintile
extract_quintile_draws <- function(prefix, quintile_num) {
  param_name <- paste0(prefix, quintile_num)
  idx <- which(design_names == param_name)
  if (length(idx) == 0) {
    return(NULL)
  }
  draws_df[[paste0("beta[", idx, "]")]]
}

# Helper function to extract diagnostics
extract_diag <- function(prefix, quintile_num, label) {
  param_name <- paste0(prefix, quintile_num)
  idx <- which(design_names == param_name)
  if (length(idx) == 0) {
    return(list(
      quintile = label,
      mean = NA_real_,
      sd = NA_real_,
      q025 = NA_real_,
      q975 = NA_real_,
      rhat = NA_real_,
      ess_bulk = NA_real_,
      ess_tail = NA_real_
    ))
  }
  row <- summ[idx, ]
  list(
    quintile = label,
    mean = row$mean,
    sd = row$sd,
    q025 = row$q5,
    q975 = row$q95,
    rhat = row$rhat,
    ess_bulk = row$ess_bulk,
    ess_tail = row$ess_tail
  )
}

# Prepare output data structure based on model type
if (MODEL_TYPE == "overall") {
  # Extract Q2-Q5 for Overall EQI
  q2_draws <- extract_quintile_draws("EQI_factor", 2)
  q3_draws <- extract_quintile_draws("EQI_factor", 3)
  q4_draws <- extract_quintile_draws("EQI_factor", 4)
  q5_draws <- extract_quintile_draws("EQI_factor", 5)

  if (is.null(q2_draws) || is.null(q3_draws) || is.null(q4_draws) || is.null(q5_draws)) {
    stop("Could not find all EQI quintile parameters in design matrix")
  }

  n_draws <- length(q2_draws)
  message("Extracted ", n_draws, " draws per quintile (Overall EQI)")

  # Wide format
  draws_wide <- data.frame(
    draw_id = 1:n_draws,
    Q2 = q2_draws,
    Q3 = q3_draws,
    Q4 = q4_draws,
    Q5 = q5_draws
  )

  # Long format
  draws_long <- data.frame(
    draw_id = rep(1:n_draws, 4),
    quintile = factor(rep(c("Q2", "Q3", "Q4", "Q5"), each = n_draws),
      levels = c("Q2", "Q3", "Q4", "Q5")
    ),
    effect = c(q2_draws, q3_draws, q4_draws, q5_draws)
  )

  # Summary
  summary_df <- do.call(rbind, lapply(
    list(
      extract_diag("EQI_factor", 2, "Q2"),
      extract_diag("EQI_factor", 3, "Q3"),
      extract_diag("EQI_factor", 4, "Q4"),
      extract_diag("EQI_factor", 5, "Q5")
    ),
    as.data.frame
  ))

  # Check convergence
  max_rhat <- max(summary_df$rhat, na.rm = TRUE)
  min_ess <- min(summary_df$ess_bulk, na.rm = TRUE)

  message("Summary statistics (Overall EQI):")
  print(summary_df)

  ridge_data <- list(
    metadata = list(
      cancer_type = CANCER_TYPE,
      eqi_period = EQI_PERIOD,
      aamr_period = AAMR_PERIOD,
      lag = LAG,
      model = "Overall",
      n_obs = nrow(scen_dt),
      n_states = length(states),
      n_draws = n_draws,
      iter = opt$iter,
      warmup = opt$warmup,
      chains = opt$chains,
      timestamp = Sys.time(),
      convergence = list(
        max_rhat = max_rhat,
        min_ess_bulk = min_ess
      )
    ),
    draws_wide = draws_wide,
    draws_long = draws_long,
    summary = summary_df,
    design_matrix_names = design_names
  )
} else if (MODEL_TYPE == "multi") {
  # Extract Q2-Q5 for each of the 5 domains
  domains <- c("Air", "Water", "Land", "Built", "Social")
  prefixes <- c(
    "EQI_Air_factor", "EQI_Water_factor", "EQI_Land_factor",
    "EQI_Built_factor", "EQI_Social_factor"
  )

  all_draws_wide <- list()
  all_draws_long <- list()
  all_summaries <- list()

  for (i in seq_along(domains)) {
    domain <- domains[i]
    prefix <- prefixes[i]

    q2 <- extract_quintile_draws(prefix, 2)
    q3 <- extract_quintile_draws(prefix, 3)
    q4 <- extract_quintile_draws(prefix, 4)
    q5 <- extract_quintile_draws(prefix, 5)

    if (is.null(q2) || is.null(q3) || is.null(q4) || is.null(q5)) {
      warning("Could not find all quintile parameters for domain: ", domain)
      next
    }

    n_draws <- length(q2)

    # Wide format
    draws_wide <- data.frame(
      draw_id = 1:n_draws,
      Q2 = q2,
      Q3 = q3,
      Q4 = q4,
      Q5 = q5
    )
    all_draws_wide[[domain]] <- draws_wide

    # Long format
    draws_long <- data.frame(
      draw_id = rep(1:n_draws, 4),
      domain = domain,
      quintile = factor(rep(c("Q2", "Q3", "Q4", "Q5"), each = n_draws),
        levels = c("Q2", "Q3", "Q4", "Q5")
      ),
      effect = c(q2, q3, q4, q5)
    )
    all_draws_long[[domain]] <- draws_long

    # Summary
    summary_df <- do.call(rbind, lapply(
      list(
        extract_diag(prefix, 2, "Q2"),
        extract_diag(prefix, 3, "Q3"),
        extract_diag(prefix, 4, "Q4"),
        extract_diag(prefix, 5, "Q5")
      ),
      as.data.frame
    ))
    summary_df$domain <- domain
    all_summaries[[domain]] <- summary_df
  }

  # Combine all domains
  combined_draws_long <- do.call(rbind, all_draws_long)
  combined_summary <- do.call(rbind, all_summaries)

  message("Extracted ", n_draws, " draws per quintile per domain (5 domains)")
  message("Summary statistics (Multi-domain):")
  print(combined_summary)

  # Check convergence
  max_rhat <- max(combined_summary$rhat, na.rm = TRUE)
  min_ess <- min(combined_summary$ess_bulk, na.rm = TRUE)

  ridge_data <- list(
    metadata = list(
      cancer_type = CANCER_TYPE,
      eqi_period = EQI_PERIOD,
      aamr_period = AAMR_PERIOD,
      lag = LAG,
      model = "MultiDomain",
      domains = domains,
      n_obs = nrow(scen_dt),
      n_states = length(states),
      n_draws = n_draws,
      iter = opt$iter,
      warmup = opt$warmup,
      chains = opt$chains,
      timestamp = Sys.time(),
      convergence = list(
        max_rhat = max_rhat,
        min_ess_bulk = min_ess
      )
    ),
    draws_wide = all_draws_wide, # Named list by domain
    draws_long = combined_draws_long, # Combined data frame
    summary = combined_summary, # Combined summary with domain column
    design_matrix_names = design_names
  )
}

# ============================================================================
# Convergence warnings
# ============================================================================
if (max_rhat > 1.05) {
  warning(
    "Some chains did not converge well (max R-hat = ",
    sprintf("%.3f", max_rhat), ")"
  )
}
if (min_ess < 400) {
  warning(
    "Low effective sample size (min ESS_bulk = ",
    round(min_ess), ")"
  )
}

# ============================================================================
# Save output
# ============================================================================
message("Preparing output...")

out_dir <- file.path(project_root, opt$`output-dir`)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  message("Created output directory: ", out_dir)
}

# Generate filename
eqi_str <- gsub("-", "_", EQI_PERIOD)
aamr_str <- gsub("-", "_", AAMR_PERIOD)
model_str <- ifelse(MODEL_TYPE == "overall", "Overall", "MultiDomain")
out_file <- file.path(out_dir, sprintf(
  "%s_%s_%s_Lag%d_%s.rds",
  CANCER_TYPE, eqi_str, aamr_str, LAG, model_str
))

saveRDS(ridge_data, out_file)

message("========================================")
message("✓ SUCCESS")
message("========================================")
message("Output saved to: ", out_file)
message("File size: ", sprintf("%.2f MB", file.size(out_file) / 1024^2))
message("")
message("Convergence diagnostics:")
message("  Max R-hat:       ", sprintf("%.4f", max_rhat))
message("  Min ESS (bulk):  ", sprintf("%.0f", min_ess))
message("")
message("To load and plot:")
message("  data <- readRDS('", basename(out_file), "')")
if (MODEL_TYPE == "overall") {
  message("  # For Overall model:")
  message("  library(ggridges)")
  message("  ggplot(data$draws_long, aes(x=effect, y=quintile, fill=quintile)) +")
  message("    geom_density_ridges(alpha=0.7)")
} else {
  message("  # For Multi-domain model:")
  message("  library(ggridges)")
  message("  ggplot(data$draws_long, aes(x=effect, y=quintile, fill=domain)) +")
  message("    geom_density_ridges(alpha=0.7) +")
  message("    facet_wrap(~domain)")
}
message("========================================")
