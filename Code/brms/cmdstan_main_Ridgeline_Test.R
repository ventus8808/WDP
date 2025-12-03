#!/usr/bin/env Rscript
# Ridgeline plot posterior extraction for C00-C97, lag=5, Overall EQI only
# Extracts full MCMC draws for Q2-Q5 quintiles to enable ridgeline visualization
# Output: Result/Ridgeline/C00-C97_Ridge_Test.rds

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(cmdstanr)
  library(posterior)
})

utils::globalVariables(c('EQI','Smoking_Rate','State_FIPS','RUCC'))

# ============================================================================
# Configuration (hardcoded for test)
# ============================================================================
CANCER_TYPE <- "C00_C97"
EQI_PERIOD <- "2000-2005"
AAMR_PERIOD <- "2006-2010"
LAG <- 5
ITER <- 800
WARMUP <- 400
CHAINS <- 4
ADAPT_DELTA <- 0.95
MAX_TREEDEPTH <- 12
SEED <- 1234
set.seed(SEED)

cores_avail <- parallel::detectCores(logical=TRUE)
cores_used <- max(1, floor(cores_avail * 0.8))
options(mc.cores = cores_used)

message("========================================")
message("Ridgeline Test: Posterior Extraction")
message("========================================")
message("Cancer Type: ", CANCER_TYPE)
message("EQI Period:  ", EQI_PERIOD)
message("AAMR Period: ", AAMR_PERIOD)
message("Lag:         ", LAG)
message("Iterations:  ", ITER, " (warmup: ", WARMUP, ")")
message("Chains:      ", CHAINS)
message("Cores used:  ", cores_used, " / ", cores_avail)
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
mod <- cmdstan_model(stan_file)
message("✓ Model compiled")

# ============================================================================
# Load and filter data
# ============================================================================
project_root <- normalizePath(".")
data_path <- file.path(project_root,
                       "Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Cluster_Climate.csv")

if (!file.exists(data_path)) {
  stop("Data file not found: ", data_path)
}

message("Loading data: ", basename(data_path))
dt <- fread(data_path)

# Check required columns
req_cols <- c("COUNTY_FIPS", "EQI_Period", "Time_Period", "Lag_Years",
              "Cancer_Type", "AAMR_Lower", "AAMR_Upper", "Smoking_Rate",
              "RUCC", "EQI")
miss_cols <- setdiff(req_cols, names(dt))
if (length(miss_cols)) {
  stop("Missing required columns: ", paste(miss_cols, collapse=", "))
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
# Build design matrix (Overall EQI model only)
# ============================================================================
message("Building design matrix...")

# Create EQI factor with treatment contrasts (Q1 = reference)
scen_dt <- scen_dt %>%
  mutate(EQI_factor = factor(EQI, levels = 1:5))

# Keep only complete cases
scen_dt <- scen_dt[complete.cases(scen_dt[, c("Smoking_Rate", "EQI_factor",
                                                "AAMR_Lower", "AAMR_Upper",
                                                "cens", "State_FIPS")]), ]

# Build model matrix
mm <- model.matrix(~ Smoking_Rate + EQI_factor, scen_dt,
                   contrasts.arg = list(EQI_factor = contr.treatment(5)))
colnames(mm) <- make.names(colnames(mm))
design_names <- colnames(mm)

message("Design matrix:")
message("  Rows:    ", nrow(mm))
message("  Columns: ", ncol(mm))
message("  Names:   ", paste(design_names, collapse=", "))

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
message("  This may take several minutes...")

fit <- mod$sample(
  data = stan_data,
  chains = CHAINS,
  iter_sampling = ITER - WARMUP,
  iter_warmup = WARMUP,
  adapt_delta = ADAPT_DELTA,
  max_treedepth = MAX_TREEDEPTH,
  parallel_chains = min(CHAINS, cores_used),
  refresh = 100,
  seed = SEED,
  init = rep(list(init_fun()), CHAINS)
)

message("✓ Sampling complete")

# ============================================================================
# Extract posterior draws
# ============================================================================
message("Extracting posterior draws...")

# Extract beta parameters
draws_df <- as_draws_df(fit$draws("beta"))

# Identify column indices for Q2-Q5
# design_names should be: "(Intercept)", "Smoking_Rate", "EQI_factor2", "EQI_factor3", "EQI_factor4", "EQI_factor5"
beta_cols <- grep("^beta\\[", colnames(draws_df), value = TRUE)
n_beta <- length(beta_cols)

if (n_beta != length(design_names)) {
  warning("Beta columns (", n_beta, ") don't match design matrix columns (",
          length(design_names), ")")
}

# Map beta columns to design names
beta_mapping <- data.frame(
  beta_col = beta_cols,
  param_name = design_names,
  stringsAsFactors = FALSE
)

message("Beta parameter mapping:")
print(beta_mapping)

# Find indices for EQI quintiles
q2_idx <- which(design_names == "EQI_factor2")
q3_idx <- which(design_names == "EQI_factor3")
q4_idx <- which(design_names == "EQI_factor4")
q5_idx <- which(design_names == "EQI_factor5")

if (length(q2_idx) == 0 || length(q3_idx) == 0 ||
    length(q4_idx) == 0 || length(q5_idx) == 0) {
  stop("Could not find all EQI quintile parameters in design matrix")
}

# Extract draws for each quintile
q2_draws <- draws_df[[paste0("beta[", q2_idx, "]")]]
q3_draws <- draws_df[[paste0("beta[", q3_idx, "]")]]
q4_draws <- draws_df[[paste0("beta[", q4_idx, "]")]]
q5_draws <- draws_df[[paste0("beta[", q5_idx, "]")]]

n_draws <- length(q2_draws)
message("Extracted ", n_draws, " draws per quintile")

# ============================================================================
# Compute summary statistics
# ============================================================================
message("Computing summary statistics...")

summ <- posterior::summarize_draws(fit$draws("beta"))

# Extract diagnostics for each quintile
extract_diag <- function(idx, label) {
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

diag_q2 <- extract_diag(q2_idx, "Q2")
diag_q3 <- extract_diag(q3_idx, "Q3")
diag_q4 <- extract_diag(q4_idx, "Q4")
diag_q5 <- extract_diag(q5_idx, "Q5")

summary_df <- do.call(rbind, lapply(
  list(diag_q2, diag_q3, diag_q4, diag_q5),
  as.data.frame
))

message("Summary statistics:")
print(summary_df)

# Check convergence
max_rhat <- max(summary_df$rhat, na.rm = TRUE)
min_ess <- min(summary_df$ess_bulk, na.rm = TRUE)

if (max_rhat > 1.05) {
  warning("Some chains did not converge well (max R-hat = ",
          sprintf("%.3f", max_rhat), ")")
}
if (min_ess < 400) {
  warning("Low effective sample size (min ESS_bulk = ",
          round(min_ess), ")")
}

# ============================================================================
# Prepare output data structure
# ============================================================================
message("Preparing output data structure...")

# Wide format (one row per draw)
draws_wide <- data.frame(
  draw_id = 1:n_draws,
  Q2 = q2_draws,
  Q3 = q3_draws,
  Q4 = q4_draws,
  Q5 = q5_draws
)

# Long format (convenient for ggridges)
draws_long <- data.frame(
  draw_id = rep(1:n_draws, 4),
  quintile = factor(rep(c("Q2", "Q3", "Q4", "Q5"), each = n_draws),
                    levels = c("Q2", "Q3", "Q4", "Q5")),
  effect = c(q2_draws, q3_draws, q4_draws, q5_draws)
)

# Full output list
ridge_data <- list(
  metadata = list(
    cancer_type = CANCER_TYPE,
    eqi_period = EQI_PERIOD,
    aamr_period = AAMR_PERIOD,
    lag = LAG,
    n_obs = nrow(scen_dt),
    n_states = length(states),
    n_draws = n_draws,
    iter = ITER,
    warmup = WARMUP,
    chains = CHAINS,
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

# ============================================================================
# Save output
# ============================================================================
out_dir <- file.path(project_root, "Result/Ridgeline")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  message("Created output directory: ", out_dir)
}

out_file <- file.path(out_dir, "C00_C97_Ridge_Test.rds")
saveRDS(ridge_data, out_file)

message("========================================")
message("✓ SUCCESS")
message("========================================")
message("Output saved to: ", out_file)
message("File size: ", sprintf("%.2f MB", file.size(out_file) / 1024^2))
message("")
message("To load and plot:")
message("  data <- readRDS('", out_file, "')")
message("  library(ggridges)")
message("  ggplot(data$draws_long, aes(x=effect, y=quintile, fill=quintile)) +")
message("    geom_density_ridges(alpha=0.7)")
message("========================================")
