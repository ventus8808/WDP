#!/usr/bin/env Rscript
# Pure cmdstanr interval-censored mixed model pipeline (no brms / no rstan)
# SVI trajectory exposure version of brms_Main.R.
#   - Exposure: SVI trajectory category (A/B/C/D), A = reference.
#   - Output:   B/C/D mortality rate differences (MRD) vs A. Because the model is
#               linear on AAMR, each exposure coefficient is the adjusted AAMR
#               difference (= MRD) relative to category A.
#   - SVI is a single static 2000-2022 trajectory class per county, so there is
#     no EQI-style period; the four AAMR periods are treated as four lags
#     (2006-2010 -> 2021-2024, lag = 5/10/15/20 years from baseline).
#   - Multi-domain model removed (SVI is a single overall index).
# Interval likelihood: exact rows (cens=0) use point normal density; interval
# rows (cens=2) use the CDF difference.

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
utils::globalVariables(c(
  "SVI", "State_FIPS",
  "Smoking_rate", "Physical_Activities_rate", "Obesity_rate"
))

option_list <- list(
  make_option(c("--data"), type = "character", default = "Data/Processed/df_SVI.csv"),
  make_option(c("--output-dir"), type = "character", default = "Result/brms_SVI_Main"),
  make_option(c("--outcomes"), type = "character", default = NA),
  make_option(c("--chains"), type = "integer", default = 6),
  make_option(c("--iter"), type = "integer", default = 1800),
  make_option(c("--warmup"), type = "integer", default = 1000),
  make_option(c("--adapt-delta"), type = "double", default = 0.95),
  make_option(c("--max-treedepth"), type = "integer", default = 12),
  make_option(c("--min-n"), type = "integer", default = 50),
  make_option(c("--seed"), type = "integer", default = 1234),
  make_option(c("--test"), action = "store_true", default = FALSE)
)
opt <- parse_args(OptionParser(option_list = option_list))

cores_avail <- parallel::detectCores(logical = TRUE)
slurm_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA)))
cores_used <- opt$chains
options(mc.cores = cores_used)

message("--- CPU Resource Report ---")
message("Environment: ", if (!is.na(slurm_cpus)) "Slurm (HPC)" else "Local Machine")
message("Total Cores Available: ", cores_avail)
message("Setting mc.cores to:   ", cores_used)
message("---------------------------")

set.seed(opt$seed)

# ── Stan model ─────────────────────────────────────────────────────────────────
stan_code <- "data {\n  int<lower=1> N;\n  int<lower=1> S;\n  array[N] int<lower=1,upper=S> state;\n  vector[N] y_lower;\n  vector[N] y_upper;\n  array[N] int<lower=0,upper=2> cens;\n  int<lower=1> K;\n  matrix[N,K] X;\n} \nparameters {\n  vector[K] beta;\n  vector[S] z_u;\n  real<lower=0> sigma;\n  real<lower=0> sigma_u;\n} \ntransformed parameters {\n  vector[S] u = sigma_u * z_u;\n} \nmodel {\n  beta ~ normal(0,5);\n  z_u ~ normal(0,1);\n  sigma ~ exponential(1);\n  sigma_u ~ exponential(1);\n  for (i in 1:N) {\n    real mu = X[i] * beta + u[state[i]];\n    if (cens[i]==0) {\n      target += normal_lpdf(y_lower[i] | mu, sigma);\n    } else {\n      real p_up = normal_cdf(y_upper[i] | mu, sigma);\n      real p_lo = normal_cdf(y_lower[i] | mu, sigma);\n      real diff = fmax(p_up - p_lo, 1e-12);\n      target += log(diff);\n    }\n  }\n}"
stan_file <- file.path(tempdir(), "interval_mixed_model.stan")
writeLines(stan_code, stan_file)
mod <- cmdstan_model(stan_file)

# Load data
project_root <- normalizePath(".")
path <- file.path(project_root, opt$data)
if (!file.exists(path)) stop("Data not found: ", path)
dt <- fread(path)

req <- c(
  "COUNTY_FIPS", "Time_Period", "Outcome", "AAMR_Lower", "AAMR_Upper",
  "SVI", "Smoking_rate", "Physical_Activities_rate", "Obesity_rate"
)
miss <- setdiff(req, names(dt))
if (length(miss)) stop("Missing cols: ", paste(miss, collapse = ","))

if (!"State_FIPS" %in% names(dt)) dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]

# interval censoring code
dt <- dt[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
dt[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]

# Four lags: SVI baseline trajectory -> each AAMR period (5/10/15/20 yr)
scenario_list <- list(
  list(key = "SVI_AAMR2006_2010", aamr = "2006-2010", lag = 5),
  list(key = "SVI_AAMR2011_2015", aamr = "2011-2015", lag = 10),
  list(key = "SVI_AAMR2016_2020", aamr = "2016-2020", lag = 15),
  list(key = "SVI_AAMR2021_2024", aamr = "2021-2024", lag = 20)
)

all_outcomes <- sort(unique(dt$Outcome))
selected <- if (is.na(opt$`outcomes`)) {
  all_outcomes
} else {
  reqc <- str_split(opt$`outcomes`, ",", simplify = TRUE) |>
    as.vector() |>
    str_trim()
  inv <- setdiff(reqc, all_outcomes)
  if (length(inv)) stop("Invalid outcomes: ", paste(inv, collapse = ","))
  reqc
}
message("Outcomes to analyze: ", paste(selected, collapse = ","))

out_dir <- file.path(project_root, opt$`output-dir`)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

format_cell <- function(draws) {
  if (length(draws) == 0) {
    return("")
  }
  ci <- quantile(draws, c(0.025, 0.975), na.rm = TRUE)
  sprintf("%0.2f(%0.2f,%0.2f)", mean(draws), ci[1], ci[2])
}
append_rows <- function(path, df) {
  if (!file.exists(path)) write_csv(df, path) else suppressWarnings(write.table(df, path, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE))
}

# Posterior diagnostics helpers
# p_posterior is the two-sided posterior tail-area probability:
#   p = 2 * min(Pr(beta > 0), Pr(beta < 0)) estimated from the proportion of draws.
# R-hat, ESS (bulk and tail) come from split-chain draws via posterior::summarize_draws.
compute_p <- function(draws) {
  if (length(draws) == 0) {
    return(NA_character_)
  }
  pos <- sum(draws > 0, na.rm = TRUE)
  neg <- sum(draws < 0, na.rm = TRUE)
  n <- pos + neg
  if (n == 0) {
    return(NA_character_)
  }
  p_pos <- (pos + 0.5) / (n + 1)
  p_neg <- (neg + 0.5) / (n + 1)
  p <- 2 * min(p_pos, p_neg)
  sprintf("%.4f", p)
}

# Extract per-category (B/C/D vs A) MRD estimates from treatment contrasts.
extract_svi <- function(draw_df, names_vec) {
  out <- list(A = "0.00", B = "", C = "", D = "")
  for (lv in c("B", "C", "D")) {
    nm <- paste0("SVI_factor", lv)
    idx <- match(nm, names_vec)
    out[[lv]] <- if (is.na(idx)) "" else format_cell(draw_df[[paste0("beta[", idx, "]")]])
  }
  out
}

extract_svi_metrics <- function(draw_df, names_vec, summ_df) {
  out <- list(
    B_p = NA_character_, C_p = NA_character_, D_p = NA_character_,
    B_rhat = NA_real_, C_rhat = NA_real_, D_rhat = NA_real_,
    B_ess_bulk = NA_real_, C_ess_bulk = NA_real_, D_ess_bulk = NA_real_,
    B_ess_tail = NA_real_, C_ess_tail = NA_real_, D_ess_tail = NA_real_
  )
  for (lv in c("B", "C", "D")) {
    nm <- paste0("SVI_factor", lv)
    idx <- match(nm, names_vec)
    if (!is.na(idx)) {
      col <- paste0("beta[", idx, "]")
      out[[paste0(lv, "_p")]] <- compute_p(draw_df[[col]])
      sr <- summ_df[summ_df$variable == col, , drop = FALSE]
      if (nrow(sr)) {
        out[[paste0(lv, "_rhat")]] <- sr$rhat
        out[[paste0(lv, "_ess_bulk")]] <- sr$ess_bulk
        out[[paste0(lv, "_ess_tail")]] <- sr$ess_tail
      }
    }
  }
  out
}

extract_covariate <- function(draw_df, names_vec, col_name, summ_df) {
  idx <- match(col_name, names_vec)
  if (is.na(idx)) {
    return(list(est = "", p = NA_character_, rhat = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_))
  }
  col <- paste0("beta[", idx, "]")
  draws_col <- draw_df[[col]]
  sr <- summ_df[summ_df$variable == col, , drop = FALSE]
  list(
    est      = format_cell(draws_col),
    p        = compute_p(draws_col),
    rhat     = if (nrow(sr)) sr$rhat else NA_real_,
    ess_bulk = if (nrow(sr)) sr$ess_bulk else NA_real_,
    ess_tail = if (nrow(sr)) sr$ess_tail else NA_real_
  )
}

build_design_overall <- function(d) {
  # SVI is unordered -> default treatment contrasts (A reference) -> names SVI_factorB/C/D
  d <- d %>% mutate(SVI_factor = factor(SVI, levels = c("A", "B", "C", "D")))
  d <- d[complete.cases(d[, c(
    "SVI_factor", "AAMR_Lower", "AAMR_Upper", "cens", "State_FIPS",
    "Smoking_rate", "Physical_Activities_rate", "Obesity_rate"
  )]), ]
  mm <- model.matrix(~ Smoking_rate + Physical_Activities_rate + Obesity_rate + SVI_factor, d,
    contrasts.arg = list(SVI_factor = contr.treatment(c("A", "B", "C", "D")))
  )
  colnames(mm) <- make.names(colnames(mm))
  list(X = mm, names = colnames(mm), df = d)
}

for (outcome in selected) {
  message("===== Outcome: ", outcome, " =====")
  outfile <- file.path(out_dir, paste0(outcome, "_SVI.csv"))
  for (sc in scenario_list) {
    scen_key <- sc$key
    aamr_p <- sc$aamr
    lagv <- sc$lag
    scen_dt <- dt[Time_Period == aamr_p & Outcome == outcome]
    if (nrow(scen_dt) < opt$`min-n`) {
      message("[Skip] Scenario ", scen_key, " n=", nrow(scen_dt))
      next
    }
    aamr_out <- gsub("-", "_", aamr_p)

    # Overall SVI + SM + PA + OB model
    des_overall <- build_design_overall(scen_dt)
    if (nrow(des_overall$df) < opt$`min-n`) {
      message("[Skip] ", scen_key, " SVI+SM+PA+OB n=", nrow(des_overall$df))
      next
    }
    states_o <- sort(unique(des_overall$df$State_FIPS))
    state_index_o <- match(des_overall$df$State_FIPS, states_o)
    data_list <- list(
      N = nrow(des_overall$df), S = length(states_o), state = state_index_o,
      y_lower = des_overall$df$AAMR_Lower, y_upper = des_overall$df$AAMR_Upper, cens = des_overall$df$cens,
      K = ncol(des_overall$X), X = des_overall$X
    )
    init_fun <- function() list(beta = rep(0, data_list$K), z_u = rep(0, data_list$S), sigma = 50, sigma_u = 10)
    fit_overall <- try(mod$sample(
      data = data_list, chains = opt$chains, iter_sampling = opt$iter - opt$warmup, iter_warmup = opt$warmup,
      adapt_delta = opt$`adapt-delta`, max_treedepth = opt$`max-treedepth`, parallel_chains = opt$chains, refresh = 0, seed = opt$seed,
      init = rep(list(init_fun()), opt$chains)
    ), silent = TRUE)
    if (inherits(fit_overall, "try-error")) {
      message("[Fail] SVI+SM+PA+OB model ", scen_key)
      next
    }
    draws <- as_draws_df(fit_overall$draws("beta"))
    colnames(draws) <- paste0("beta[", seq_len(ncol(draws)), "]")
    sv <- extract_svi(draws, des_overall$names)
    summ_over <- posterior::summarize_draws(fit_overall$draws("beta"))
    met <- extract_svi_metrics(draws, des_overall$names, summ_over)
    sm_o <- extract_covariate(draws, des_overall$names, "Smoking_rate", summ_over)
    pa_o <- extract_covariate(draws, des_overall$names, "Physical_Activities_rate", summ_over)
    ob_o <- extract_covariate(draws, des_overall$names, "Obesity_rate", summ_over)
    row_over <- tibble(
      ICD_Code = outcome, AAMR_Period = aamr_out, Lag = lagv, Model = "SVI+SM+PA+OB",
      A = sv$A, B = sv$B, C = sv$C, D = sv$D,
      B_p = met$B_p, C_p = met$C_p, D_p = met$D_p,
      B_rhat = sprintf("%.4f", met$B_rhat), C_rhat = sprintf("%.4f", met$C_rhat), D_rhat = sprintf("%.4f", met$D_rhat),
      B_ess_bulk = as.integer(round(met$B_ess_bulk)), C_ess_bulk = as.integer(round(met$C_ess_bulk)), D_ess_bulk = as.integer(round(met$D_ess_bulk)),
      B_ess_tail = as.integer(round(met$B_ess_tail)), C_ess_tail = as.integer(round(met$C_ess_tail)), D_ess_tail = as.integer(round(met$D_ess_tail)),
      SM = sm_o$est, SM_p = sm_o$p, SM_rhat = sprintf("%.4f", sm_o$rhat), SM_ess_bulk = as.integer(round(sm_o$ess_bulk)), SM_ess_tail = as.integer(round(sm_o$ess_tail)),
      PA = pa_o$est, PA_p = pa_o$p, PA_rhat = sprintf("%.4f", pa_o$rhat), PA_ess_bulk = as.integer(round(pa_o$ess_bulk)), PA_ess_tail = as.integer(round(pa_o$ess_tail)),
      OB = ob_o$est, OB_p = ob_o$p, OB_rhat = sprintf("%.4f", ob_o$rhat), OB_ess_bulk = as.integer(round(ob_o$ess_bulk)), OB_ess_tail = as.integer(round(ob_o$ess_tail))
    )
    append_rows(outfile, row_over)
    message("[OK] ", scen_key, " SVI+SM+PA+OB")
  }
  message("===== Completed: ", outcome, " =====")
}
message("All requested analyses complete. Output directory: ", out_dir)
