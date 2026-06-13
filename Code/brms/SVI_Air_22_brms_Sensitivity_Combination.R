#!/usr/bin/env Rscript
# EQI-Air × SVI 2x2 joint sensitivity combination analysis.
# Exposure: EQI_Air_SVI_2 (4 joint levels a1s1/a1s2/a2s1/a2s2), a1s1 = reference;
#           output a1s2 / a2s1 / a2s2 mortality rate differences (MRD) vs a1s1.
# All 2^6 = 64 subsets of 6 control covariates, single representative scenario
# (AAMR 2011-2015, lag 10). Same df as the 3x3 script (Data/Processed/df_SVI.csv).
# Output format: ICD_Code, AAMR_Period, Lag, Model, joint levels + diagnostics.

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

candidates <- list.dirs(path.expand("~/.cmdstan"), recursive = FALSE)
if (length(candidates) == 0) stop("未找到 cmdstan，请先运行 install_cmdstan()")
cmdstanr::set_cmdstan_path(tail(sort(candidates), 1))

utils::globalVariables(c(
  "EQI_Air_SVI_2", "Smoking_rate", "Physical_Activities_rate", "Obesity_rate",
  "Uninsured_rate", "Physician_Density_per100k", "Diabetes_Prevalence_rate", "State_FIPS"
))

# ── Joint exposure configuration (2x2) ──────────────────────────────────────────
JOINT_COL    <- "EQI_Air_SVI_2"
JOINT_LEVELS <- c("a1s1", "a1s2", "a2s1", "a2s2")   # a1s1 = reference (first)
REF_LEVEL    <- "a1s1"
NON_REF      <- setdiff(JOINT_LEVELS, REF_LEVEL)

option_list <- list(
  make_option(c("--data"), type = "character", default = "Data/Processed/df_SVI.csv"),
  make_option(c("--output-dir"), type = "character", default = "Result/brms_SVI_Air_22_Sensitivity_Combination"),
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
  "COUNTY_FIPS", "Time_Period", "Outcome", "AAMR_Lower", "AAMR_Upper", JOINT_COL,
  "Smoking_rate", "Physical_Activities_rate", "Obesity_rate", "Uninsured_rate", "Physician_Density_per100k", "Diabetes_Prevalence_rate"
)
miss <- setdiff(req, names(dt))
if (length(miss)) stop("Missing cols: ", paste(miss, collapse = ","))

if (!"State_FIPS" %in% names(dt)) dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]

dt <- dt[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
dt[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]

# ICD code → short disease name for output file naming (overall outcomes only)
icd_to_name <- c(
  "I00_I99" = "CVD",
  "J40_J47_J60_J70_J84_D86_C34" = "CRD",
  "K70_K76_C22" = "CLD",
  "N00_N29_C64_C65" = "CKD",
  "X60_X84_Y87.0" = "Suicide",
  "G20_G30_G12.2_F01_F03" = "NDD",
  "C00_C97" = "Cancer"
)

disease_label <- function(icd) {
  nm <- icd_to_name[icd]
  if (is.na(nm)) icd else unname(nm)
}

# Single representative scenario: joint baseline exposure -> AAMR 2011-2015 (lag 10)
scenario_list <- list(
  list(key = "AirSVI_AAMR2011_2015", aamr = "2011-2015", lag = 10)
)

# 6 control covariates and their abbreviations for model labels
all_covariates <- c(
  "Smoking_rate", "Physical_Activities_rate", "Obesity_rate",
  "Uninsured_rate", "Physician_Density_per100k", "Diabetes_Prevalence_rate"
)
covar_abbrev <- c(
  Smoking_rate = "SM", Physical_Activities_rate = "PA", Obesity_rate = "OB",
  Uninsured_rate = "UN", Physician_Density_per100k = "PD", Diabetes_Prevalence_rate = "DB"
)

# All 2^6 = 64 covariate subsets (k = 0 to 6)
combo_list <- list(character(0)) # k=0: joint only
for (k in seq_along(all_covariates)) {
  combo_list <- c(combo_list, combn(all_covariates, k, simplify = FALSE))
}
message(
  "Total covariate combinations: ", length(combo_list),
  " (k=0:1, k=1:6, k=2:15, k=3:20, k=4:15, k=5:6, k=6:1)"
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

# p_posterior: two-sided posterior tail-area probability with Jeffreys correction.
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
  sprintf("%.4f", 2 * min(p_pos, p_neg))
}

# Per joint-level (vs reference) MRD estimates from treatment contrasts.
extract_joint <- function(draw_df, names_vec) {
  out <- list()
  for (lv in JOINT_LEVELS) {
    if (lv == REF_LEVEL) {
      out[[lv]] <- "0.00"
      next
    }
    nm <- make.names(paste0("Joint", lv))
    idx <- match(nm, names_vec)
    out[[lv]] <- if (is.na(idx)) "" else format_cell(draw_df[[paste0("beta[", idx, "]")]])
  }
  out
}

# Per joint-level p / R-hat / ESS diagnostics (flat, formatted, ready to splice).
extract_joint_metrics <- function(draw_df, names_vec, summ_df) {
  out <- list()
  for (lv in NON_REF) {
    p <- NA_character_; rhat <- NA_real_; eb <- NA_real_; et <- NA_real_
    nm <- make.names(paste0("Joint", lv))
    idx <- match(nm, names_vec)
    if (!is.na(idx)) {
      col <- paste0("beta[", idx, "]")
      p <- compute_p(draw_df[[col]])
      sr <- summ_df[summ_df$variable == col, , drop = FALSE]
      if (nrow(sr)) {
        rhat <- sr$rhat; eb <- sr$ess_bulk; et <- sr$ess_tail
      }
    }
    out[[paste0(lv, "_p")]]        <- p
    out[[paste0(lv, "_rhat")]]     <- sprintf("%.4f", rhat)
    out[[paste0(lv, "_ess_bulk")]] <- as.integer(round(eb))
    out[[paste0(lv, "_ess_tail")]] <- as.integer(round(et))
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

covar_block <- function(draws, names_vec, summ) {
  out <- list()
  for (cn in all_covariates) {
    ab <- unname(covar_abbrev[cn])
    cc <- extract_covariate(draws, names_vec, cn, summ)
    out[[ab]]               <- cc$est
    out[[paste0(ab, "_p")]] <- cc$p
    out[[paste0(ab, "_rhat")]]     <- sprintf("%.4f", cc$rhat)
    out[[paste0(ab, "_ess_bulk")]] <- as.integer(round(cc$ess_bulk))
    out[[paste0(ab, "_ess_tail")]] <- as.integer(round(cc$ess_tail))
  }
  out
}

# Build design: intercept + optional covariates + Joint factor (treatment, ref a1s1).
build_design_combo <- function(d, covariates) {
  d <- d %>% mutate(Joint = factor(.data[[JOINT_COL]], levels = JOINT_LEVELS))
  needed <- c("Joint", "AAMR_Lower", "AAMR_Upper", "cens", "State_FIPS", covariates)
  d <- d[complete.cases(d[, ..needed]), ]
  form <- if (length(covariates) == 0) {
    ~Joint
  } else {
    as.formula(paste("~", paste(c(covariates, "Joint"), collapse = " + ")))
  }
  mm <- model.matrix(form, d, contrasts.arg = list(Joint = contr.treatment(JOINT_LEVELS)))
  colnames(mm) <- make.names(colnames(mm))
  list(X = mm, names = colnames(mm), df = d)
}

for (outcome in selected) {
  dlabel <- disease_label(outcome)
  message("===== Disease: ", dlabel, " (", outcome, ") =====")
  outfile <- file.path(out_dir, paste0(dlabel, "_SVI_Air_22_Sensitivity_Combination.csv"))
  for (sc in scenario_list) {
    scen_key <- sc$key
    aamr_p <- sc$aamr
    lagv <- sc$lag
    scen_dt <- dt[Time_Period == aamr_p & Outcome == outcome]
    if (nrow(scen_dt) < opt$`min-n`) {
      message("[Skip] Scenario ", scen_key, " overall n=", nrow(scen_dt))
      next
    }
    aamr_out <- gsub("-", "_", aamr_p)

    for (combo in combo_list) {
      model_label <- if (length(combo) == 0) "AirSVI" else paste(c("AirSVI", covar_abbrev[combo]), collapse = "+")

      des <- build_design_combo(scen_dt, combo)
      if (nrow(des$df) < opt$`min-n`) {
        message("[Skip] ", scen_key, " ", model_label, " n=", nrow(des$df))
        next
      }

      states_c <- sort(unique(des$df$State_FIPS))
      state_index_c <- match(des$df$State_FIPS, states_c)
      data_list <- list(
        N = nrow(des$df), S = length(states_c), state = state_index_c,
        y_lower = des$df$AAMR_Lower, y_upper = des$df$AAMR_Upper, cens = des$df$cens,
        K = ncol(des$X), X = des$X
      )
      init_fn <- function() list(beta = rep(0, data_list$K), z_u = rep(0, data_list$S), sigma = 50, sigma_u = 10)

      fit <- try(mod$sample(
        data = data_list, chains = opt$chains, iter_sampling = opt$iter - opt$warmup, iter_warmup = opt$warmup,
        adapt_delta = opt$`adapt-delta`, max_treedepth = opt$`max-treedepth`,
        parallel_chains = opt$chains, refresh = 0, seed = opt$seed,
        init = rep(list(init_fn()), opt$chains)
      ), silent = TRUE)

      if (inherits(fit, "try-error")) {
        message("[Fail] ", scen_key, " ", model_label)
        next
      }

      draws <- as_draws_df(fit$draws("beta"))
      colnames(draws) <- paste0("beta[", seq_len(ncol(draws)), "]")
      summ <- posterior::summarize_draws(fit$draws("beta"))

      row_list <- c(
        list(ICD_Code = outcome, AAMR_Period = aamr_out, Lag = lagv, Model = model_label),
        extract_joint(draws, des$names),
        extract_joint_metrics(draws, des$names, summ),
        covar_block(draws, des$names, summ)
      )
      row <- do.call(tibble, row_list)
      append_rows(outfile, row)
      message("[OK] ", scen_key, " ", model_label)
    }
  }
  message("===== Completed: ", dlabel, " =====")
}
message("All requested analyses complete. Output directory: ", out_dir)
