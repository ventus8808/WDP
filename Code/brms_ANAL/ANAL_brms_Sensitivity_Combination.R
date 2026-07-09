#!/usr/bin/env Rscript
# ANAL sensitivity-combination analysis for overall disease outcomes.
# Design:
#   - Outcomes: 7 overall disease groups only.
#   - Exposures: four ANAL metrics, one at a time:
#       popw_mean_rad, mean_rad, sol, lit_area_km2
#   - Exposure coding: period-specific quintiles within ANAL_Period, Q1 reference.
#   - Pooled model: default pooled 10-year lag; Matched_Pair is included as a
#     fixed effect so the pooled estimate is adjusted for exposure/outcome window.
#   - Covariates: baseline exposure-only, then all +1, +2, +3 combinations from
#     the same six controls used by brms_Sensitivity_Combination.R.
#
# Output: one CSV per disease x exposure metric.

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
  "State_FIPS", "Matched_Pair", "Exposure_factor", "Matched_Pair_factor",
  "Smoking_rate", "Physical_Activities_rate", "Obesity_rate",
  "Uninsured_rate", "Physician_Density_per100k", "Diabetes_Prevalence_rate"
))

option_list <- list(
  make_option(c("--data"), type = "character", default = "Data/Processed/df_ANAL.csv"),
  make_option(c("--output-dir"), type = "character", default = "Result/brms_ANAL_Sensitivity_Combination"),
  make_option(c("--outcomes"), type = "character", default = NA,
              help = "Comma-separated ICD outcomes. Default: 7 overall disease groups."),
  make_option(c("--exposures"), type = "character", default = "popw_mean_rad,mean_rad,sol,lit_area_km2",
              help = "Comma-separated ANAL metrics."),
  make_option(c("--lag"), type = "integer", default = 10,
              help = "Pooled lag to run. Default: 10."),
  make_option(c("--max-covariates"), type = "integer", default = 3,
              help = "Maximum covariates in a combination. Default: 3."),
  make_option(c("--chains"), type = "integer", default = 6),
  make_option(c("--iter"), type = "integer", default = 1800),
  make_option(c("--warmup"), type = "integer", default = 1000),
  make_option(c("--adapt-delta"), type = "double", default = 0.95),
  make_option(c("--max-treedepth"), type = "integer", default = 12),
  make_option(c("--min-n"), type = "integer", default = 50),
  make_option(c("--seed"), type = "integer", default = 1234),
  make_option(c("--overwrite"), action = "store_true", default = FALSE),
  make_option(c("--test"), action = "store_true", default = FALSE)
)
opt <- parse_args(OptionParser(option_list = option_list))

if (opt$`max-covariates` < 0 || opt$`max-covariates` > 6) {
  stop("--max-covariates must be between 0 and 6")
}

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
stan_file <- file.path(tempdir(), "anal_interval_mixed_model.stan")
writeLines(stan_code, stan_file)
mod <- cmdstan_model(stan_file)

# ── Load data ──────────────────────────────────────────────────────────────────
project_root <- normalizePath(".")
path <- file.path(project_root, opt$data)
if (!file.exists(path)) stop("Data not found: ", path)
dt <- fread(path)

exposure_metrics <- str_split(opt$exposures, ",", simplify = TRUE) |>
  as.vector() |>
  str_trim()

all_covariates <- c(
  "Smoking_rate", "Physical_Activities_rate", "Obesity_rate",
  "Uninsured_rate", "Physician_Density_per100k", "Diabetes_Prevalence_rate"
)
covar_abbrev <- c(
  Smoking_rate = "SM", Physical_Activities_rate = "PA", Obesity_rate = "OB",
  Uninsured_rate = "UN", Physician_Density_per100k = "PD", Diabetes_Prevalence_rate = "DB"
)

req <- c(
  "COUNTY_FIPS", "ANAL_Period", "AAMR_Period", "Time_Period", "Lag_Years",
  "Matched_Pair", "Outcome", "AAMR_Lower", "AAMR_Upper",
  exposure_metrics, all_covariates
)
miss <- setdiff(req, names(dt))
if (length(miss)) stop("Missing cols: ", paste(miss, collapse = ","))

if (!"State_FIPS" %in% names(dt)) {
  dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]
}
dt[, State_FIPS := as.character(State_FIPS)]

dt <- dt[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
dt[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]
dt <- dt[Lag_Years == opt$lag]
if (nrow(dt) == 0) stop("No rows after filtering Lag_Years == ", opt$lag)

# ── Overall disease outcomes ───────────────────────────────────────────────────
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

overall_outcomes <- names(icd_to_name)
selected <- if (is.na(opt$outcomes)) {
  overall_outcomes
} else {
  str_split(opt$outcomes, ",", simplify = TRUE) |>
    as.vector() |>
    str_trim()
}
missing_outcomes <- setdiff(selected, unique(dt$Outcome))
if (length(missing_outcomes)) {
  stop("Requested outcomes not present after lag filter: ", paste(missing_outcomes, collapse = ","))
}

if (opt$test) {
  selected <- selected[1]
  exposure_metrics <- exposure_metrics[1]
  opt$chains <- min(opt$chains, 2)
  opt$iter <- min(opt$iter, 400)
  opt$warmup <- min(opt$warmup, 200)
  message("Test mode: first outcome/exposure only with reduced iterations.")
}

# ── Covariate combinations: k = 0, 1, 2, 3 by default ─────────────────────────
combo_list <- list(character(0))
if (opt$`max-covariates` >= 1) {
  for (k in seq_len(opt$`max-covariates`)) {
    combo_list <- c(combo_list, combn(all_covariates, k, simplify = FALSE))
  }
}
message(
  "Total covariate combinations: ", length(combo_list),
  " (k=0 through k=", opt$`max-covariates`, ")"
)
message("Outcomes: ", paste(selected, collapse = ","))
message("Exposures: ", paste(exposure_metrics, collapse = ","))
message("Pooled lag: ", opt$lag, " years")

out_dir <- file.path(project_root, opt$`output-dir`)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

# ── Helpers ────────────────────────────────────────────────────────────────────
add_period_quintile <- function(d, metric) {
  q_col <- paste0(metric, "_Quintile")
  if (metric == "popw_mean_rad" && "ANAL_Quintile" %in% names(d)) {
    d[, (q_col) := as.integer(ANAL_Quintile)]
    return(q_col)
  }

  lookup <- unique(d[, .(COUNTY_FIPS, ANAL_Period, exposure_value = get(metric))])
  lookup[!is.na(exposure_value),
         q_tmp := pmin(5L, as.integer(floor((frank(exposure_value, ties.method = "first") - 1) * 5 / .N) + 1L)),
         by = ANAL_Period]
  lookup[, q_tmp := as.integer(q_tmp)]
  setnames(lookup, "q_tmp", q_col)
  lookup <- lookup[, .(COUNTY_FIPS, ANAL_Period, q_value = get(q_col))]
  setnames(lookup, "q_value", q_col)
  d[lookup,
    (q_col) := get(paste0("i.", q_col)),
    on = .(COUNTY_FIPS, ANAL_Period)]
  q_col
}

format_cell <- function(draws) {
  if (length(draws) == 0) return("")
  ci <- quantile(draws, c(0.025, 0.975), na.rm = TRUE)
  sprintf("%0.2f(%0.2f,%0.2f)", mean(draws), ci[1], ci[2])
}

append_rows <- function(path, df) {
  if (!file.exists(path)) {
    write_csv(df, path)
  } else {
    suppressWarnings(write.table(df, path, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE))
  }
}

compute_p <- function(draws) {
  if (length(draws) == 0) return(NA_character_)
  pos <- sum(draws > 0, na.rm = TRUE)
  neg <- sum(draws < 0, na.rm = TRUE)
  n <- pos + neg
  if (n == 0) return(NA_character_)
  p_pos <- (pos + 0.5) / (n + 1)
  p_neg <- (neg + 0.5) / (n + 1)
  sprintf("%.4f", 2 * min(p_pos, p_neg))
}

extract_quintiles <- function(draw_df, names_vec, prefix) {
  out <- list(Q1 = "0.00", Q2 = "", Q3 = "", Q4 = "", Q5 = "")
  if (any(grepl(paste0(prefix, "\\.L"), names_vec))) return(out)
  for (q in 2:5) {
    nm <- paste0(prefix, q)
    idx <- match(nm, names_vec)
    out[[paste0("Q", q)]] <- if (is.na(idx)) "" else format_cell(draw_df[[paste0("beta[", idx, "]")]])
  }
  out
}

extract_quintile_metrics <- function(draw_df, names_vec, prefix, summ_df) {
  out <- list(
    Q2_p = NA_real_, Q3_p = NA_real_, Q4_p = NA_real_, Q5_p = NA_real_,
    Q2_rhat = NA_real_, Q3_rhat = NA_real_, Q4_rhat = NA_real_, Q5_rhat = NA_real_,
    Q2_ess_bulk = NA_real_, Q3_ess_bulk = NA_real_, Q4_ess_bulk = NA_real_, Q5_ess_bulk = NA_real_,
    Q2_ess_tail = NA_real_, Q3_ess_tail = NA_real_, Q4_ess_tail = NA_real_, Q5_ess_tail = NA_real_
  )
  if (any(grepl(paste0(prefix, "\\.L"), names_vec))) return(out)
  for (q in 2:5) {
    nm <- paste0(prefix, q)
    idx <- match(nm, names_vec)
    if (!is.na(idx)) {
      col <- paste0("beta[", idx, "]")
      out[[paste0("Q", q, "_p")]] <- compute_p(draw_df[[col]])
      sr <- summ_df[summ_df$variable == col, , drop = FALSE]
      if (nrow(sr)) {
        out[[paste0("Q", q, "_rhat")]] <- sr$rhat
        out[[paste0("Q", q, "_ess_bulk")]] <- sr$ess_bulk
        out[[paste0("Q", q, "_ess_tail")]] <- sr$ess_tail
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
  sr <- summ_df[summ_df$variable == col, , drop = FALSE]
  list(
    est = format_cell(draw_df[[col]]),
    p = compute_p(draw_df[[col]]),
    rhat = if (nrow(sr)) sr$rhat else NA_real_,
    ess_bulk = if (nrow(sr)) sr$ess_bulk else NA_real_,
    ess_tail = if (nrow(sr)) sr$ess_tail else NA_real_
  )
}

build_design_combo <- function(d, q_col, covariates) {
  d <- copy(d)
  d[, Exposure_factor := factor(get(q_col), levels = 1:5)]
  d[, Matched_Pair_factor := factor(Matched_Pair)]
  needed <- c("Exposure_factor", "Matched_Pair_factor", "AAMR_Lower", "AAMR_Upper", "cens", "State_FIPS", covariates)
  d <- d[complete.cases(d[, ..needed])]

  rhs <- c(covariates, "Matched_Pair_factor", "Exposure_factor")
  form <- as.formula(paste("~", paste(rhs, collapse = " + ")))
  mm <- model.matrix(
    form,
    d,
    contrasts.arg = list(
      Exposure_factor = contr.treatment(5),
      Matched_Pair_factor = contr.treatment(nlevels(d$Matched_Pair_factor))
    )
  )
  colnames(mm) <- make.names(colnames(mm))
  list(X = mm, names = colnames(mm), df = d)
}

# ── Main loop ──────────────────────────────────────────────────────────────────
for (metric in exposure_metrics) {
  if (!metric %in% names(dt)) stop("Exposure metric not found: ", metric)
  q_col <- add_period_quintile(dt, metric)
  exposure_label <- safe_name(metric)

  for (outcome in selected) {
    dlabel <- disease_label(outcome)
    outfile <- file.path(out_dir, paste0(dlabel, "_", exposure_label, "_ANAL_Sensitivity_Combination.csv"))
    if (opt$overwrite && file.exists(outfile)) {
      unlink(outfile)
    }
    message("===== Disease: ", dlabel, " (", outcome, "), Exposure: ", metric, " =====")

    scen_dt <- dt[Outcome == outcome]
    if (nrow(scen_dt) < opt$`min-n`) {
      message("[Skip] ", outcome, " n=", nrow(scen_dt))
      next
    }

    pair_list <- paste(sort(unique(scen_dt$Matched_Pair)), collapse = ";")
    aamr_list <- paste(sort(unique(scen_dt$AAMR_Period)), collapse = ";")
    anal_list <- paste(sort(unique(scen_dt$ANAL_Period)), collapse = ";")

    for (combo in combo_list) {
      model_label <- if (length(combo) == 0) {
        "Baseline"
      } else {
        paste0("+", length(combo), "_", paste(covar_abbrev[combo], collapse = "+"))
      }
      formula_label <- paste(c(metric, "Matched_Pair", covar_abbrev[combo]), collapse = "+")

      des <- build_design_combo(scen_dt, q_col, combo)
      if (nrow(des$df) < opt$`min-n`) {
        message("[Skip] ", dlabel, " ", metric, " ", model_label, " n=", nrow(des$df))
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
        data = data_list,
        chains = opt$chains,
        iter_sampling = opt$iter - opt$warmup,
        iter_warmup = opt$warmup,
        adapt_delta = opt$`adapt-delta`,
        max_treedepth = opt$`max-treedepth`,
        parallel_chains = opt$chains,
        refresh = 0,
        seed = opt$seed,
        init = rep(list(init_fn()), opt$chains)
      ), silent = TRUE)

      if (inherits(fit, "try-error")) {
        message("[Fail] ", dlabel, " ", metric, " ", model_label)
        next
      }

      draws <- as_draws_df(fit$draws("beta"))
      colnames(draws) <- paste0("beta[", seq_len(ncol(draws)), "]")
      summ <- posterior::summarize_draws(fit$draws("beta"))

      q_vals <- extract_quintiles(draws, des$names, "Exposure_factor")
      met <- extract_quintile_metrics(draws, des$names, "Exposure_factor", summ)
      sm_c <- extract_covariate(draws, des$names, "Smoking_rate", summ)
      pa_c <- extract_covariate(draws, des$names, "Physical_Activities_rate", summ)
      ob_c <- extract_covariate(draws, des$names, "Obesity_rate", summ)
      un_c <- extract_covariate(draws, des$names, "Uninsured_rate", summ)
      pd_c <- extract_covariate(draws, des$names, "Physician_Density_per100k", summ)
      db_c <- extract_covariate(draws, des$names, "Diabetes_Prevalence_rate", summ)

      row <- tibble(
        ICD_Code = outcome,
        Disease = dlabel,
        Exposure = metric,
        Exposure_Quintile_Column = q_col,
        ANAL_Period = "Pooled",
        AAMR_Period = "Pooled",
        ANAL_Periods = anal_list,
        AAMR_Periods = aamr_list,
        Matched_Pairs = pair_list,
        Lag = opt$lag,
        Model = model_label,
        Formula_Model = formula_label,
        Covariate_Count = length(combo),
        N = nrow(des$df),
        States = length(states_c),
        Pairs = length(unique(des$df$Matched_Pair)),
        Q1 = q_vals$Q1, Q2 = q_vals$Q2, Q3 = q_vals$Q3, Q4 = q_vals$Q4, Q5 = q_vals$Q5,
        Q2_p = met$Q2_p, Q3_p = met$Q3_p, Q4_p = met$Q4_p, Q5_p = met$Q5_p,
        Q2_rhat = sprintf("%.4f", met$Q2_rhat), Q3_rhat = sprintf("%.4f", met$Q3_rhat),
        Q4_rhat = sprintf("%.4f", met$Q4_rhat), Q5_rhat = sprintf("%.4f", met$Q5_rhat),
        Q2_ess_bulk = as.integer(round(met$Q2_ess_bulk)), Q3_ess_bulk = as.integer(round(met$Q3_ess_bulk)),
        Q4_ess_bulk = as.integer(round(met$Q4_ess_bulk)), Q5_ess_bulk = as.integer(round(met$Q5_ess_bulk)),
        Q2_ess_tail = as.integer(round(met$Q2_ess_tail)), Q3_ess_tail = as.integer(round(met$Q3_ess_tail)),
        Q4_ess_tail = as.integer(round(met$Q4_ess_tail)), Q5_ess_tail = as.integer(round(met$Q5_ess_tail)),
        SM = sm_c$est, SM_p = sm_c$p, SM_rhat = sprintf("%.4f", sm_c$rhat), SM_ess_bulk = as.integer(round(sm_c$ess_bulk)), SM_ess_tail = as.integer(round(sm_c$ess_tail)),
        PA = pa_c$est, PA_p = pa_c$p, PA_rhat = sprintf("%.4f", pa_c$rhat), PA_ess_bulk = as.integer(round(pa_c$ess_bulk)), PA_ess_tail = as.integer(round(pa_c$ess_tail)),
        OB = ob_c$est, OB_p = ob_c$p, OB_rhat = sprintf("%.4f", ob_c$rhat), OB_ess_bulk = as.integer(round(ob_c$ess_bulk)), OB_ess_tail = as.integer(round(ob_c$ess_tail)),
        UN = un_c$est, UN_p = un_c$p, UN_rhat = sprintf("%.4f", un_c$rhat), UN_ess_bulk = as.integer(round(un_c$ess_bulk)), UN_ess_tail = as.integer(round(un_c$ess_tail)),
        PD = pd_c$est, PD_p = pd_c$p, PD_rhat = sprintf("%.4f", pd_c$rhat), PD_ess_bulk = as.integer(round(pd_c$ess_bulk)), PD_ess_tail = as.integer(round(pd_c$ess_tail)),
        DB = db_c$est, DB_p = db_c$p, DB_rhat = sprintf("%.4f", db_c$rhat), DB_ess_bulk = as.integer(round(db_c$ess_bulk)), DB_ess_tail = as.integer(round(db_c$ess_tail))
      )
      append_rows(outfile, row)
      message("[OK] ", dlabel, " ", metric, " ", model_label)
    }
    message("===== Completed: ", dlabel, " / ", metric, " =====")
  }
}

message("All requested ANAL analyses complete. Output directory: ", out_dir)
