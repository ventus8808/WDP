#!/usr/bin/env Rscript
# Interval-censored mixed model pipeline stratified by demographic factors.
# Stratification variables: Stratum (Male, Female, White, Black, Others)
# EQI 2000-2005 only (4 lag scenarios: 5, 10, 15, 20 years).

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
  "EQI", "EQI_Air", "EQI_Water", "EQI_Land", "EQI_Built", "EQI_Social", "State_FIPS",
  "Smoking_rate", "Physical_Activities_rate", "Obesity_rate"
))

option_list <- list(
  make_option(c("--data"), type = "character", default = "Data/Processed/df_Stratified.csv"),
  make_option(c("--output-dir"), type = "character", default = "Result/brms_Stratified_Demo"),
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

# 种子设置保留
set.seed(opt$seed)

# ── Stan model ─────────────────────────────────────────────────────────────────
stan_code <- "data {\n  int<lower=1> N;\n  int<lower=1> S;\n  array[N] int<lower=1,upper=S> state;\n  vector[N] y_lower;\n  vector[N] y_upper;\n  array[N] int<lower=0,upper=2> cens;\n  int<lower=1> K;\n  matrix[N,K] X;\n} \nparameters {\n  vector[K] beta;\n  vector[S] z_u;\n  real<lower=0> sigma;\n  real<lower=0> sigma_u;\n} \ntransformed parameters {\n  vector[S] u = sigma_u * z_u;\n} \nmodel {\n  beta ~ normal(0,5);\n  z_u ~ normal(0,1);\n  sigma ~ exponential(1);\n  sigma_u ~ exponential(1);\n  for (i in 1:N) {\n    real mu = X[i] * beta + u[state[i]];\n    if (cens[i]==0) {\n      target += normal_lpdf(y_lower[i] | mu, sigma);\n    } else {\n      real p_up = normal_cdf(y_upper[i] | mu, sigma);\n      real p_lo = normal_cdf(y_lower[i] | mu, sigma);\n      real diff = fmax(p_up - p_lo, 1e-12);\n      target += log(diff);\n    }\n  }\n}"
stan_file <- file.path(tempdir(), "interval_mixed_model.stan")
writeLines(stan_code, stan_file)
mod <- cmdstan_model(stan_file)

project_root <- normalizePath(".")
path <- file.path(project_root, opt$data)
if (!file.exists(path)) stop("Data not found: ", path)
dt <- fread(path)

req <- c(
  "COUNTY_FIPS", "EQI_Period", "Time_Period", "Lag_Years", "Outcome", "Stratum",
  "AAMR_Lower", "AAMR_Upper",
  "EQI", "EQI_Air", "EQI_Water", "EQI_Land", "EQI_Built", "EQI_Social",
  "Smoking_rate", "Physical_Activities_rate", "Obesity_rate"
)
miss <- setdiff(req, names(dt))
if (length(miss)) stop("Missing cols: ", paste(miss, collapse = ","))

if (!"State_FIPS" %in% names(dt)) dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]

dt <- dt[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
dt[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]

# EQI 2000-2005 only, 4 lag scenarios
scenario_list <- list(
  list(key = "EQI0005_AAMR2006_2010", eqi = "2000-2005", aamr = "2006-2010", lag = 5),
  list(key = "EQI0005_AAMR2011_2015", eqi = "2000-2005", aamr = "2011-2015", lag = 10),
  list(key = "EQI0005_AAMR2016_2020", eqi = "2000-2005", aamr = "2016-2020", lag = 15),
  list(key = "EQI0005_AAMR2021_2024", eqi = "2000-2005", aamr = "2021-2024", lag = 20)
)

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

extract_quintile_metrics <- function(draw_df, names_vec, prefix, summ_df) {
  out <- list(
    Q2_p = NA_real_, Q3_p = NA_real_, Q4_p = NA_real_, Q5_p = NA_real_,
    Q2_rhat = NA_real_, Q3_rhat = NA_real_, Q4_rhat = NA_real_, Q5_rhat = NA_real_,
    Q2_ess_bulk = NA_real_, Q3_ess_bulk = NA_real_, Q4_ess_bulk = NA_real_, Q5_ess_bulk = NA_real_,
    Q2_ess_tail = NA_real_, Q3_ess_tail = NA_real_, Q4_ess_tail = NA_real_, Q5_ess_tail = NA_real_
  )
  if (any(grepl(paste0(prefix, "\\.L"), names_vec))) {
    return(out)
  }
  for (q in 2:5) {
    nm <- paste0(prefix, q)
    idx <- match(nm, names_vec)
    if (!is.na(idx)) {
      col <- paste0("beta[", idx, "]")
      draws_col <- draw_df[[col]]
      out[[paste0("Q", q, "_p")]] <- compute_p(draws_col)
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

extract_quintiles <- function(draw_df, names_vec, prefix) {
  out <- list(Q1 = "0.00", Q2 = "", Q3 = "", Q4 = "", Q5 = "")
  if (any(grepl(paste0(prefix, "\\.L"), names_vec))) {
    return(out)
  }
  for (q in 2:5) {
    nm <- paste0(prefix, q)
    idx <- match(nm, names_vec)
    out[[paste0("Q", q)]] <- if (is.na(idx)) "" else format_cell(draw_df[[paste0("beta[", idx, "]")]])
  }
  out
}

compute_mu_q1 <- function(draw_df, names_vec, d) {
  sm_idx <- match("Smoking_rate", names_vec)
  pa_idx <- match("Physical_Activities_rate", names_vec)
  ob_idx <- match("Obesity_rate", names_vec)
  if (any(is.na(c(sm_idx, pa_idx, ob_idx)))) {
    return(NULL)
  }
  draw_df[["beta[1]"]] +
    draw_df[[paste0("beta[", sm_idx, "]")]] * mean(d$Smoking_rate, na.rm = TRUE) +
    draw_df[[paste0("beta[", pa_idx, "]")]] * mean(d$Physical_Activities_rate, na.rm = TRUE) +
    draw_df[[paste0("beta[", ob_idx, "]")]] * mean(d$Obesity_rate, na.rm = TRUE)
}

extract_mrr <- function(draw_df, names_vec, prefix, d) {
  out_q <- list(Q1 = "", Q2 = "", Q3 = "", Q4 = "", Q5 = "")
  out_p <- list(Q2_p = NA_character_, Q3_p = NA_character_, Q4_p = NA_character_, Q5_p = NA_character_)
  if (any(grepl(paste0(prefix, "\\.L"), names_vec))) {
    return(list(quintiles = out_q, pvals = out_p))
  }

  mu_q1 <- compute_mu_q1(draw_df, names_vec, d)
  if (is.null(mu_q1)) {
    return(list(quintiles = out_q, pvals = out_p))
  }

  q_draws <- vector("list", 5)
  q_draws[[1]] <- mu_q1 / mu_q1
  out_q$Q1 <- format_cell(q_draws[[1]])
  for (q in 2:5) {
    idx <- match(paste0(prefix, q), names_vec)
    if (is.na(idx)) {
      next
    }
    q_draws[[q]] <- (mu_q1 + draw_df[[paste0("beta[", idx, "]")]]) / mu_q1
    out_q[[paste0("Q", q)]] <- format_cell(q_draws[[q]])
    out_p[[paste0("Q", q, "_p")]] <- compute_p(q_draws[[q]] - 1)
  }
  list(quintiles = out_q, pvals = out_p)
}

build_design_overall <- function(d) {
  d <- d %>% mutate(EQI_factor = factor(EQI, levels = 1:5))
  d <- d[complete.cases(d[, c(
    "EQI_factor", "AAMR_Lower", "AAMR_Upper", "cens", "State_FIPS",
    "Smoking_rate", "Physical_Activities_rate", "Obesity_rate"
  )]), ]
  mm <- model.matrix(~ Smoking_rate + Physical_Activities_rate + Obesity_rate + EQI_factor, d,
    contrasts.arg = list(EQI_factor = contr.treatment(5))
  )
  colnames(mm) <- make.names(colnames(mm))
  list(X = mm, names = colnames(mm), df = d)
}

build_design_multi <- function(d) {
  d <- d %>% mutate(
    EQI_Air_factor    = factor(EQI_Air, levels = 1:5),
    EQI_Water_factor  = factor(EQI_Water, levels = 1:5),
    EQI_Land_factor   = factor(EQI_Land, levels = 1:5),
    EQI_Built_factor  = factor(EQI_Built, levels = 1:5),
    EQI_Social_factor = factor(EQI_Social, levels = 1:5)
  )
  d <- d[complete.cases(d[, c(
    "EQI_Air_factor", "EQI_Water_factor", "EQI_Land_factor", "EQI_Built_factor", "EQI_Social_factor",
    "AAMR_Lower", "AAMR_Upper", "cens", "State_FIPS",
    "Smoking_rate", "Physical_Activities_rate", "Obesity_rate"
  )]), ]
  form <- as.formula("~ Smoking_rate + Physical_Activities_rate + Obesity_rate + EQI_Air_factor + EQI_Water_factor + EQI_Land_factor + EQI_Built_factor + EQI_Social_factor")
  mm <- model.matrix(form, d,
    contrasts.arg = list(
      EQI_Air_factor    = contr.treatment(5), EQI_Water_factor = contr.treatment(5),
      EQI_Land_factor   = contr.treatment(5), EQI_Built_factor = contr.treatment(5),
      EQI_Social_factor = contr.treatment(5)
    )
  )
  colnames(mm) <- make.names(colnames(mm))
  list(X = mm, names = colnames(mm), df = d)
}

for (outcome in selected) {
  dlabel <- disease_label(outcome)
  message("===== Outcome: ", dlabel, " (", outcome, ") =====")
  outfile_mrd <- file.path(out_dir, paste0(dlabel, "_Stratified_Demo_MRD.csv"))
  outfile_mrr <- file.path(out_dir, paste0(dlabel, "_Stratified_Demo_MRR.csv"))

  strata_vals <- sort(unique(na.omit(dt[Outcome == outcome]$Stratum)))
  message("  Strata: ", paste(strata_vals, collapse = ", "))

  for (stratum in strata_vals) {
    strat_dt <- dt[Outcome == outcome & Stratum == stratum]
    message("  Stratum: ", stratum, " (n=", nrow(strat_dt), ")")

    for (sc in scenario_list) {
      scen_key <- sc$key
      eqi_p <- sc$eqi
      aamr_p <- sc$aamr
      lagv <- sc$lag
      scen_dt <- strat_dt[EQI_Period == eqi_p & Time_Period == aamr_p]
      if (nrow(scen_dt) < opt$`min-n`) {
        message("[Skip] ", stratum, " ", scen_key, " n=", nrow(scen_dt))
        next
      }
      eqi_out <- gsub("-", "_", eqi_p)
      aamr_out <- gsub("-", "_", aamr_p)

      # Overall EQI + SM + PA + OB
      des_overall <- build_design_overall(scen_dt)
      if (nrow(des_overall$df) < opt$`min-n`) {
        message("[Skip] ", stratum, " ", scen_key, " EQI n=", nrow(des_overall$df))
      } else {
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
          adapt_delta = opt$`adapt-delta`, max_treedepth = opt$`max-treedepth`,
          parallel_chains = opt$chains, refresh = 0, seed = opt$seed,
          init = rep(list(init_fun()), opt$chains)
        ), silent = TRUE)
        if (inherits(fit_overall, "try-error")) {
          message("[Fail] ", stratum, " ", scen_key, " EQI+SM+PA+OB")
        } else {
          draws <- as_draws_df(fit_overall$draws("beta"))
          colnames(draws) <- paste0("beta[", seq_len(ncol(draws)), "]")
          q_over <- extract_quintiles(draws, des_overall$names, "EQI_factor")
          mrr_over <- extract_mrr(draws, des_overall$names, "EQI_factor", des_overall$df)
          summ_over <- posterior::summarize_draws(fit_overall$draws("beta"))
          met_over <- extract_quintile_metrics(draws, des_overall$names, "EQI_factor", summ_over)
          row_over_mrd <- tibble(
            ICD_Code = outcome, EQI_Period = eqi_out, AAMR_Period = aamr_out, Lag = lagv,
            Model = paste0(stratum, "_EQI+SM+PA+OB"),
            Q1 = q_over$Q1, Q2 = q_over$Q2, Q3 = q_over$Q3, Q4 = q_over$Q4, Q5 = q_over$Q5,
            Q2_p = met_over$Q2_p, Q3_p = met_over$Q3_p, Q4_p = met_over$Q4_p, Q5_p = met_over$Q5_p,
            Q2_rhat = sprintf("%.4f", met_over$Q2_rhat), Q3_rhat = sprintf("%.4f", met_over$Q3_rhat), Q4_rhat = sprintf("%.4f", met_over$Q4_rhat), Q5_rhat = sprintf("%.4f", met_over$Q5_rhat),
            Q2_ess_bulk = as.integer(round(met_over$Q2_ess_bulk)), Q3_ess_bulk = as.integer(round(met_over$Q3_ess_bulk)), Q4_ess_bulk = as.integer(round(met_over$Q4_ess_bulk)), Q5_ess_bulk = as.integer(round(met_over$Q5_ess_bulk)),
            Q2_ess_tail = as.integer(round(met_over$Q2_ess_tail)), Q3_ess_tail = as.integer(round(met_over$Q3_ess_tail)), Q4_ess_tail = as.integer(round(met_over$Q4_ess_tail)), Q5_ess_tail = as.integer(round(met_over$Q5_ess_tail))
          )
          row_over_mrr <- tibble(
            ICD_Code = outcome, EQI_Period = eqi_out, AAMR_Period = aamr_out, Lag = lagv,
            Model = paste0(stratum, "_EQI+SM+PA+OB"),
            Q1 = mrr_over$quintiles$Q1, Q2 = mrr_over$quintiles$Q2, Q3 = mrr_over$quintiles$Q3, Q4 = mrr_over$quintiles$Q4, Q5 = mrr_over$quintiles$Q5,
            Q2_p = mrr_over$pvals$Q2_p, Q3_p = mrr_over$pvals$Q3_p, Q4_p = mrr_over$pvals$Q4_p, Q5_p = mrr_over$pvals$Q5_p,
            Q2_rhat = sprintf("%.4f", met_over$Q2_rhat), Q3_rhat = sprintf("%.4f", met_over$Q3_rhat), Q4_rhat = sprintf("%.4f", met_over$Q4_rhat), Q5_rhat = sprintf("%.4f", met_over$Q5_rhat),
            Q2_ess_bulk = as.integer(round(met_over$Q2_ess_bulk)), Q3_ess_bulk = as.integer(round(met_over$Q3_ess_bulk)), Q4_ess_bulk = as.integer(round(met_over$Q4_ess_bulk)), Q5_ess_bulk = as.integer(round(met_over$Q5_ess_bulk)),
            Q2_ess_tail = as.integer(round(met_over$Q2_ess_tail)), Q3_ess_tail = as.integer(round(met_over$Q3_ess_tail)), Q4_ess_tail = as.integer(round(met_over$Q4_ess_tail)), Q5_ess_tail = as.integer(round(met_over$Q5_ess_tail))
          )
          append_rows(outfile_mrd, row_over_mrd)
          append_rows(outfile_mrr, row_over_mrr)
          message("[OK] ", stratum, " ", scen_key, " EQI+SM+PA+OB")
        }
      }

      # Multi-domain EQI sub-indices + SM + PA + OB
      des_multi <- build_design_multi(scen_dt)
      if (nrow(des_multi$df) < opt$`min-n`) {
        message("[Skip] ", stratum, " ", scen_key, " Multi-domain n=", nrow(des_multi$df))
      } else {
        states_m <- sort(unique(des_multi$df$State_FIPS))
        state_index_m <- match(des_multi$df$State_FIPS, states_m)
        data_list2 <- list(
          N = nrow(des_multi$df), S = length(states_m), state = state_index_m,
          y_lower = des_multi$df$AAMR_Lower, y_upper = des_multi$df$AAMR_Upper, cens = des_multi$df$cens,
          K = ncol(des_multi$X), X = des_multi$X
        )
        init_fun2 <- function() list(beta = rep(0, data_list2$K), z_u = rep(0, data_list2$S), sigma = 50, sigma_u = 10)
        fit_multi <- try(mod$sample(
          data = data_list2, chains = opt$chains, iter_sampling = opt$iter - opt$warmup, iter_warmup = opt$warmup,
          adapt_delta = opt$`adapt-delta`, max_treedepth = opt$`max-treedepth`,
          parallel_chains = opt$chains, refresh = 0, seed = opt$seed,
          init = rep(list(init_fun2()), opt$chains)
        ), silent = TRUE)
        if (inherits(fit_multi, "try-error")) {
          message("[Fail] ", stratum, " ", scen_key, " Multi-domain")
        } else {
          draws_m <- as_draws_df(fit_multi$draws("beta"))
          colnames(draws_m) <- paste0("beta[", seq_len(ncol(draws_m)), "]")
          domain_prefix <- c("EQI_Air_factor", "EQI_Water_factor", "EQI_Land_factor", "EQI_Built_factor", "EQI_Social_factor")
          summ_m <- posterior::summarize_draws(fit_multi$draws("beta"))
          for (dom in domain_prefix) {
            qd <- extract_quintiles(draws_m, des_multi$names, dom)
            md_mrr <- extract_mrr(draws_m, des_multi$names, dom, des_multi$df)
            md <- extract_quintile_metrics(draws_m, des_multi$names, dom, summ_m)
            row_dom_mrd <- tibble(
              ICD_Code = outcome, EQI_Period = eqi_out, AAMR_Period = aamr_out, Lag = lagv,
              Model = paste0(stratum, "_", sub("_factor", "", dom), "+SM+PA+OB"),
              Q1 = qd$Q1, Q2 = qd$Q2, Q3 = qd$Q3, Q4 = qd$Q4, Q5 = qd$Q5,
              Q2_p = md$Q2_p, Q3_p = md$Q3_p, Q4_p = md$Q4_p, Q5_p = md$Q5_p,
              Q2_rhat = sprintf("%.4f", md$Q2_rhat), Q3_rhat = sprintf("%.4f", md$Q3_rhat), Q4_rhat = sprintf("%.4f", md$Q4_rhat), Q5_rhat = sprintf("%.4f", md$Q5_rhat),
              Q2_ess_bulk = as.integer(round(md$Q2_ess_bulk)), Q3_ess_bulk = as.integer(round(md$Q3_ess_bulk)), Q4_ess_bulk = as.integer(round(md$Q4_ess_bulk)), Q5_ess_bulk = as.integer(round(md$Q5_ess_bulk)),
              Q2_ess_tail = as.integer(round(md$Q2_ess_tail)), Q3_ess_tail = as.integer(round(md$Q3_ess_tail)), Q4_ess_tail = as.integer(round(md$Q4_ess_tail)), Q5_ess_tail = as.integer(round(md$Q5_ess_tail))
            )
            row_dom_mrr <- tibble(
              ICD_Code = outcome, EQI_Period = eqi_out, AAMR_Period = aamr_out, Lag = lagv,
              Model = paste0(stratum, "_", sub("_factor", "", dom), "+SM+PA+OB"),
              Q1 = md_mrr$quintiles$Q1, Q2 = md_mrr$quintiles$Q2, Q3 = md_mrr$quintiles$Q3, Q4 = md_mrr$quintiles$Q4, Q5 = md_mrr$quintiles$Q5,
              Q2_p = md_mrr$pvals$Q2_p, Q3_p = md_mrr$pvals$Q3_p, Q4_p = md_mrr$pvals$Q4_p, Q5_p = md_mrr$pvals$Q5_p,
              Q2_rhat = sprintf("%.4f", md$Q2_rhat), Q3_rhat = sprintf("%.4f", md$Q3_rhat), Q4_rhat = sprintf("%.4f", md$Q4_rhat), Q5_rhat = sprintf("%.4f", md$Q5_rhat),
              Q2_ess_bulk = as.integer(round(md$Q2_ess_bulk)), Q3_ess_bulk = as.integer(round(md$Q3_ess_bulk)), Q4_ess_bulk = as.integer(round(md$Q4_ess_bulk)), Q5_ess_bulk = as.integer(round(md$Q5_ess_bulk)),
              Q2_ess_tail = as.integer(round(md$Q2_ess_tail)), Q3_ess_tail = as.integer(round(md$Q3_ess_tail)), Q4_ess_tail = as.integer(round(md$Q4_ess_tail)), Q5_ess_tail = as.integer(round(md$Q5_ess_tail))
            )
            append_rows(outfile_mrd, row_dom_mrd)
            append_rows(outfile_mrr, row_dom_mrr)
          }
          message("[OK] ", stratum, " ", scen_key, " Multi-domain+SM+PA+OB")
        }
      }
    }
  }
  message("===== Completed: ", outcome, " =====")
}
message("All requested analyses complete. Output directory: ", out_dir)
