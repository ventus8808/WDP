#!/usr/bin/env Rscript
# MRR pipeline — EQI+SM+PA+OB, EQI 2000-2005 only (lags 5, 10, 15, 20).
#
# Reference designs:
#   SameRef — universal Lag5 Q1 mu denominator.
#              Q1@Lag5 = 1.0 by construction; Q1@other lags captures secular trend.
#   LagRef  — per-lag Q1 mu denominator (within-lag relative comparison).
#              Q1 = 1.0 at every lag.
#
# Output files per disease (Result/brms_MRR/):
#   {dlabel}_MRD.csv                   — wide, MRD (beta coefficients), one row per lag
#   {dlabel}_MRR_SameRef.csv           — wide, MRR vs universal Lag5 Q1, one row per lag
#   {dlabel}_MRR_LagRef.csv            — wide, MRR vs per-lag Q1, one row per lag
#   {dlabel}_lag_test_MRD.csv          — pairwise Q5 lag comparison (MRD scale)
#   {dlabel}_lag_test_MRR_SameRef.csv  — pairwise Q5 lag comparison (SameRef MRR)
#   {dlabel}_lag_test_MRR_LagRef.csv   — pairwise Q5 lag comparison (LagRef MRR)

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
  "EQI", "State_FIPS",
  "Smoking_rate", "Physical_Activities_rate", "Obesity_rate"
))

option_list <- list(
  make_option(c("--data"), type = "character", default = "Data/Processed/df.csv"),
  make_option(c("--output-dir"), type = "character", default = "Result/brms_Main_MRR"),
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
stan_file <- file.path(tempdir(), "interval_mixed_mrr.stan")
writeLines(stan_code, stan_file)
mod <- cmdstan_model(stan_file)

project_root <- normalizePath(".")
path <- file.path(project_root, opt$data)
if (!file.exists(path)) stop("Data not found: ", path)
dt <- fread(path)

req <- c(
  "COUNTY_FIPS", "EQI_Period", "Time_Period", "Lag_Years", "Outcome",
  "AAMR_Lower", "AAMR_Upper",
  "EQI", "Smoking_rate", "Physical_Activities_rate", "Obesity_rate"
)
miss <- setdiff(req, names(dt))
if (length(miss)) stop("Missing cols: ", paste(miss, collapse = ","))

if (!"State_FIPS" %in% names(dt)) dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]
dt <- dt[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
dt[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]

# EQI 2000-2005 only — 4 lag scenarios
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
  intersect(names(icd_to_name), all_outcomes)
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

# ── Helper functions ───────────────────────────────────────────────────────────
format_cell <- function(draws) {
  if (length(draws) == 0) {
    return("")
  }
  ci <- quantile(draws, c(0.025, 0.975), na.rm = TRUE)
  sprintf("%0.2f(%0.2f,%0.2f)", mean(draws), ci[1], ci[2])
}

format_mrr_cell <- function(draws) {
  if (length(draws) == 0) {
    return("")
  }
  ci <- quantile(draws, c(0.025, 0.975), na.rm = TRUE)
  sprintf("%0.4f(%0.4f,%0.4f)", mean(draws), ci[1], ci[2])
}

append_rows <- function(path, df) {
  if (!file.exists(path)) {
    write_csv(df, path)
  } else {
    suppressWarnings(write.table(df, path, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE))
  }
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
  sprintf("%.4f", 2 * min((pos + 0.5) / (n + 1), (neg + 0.5) / (n + 1)))
}

extract_quintile_metrics <- function(draw_df, names_vec, prefix, summ_df) {
  out <- list(
    Q2_p = NA_character_, Q3_p = NA_character_, Q4_p = NA_character_, Q5_p = NA_character_,
    Q2_rhat = NA_real_, Q3_rhat = NA_real_, Q4_rhat = NA_real_, Q5_rhat = NA_real_,
    Q2_ess_bulk = NA_real_, Q3_ess_bulk = NA_real_, Q4_ess_bulk = NA_real_, Q5_ess_bulk = NA_real_,
    Q2_ess_tail = NA_real_, Q3_ess_tail = NA_real_, Q4_ess_tail = NA_real_, Q5_ess_tail = NA_real_
  )
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

extract_quintiles <- function(draw_df, names_vec, prefix) {
  out <- list(Q1 = "0.00", Q2 = "", Q3 = "", Q4 = "", Q5 = "")
  for (q in 2:5) {
    nm <- paste0(prefix, q)
    idx <- match(nm, names_vec)
    out[[paste0("Q", q)]] <- if (is.na(idx)) "" else format_cell(draw_df[[paste0("beta[", idx, "]")]])
  }
  out
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

# Q1 mu at covariate means: intercept + beta_SM*mean_SM + beta_PA*mean_PA + beta_OB*mean_OB
compute_mu_Q1 <- function(draws, names_vec, d) {
  sm_idx <- match("Smoking_rate", names_vec)
  pa_idx <- match("Physical_Activities_rate", names_vec)
  ob_idx <- match("Obesity_rate", names_vec)
  draws[["beta[1]"]] +
    draws[[paste0("beta[", sm_idx, "]")]] * mean(d$Smoking_rate, na.rm = TRUE) +
    draws[[paste0("beta[", pa_idx, "]")]] * mean(d$Physical_Activities_rate, na.rm = TRUE) +
    draws[[paste0("beta[", ob_idx, "]")]] * mean(d$Obesity_rate, na.rm = TRUE)
}

# Compute MRR draws for Q1-Q5. ref_draws is the denominator (either Lag5 Q1 or this lag's Q1).
# Returns a list of 5 draw vectors. NULL if a quintile's beta is missing.
compute_mrr_draws <- function(draws, names_vec, mu_Q1, ref_draws) {
  lapply(1:5, function(q) {
    if (q == 1) {
      mu_Q1 / ref_draws
    } else {
      idx <- match(paste0("EQI_factor", q), names_vec)
      if (is.na(idx)) {
        return(NULL)
      }
      (mu_Q1 + draws[[paste0("beta[", idx, "]")]]) / ref_draws
    }
  })
}

# Build a wide MRR row — same column layout as MRD but on MRR scale (4 decimal places).
# Includes Q1_p (tests MRR vs 1.0; NA when Q1 is exactly 1.0 by construction).
# rhat/ESS come from the underlying beta draws (same parameters as MRD).
build_mrr_row <- function(mrr_qdraws, q_metrics, outcome, eqi_out, aamr_out, lagv) {
  fmt <- function(d) if (is.null(d)) "" else format_mrr_cell(d)
  p1 <- function(d) if (is.null(d)) NA_character_ else compute_p(d - 1)
  tibble(
    ICD_Code = outcome,
    EQI_Period = eqi_out,
    AAMR_Period = aamr_out,
    Lag = lagv,
    Q1 = fmt(mrr_qdraws[[1]]), Q2 = fmt(mrr_qdraws[[2]]),
    Q3 = fmt(mrr_qdraws[[3]]), Q4 = fmt(mrr_qdraws[[4]]),
    Q5 = fmt(mrr_qdraws[[5]]),
    Q1_p = p1(mrr_qdraws[[1]]), Q2_p = p1(mrr_qdraws[[2]]),
    Q3_p = p1(mrr_qdraws[[3]]), Q4_p = p1(mrr_qdraws[[4]]),
    Q5_p = p1(mrr_qdraws[[5]]),
    Q2_rhat = sprintf("%.4f", q_metrics$Q2_rhat),
    Q3_rhat = sprintf("%.4f", q_metrics$Q3_rhat),
    Q4_rhat = sprintf("%.4f", q_metrics$Q4_rhat),
    Q5_rhat = sprintf("%.4f", q_metrics$Q5_rhat),
    Q2_ess_bulk = as.integer(round(q_metrics$Q2_ess_bulk)),
    Q3_ess_bulk = as.integer(round(q_metrics$Q3_ess_bulk)),
    Q4_ess_bulk = as.integer(round(q_metrics$Q4_ess_bulk)),
    Q5_ess_bulk = as.integer(round(q_metrics$Q5_ess_bulk)),
    Q2_ess_tail = as.integer(round(q_metrics$Q2_ess_tail)),
    Q3_ess_tail = as.integer(round(q_metrics$Q3_ess_tail)),
    Q4_ess_tail = as.integer(round(q_metrics$Q4_ess_tail)),
    Q5_ess_tail = as.integer(round(q_metrics$Q5_ess_tail))
  )
}

# Pairwise lag test on Q5 draws stored per lag. Two-sided Jeffreys p-value.
# store: named list keyed by lag (character), values = Q5 draw vectors.
run_lag_test <- function(store, outcome, out_file) {
  lag_keys <- as.character(sort(as.integer(names(store))))
  if (length(lag_keys) < 2) {
    return(invisible(NULL))
  }
  pairs <- combn(lag_keys, 2, simplify = FALSE)
  rows <- lapply(pairs, function(p) {
    la <- p[1]
    lb <- p[2]
    diff_draws <- store[[la]] - store[[lb]]
    tibble(
      ICD_Code      = outcome,
      comparison    = paste0("Lag", la, "_vs_Lag", lb),
      Q5_diff_mean  = round(mean(diff_draws, na.rm = TRUE), 4),
      Q5_diff_lower = round(quantile(diff_draws, 0.025, na.rm = TRUE), 4),
      Q5_diff_upper = round(quantile(diff_draws, 0.975, na.rm = TRUE), 4),
      p             = compute_p(diff_draws)
    )
  })
  result <- bind_rows(Filter(Negate(is.null), rows))
  if (nrow(result) > 0) {
    append_rows(out_file, result)
    message("[LAG TEST] ", basename(out_file), " — ", nrow(result), " comparisons written")
  }
}

# ── Main loop ──────────────────────────────────────────────────────────────────
for (outcome in selected) {
  dlabel <- disease_label(outcome)
  message("===== Outcome: ", dlabel, " (", outcome, ") =====")

  mrd_file <- file.path(out_dir, paste0(dlabel, "_MRD.csv"))
  sameref_file <- file.path(out_dir, paste0(dlabel, "_MRR_SameRef.csv"))
  lagref_file <- file.path(out_dir, paste0(dlabel, "_MRR_LagRef.csv"))
  lag_mrd_file <- file.path(out_dir, paste0(dlabel, "_lag_test_MRD.csv"))
  lag_sr_file <- file.path(out_dir, paste0(dlabel, "_lag_test_MRR_SameRef.csv"))
  lag_lr_file <- file.path(out_dir, paste0(dlabel, "_lag_test_MRR_LagRef.csv"))

  lag5_ref_draws <- NULL # Lag5 Q1 mu — universal SameRef denominator
  lag_q5_mrd <- list() # Q5 beta draws per lag (for MRD lag test)
  lag_q5_sameref <- list() # Q5 SameRef MRR draws per lag
  lag_q5_lagref <- list() # Q5 LagRef MRR draws per lag

  for (sc in scenario_list) {
    scen_key <- sc$key
    eqi_p <- sc$eqi
    aamr_p <- sc$aamr
    lagv <- sc$lag
    scen_dt <- dt[EQI_Period == eqi_p & Time_Period == aamr_p & Outcome == outcome]
    if (nrow(scen_dt) < opt$`min-n`) {
      message("[Skip] ", scen_key, " n=", nrow(scen_dt))
      next
    }
    eqi_out <- gsub("-", "_", eqi_p)
    aamr_out <- gsub("-", "_", aamr_p)

    des <- build_design_overall(scen_dt)
    if (nrow(des$df) < opt$`min-n`) {
      message("[Skip] ", scen_key, " after design n=", nrow(des$df))
      next
    }
    states <- sort(unique(des$df$State_FIPS))
    si <- match(des$df$State_FIPS, states)
    dl <- list(
      N = nrow(des$df), S = length(states), state = si,
      y_lower = des$df$AAMR_Lower, y_upper = des$df$AAMR_Upper, cens = des$df$cens,
      K = ncol(des$X), X = des$X
    )
    init_fn <- function() list(beta = rep(0, dl$K), z_u = rep(0, dl$S), sigma = 50, sigma_u = 10)

    fit <- try(mod$sample(
      data = dl, chains = opt$chains,
      iter_sampling = opt$iter - opt$warmup, iter_warmup = opt$warmup,
      adapt_delta = opt$`adapt-delta`, max_treedepth = opt$`max-treedepth`,
      parallel_chains = opt$chains, refresh = 0, seed = opt$seed,
      init = rep(list(init_fn()), opt$chains)
    ), silent = TRUE)
    if (inherits(fit, "try-error")) {
      message("[Fail] ", scen_key)
      next
    }

    draws <- as_draws_df(fit$draws("beta"))
    colnames(draws) <- paste0("beta[", seq_len(ncol(draws)), "]")
    summ <- posterior::summarize_draws(fit$draws("beta"))

    mu_Q1_current <- compute_mu_Q1(draws, des$names, des$df)
    if (lagv == 5) lag5_ref_draws <- mu_Q1_current

    # ── MRD output ─────────────────────────────────────────────────────────────
    q_over <- extract_quintiles(draws, des$names, "EQI_factor")
    met <- extract_quintile_metrics(draws, des$names, "EQI_factor", summ)
    sm_o <- extract_covariate(draws, des$names, "Smoking_rate", summ)
    pa_o <- extract_covariate(draws, des$names, "Physical_Activities_rate", summ)
    ob_o <- extract_covariate(draws, des$names, "Obesity_rate", summ)
    row_mrd <- tibble(
      ICD_Code = outcome, EQI_Period = eqi_out, AAMR_Period = aamr_out, Lag = lagv, Model = "EQI+SM+PA+OB",
      Q1 = q_over$Q1, Q2 = q_over$Q2, Q3 = q_over$Q3, Q4 = q_over$Q4, Q5 = q_over$Q5,
      Q2_p = met$Q2_p, Q3_p = met$Q3_p, Q4_p = met$Q4_p, Q5_p = met$Q5_p,
      Q2_rhat = sprintf("%.4f", met$Q2_rhat), Q3_rhat = sprintf("%.4f", met$Q3_rhat),
      Q4_rhat = sprintf("%.4f", met$Q4_rhat), Q5_rhat = sprintf("%.4f", met$Q5_rhat),
      Q2_ess_bulk = as.integer(round(met$Q2_ess_bulk)), Q3_ess_bulk = as.integer(round(met$Q3_ess_bulk)),
      Q4_ess_bulk = as.integer(round(met$Q4_ess_bulk)), Q5_ess_bulk = as.integer(round(met$Q5_ess_bulk)),
      Q2_ess_tail = as.integer(round(met$Q2_ess_tail)), Q3_ess_tail = as.integer(round(met$Q3_ess_tail)),
      Q4_ess_tail = as.integer(round(met$Q4_ess_tail)), Q5_ess_tail = as.integer(round(met$Q5_ess_tail)),
      SM = sm_o$est, SM_p = sm_o$p,
      SM_rhat = sprintf("%.4f", sm_o$rhat),
      SM_ess_bulk = as.integer(round(sm_o$ess_bulk)), SM_ess_tail = as.integer(round(sm_o$ess_tail)),
      PA = pa_o$est, PA_p = pa_o$p,
      PA_rhat = sprintf("%.4f", pa_o$rhat),
      PA_ess_bulk = as.integer(round(pa_o$ess_bulk)), PA_ess_tail = as.integer(round(pa_o$ess_tail)),
      OB = ob_o$est, OB_p = ob_o$p,
      OB_rhat = sprintf("%.4f", ob_o$rhat),
      OB_ess_bulk = as.integer(round(ob_o$ess_bulk)), OB_ess_tail = as.integer(round(ob_o$ess_tail))
    )
    append_rows(mrd_file, row_mrd)
    message("[OK] ", scen_key, " MRD")

    # Q5 beta draws for MRD lag test
    q5_idx <- match("EQI_factor5", des$names)
    if (!is.na(q5_idx)) lag_q5_mrd[[as.character(lagv)]] <- draws[[paste0("beta[", q5_idx, "]")]]

    # ── MRR_SameRef: divide all quintile mu by universal Lag5 Q1 ──────────────
    if (!is.null(lag5_ref_draws)) {
      sr_qdraws <- compute_mrr_draws(draws, des$names, mu_Q1_current, lag5_ref_draws)
      append_rows(sameref_file, build_mrr_row(sr_qdraws, met, outcome, eqi_out, aamr_out, lagv))
      message("[OK] ", scen_key, " MRR_SameRef")
      if (!is.null(sr_qdraws[[5]])) lag_q5_sameref[[as.character(lagv)]] <- sr_qdraws[[5]]
    } else {
      message("[WARN] lag5_ref not yet set — MRR_SameRef skipped for ", scen_key)
    }

    # ── MRR_LagRef: divide by this lag's own Q1. Q1 = mu_Q1/mu_Q1 = 1.0 ─────
    # p-value for Q1 will be NA (all draws exactly 1.0, centered diff = 0).
    lr_qdraws <- compute_mrr_draws(draws, des$names, mu_Q1_current, mu_Q1_current)
    append_rows(lagref_file, build_mrr_row(lr_qdraws, met, outcome, eqi_out, aamr_out, lagv))
    message("[OK] ", scen_key, " MRR_LagRef")
    if (!is.null(lr_qdraws[[5]])) lag_q5_lagref[[as.character(lagv)]] <- lr_qdraws[[5]]
  }

  # ── Pairwise lag tests (all pairs) ────────────────────────────────────────
  run_lag_test(lag_q5_mrd, outcome, lag_mrd_file)
  run_lag_test(lag_q5_sameref, outcome, lag_sr_file)
  run_lag_test(lag_q5_lagref, outcome, lag_lr_file)

  message("===== Completed: ", dlabel, " =====")
}
message("All analyses complete. Output directory: ", out_dir)
