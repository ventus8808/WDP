#!/usr/bin/env Rscript
# Sensitivity analysis for the main EQI+SM+PA+OB model.
# Base covariates (SM, PA, OB) are always included.
# Sensitivity covariates (UN, PD, DB, FC) are added:
#   - individually: one at a time (+UN, +PD, +DB, +FC)
#   - simultaneously: all four together (+UN+PD+DB+FC)
# Runs both Overall EQI and multi-domain models for each variant.
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
  "Smoking_rate", "Physical_Activities_rate", "Obesity_rate",
  "Uninsured_rate", "Physician_Density_per100k", "Diabetes_Prevalence_rate", "Forest_Coverage"
))

option_list <- list(
  make_option(c("--data"),          type = "character", default = "Data/Processed/df.csv",      help = "Input data"),
  make_option(c("--output-dir"),    type = "character", default = "Result/brms_Sensitivity",    help = "Output directory"),
  make_option(c("--outcomes"),      type = "character", default = NA,                            help = "Comma separated ICD codes"),
  make_option(c("--chains"),        type = "integer",   default = 4),
  make_option(c("--iter"),          type = "integer",   default = 2000),
  make_option(c("--warmup"),        type = "integer",   default = 1000),
  make_option(c("--adapt-delta"),   type = "double",    default = 0.95),
  make_option(c("--max-treedepth"), type = "integer",   default = 12),
  make_option(c("--min-n"),         type = "integer",   default = 50),
  make_option(c("--seed"),          type = "integer",   default = 1234),
  make_option(c("--test"),          action = "store_true", default = FALSE)
)
opt <- parse_args(OptionParser(option_list = option_list))
if (opt$test) {
  opt$iter   <- min(opt$iter,   800)
  opt$warmup <- min(opt$warmup, 300)
  message("[TEST MODE] iter=", opt$iter, " warmup=", opt$warmup)
}
set.seed(opt$seed)

cores_avail <- parallel::detectCores(logical = TRUE)
cores_used  <- max(1, floor(cores_avail * 0.8))
options(mc.cores = cores_used)
message("Detected cores: ", cores_avail, " | Using: ", cores_used)

stan_code <- "data {\n  int<lower=1> N;\n  int<lower=1> S;\n  array[N] int<lower=1,upper=S> state;\n  vector[N] y_lower;\n  vector[N] y_upper;\n  array[N] int<lower=0,upper=2> cens;\n  int<lower=1> K;\n  matrix[N,K] X;\n} \nparameters {\n  vector[K] beta;\n  vector[S] z_u;\n  real<lower=0> sigma;\n  real<lower=0> sigma_u;\n} \ntransformed parameters {\n  vector[S] u = sigma_u * z_u;\n} \nmodel {\n  beta ~ normal(0,5);\n  z_u ~ normal(0,1);\n  sigma ~ exponential(1);\n  sigma_u ~ exponential(1);\n  for (i in 1:N) {\n    real mu = X[i] * beta + u[state[i]];\n    if (cens[i]==0) {\n      target += normal_lpdf(y_lower[i] | mu, sigma);\n    } else {\n      real p_up = normal_cdf(y_upper[i] | mu, sigma);\n      real p_lo = normal_cdf(y_lower[i] | mu, sigma);\n      real diff = fmax(p_up - p_lo, 1e-12);\n      target += log(diff);\n    }\n  }\n}"
stan_file <- file.path(tempdir(), "interval_mixed_model.stan")
writeLines(stan_code, stan_file)
mod <- cmdstan_model(stan_file)

project_root <- normalizePath(".")
path <- file.path(project_root, opt$data)
if (!file.exists(path)) stop("Data not found: ", path)
dt <- fread(path)

req <- c(
  "COUNTY_FIPS", "EQI_Period", "Time_Period", "Lag_Years", "Outcome",
  "AAMR_Lower", "AAMR_Upper",
  "EQI", "EQI_Air", "EQI_Water", "EQI_Land", "EQI_Built", "EQI_Social",
  "Smoking_rate", "Physical_Activities_rate", "Obesity_rate",
  "Uninsured_rate", "Physician_Density_per100k", "Diabetes_Prevalence_rate", "Forest_Coverage"
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
  "I00_I99"                      = "CVD",
  "J40_J47_J60_J70_J84_D86_C34" = "CRD",
  "K70_K76_C22"                  = "CLD",
  "N00_N29_C64_C65"              = "CKD",
  "X60_X84_Y87.0"                = "Suicide",
  "G20_G30_G12.2_F01_F03"        = "NDD",
  "C00_C97"                      = "Cancer"
)
disease_label <- function(icd) { nm <- icd_to_name[icd]; if (is.na(nm)) icd else unname(nm) }

all_outcomes <- sort(unique(dt$Outcome))
selected <- if (is.na(opt$`outcomes`)) {
  intersect(names(icd_to_name), all_outcomes)
} else {
  reqc <- str_split(opt$`outcomes`, ",", simplify = TRUE) |> as.vector() |> str_trim()
  inv  <- setdiff(reqc, all_outcomes)
  if (length(inv)) stop("Invalid outcomes: ", paste(inv, collapse = ","))
  reqc
}
message("Outcomes to analyze: ", paste(selected, collapse = ","))

out_dir <- file.path(project_root, opt$`output-dir`)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── Sensitivity covariate specs ────────────────────────────────────────────────
# Base covariates (always in model): SM, PA, OB
base_covs <- c("Smoking_rate", "Physical_Activities_rate", "Obesity_rate")

# Sensitivity covariates and their abbreviations
sens_abbrev <- c(
  Uninsured_rate             = "UN",
  Physician_Density_per100k  = "PD",
  Diabetes_Prevalence_rate   = "DB",
  Forest_Coverage            = "FC"
)
sens_vars <- names(sens_abbrev)

# Model variants: individual additions + all simultaneous
sens_specs <- c(
  lapply(sens_vars, function(v) list(extra = v,        label = paste0("+", sens_abbrev[v]))),
  list(list(extra = sens_vars, label = paste0("+", paste(sens_abbrev, collapse = "+"))))
)

# ── Helper functions ───────────────────────────────────────────────────────────
format_cell <- function(draws) {
  if (length(draws) == 0) return("")
  ci <- quantile(draws, c(0.025, 0.975), na.rm = TRUE)
  sprintf("%0.2f(%0.2f,%0.2f)", mean(draws), ci[1], ci[2])
}
append_rows <- function(path, df) {
  if (!file.exists(path)) write_csv(df, path) else suppressWarnings(write.table(df, path, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE))
}

compute_p <- function(draws) {
  if (length(draws) == 0) return(NA_character_)
  pos <- sum(draws > 0, na.rm = TRUE)
  neg <- sum(draws < 0, na.rm = TRUE)
  n   <- pos + neg
  if (n == 0) return(NA_character_)
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
  if (any(grepl(paste0(prefix, "\\.L"), names_vec))) return(out)
  for (q in 2:5) {
    nm  <- paste0(prefix, q)
    idx <- match(nm, names_vec)
    if (!is.na(idx)) {
      col       <- paste0("beta[", idx, "]")
      draws_col <- draw_df[[col]]
      out[[paste0("Q", q, "_p")]] <- compute_p(draws_col)
      sr <- summ_df[summ_df$variable == col, , drop = FALSE]
      if (nrow(sr)) {
        out[[paste0("Q", q, "_rhat")]]     <- sr$rhat
        out[[paste0("Q", q, "_ess_bulk")]] <- sr$ess_bulk
        out[[paste0("Q", q, "_ess_tail")]] <- sr$ess_tail
      }
    }
  }
  out
}

extract_quintiles <- function(draw_df, names_vec, prefix) {
  out <- list(Q1 = "0.00", Q2 = "", Q3 = "", Q4 = "", Q5 = "")
  if (any(grepl(paste0(prefix, "\\.L"), names_vec))) return(out)
  for (q in 2:5) {
    nm  <- paste0(prefix, q)
    idx <- match(nm, names_vec)
    out[[paste0("Q", q)]] <- if (is.na(idx)) "" else format_cell(draw_df[[paste0("beta[", idx, "]")]])
  }
  out
}

extract_covariate <- function(draw_df, names_vec, col_name, summ_df) {
  idx <- match(col_name, names_vec)
  if (is.na(idx)) return(list(est = "", p = NA_character_, rhat = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_))
  col      <- paste0("beta[", idx, "]")
  draws_col <- draw_df[[col]]
  sr       <- summ_df[summ_df$variable == col, , drop = FALSE]
  list(
    est      = format_cell(draws_col),
    p        = compute_p(draws_col),
    rhat     = if (nrow(sr)) sr$rhat      else NA_real_,
    ess_bulk = if (nrow(sr)) sr$ess_bulk  else NA_real_,
    ess_tail = if (nrow(sr)) sr$ess_tail  else NA_real_
  )
}

# Build design matrix: base covariates (SM+PA+OB) + extra sensitivity covariates + EQI_factor
build_design_overall <- function(d, extra_covs = character(0)) {
  d <- d %>% mutate(EQI_factor = factor(EQI, levels = 1:5))
  all_covs <- c(base_covs, extra_covs)
  needed   <- c("EQI_factor", "AAMR_Lower", "AAMR_Upper", "cens", "State_FIPS", all_covs)
  d <- d[complete.cases(d[, ..needed]), ]
  form <- as.formula(paste("~", paste(c(all_covs, "EQI_factor"), collapse = " + ")))
  mm <- model.matrix(form, d, contrasts.arg = list(EQI_factor = contr.treatment(5)))
  colnames(mm) <- make.names(colnames(mm))
  list(X = mm, names = colnames(mm), df = d)
}

build_design_multi <- function(d, extra_covs = character(0)) {
  d <- d %>% mutate(
    EQI_Air_factor    = factor(EQI_Air,    levels = 1:5),
    EQI_Water_factor  = factor(EQI_Water,  levels = 1:5),
    EQI_Land_factor   = factor(EQI_Land,   levels = 1:5),
    EQI_Built_factor  = factor(EQI_Built,  levels = 1:5),
    EQI_Social_factor = factor(EQI_Social, levels = 1:5)
  )
  all_covs <- c(base_covs, extra_covs)
  needed   <- c("EQI_Air_factor", "EQI_Water_factor", "EQI_Land_factor", "EQI_Built_factor", "EQI_Social_factor",
                "AAMR_Lower", "AAMR_Upper", "cens", "State_FIPS", all_covs)
  d <- d[complete.cases(d[, ..needed]), ]
  form <- as.formula(paste("~", paste(c(all_covs,
    "EQI_Air_factor", "EQI_Water_factor", "EQI_Land_factor", "EQI_Built_factor", "EQI_Social_factor"),
    collapse = " + ")))
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

run_stan <- function(des) {
  states <- sort(unique(des$df$State_FIPS))
  si     <- match(des$df$State_FIPS, states)
  dl <- list(
    N = nrow(des$df), S = length(states), state = si,
    y_lower = des$df$AAMR_Lower, y_upper = des$df$AAMR_Upper, cens = des$df$cens,
    K = ncol(des$X), X = des$X
  )
  init_fn <- function() list(beta = rep(0, dl$K), z_u = rep(0, dl$S), sigma = 50, sigma_u = 10)
  fit <- try(mod$sample(
    data = dl, chains = opt$chains, iter_sampling = opt$iter - opt$warmup, iter_warmup = opt$warmup,
    adapt_delta = opt$`adapt-delta`, max_treedepth = opt$`max-treedepth`,
    parallel_chains = min(opt$chains, cores_used), refresh = 0, seed = opt$seed,
    init = rep(list(init_fn()), opt$chains)
  ), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  draws <- as_draws_df(fit$draws("beta"))
  colnames(draws) <- paste0("beta[", seq_len(ncol(draws)), "]")
  list(draws = draws, summ = posterior::summarize_draws(fit$draws("beta")))
}

build_row <- function(outcome, eqi_out, aamr_out, lagv, model_label,
                       res, des, eqi_prefix, extra_covs) {
  q   <- extract_quintiles(res$draws, des$names, eqi_prefix)
  met <- extract_quintile_metrics(res$draws, des$names, eqi_prefix, res$summ)

  all_extra <- c(base_covs, extra_covs)
  abbrevs   <- c(SM = "Smoking_rate", PA = "Physical_Activities_rate", OB = "Obesity_rate",
                 UN = "Uninsured_rate", PD = "Physician_Density_per100k",
                 DB = "Diabetes_Prevalence_rate", FC = "Forest_Coverage")

  cov_cols <- list()
  for (ab in names(abbrevs)) {
    cv <- extract_covariate(res$draws, des$names, abbrevs[[ab]], res$summ)
    cov_cols[[ab]]                <- cv$est
    cov_cols[[paste0(ab, "_p")]]  <- cv$p
    cov_cols[[paste0(ab, "_rhat")]]     <- sprintf("%.4f", cv$rhat)
    cov_cols[[paste0(ab, "_ess_bulk")]] <- as.integer(round(cv$ess_bulk))
    cov_cols[[paste0(ab, "_ess_tail")]] <- as.integer(round(cv$ess_tail))
  }

  bind_cols(
    tibble(
      ICD_Code = outcome, EQI_Period = eqi_out, AAMR_Period = aamr_out, Lag = lagv,
      Model = model_label,
      Q1 = q$Q1, Q2 = q$Q2, Q3 = q$Q3, Q4 = q$Q4, Q5 = q$Q5,
      Q2_p = met$Q2_p, Q3_p = met$Q3_p, Q4_p = met$Q4_p, Q5_p = met$Q5_p,
      Q2_rhat = sprintf("%.4f", met$Q2_rhat), Q3_rhat = sprintf("%.4f", met$Q3_rhat),
      Q4_rhat = sprintf("%.4f", met$Q4_rhat), Q5_rhat = sprintf("%.4f", met$Q5_rhat),
      Q2_ess_bulk = as.integer(round(met$Q2_ess_bulk)), Q3_ess_bulk = as.integer(round(met$Q3_ess_bulk)),
      Q4_ess_bulk = as.integer(round(met$Q4_ess_bulk)), Q5_ess_bulk = as.integer(round(met$Q5_ess_bulk)),
      Q2_ess_tail = as.integer(round(met$Q2_ess_tail)), Q3_ess_tail = as.integer(round(met$Q3_ess_tail)),
      Q4_ess_tail = as.integer(round(met$Q4_ess_tail)), Q5_ess_tail = as.integer(round(met$Q5_ess_tail))
    ),
    as_tibble(cov_cols)
  )
}

# ── Main loop ──────────────────────────────────────────────────────────────────
for (outcome in selected) {
  dlabel  <- disease_label(outcome)
  message("===== Outcome: ", dlabel, " (", outcome, ") =====")
  outfile <- file.path(out_dir, paste0(dlabel, "_Sensitivity.csv"))

  for (sc in scenario_list) {
    scen_key <- sc$key
    eqi_p    <- sc$eqi
    aamr_p   <- sc$aamr
    lagv     <- sc$lag
    scen_dt  <- dt[EQI_Period == eqi_p & Time_Period == aamr_p & Outcome == outcome]
    if (nrow(scen_dt) < opt$`min-n`) {
      message("[Skip] ", scen_key, " n=", nrow(scen_dt))
      next
    }
    eqi_out  <- gsub("-", "_", eqi_p)
    aamr_out <- gsub("-", "_", aamr_p)

    for (spec in sens_specs) {
      extra <- spec$extra
      slabel <- spec$label

      # ── Overall EQI + SM+PA+OB + extra ──────────────────────────────────────
      des_o <- build_design_overall(scen_dt, extra)
      if (nrow(des_o$df) < opt$`min-n`) {
        message("[Skip] ", scen_key, " EQI", slabel, " n=", nrow(des_o$df))
      } else {
        res_o <- run_stan(des_o)
        if (is.null(res_o)) {
          message("[Fail] ", scen_key, " EQI", slabel)
        } else {
          model_label <- paste0("EQI+SM+PA+OB", slabel)
          row_o <- build_row(outcome, eqi_out, aamr_out, lagv, model_label,
                             res_o, des_o, "EQI_factor", extra)
          append_rows(outfile, row_o)
          message("[OK] ", scen_key, " ", model_label)
        }
      }

      # ── Multi-domain + SM+PA+OB + extra ─────────────────────────────────────
      des_m <- build_design_multi(scen_dt, extra)
      if (nrow(des_m$df) < opt$`min-n`) {
        message("[Skip] ", scen_key, " Multi", slabel, " n=", nrow(des_m$df))
      } else {
        res_m <- run_stan(des_m)
        if (is.null(res_m)) {
          message("[Fail] ", scen_key, " Multi", slabel)
        } else {
          domain_prefix <- c("EQI_Air_factor", "EQI_Water_factor", "EQI_Land_factor",
                             "EQI_Built_factor", "EQI_Social_factor")
          for (dom in domain_prefix) {
            dom_name    <- sub("_factor", "", dom)
            model_label <- paste0(dom_name, "+SM+PA+OB", slabel)
            row_d <- build_row(outcome, eqi_out, aamr_out, lagv, model_label,
                               res_m, des_m, dom, extra)
            append_rows(outfile, row_d)
          }
          message("[OK] ", scen_key, " Multi+SM+PA+OB", slabel)
        }
      }
    }
  }
  message("===== Completed: ", dlabel, " =====")
}
message("All requested analyses complete. Output directory: ", out_dir)
