#!/usr/bin/env Rscript
# Interval-censored EQI ~ AAMR brms pipeline (incremental per disease)
# Requirements (conda env 'brms'): r-base, brms, tidyverse, data.table, yaml, optparse
# Tries cmdstanr backend if available; otherwise falls back to rstan.

use_cmdstan <- FALSE  # will set TRUE if cmdstanr + CmdStan available

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(readr)
  library(posterior)
  # Attempt to load cmdstanr before brms so backend can be selected cleanly
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    library(cmdstanr)
    ver <- try(cmdstanr::cmdstan_version(error_on_NA = FALSE), silent = TRUE)
    if (inherits(ver, "try-error") || is.na(ver)) {
      message("CmdStan not found; attempting installation (this may take a few minutes)...")
      try(cmdstanr::install_cmdstan(cores = parallel::detectCores(), quiet = TRUE), silent = TRUE)
      ver <- try(cmdstanr::cmdstan_version(error_on_NA = FALSE), silent = TRUE)
    }
    if (!inherits(ver, "try-error") && !is.na(ver)) {
      use_cmdstan <- TRUE
      options(brms.backend = "cmdstanr")
      message("Using cmdstanr backend (CmdStan v", ver, ")")
    } else {
      message("CmdStan installation attempt did not succeed; will fallback to rstan (may error if toolchain mismatched).")
    }
  } else {
    message("Package cmdstanr not available at script start; install r-cmdstanr for stable interval modeling.")
  }
  library(brms)
})

# --- Command line options ---------------------------------------------------
option_list <- list(
  make_option(c("--data"), type="character", default="Data/Processed/df_EQI_AAMR/EQI_AAMR_Interval.csv", help="Input interval data CSV"),
  make_option(c("--output-dir"), type="character", default="Result/brms", help="Output directory"),
  make_option(c("--cancer-types"), type="character", default=NA, help="Comma separated ICD codes (default: all)"),
  make_option(c("--chains"), type="integer", default=4, help="MCMC chains"),
  make_option(c("--iter"), type="integer", default=2000, help="Total iterations"),
  make_option(c("--warmup"), type="integer", default=1000, help="Warmup iterations"),
  make_option(c("--adapt-delta"), type="double", default=0.95, help="adapt_delta for HMC"),
  make_option(c("--max-treedepth"), type="integer", default=12, help="max_treedepth for HMC"),
  make_option(c("--min-n"), type="integer", default=50, help="Minimum sample size per overall model"),
  make_option(c("--min-n-rucc"), type="integer", default=30, help="Minimum sample size per RUCC layer model"),
  make_option(c("--test"), action="store_true", default=FALSE, help="Test mode: fewer iterations (iter=1000,warmup=500)"),
  make_option(c("--seed"), type="integer", default=1234, help="Random seed")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (opt$test) {
  opt$iter <- min(opt$iter, 1000)
  opt$warmup <- min(opt$warmup, 500)
  message("[TEST MODE] iter=", opt$iter, ", warmup=", opt$warmup)
}

set.seed(opt$seed)

# --- Detect backend cores (80% rule) ---------------------------------------
cores_avail <- parallel::detectCores(logical = TRUE)
cores_used <- max(1, floor(cores_avail * 0.8))
options(mc.cores = cores_used)
message("Detected cores: ", cores_avail, " | Using: ", cores_used)

## (Backend selection already attempted above prior to loading brms.)

# --- Load data --------------------------------------------------------------
project_root <- normalizePath(".")
data_path <- file.path(project_root, opt$data)
if (!file.exists(data_path)) {
  stop("Input data not found: ", data_path)
}
dt <- fread(data_path)

# --- Basic validations ------------------------------------------------------
required_cols <- c("COUNTY_FIPS","EQI_Period","Time_Period","Lag_Years","Cancer_Type","AAMR_lower","AAMR_upper","Smoking_Rate","RUCC",
                   "EQI","EQI_Air","EQI_Water","EQI_Land","EQI_Built","EQI_Social")
missing_req <- setdiff(required_cols, names(dt))
if (length(missing_req)) stop("Missing required columns: ", paste(missing_req, collapse=","))

# Harmonize state identifier
if (!"State_FIPS" %in% names(dt)) {
  if ("State" %in% names(dt)) {
    dt[, State_FIPS := as.character(State)]
  } else {
    dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]
  }
}

# Interval censoring codes (0 = exact, 2 = interval)
dt <- dt[!is.na(AAMR_lower) & !is.na(AAMR_upper)]
dt[, cens := ifelse(AAMR_lower == AAMR_upper, 0, 2)]

# Factor conversions (ordered quintiles 1..5)
quintile_vars <- c("EQI","EQI_Air","EQI_Water","EQI_Land","EQI_Built","EQI_Social")
for (v in quintile_vars) {
  dt[[paste0(v,"_factor")]] <- factor(dt[[v]], levels = 1:5, ordered = TRUE)
}

# RUCC limit to 1..4 (ignore others silently)
dt <- dt[RUCC %in% 1:4 | is.na(RUCC)]

# Scenario mapping -----------------------------------------------------------
scenario_list <- list(
  list(key="EQI0005_AAMR2006_2010", eqi="2000-2005", aamr="2006-2010", lag=5),
  list(key="EQI0005_AAMR2011_2015", eqi="2000-2005", aamr="2011-2015", lag=10),
  list(key="EQI0610_AAMR2011_2015", eqi="2006-2010", aamr="2011-2015", lag=5),
  list(key="EQI0610_AAMR2016_2020", eqi="2006-2010", aamr="2016-2020", lag=10)
)

# Cancer type selection ------------------------------------------------------
all_cancers <- sort(unique(dt$Cancer_Type))
if (is.na(opt$`cancer-types`)) {
  selected_cancers <- all_cancers
} else {
  requested <- str_split(opt$`cancer-types`, ",", simplify = TRUE) |> as.vector() |> str_trim()
  invalid <- setdiff(requested, all_cancers)
  if (length(invalid)) stop("Invalid cancer types: ", paste(invalid, collapse=","))
  selected_cancers <- requested
}
message("Cancer types to analyze: ", paste(selected_cancers, collapse=", "))

# Output setup ---------------------------------------------------------------
out_dir <- file.path(project_root, opt$`output-dir`)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Significance markers -------------------------------------------------------
sig_mark <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("")
}

format_cell <- function(draws) {
  if (is.null(draws) || length(draws) == 0) return("")
  ci <- quantile(draws, c(0.025, 0.975), na.rm=TRUE)
  p <- 2 * min(mean(draws > 0, na.rm=TRUE), mean(draws < 0, na.rm=TRUE))
  sprintf("%0.2f(%0.2f, %0.2f)%s", mean(draws, na.rm=TRUE), ci[1], ci[2], sig_mark(p))
}

# Extract posterior draws for factor levels Q2..Q5 given base variable name pattern
extract_quintiles <- function(post, base_name) {
  out <- list(Q1 = "0.00")
  for (q in 2:5) {
    nm <- paste0("b_", base_name, q)  # e.g. b_EQI_factor2
    if (!(nm %in% names(post))) {
      out[[paste0("Q", q)]] <- ""
    } else {
      out[[paste0("Q", q)]] <- format_cell(post[[nm]])
    }
  }
  out
}

# Append rows to CSV incrementally
append_rows <- function(path, rows_df) {
  if (!file.exists(path)) {
    write_csv(rows_df, path)
  } else {
    suppressWarnings(write.table(rows_df, path, sep=",", col.names=FALSE, row.names=FALSE, append=TRUE))
  }
}

# Fit a single model safely --------------------------------------------------
fit_model <- function(data_model, formula, iter, warmup, chains, adapt_delta, max_treedepth) {
  priors <- c(
    set_prior("student_t(3,0,10)", class="Intercept"),
    set_prior("normal(0,5)", class="b"),
    set_prior("exponential(1)", class="sd"),
    set_prior("exponential(1)", class="sigma")
  )
  ctrl <- list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  # Interval censoring: y | cens(cens) + y2(upper)
  # Build bf object
  bform <- brms::bf(formula)
  backend_sel <- if (use_cmdstan) "cmdstanr" else "rstan"
  result <- try(brm(
    bform,
    data = data_model,
    family = gaussian(),
    prior = priors,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = chains,
    control = ctrl,
    silent = 2,
    refresh = 0,
    backend = backend_sel
  ), silent = TRUE)
  if (inherits(result, "try-error")) {
    msg <- as.character(result)
    message("[Model Compile/Error] ", substr(msg,1,240))
    return(NULL)
  }
  result
}

# Main loop ------------------------------------------------------------------
for (cancer in selected_cancers) {
  message("===== Disease: ", cancer, " =====")
  out_file <- file.path(out_dir, paste0(cancer, "_brms.csv"))
  for (sc in scenario_list) {
    scen_key <- sc$key
    eqi_p <- sc$eqi
    aamr_p <- sc$aamr
    lag_val <- sc$lag
    scen_dt <- dt[EQI_Period == eqi_p & Time_Period == aamr_p & Cancer_Type == cancer]
    if (nrow(scen_dt) < opt$`min-n`) {
      message("[Skip] Scenario ", scen_key, " overall insufficient n=", nrow(scen_dt))
      next
    }
    # Prepare common columns for output rows
    eqi_period_out <- gsub('-', '_', eqi_p)
    aamr_period_out <- gsub('-', '_', aamr_p)

    # Layers: Overall + RUCC 1..4
    layers <- c("Overall", paste0("RUCC", 1:4))
    for (lay in layers) {
      if (lay == "Overall") {
        layer_dt <- scen_dt
        layer_tag <- ""
        min_required <- opt$`min-n`
      } else {
        rucc_val <- as.integer(sub("RUCC", "", lay))
        layer_dt <- scen_dt[RUCC == rucc_val]
        layer_tag <- paste0("RUCC", rucc_val, "_")
        min_required <- opt$`min-n-rucc`
      }
      if (nrow(layer_dt) < min_required) {
        message("[Skip] ", scen_key, " ", lay, " n=", nrow(layer_dt))
        next
      }

      # Build model dataset with lower, upper and censoring code
      model_dt <- layer_dt %>% mutate(
        AAMR_lower = AAMR_lower,
        AAMR_upper = AAMR_upper,
        cens = cens
      )

      # Overall EQI model ----------------------------------------------------
      model_dt$EQI_factor <- factor(model_dt$EQI, levels=1:5, ordered=TRUE)
    overall_formula <- "AAMR_lower | cens(cens, y2 = AAMR_upper) ~ Smoking_Rate + EQI_factor + (1|State_FIPS)"
      fit_overall <- fit_model(model_dt, overall_formula, opt$iter, opt$warmup, opt$chains, opt$`adapt-delta`, opt$`max-treedepth`)
      if (is.null(fit_overall)) {
        message("[Fail] Overall EQI model ", scen_key, " ", lay)
      } else {
  post_overall <- posterior::as_draws_df(fit_overall)
        q_overall <- extract_quintiles(post_overall, "EQI_factor")
        row_overall <- tibble(
          ICD_Code = cancer,
          EQI_Period = eqi_period_out,
          AAMR_Period = aamr_period_out,
          Lag = lag_val,
          Model = paste0(layer_tag, "EQI"),
          Q1 = q_overall$Q1, Q2 = q_overall$Q2, Q3 = q_overall$Q3, Q4 = q_overall$Q4, Q5 = q_overall$Q5
        )
        append_rows(out_file, row_overall)
        message("[OK] ", scen_key, " ", lay, " Model=", layer_tag, "EQI")
      }

      # Multi-domain model ---------------------------------------------------
      for (dom in quintile_vars[-1]) { # ensure factors exist
        model_dt[[paste0(dom, "_factor")]] <- factor(model_dt[[dom]], levels=1:5, ordered=TRUE)
      }
      multi_formula <- paste(
        "AAMR_lower | cens(cens, y2 = AAMR_upper) ~ Smoking_Rate +",
        paste(paste0(quintile_vars[-1], "_factor"), collapse=" + "),
        "+ (1|State_FIPS)"
      )
      fit_multi <- fit_model(model_dt, multi_formula, opt$iter, opt$warmup, opt$chains, opt$`adapt-delta`, opt$`max-treedepth`)
      if (is.null(fit_multi)) {
        message("[Fail] Multi-domain model ", scen_key, " ", lay)
        next
      }
  post_multi <- posterior::as_draws_df(fit_multi)
      # Produce rows for 5 domains (Air, Water, Land, Built, Social)
      domain_labels <- c("EQI_Air","EQI_Water","EQI_Land","EQI_Built","EQI_Social")
      for (dom in domain_labels) {
        base <- paste0(dom, "_factor")
        q_dom <- extract_quintiles(post_multi, base)
        row_dom <- tibble(
          ICD_Code = cancer,
          EQI_Period = eqi_period_out,
          AAMR_Period = aamr_period_out,
          Lag = lag_val,
          Model = paste0(layer_tag, dom),
          Q1 = q_dom$Q1, Q2 = q_dom$Q2, Q3 = q_dom$Q3, Q4 = q_dom$Q4, Q5 = q_dom$Q5
        )
        append_rows(out_file, row_dom)
      }
      message("[OK] ", scen_key, " ", lay, " Multi-domain (5 rows)")
    }
  }
  message("===== Completed: ", cancer, " =====")
}

message("All requested analyses complete. Output directory: ", out_dir)
