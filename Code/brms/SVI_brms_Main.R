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
#   - DISEASE-SPECIFIC adjustment: each disease category uses its own covariate
#     set (see DISEASE_COVSET below). The outcome's category is resolved from
#     config.yaml (diseases -> overall/subtypes icd_code), so subtypes inherit
#     their family's covariate set automatically.
#         cancer / liver(CLD) / cvd / kidney(CKD) : SM + UN
#         ndd                       : PA + UN
#         suicide                   : UN + PD + DB
#         respiratory(CRD)          : UN
#     Cancer takes precedence: any (sub)type that is a cancer (belongs to the
#     cancer family in config) uses the cancer set (SM + UN), even when it is
#     also listed under an organ system (e.g. C34 lung cancer, C64_C65 RCC).
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
  "Smoking_rate", "Physical_Activities_rate", "Obesity_rate",
  "Uninsured_rate", "Physician_Density_per100k", "Diabetes_Prevalence_rate"
))

# ── Disease-specific covariate adjustment sets ──────────────────────────────────
# Covariate column names keyed by config.yaml disease-category key.
# respiratory (CRD) is intentionally absent: no set provided -> skipped + warned.
DISEASE_COVSET <- list(
  cancer      = c("Smoking_rate", "Uninsured_rate"),
  liver       = c("Smoking_rate", "Uninsured_rate"),
  cvd         = c("Smoking_rate", "Uninsured_rate"),
  ndd         = c("Physical_Activities_rate", "Uninsured_rate"),
  suicide     = c("Uninsured_rate", "Physician_Density_per100k", "Diabetes_Prevalence_rate"),
  kidney      = c("Smoking_rate", "Uninsured_rate"),
  respiratory = c("Uninsured_rate")
)
covar_abbrev <- c(
  Smoking_rate = "SM", Physical_Activities_rate = "PA", Obesity_rate = "OB",
  Uninsured_rate = "UN", Physician_Density_per100k = "PD", Diabetes_Prevalence_rate = "DB"
)
ALL_COVS <- names(covar_abbrev)   # fixed output column order

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
  "COUNTY_FIPS", "Time_Period", "Outcome", "AAMR_Lower", "AAMR_Upper", "SVI",
  ALL_COVS
)
miss <- setdiff(req, names(dt))
if (length(miss)) stop("Missing cols: ", paste(miss, collapse = ","))

if (!"State_FIPS" %in% names(dt)) dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]

# interval censoring code
dt <- dt[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
dt[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]

# ── ICD code -> disease category map ────────────────────────────────────────────
# Mirrors config.yaml `diseases` (overall + subtypes icd_code) but hardcoded so the
# pipeline needs no extra R package (the `yaml` package is absent on the cluster).
# If config.yaml disease definitions change, update this list to match.
DISEASE_CODES <- list(
  liver       = c("K70_K76_C22", "K70_K76", "K70", "K71", "K73", "K74",
                  "K71_K73_K74", "K76", "K76.7", "C22"),
  respiratory = c("J40_J47_J60_J70_J84_D86_C34", "J40_J47_J60_J70_J84_D86",
                  "J43_J44", "J45", "J84_D86", "J60_J66", "C34"),
  kidney      = c("N00_N29_C64_C65", "N00_N29", "N18_N19", "N00_N15", "C64_C65"),
  cvd         = c("I00_I99", "I20_I25", "I60_I69", "I10_I15", "I50"),
  suicide     = c("X60_X84_Y87.0", "X60_X69", "X70_X84", "Y87.0"),
  ndd         = c("G20_G30_G12.2_F01_F03", "G30_F01_F03", "G20", "G10", "G12.2"),
  cancer      = c("C00_C97", "C18_C21", "C22", "C25", "C34", "C50", "C56",
                  "C61", "C64_C65", "C82_C85", "C91_C95")
)
icd_to_cats <- list()   # icd_code -> character vector of category keys
for (cat in names(DISEASE_CODES)) {
  for (code in DISEASE_CODES[[cat]]) {
    icd_to_cats[[code]] <- unique(c(icd_to_cats[[code]], cat))
  }
}

# Resolve the covariate set for an outcome (ICD code), or report why it can't be.
resolve_covset <- function(outcome) {
  cats <- icd_to_cats[[outcome]]
  if (is.null(cats)) {
    return(list(ok = FALSE, reason = paste0("no category in config for '", outcome, "'")))
  }
  # Cancer takes precedence: any cancer (sub)type uses the cancer covariate set.
  if ("cancer" %in% cats && "cancer" %in% names(DISEASE_COVSET)) {
    return(list(ok = TRUE, covset = DISEASE_COVSET[["cancer"]], cat = "cancer"))
  }
  defined <- cats[cats %in% names(DISEASE_COVSET)]
  if (length(defined) == 0) {
    return(list(ok = FALSE, reason = paste0("no covariate set for category(ies): ",
                                            paste(cats, collapse = "/"))))
  }
  uniq <- unique(lapply(defined, function(k) sort(DISEASE_COVSET[[k]])))
  if (length(uniq) > 1) {
    return(list(ok = FALSE, reason = paste0("ambiguous covariate sets across categories: ",
                                            paste(defined, collapse = "/"))))
  }
  list(ok = TRUE, covset = DISEASE_COVSET[[defined[1]]], cat = defined[1])
}

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

# Per-category (B/C/D) p / R-hat / ESS, flat & formatted (ready to splice into a row).
svi_metrics_block <- function(draw_df, names_vec, summ_df) {
  out <- list()
  for (lv in c("B", "C", "D")) {
    p <- NA_character_; rhat <- NA_real_; eb <- NA_real_; et <- NA_real_
    nm <- paste0("SVI_factor", lv)
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

# All 6 covariate blocks (only the adjusted ones are populated; others blank).
covar_block <- function(draw_df, names_vec, summ_df) {
  out <- list()
  for (cn in ALL_COVS) {
    ab <- unname(covar_abbrev[cn])
    cc <- extract_covariate(draw_df, names_vec, cn, summ_df)
    out[[ab]]                      <- cc$est
    out[[paste0(ab, "_p")]]        <- cc$p
    out[[paste0(ab, "_rhat")]]     <- sprintf("%.4f", cc$rhat)
    out[[paste0(ab, "_ess_bulk")]] <- as.integer(round(cc$ess_bulk))
    out[[paste0(ab, "_ess_tail")]] <- as.integer(round(cc$ess_tail))
  }
  out
}

# Build design: intercept + disease-specific covariates + SVI_factor (treatment, A ref).
build_design <- function(d, covariates) {
  d <- as.data.frame(d)
  d$SVI_factor <- factor(d$SVI, levels = c("A", "B", "C", "D"))
  needed <- c("SVI_factor", "AAMR_Lower", "AAMR_Upper", "cens", "State_FIPS", covariates)
  d <- d[stats::complete.cases(d[, needed, drop = FALSE]), , drop = FALSE]
  form <- if (length(covariates) == 0) {
    ~SVI_factor
  } else {
    as.formula(paste("~", paste(c(covariates, "SVI_factor"), collapse = " + ")))
  }
  mm <- model.matrix(form, d, contrasts.arg = list(SVI_factor = contr.treatment(c("A", "B", "C", "D"))))
  colnames(mm) <- make.names(colnames(mm))
  list(X = mm, names = colnames(mm), df = d)
}

for (outcome in selected) {
  rc <- resolve_covset(outcome)
  if (!isTRUE(rc$ok)) {
    message("[Skip] Outcome ", outcome, " -> ", rc$reason)
    next
  }
  covset <- rc$covset
  model_label <- paste(c("SVI", covar_abbrev[covset]), collapse = "+")
  message("===== Outcome: ", outcome, "  [", rc$cat, "]  Model: ", model_label, " =====")
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

    des <- build_design(scen_dt, covset)
    if (nrow(des$df) < opt$`min-n`) {
      message("[Skip] ", scen_key, " ", model_label, " n=", nrow(des$df))
      next
    }
    states_o <- sort(unique(des$df$State_FIPS))
    state_index_o <- match(des$df$State_FIPS, states_o)
    data_list <- list(
      N = nrow(des$df), S = length(states_o), state = state_index_o,
      y_lower = des$df$AAMR_Lower, y_upper = des$df$AAMR_Upper, cens = des$df$cens,
      K = ncol(des$X), X = des$X
    )
    init_fun <- function() list(beta = rep(0, data_list$K), z_u = rep(0, data_list$S), sigma = 50, sigma_u = 10)
    fit_overall <- try(mod$sample(
      data = data_list, chains = opt$chains, iter_sampling = opt$iter - opt$warmup, iter_warmup = opt$warmup,
      adapt_delta = opt$`adapt-delta`, max_treedepth = opt$`max-treedepth`, parallel_chains = opt$chains, refresh = 0, seed = opt$seed,
      init = rep(list(init_fun()), opt$chains)
    ), silent = TRUE)
    if (inherits(fit_overall, "try-error")) {
      message("[Fail] ", model_label, " model ", scen_key)
      next
    }
    draws <- as_draws_df(fit_overall$draws("beta"))
    colnames(draws) <- paste0("beta[", seq_len(ncol(draws)), "]")
    summ_over <- posterior::summarize_draws(fit_overall$draws("beta"))

    row_list <- c(
      list(ICD_Code = outcome, AAMR_Period = aamr_out, Lag = lagv, Model = model_label),
      extract_svi(draws, des$names),
      svi_metrics_block(draws, des$names, summ_over),
      covar_block(draws, des$names, summ_over)
    )
    row <- do.call(tibble, row_list)
    append_rows(outfile, row)
    message("[OK] ", scen_key, " ", model_label)
  }
  message("===== Completed: ", outcome, " =====")
}
message("All requested analyses complete. Output directory: ", out_dir)
