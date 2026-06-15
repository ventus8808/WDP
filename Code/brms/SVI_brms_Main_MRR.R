#!/usr/bin/env Rscript
# SVI MRR pipeline — interval-censored mixed model, SVI trajectory exposure.
#   - Exposure: SVI category A/B/C/D (A = reference).
#   - Disease-specific covariate adjustment (see DISEASE_COVSET); the outcome's
#     category is resolved from a hardcoded ICD->category map mirroring config.yaml
#     (no yaml package needed). Cancer takes precedence for cancer (sub)types.
#   - SVI is static -> four AAMR periods treated as four lags (5/10/15/20 yr).
#
# Reference designs (denominator = adjusted category-A mean rate, mu_A):
#   SameRef — universal Lag5 mu_A denominator.
#              A@Lag5 = 1.0 by construction; A@other lags captures secular trend.
#   LagRef  — per-lag mu_A denominator (within-lag relative comparison); A = 1.0.
#
# Output files per disease:
#   {dlabel}_MRD.csv                   — wide, MRD (deaths/100k) vs A, one row per lag
#   {dlabel}_MRR_SameRef.csv           — wide, MRR vs universal Lag5 A, one row per lag
#   {dlabel}_MRR_LagRef.csv            — wide, MRR vs per-lag A, one row per lag
#   {dlabel}_lag_test_MRD.csv          — pairwise category-D lag comparison (MRD)
#   {dlabel}_lag_test_MRR_SameRef.csv  — pairwise category-D lag comparison (SameRef MRR)
#   {dlabel}_lag_test_MRR_LagRef.csv   — pairwise category-D lag comparison (LagRef MRR)

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

# ── Disease-specific covariate adjustment sets (config.yaml category keys) ───────
DISEASE_COVSET <- list(
  cancer      = c("Smoking_rate", "Uninsured_rate"),
  liver       = c("Smoking_rate", "Uninsured_rate"),
  cvd         = c("Smoking_rate", "Uninsured_rate"),
  kidney      = c("Smoking_rate", "Uninsured_rate"),
  ndd         = c("Physical_Activities_rate", "Uninsured_rate"),
  suicide     = c("Uninsured_rate", "Physician_Density_per100k", "Diabetes_Prevalence_rate"),
  respiratory = c("Uninsured_rate")
)
covar_abbrev <- c(
  Smoking_rate = "SM", Physical_Activities_rate = "PA", Obesity_rate = "OB",
  Uninsured_rate = "UN", Physician_Density_per100k = "PD", Diabetes_Prevalence_rate = "DB"
)
ALL_COVS <- names(covar_abbrev)

# ICD code -> category (mirrors config.yaml diseases; hardcoded, no yaml dependency)
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
icd_to_cats <- list()
for (cat in names(DISEASE_CODES)) {
  for (code in DISEASE_CODES[[cat]]) {
    icd_to_cats[[code]] <- unique(c(icd_to_cats[[code]], cat))
  }
}
resolve_covset <- function(outcome) {
  cats <- icd_to_cats[[outcome]]
  if (is.null(cats)) {
    return(list(ok = FALSE, reason = paste0("no category for '", outcome, "'")))
  }
  if ("cancer" %in% cats && "cancer" %in% names(DISEASE_COVSET)) {
    return(list(ok = TRUE, covset = DISEASE_COVSET[["cancer"]], cat = "cancer"))
  }
  defined <- cats[cats %in% names(DISEASE_COVSET)]
  if (length(defined) == 0) {
    return(list(ok = FALSE, reason = paste0("no covariate set for: ", paste(cats, collapse = "/"))))
  }
  uniq <- unique(lapply(defined, function(k) sort(DISEASE_COVSET[[k]])))
  if (length(uniq) > 1) {
    return(list(ok = FALSE, reason = paste0("ambiguous covariate sets: ", paste(defined, collapse = "/"))))
  }
  list(ok = TRUE, covset = DISEASE_COVSET[[defined[1]]], cat = defined[1])
}

option_list <- list(
  make_option(c("--data"), type = "character", default = "Data/Processed/df_SVI.csv"),
  make_option(c("--output-dir"), type = "character", default = "Result/brms_SVI_Main_MRR"),
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
stan_file <- file.path(tempdir(), "interval_mixed_svi_mrr.stan")
writeLines(stan_code, stan_file)
mod <- cmdstan_model(stan_file)

project_root <- normalizePath(".")
path <- file.path(project_root, opt$data)
if (!file.exists(path)) stop("Data not found: ", path)
dt <- fread(path)

req <- c(
  "COUNTY_FIPS", "Time_Period", "Outcome", "AAMR_Lower", "AAMR_Upper", "SVI", ALL_COVS
)
miss <- setdiff(req, names(dt))
if (length(miss)) stop("Missing cols: ", paste(miss, collapse = ","))

if (!"State_FIPS" %in% names(dt)) dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]
dt <- dt[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
dt[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]

# SVI static -> 4 lag scenarios (filter by AAMR Time_Period only)
scenario_list <- list(
  list(key = "SVI_AAMR2006_2010", aamr = "2006-2010", lag = 5),
  list(key = "SVI_AAMR2011_2015", aamr = "2011-2015", lag = 10),
  list(key = "SVI_AAMR2016_2020", aamr = "2016-2020", lag = 15),
  list(key = "SVI_AAMR2021_2024", aamr = "2021-2024", lag = 20)
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
  if (length(draws) == 0) return("")
  ci <- quantile(draws, c(0.025, 0.975), na.rm = TRUE)
  sprintf("%0.2f(%0.2f,%0.2f)", mean(draws), ci[1], ci[2])
}
format_mrr_cell <- function(draws) {
  if (length(draws) == 0) return("")
  ci <- quantile(draws, c(0.025, 0.975), na.rm = TRUE)
  sprintf("%0.4f(%0.4f,%0.4f)", mean(draws), ci[1], ci[2])
}
append_rows <- function(path, df) {
  if (!file.exists(path)) write_csv(df, path) else suppressWarnings(write.table(df, path, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE))
}
compute_p <- function(draws) {
  if (length(draws) == 0) return(NA_character_)
  pos <- sum(draws > 0, na.rm = TRUE)
  neg <- sum(draws < 0, na.rm = TRUE)
  n <- pos + neg
  if (n == 0) return(NA_character_)
  sprintf("%.4f", 2 * min((pos + 0.5) / (n + 1), (neg + 0.5) / (n + 1)))
}

# Per-category (B/C/D vs A) p / R-hat / ESS (formatted, flat).
svi_metrics_block <- function(draw_df, names_vec, summ_df) {
  out <- list()
  for (lv in c("B", "C", "D")) {
    p <- NA_character_; rhat <- NA_real_; eb <- NA_real_; et <- NA_real_
    idx <- match(paste0("SVI_factor", lv), names_vec)
    if (!is.na(idx)) {
      col <- paste0("beta[", idx, "]")
      p <- compute_p(draw_df[[col]])
      sr <- summ_df[summ_df$variable == col, , drop = FALSE]
      if (nrow(sr)) { rhat <- sr$rhat; eb <- sr$ess_bulk; et <- sr$ess_tail }
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
  sr <- summ_df[summ_df$variable == col, , drop = FALSE]
  list(
    est = format_cell(draw_df[[col]]), p = compute_p(draw_df[[col]]),
    rhat = if (nrow(sr)) sr$rhat else NA_real_,
    ess_bulk = if (nrow(sr)) sr$ess_bulk else NA_real_,
    ess_tail = if (nrow(sr)) sr$ess_tail else NA_real_
  )
}
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
extract_svi <- function(draw_df, names_vec) {
  out <- list(A = "0.00", B = "", C = "", D = "")
  for (lv in c("B", "C", "D")) {
    idx <- match(paste0("SVI_factor", lv), names_vec)
    out[[lv]] <- if (is.na(idx)) "" else format_cell(draw_df[[paste0("beta[", idx, "]")]])
  }
  out
}

# Build design: intercept + disease covariates + SVI_factor (treatment, A ref).
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

# Adjusted category-A mean rate: intercept + sum(beta_cov * mean_cov).
compute_mu_A <- function(draws, names_vec, d, covariates) {
  mu <- draws[["beta[1]"]]
  for (cn in covariates) {
    idx <- match(cn, names_vec)
    if (!is.na(idx)) mu <- mu + draws[[paste0("beta[", idx, "]")]] * mean(d[[cn]], na.rm = TRUE)
  }
  mu
}

# MRR draws for A/B/C/D. ref_draws = denominator (Lag5 A mu, or this lag's A mu).
compute_mrr_draws <- function(draws, names_vec, mu_A, ref_draws) {
  out <- lapply(c("A", "B", "C", "D"), function(lv) {
    if (lv == "A") {
      mu_A / ref_draws
    } else {
      idx <- match(paste0("SVI_factor", lv), names_vec)
      if (is.na(idx)) NULL else (mu_A + draws[[paste0("beta[", idx, "]")]]) / ref_draws
    }
  })
  names(out) <- c("A", "B", "C", "D")
  out
}

# Wide MRR row (4-decimal MRR). *_p tests MRR vs 1.0 (NA when level is exactly 1.0).
# rhat/ESS come from the underlying beta draws (B/C/D), passed via `met`.
build_mrr_row <- function(mrr, met, outcome, aamr_out, lagv) {
  fmt <- function(d) if (is.null(d)) "" else format_mrr_cell(d)
  p1 <- function(d) if (is.null(d)) NA_character_ else compute_p(d - 1)
  tibble(
    ICD_Code = outcome, AAMR_Period = aamr_out, Lag = lagv,
    A = fmt(mrr$A), B = fmt(mrr$B), C = fmt(mrr$C), D = fmt(mrr$D),
    A_p = p1(mrr$A), B_p = p1(mrr$B), C_p = p1(mrr$C), D_p = p1(mrr$D),
    B_rhat = met$B_rhat, C_rhat = met$C_rhat, D_rhat = met$D_rhat,
    B_ess_bulk = met$B_ess_bulk, C_ess_bulk = met$C_ess_bulk, D_ess_bulk = met$D_ess_bulk,
    B_ess_tail = met$B_ess_tail, C_ess_tail = met$C_ess_tail, D_ess_tail = met$D_ess_tail
  )
}

# Pairwise lag test on stored category-D draws. store: lag(char) -> D draw vector.
run_lag_test <- function(store, outcome, out_file) {
  lag_keys <- as.character(sort(as.integer(names(store))))
  if (length(lag_keys) < 2) return(invisible(NULL))
  pairs <- combn(lag_keys, 2, simplify = FALSE)
  rows <- lapply(pairs, function(p) {
    diff_draws <- store[[p[1]]] - store[[p[2]]]
    tibble(
      ICD_Code     = outcome,
      comparison   = paste0("Lag", p[1], "_vs_Lag", p[2]),
      D_diff_mean  = round(mean(diff_draws, na.rm = TRUE), 4),
      D_diff_lower = round(quantile(diff_draws, 0.025, na.rm = TRUE), 4),
      D_diff_upper = round(quantile(diff_draws, 0.975, na.rm = TRUE), 4),
      p            = compute_p(diff_draws)
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
  rc <- resolve_covset(outcome)
  if (!isTRUE(rc$ok)) {
    message("[Skip] Outcome ", outcome, " -> ", rc$reason)
    next
  }
  covset <- rc$covset
  model_label <- paste(c("SVI", covar_abbrev[covset]), collapse = "+")
  dlabel <- disease_label(outcome)
  message("===== Outcome: ", dlabel, " (", outcome, ")  [", rc$cat, "]  Model: ", model_label, " =====")

  mrd_file     <- file.path(out_dir, paste0(dlabel, "_MRD.csv"))
  sameref_file <- file.path(out_dir, paste0(dlabel, "_MRR_SameRef.csv"))
  lagref_file  <- file.path(out_dir, paste0(dlabel, "_MRR_LagRef.csv"))
  lag_mrd_file <- file.path(out_dir, paste0(dlabel, "_lag_test_MRD.csv"))
  lag_sr_file  <- file.path(out_dir, paste0(dlabel, "_lag_test_MRR_SameRef.csv"))
  lag_lr_file  <- file.path(out_dir, paste0(dlabel, "_lag_test_MRR_LagRef.csv"))

  lag5_ref_draws <- NULL          # Lag5 mu_A — universal SameRef denominator
  lag_d_mrd <- list()             # D beta draws per lag (MRD lag test)
  lag_d_sameref <- list()         # D SameRef MRR draws per lag
  lag_d_lagref <- list()          # D LagRef MRR draws per lag

  for (sc in scenario_list) {
    scen_key <- sc$key
    aamr_p <- sc$aamr
    lagv <- sc$lag
    scen_dt <- dt[Time_Period == aamr_p & Outcome == outcome]
    if (nrow(scen_dt) < opt$`min-n`) {
      message("[Skip] ", scen_key, " n=", nrow(scen_dt))
      next
    }
    aamr_out <- gsub("-", "_", aamr_p)

    des <- build_design(scen_dt, covset)
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

    mu_A_current <- compute_mu_A(draws, des$names, des$df, covset)
    if (lagv == 5) lag5_ref_draws <- mu_A_current

    met <- svi_metrics_block(draws, des$names, summ)

    # ── MRD output ────────────────────────────────────────────────────────────
    row_mrd <- do.call(tibble, c(
      list(ICD_Code = outcome, AAMR_Period = aamr_out, Lag = lagv, Model = model_label),
      extract_svi(draws, des$names),
      met,
      covar_block(draws, des$names, summ)
    ))
    append_rows(mrd_file, row_mrd)
    message("[OK] ", scen_key, " MRD")

    d_idx <- match("SVI_factorD", des$names)
    if (!is.na(d_idx)) lag_d_mrd[[as.character(lagv)]] <- draws[[paste0("beta[", d_idx, "]")]]

    # ── MRR_SameRef: divide all category mu by universal Lag5 mu_A ─────────────
    if (!is.null(lag5_ref_draws)) {
      sr <- compute_mrr_draws(draws, des$names, mu_A_current, lag5_ref_draws)
      append_rows(sameref_file, build_mrr_row(sr, met, outcome, aamr_out, lagv))
      message("[OK] ", scen_key, " MRR_SameRef")
      if (!is.null(sr$D)) lag_d_sameref[[as.character(lagv)]] <- sr$D
    } else {
      message("[WARN] lag5_ref not yet set — MRR_SameRef skipped for ", scen_key)
    }

    # ── MRR_LagRef: divide by this lag's own mu_A (A = 1.0) ────────────────────
    lr <- compute_mrr_draws(draws, des$names, mu_A_current, mu_A_current)
    append_rows(lagref_file, build_mrr_row(lr, met, outcome, aamr_out, lagv))
    message("[OK] ", scen_key, " MRR_LagRef")
    if (!is.null(lr$D)) lag_d_lagref[[as.character(lagv)]] <- lr$D
  }

  run_lag_test(lag_d_mrd, outcome, lag_mrd_file)
  run_lag_test(lag_d_sameref, outcome, lag_sr_file)
  run_lag_test(lag_d_lagref, outcome, lag_lr_file)

  message("===== Completed: ", dlabel, " =====")
}
message("All analyses complete. Output directory: ", out_dir)
