#!/usr/bin/env Rscript
# Combined MRR pipeline: Stratified (sex/race) + Typology + LandUse
# For a given base ICD code, computes MRR (Q1-Q5 vs Lag5 Q1) for:
#   - Overall EQI and all 5 EQI domains (Air/Water/Land/Built/Social)
#   - Sex/race strata (from EQI_AAMR_Stratifed.csv)
#   - Typology strata (from EQI_AAMR_Cluster_Climate_Typology_LandUse.csv)
#
# Output format mirrors MRD: ICD_Code carries the stratum for sex/race
# (e.g. NDD_Asian), Model carries the domain (Stratified_EQI, Stratified_Air …)
#
# Outputs (all in --output-dir):
#   {cancer}_Stratified_MRR.csv — sex/race strata, all 6 models
#   {cancer}_Typology_MRR.csv   — Typology strata, all 6 models
#   {cancer}_lag_test.csv       — EQI0005 Q5 pairwise lag comparison

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
utils::globalVariables(c('EQI','EQI_Air','EQI_Water','EQI_Land','EQI_Built','EQI_Social','Smoking_Rate','State_FIPS'))

option_list <- list(
  make_option(c("--cancer-type"),    type="character", default=NA,
              help="Base ICD code (e.g. C00_C97, G20). Required."),
  make_option(c("--data-typo"),      type="character",
              default="Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Cluster_Climate_Typology_LandUse.csv",
              help="Typology/LandUse data file"),
  make_option(c("--data-strat"),     type="character",
              default="Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Stratifed.csv",
              help="Sex/race stratified data file"),
  make_option(c("--output-dir"),     type="character", default="Result/brms_Stratified_Typo_LandUse_MRR"),
  make_option(c("--chains"),         type="integer",   default=4),
  make_option(c("--iter"),           type="integer",   default=2000),
  make_option(c("--warmup"),         type="integer",   default=1000),
  make_option(c("--adapt-delta"),    type="double",    default=0.95),
  make_option(c("--max-treedepth"), type="integer",   default=12),
  make_option(c("--min-n"),          type="integer",   default=50),
  make_option(c("--seed"),           type="integer",   default=1234),
  make_option(c("--test"),           action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.na(opt$`cancer-type`)) stop("--cancer-type is required")
if (opt$test) {
  opt$iter   <- min(opt$iter,   800)
  opt$warmup <- min(opt$warmup, 300)
  message("[TEST MODE] iter=", opt$iter, " warmup=", opt$warmup)
}
set.seed(opt$seed)

cores_avail <- parallel::detectCores(logical=TRUE)
cores_used  <- max(1, floor(cores_avail * 0.8))
options(mc.cores=cores_used)
message("Detected cores: ", cores_avail, " | Using: ", cores_used)

# ---------------------------------------------------------------------------
# Stan model (interval-censored mixed model)
# ---------------------------------------------------------------------------
stan_code <- "data {\n  int<lower=1> N;\n  int<lower=1> S;\n  array[N] int<lower=1,upper=S> state;\n  vector[N] y_lower;\n  vector[N] y_upper;\n  array[N] int<lower=0,upper=2> cens;\n  int<lower=1> K;\n  matrix[N,K] X;\n} \nparameters {\n  vector[K] beta;\n  vector[S] z_u;\n  real<lower=0> sigma;\n  real<lower=0> sigma_u;\n} \ntransformed parameters {\n  vector[S] u = sigma_u * z_u;\n} \nmodel {\n  beta ~ normal(0,5);\n  z_u ~ normal(0,1);\n  sigma ~ exponential(1);\n  sigma_u ~ exponential(1);\n  for (i in 1:N) {\n    real mu = X[i] * beta + u[state[i]];\n    if (cens[i]==0) {\n      target += normal_lpdf(y_lower[i] | mu, sigma);\n    } else {\n      real p_up = normal_cdf(y_upper[i] | mu, sigma);\n      real p_lo = normal_cdf(y_lower[i] | mu, sigma);\n      real diff = fmax(p_up - p_lo, 1e-12);\n      target += log(diff);\n    }\n  }\n}"
stan_file <- file.path(tempdir(), "interval_mixed_model.stan")
writeLines(stan_code, stan_file)
mod <- cmdstan_model(stan_file)

project_root <- normalizePath(".")
base_cancer  <- opt$`cancer-type`
out_dir      <- file.path(project_root, opt$`output-dir`)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

# ---------------------------------------------------------------------------
# Scenarios
# ---------------------------------------------------------------------------
scenario_list <- list(
  list(key="EQI0005_AAMR2006_2010", eqi="2000-2005", aamr="2006-2010", lag=5),
  list(key="EQI0005_AAMR2011_2015", eqi="2000-2005", aamr="2011-2015", lag=10),
  list(key="EQI0005_AAMR2016_2020", eqi="2000-2005", aamr="2016-2020", lag=15),
  list(key="EQI0610_AAMR2011_2015", eqi="2006-2010", aamr="2011-2015", lag=5),
  list(key="EQI0610_AAMR2016_2020", eqi="2006-2010", aamr="2016-2020", lag=10)
)

# Domain models: (model suffix, EQI factor prefix in design matrix)
domain_models <- list(
  list(suffix="Air",    prefix="EQI_Air_factor"),
  list(suffix="Water",  prefix="EQI_Water_factor"),
  list(suffix="Land",   prefix="EQI_Land_factor"),
  list(suffix="Built",  prefix="EQI_Built_factor"),
  list(suffix="Social", prefix="EQI_Social_factor")
)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
format_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.0001) return("p<0.0001")
  sprintf("%.4f", p)
}

append_rows <- function(path, df) {
  if (!file.exists(path)) write_csv(df, path) else
    suppressWarnings(write.table(df, path, sep=",", col.names=FALSE, row.names=FALSE, append=TRUE))
}

build_design_overall <- function(d) {
  d <- d %>% mutate(EQI_factor = factor(EQI, levels=1:5))
  d <- d[complete.cases(d[, c("Smoking_Rate","EQI_factor","AAMR_Lower","AAMR_Upper","cens","State_FIPS")]), ]
  mm <- model.matrix(~ Smoking_Rate + EQI_factor, d,
                     contrasts.arg = list(EQI_factor=contr.treatment(5)))
  colnames(mm) <- make.names(colnames(mm))
  list(X=mm, names=colnames(mm), df=d)
}

build_design_multi <- function(d) {
  d <- d %>% mutate(
    EQI_Air_factor    = factor(EQI_Air,    levels=1:5),
    EQI_Water_factor  = factor(EQI_Water,  levels=1:5),
    EQI_Land_factor   = factor(EQI_Land,   levels=1:5),
    EQI_Built_factor  = factor(EQI_Built,  levels=1:5),
    EQI_Social_factor = factor(EQI_Social, levels=1:5)
  )
  d <- d[complete.cases(d[, c("Smoking_Rate","EQI_Air_factor","EQI_Water_factor","EQI_Land_factor",
                               "EQI_Built_factor","EQI_Social_factor",
                               "AAMR_Lower","AAMR_Upper","cens","State_FIPS")]), ]
  form <- ~ Smoking_Rate + EQI_Air_factor + EQI_Water_factor + EQI_Land_factor +
             EQI_Built_factor + EQI_Social_factor
  mm <- model.matrix(form, d,
                     contrasts.arg = list(
                       EQI_Air_factor    = contr.treatment(5),
                       EQI_Water_factor  = contr.treatment(5),
                       EQI_Land_factor   = contr.treatment(5),
                       EQI_Built_factor  = contr.treatment(5),
                       EQI_Social_factor = contr.treatment(5)))
  colnames(mm) <- make.names(colnames(mm))
  list(X=mm, names=colnames(mm), df=d)
}

# Fit Stan model; returns draws df (beta columns renamed) or NULL on failure.
fit_model <- function(des, label) {
  states  <- sort(unique(des$df$State_FIPS))
  s_idx   <- match(des$df$State_FIPS, states)
  data_list <- list(
    N       = nrow(des$df),
    S       = length(states),
    state   = s_idx,
    y_lower = des$df$AAMR_Lower,
    y_upper = des$df$AAMR_Upper,
    cens    = des$df$cens,
    K       = ncol(des$X),
    X       = des$X
  )
  init_fn <- function() list(beta=rep(0, data_list$K), z_u=rep(0, data_list$S), sigma=50, sigma_u=10)
  fit <- try(mod$sample(
    data            = data_list,
    chains          = opt$chains,
    iter_sampling   = opt$iter - opt$warmup,
    iter_warmup     = opt$warmup,
    adapt_delta     = opt$`adapt-delta`,
    max_treedepth   = opt$`max-treedepth`,
    parallel_chains = min(opt$chains, cores_used),
    refresh         = 0,
    seed            = opt$seed,
    init            = rep(list(init_fn()), opt$chains)
  ), silent=TRUE)
  if (inherits(fit, "try-error")) { message("[Fail] ", label); return(NULL) }
  draws <- as_draws_df(fit$draws("beta"))
  colnames(draws) <- paste0("beta[", seq_len(ncol(draws)), "]")
  draws
}

# Compute MRR for quintiles Q1-Q5 vs ref_draws (Lag5 Q1).
# eqi_prefix: column name prefix in names_vec, e.g. "EQI_factor" or "EQI_Air_factor"
compute_mrr <- function(draws, names_vec, layer_dt, cancer_label, eqi_out, aamr_out, lagv,
                        ref_draws, model_label, eqi_prefix="EQI_factor") {
  smoking_idx  <- match("Smoking_Rate", names_vec)
  mean_smoking <- mean(layer_dt$Smoking_Rate, na.rm=TRUE)
  mu_Q1_draws  <- draws[["beta[1]"]] + draws[[paste0("beta[", smoking_idx, "]")]] * mean_smoking
  mrr_rows <- lapply(1:5, function(q) {
    if (q == 1) {
      mu_Qq <- mu_Q1_draws
    } else {
      q_idx <- match(paste0(eqi_prefix, q), names_vec)
      if (is.na(q_idx)) return(NULL)
      mu_Qq <- mu_Q1_draws + draws[[paste0("beta[", q_idx, "]")]]
    }
    MRR_q <- mu_Qq / ref_draws
    mrr_c <- MRR_q - 1
    pos <- sum(mrr_c > 0, na.rm=TRUE); neg <- sum(mrr_c < 0, na.rm=TRUE); nn <- pos + neg
    p_raw <- if (nn == 0) NA_real_ else 2 * min((pos + 0.5)/(nn + 1), (neg + 0.5)/(nn + 1))
    tibble(
      ICD_Code    = cancer_label,
      EQI_Period  = eqi_out,
      AAMR_Period = aamr_out,
      Lag         = lagv,
      Model       = model_label,
      Quintile    = paste0("Q", q),
      MRR_mean    = round(mean(MRR_q, na.rm=TRUE), 4),
      MRR_lower   = round(quantile(MRR_q, 0.025, na.rm=TRUE), 4),
      MRR_upper   = round(quantile(MRR_q, 0.975, na.rm=TRUE), 4),
      pct_diff    = round((mean(MRR_q, na.rm=TRUE) - 1) * 100, 4),
      p           = format_p(p_raw)
    )
  })
  bind_rows(Filter(Negate(is.null), mrr_rows))
}

# Pairwise lag test on overall EQI Q5 beta draws (EQI0005 lags only).
run_lag_test <- function(lag_q5_store, cancer_label, model_label, out_dir) {
  pairs <- list(c("5","10"), c("10","15"), c("15","5"))
  rows  <- lapply(pairs, function(p) {
    la <- p[1]; lb <- p[2]
    if (!la %in% names(lag_q5_store) || !lb %in% names(lag_q5_store)) return(NULL)
    da <- lag_q5_store[[la]]; db <- lag_q5_store[[lb]]
    diff_draws <- da - db
    tibble(
      ICD_Code   = cancer_label,
      Model      = model_label,
      comparison = paste0("lag", la, "_vs_lag", lb),
      diff_mean  = round(mean(diff_draws, na.rm=TRUE), 4),
      diff_lower = round(quantile(diff_draws, 0.025, na.rm=TRUE), 4),
      diff_upper = round(quantile(diff_draws, 0.975, na.rm=TRUE), 4),
      P_a_gt_b   = format_p(mean(da > db, na.rm=TRUE))
    )
  })
  result <- bind_rows(Filter(Negate(is.null), rows))
  if (nrow(result) > 0) {
    lag_file <- file.path(out_dir, paste0(cancer_label, "_lag_test.csv"))
    append_rows(lag_file, result)
    message("[LAG TEST] ", nrow(result), " comparisons — ", cancer_label, " / ", model_label)
  }
}

# ---------------------------------------------------------------------------
# Core per-stratum runner
#
# cancer_label  : ICD_Code written to output (e.g. "NDD_Asian", "C00_C97")
# model_prefix  : prefix for the Model column:
#                 sex/race → "Stratified"   → Model = "Stratified_EQI", "Stratified_Air" …
#                 Typology → "Typology_Farming" → Model = "Typology_Farming_EQI" …
# mrr_file      : output CSV path
# ---------------------------------------------------------------------------
run_stratum_mrr <- function(dt_sub, cancer_label, model_prefix, mrr_file, out_dir) {
  req_cols <- c("COUNTY_FIPS","EQI_Period","Time_Period","Cancer_Type",
                "AAMR_Lower","AAMR_Upper","Smoking_Rate","RUCC",
                "EQI","EQI_Air","EQI_Water","EQI_Land","EQI_Built","EQI_Social",
                "State_FIPS")
  miss <- setdiff(req_cols, names(dt_sub))
  if (length(miss)) {
    message("[Skip] Missing cols for ", cancer_label, ": ", paste(miss, collapse=",")); return(invisible(NULL))
  }

  lag5_ref     <- list()   # key: eqi_p — overall EQI Q1 Lag5 draws (universal reference)
  lag_q5_store <- list()   # key: lag value — EQI0005 Q5 draws for lag test

  for (sc in scenario_list) {
    scen_key <- sc$key; eqi_p <- sc$eqi; aamr_p <- sc$aamr; lagv <- sc$lag
    scen_dt  <- dt_sub[EQI_Period == eqi_p & Time_Period == aamr_p & Cancer_Type == cancer_label]
    if (nrow(scen_dt) < opt$`min-n`) {
      message("[Skip] ", model_prefix, " ", scen_key, " n=", nrow(scen_dt)); next
    }
    eqi_out  <- gsub('-', '_', eqi_p)
    aamr_out <- gsub('-', '_', aamr_p)

    # ---- Overall EQI model ----
    des_o  <- build_design_overall(scen_dt)
    draws_o <- fit_model(des_o, paste(model_prefix, scen_key, "EQI"))
    if (!is.null(draws_o)) {
      # Cache Lag5 Q1 as universal reference for this EQI period
      smoking_idx  <- match("Smoking_Rate", des_o$names)
      mean_smoking <- mean(scen_dt$Smoking_Rate, na.rm=TRUE)
      mu_Q1 <- draws_o[["beta[1]"]] + draws_o[[paste0("beta[", smoking_idx, "]")]] * mean_smoking
      if (lagv == 5) lag5_ref[[eqi_p]] <- mu_Q1

      if (!is.null(lag5_ref[[eqi_p]])) {
        mrr_df <- compute_mrr(draws_o, des_o$names, scen_dt, cancer_label, eqi_out, aamr_out,
                              lagv, ref_draws=lag5_ref[[eqi_p]],
                              model_label=paste0(model_prefix, "_EQI"),
                              eqi_prefix="EQI_factor")
        if (nrow(mrr_df) > 0) append_rows(mrr_file, mrr_df)
        message("[OK] MRR ", model_prefix, "_EQI ", scen_key)
      } else {
        message("[WARN] lag5_ref not available for ", model_prefix, " eqi=", eqi_p, " lag=", lagv)
      }

      # Store EQI0005 Q5 draws for lag test
      if (eqi_p == "2000-2005") {
        q5_idx <- match("EQI_factor5", des_o$names)
        if (!is.na(q5_idx))
          lag_q5_store[[as.character(lagv)]] <- draws_o[[paste0("beta[", q5_idx, "]")]]
      }
    }

    # ---- Multi-domain (Air/Water/Land/Built/Social) models ----
    if (!is.null(lag5_ref[[eqi_p]])) {
      des_m  <- build_design_multi(scen_dt)
      draws_m <- fit_model(des_m, paste(model_prefix, scen_key, "Multi"))
      if (!is.null(draws_m)) {
        for (dm in domain_models) {
          mrr_df <- compute_mrr(draws_m, des_m$names, scen_dt, cancer_label, eqi_out, aamr_out,
                                lagv, ref_draws=lag5_ref[[eqi_p]],
                                model_label=paste0(model_prefix, "_", dm$suffix),
                                eqi_prefix=dm$prefix)
          if (nrow(mrr_df) > 0) append_rows(mrr_file, mrr_df)
        }
        message("[OK] MRR ", model_prefix, " domains ", scen_key)
      }
    }
  }

  # Pairwise lag test on overall EQI Q5 (EQI0005 only)
  run_lag_test(lag_q5_store, cancer_label, paste0(model_prefix, "_EQI"), out_dir)
  invisible(NULL)
}

# ===========================================================================
# Section 1: Sex/race stratified MRR
# ICD_Code = cancer_full (e.g. NDD_Asian); Model = Stratified_EQI / Stratified_Air …
# ===========================================================================
message("\n===== Section 1: Stratified (sex/race) — Base disease: ", base_cancer, " =====")

strat_path <- file.path(project_root, opt$`data-strat`)
if (!file.exists(strat_path)) {
  message("[WARN] Stratified data not found: ", strat_path, " — skipping section 1")
} else {
  dt_strat <- fread(strat_path)
  if (!"State_FIPS" %in% names(dt_strat))
    dt_strat[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]
  dt_strat <- dt_strat[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
  dt_strat[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]
  dt_strat <- dt_strat[RUCC %in% 1:4 | is.na(RUCC)]

  # Discover all sex/race variants: e.g. C00_C97_Male, NDD_Female
  all_types     <- unique(dt_strat$Cancer_Type)
  strat_cancers <- all_types[startsWith(all_types, paste0(base_cancer, "_"))]

  if (length(strat_cancers) == 0) {
    message("[INFO] No sex/race strata found for ", base_cancer, " — skipping section 1")
  } else {
    message("Found strata: ", paste(strat_cancers, collapse=", "))
    mrr_file_strat <- file.path(out_dir, paste0(base_cancer, "_Stratified_MRR.csv"))

    for (cancer_full in strat_cancers) {
      message("  Processing: ", cancer_full)
      dt_sub <- dt_strat[Cancer_Type == cancer_full]
      # model_prefix = "Stratified"; stratum is encoded in cancer_full (ICD_Code column)
      run_stratum_mrr(dt_sub, cancer_full, "Stratified", mrr_file_strat, out_dir)
    }
    message("[Done] Stratified MRR for ", base_cancer)
  }
}

# ===========================================================================
# Section 2: Typology MRR
# ICD_Code = base_cancer; Model = Typology_Farming_EQI / Typology_Farming_Air …
# ===========================================================================
message("\n===== Section 2: Typology — Base disease: ", base_cancer, " =====")

typo_path <- file.path(project_root, opt$`data-typo`)
if (!file.exists(typo_path)) {
  message("[WARN] Typology/LandUse data not found: ", typo_path, " — skipping section 2")
} else {
  dt_typo <- fread(typo_path)
  if (!"State_FIPS" %in% names(dt_typo))
    dt_typo[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]
  dt_typo <- dt_typo[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
  dt_typo[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]
  dt_typo <- dt_typo[RUCC %in% 1:4 | is.na(RUCC)]

  if (!base_cancer %in% unique(dt_typo$Cancer_Type)) {
    message("[WARN] ", base_cancer, " not found in Typology data — skipping section 2")
  } else {
    if (!"econdep" %in% names(dt_typo)) {
      message("[Skip] Typology — column 'econdep' not in data")
    } else {
      typo_values <- 1:6
      typo_labels <- c("Farming","Mining","Manufacturing","Government","Services","Nonspecialized")
      mrr_file_typo <- file.path(out_dir, paste0(base_cancer, "_Typology_MRR.csv"))

      for (i in seq_along(typo_values)) {
        sv    <- typo_values[i]
        label <- typo_labels[i]
        dt_sub <- dt_typo[Cancer_Type == base_cancer & econdep == sv]
        message("  Typology stratum: ", label, " (n=", nrow(dt_sub), ")")
        if (nrow(dt_sub) < opt$`min-n`) { message("  [Skip] too few rows"); next }
        run_stratum_mrr(dt_sub, base_cancer, paste0("Typology_", label), mrr_file_typo, out_dir)
      }
      message("[Done] Typology MRR for ", base_cancer)
    }
  }
}

message("\nAll analyses complete for: ", base_cancer)
message("Output directory: ", out_dir)
