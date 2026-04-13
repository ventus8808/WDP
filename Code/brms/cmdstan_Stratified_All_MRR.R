#!/usr/bin/env Rscript
# MRR pipeline for all stratified analyses — C00_C97 and G20_G30_G12.2_F01_F03 only.
#
# Re-runs Stan on raw data subsets to compute MRR (Q1-Q5 vs Lag5 Q1 of same stratum).
# --section controls which stratification this task handles (one task per Slurm job):
#   Male | Female | White | Black | Asian | Indian  → sex/race stratum
#   Typology | RUCC | Koppen | CensusRegion         → full section (all strata inside)
#
# Output directory: Result/Stratified_MRR/
# Output files:
#   {cancer}_{sex/race}_MRR.csv   — e.g. NDD_Male_MRR.csv
#   {cancer}_Typology_MRR.csv
#   {cancer}_RUCC_MRR.csv
#   {cancer}_Koppen_MRR.csv
#   {cancer}_CensusRegion_MRR.csv

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(readr)
  library(cmdstanr)
  library(posterior)
})
utils::globalVariables(c('EQI','EQI_Air','EQI_Water','EQI_Land','EQI_Built','EQI_Social',
                          'Smoking_Rate','State_FIPS','AAMR_Lower','AAMR_Upper','cens','RUCC'))

valid_sections <- c("Male","Female","White","Black","Asian","Indian",
                    "Typology","RUCC","Koppen","CensusRegion")

option_list <- list(
  make_option(c("--cancer-type"), type="character", default=NA,
              help="Base ICD code: C00_C97 or G20_G30_G12.2_F01_F03. Required."),
  make_option(c("--section"),     type="character", default=NA,
              help=paste("Which stratification to run. One of:", paste(valid_sections, collapse=", "))),
  make_option(c("--data-strat"), type="character",
              default="Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Stratifed.csv",
              help="Sex/race stratified data"),
  make_option(c("--data-main"), type="character",
              default="Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Cluster_Climate_Typology_LandUse.csv",
              help="Main data (Typology, RUCC, Köppen, Census Region)"),
  make_option(c("--output-dir"), type="character", default="Result/Stratified_MRR"),
  make_option(c("--chains"),      type="integer", default=4),
  make_option(c("--iter"),        type="integer", default=2000),
  make_option(c("--warmup"),      type="integer", default=1000),
  make_option(c("--adapt-delta"), type="double",  default=0.95),
  make_option(c("--max-treedepth"), type="integer", default=12),
  make_option(c("--min-n"),       type="integer", default=50),
  make_option(c("--seed"),        type="integer", default=1234),
  make_option(c("--test"),        action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.na(opt$`cancer-type`)) stop("--cancer-type is required")
if (!opt$`cancer-type` %in% c("C00_C97", "G20_G30_G12.2_F01_F03"))
  stop("--cancer-type must be C00_C97 or G20_G30_G12.2_F01_F03")
if (is.na(opt$section) || !opt$section %in% valid_sections)
  stop("--section is required. Must be one of: ", paste(valid_sections, collapse=", "))
if (opt$test) {
  opt$iter   <- min(opt$iter,   800)
  opt$warmup <- min(opt$warmup, 300)
  message("[TEST MODE] iter=", opt$iter, " warmup=", opt$warmup)
}
set.seed(opt$seed)

cores_avail <- parallel::detectCores(logical=TRUE)
cores_used  <- max(1, floor(cores_avail * 0.8))
options(mc.cores=cores_used)
message("Cores: ", cores_avail, " detected | ", cores_used, " used")

# ---------------------------------------------------------------------------
# Stan model
# ---------------------------------------------------------------------------
stan_code <- "data {\n  int<lower=1> N;\n  int<lower=1> S;\n  array[N] int<lower=1,upper=S> state;\n  vector[N] y_lower;\n  vector[N] y_upper;\n  array[N] int<lower=0,upper=2> cens;\n  int<lower=1> K;\n  matrix[N,K] X;\n}\nparameters {\n  vector[K] beta;\n  vector[S] z_u;\n  real<lower=0> sigma;\n  real<lower=0> sigma_u;\n}\ntransformed parameters {\n  vector[S] u = sigma_u * z_u;\n}\nmodel {\n  beta ~ normal(0,5);\n  z_u  ~ normal(0,1);\n  sigma   ~ exponential(1);\n  sigma_u ~ exponential(1);\n  for (i in 1:N) {\n    real mu = X[i] * beta + u[state[i]];\n    if (cens[i]==0) {\n      target += normal_lpdf(y_lower[i] | mu, sigma);\n    } else {\n      real p_up = normal_cdf(y_upper[i] | mu, sigma);\n      real p_lo = normal_cdf(y_lower[i] | mu, sigma);\n      real diff = fmax(p_up - p_lo, 1e-12);\n      target += log(diff);\n    }\n  }\n}"
stan_file <- file.path(tempdir(), "interval_mixed.stan")
writeLines(stan_code, stan_file)
mod <- cmdstan_model(stan_file)

project_root <- normalizePath(".")
base_cancer  <- opt$`cancer-type`
out_dir      <- file.path(project_root, opt$`output-dir`)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

# ---------------------------------------------------------------------------
# Scenarios: EQI0005 (lag 5/10/15) + EQI0610 (lag 5/10)
# ---------------------------------------------------------------------------
scenario_list <- list(
  list(key="EQI0005_lag05",  eqi="2000-2005", aamr="2006-2010", lag=5),
  list(key="EQI0005_lag10",  eqi="2000-2005", aamr="2011-2015", lag=10),
  list(key="EQI0005_lag15",  eqi="2000-2005", aamr="2016-2020", lag=15),
  list(key="EQI0610_lag05",  eqi="2006-2010", aamr="2011-2015", lag=5),
  list(key="EQI0610_lag10",  eqi="2006-2010", aamr="2016-2020", lag=10)
)

# Domain models: (display suffix for Model column, EQI factor prefix in design matrix)
domain_models <- list(
  list(suffix="EQI_Air",    prefix="EQI_Air_factor"),
  list(suffix="EQI_Water",  prefix="EQI_Water_factor"),
  list(suffix="EQI_Land",   prefix="EQI_Land_factor"),
  list(suffix="EQI_Built",  prefix="EQI_Built_factor"),
  list(suffix="EQI_Social", prefix="EQI_Social_factor")
)

# ---------------------------------------------------------------------------
# Helpers
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

prep_data <- function(dt) {
  if (!"State_FIPS" %in% names(dt))
    dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]
  dt <- dt[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
  dt[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]
  dt
}

build_design_overall <- function(d) {
  d <- d %>% mutate(EQI_factor = factor(EQI, levels=1:5))
  d <- d[complete.cases(d[, c("Smoking_Rate","EQI_factor","AAMR_Lower","AAMR_Upper","cens","State_FIPS")]), ]
  mm <- model.matrix(~ Smoking_Rate + EQI_factor, d,
                     contrasts.arg=list(EQI_factor=contr.treatment(5)))
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
  mm <- model.matrix(
    ~ Smoking_Rate + EQI_Air_factor + EQI_Water_factor + EQI_Land_factor +
        EQI_Built_factor + EQI_Social_factor, d,
    contrasts.arg=list(
      EQI_Air_factor    = contr.treatment(5), EQI_Water_factor  = contr.treatment(5),
      EQI_Land_factor   = contr.treatment(5), EQI_Built_factor  = contr.treatment(5),
      EQI_Social_factor = contr.treatment(5)))
  colnames(mm) <- make.names(colnames(mm))
  list(X=mm, names=colnames(mm), df=d)
}

fit_stan <- function(des, label) {
  states <- sort(unique(des$df$State_FIPS))
  s_idx  <- match(des$df$State_FIPS, states)
  dl <- list(N=nrow(des$df), S=length(states), state=s_idx,
             y_lower=des$df$AAMR_Lower, y_upper=des$df$AAMR_Upper, cens=des$df$cens,
             K=ncol(des$X), X=des$X)
  init_fn <- function() list(beta=rep(0, dl$K), z_u=rep(0, dl$S), sigma=50, sigma_u=10)
  fit <- try(mod$sample(
    data=dl, chains=opt$chains,
    iter_sampling=opt$iter - opt$warmup, iter_warmup=opt$warmup,
    adapt_delta=opt$`adapt-delta`, max_treedepth=opt$`max-treedepth`,
    parallel_chains=min(opt$chains, cores_used),
    refresh=0, seed=opt$seed, init=rep(list(init_fn()), opt$chains)
  ), silent=TRUE)
  if (inherits(fit, "try-error")) { message("[Fail] ", label); return(NULL) }
  draws <- as_draws_df(fit$draws("beta"))
  colnames(draws) <- paste0("beta[", seq_len(ncol(draws)), "]")
  draws
}

compute_mrr <- function(draws, names_vec, scen_dt, icd_label, eqi_out, aamr_out, lagv,
                        ref_draws, model_label, eqi_prefix="EQI_factor") {
  sm_idx       <- match("Smoking_Rate", names_vec)
  mean_smoking <- mean(scen_dt$Smoking_Rate, na.rm=TRUE)
  mu_Q1        <- draws[["beta[1]"]] + draws[[paste0("beta[", sm_idx, "]")]] * mean_smoking
  rows <- lapply(1:5, function(q) {
    mu_Qq <- if (q == 1) mu_Q1 else {
      qi <- match(paste0(eqi_prefix, q), names_vec)
      if (is.na(qi)) return(NULL)
      mu_Q1 + draws[[paste0("beta[", qi, "]")]]
    }
    mrr   <- mu_Qq / ref_draws
    delta <- mrr - 1
    pos <- sum(delta > 0, na.rm=TRUE); neg <- sum(delta < 0, na.rm=TRUE); nn <- pos + neg
    p_raw <- if (nn == 0) NA_real_ else 2 * min((pos+0.5)/(nn+1), (neg+0.5)/(nn+1))
    tibble(
      ICD_Code    = icd_label,
      EQI_Period  = eqi_out,
      AAMR_Period = aamr_out,
      Lag         = lagv,
      Model       = model_label,
      Quintile    = paste0("Q", q),
      MRR_mean    = round(mean(mrr, na.rm=TRUE), 4),
      MRR_lower   = round(quantile(mrr, 0.025, na.rm=TRUE), 4),
      MRR_upper   = round(quantile(mrr, 0.975, na.rm=TRUE), 4),
      pct_diff    = round((mean(mrr, na.rm=TRUE) - 1) * 100, 4),
      p           = format_p(p_raw)
    )
  })
  bind_rows(Filter(Negate(is.null), rows))
}

# ---------------------------------------------------------------------------
# Core runner: processes all scenarios for one (data subset, icd_label, model_prefix, output_file)
# model_prefix_overall: e.g. "Stratified_EQI", "Typology_Farming_EQI", "RUCC1_EQI", "koppen_major_B_EQI"
# model_prefix_domain:  e.g. "Stratified_", "Typology_Farming_", "RUCC1_", "koppen_major_B_"
#                        → domain suffix appended: "EQI_Air", "EQI_Water" etc.
# ---------------------------------------------------------------------------
run_strat_mrr <- function(dt_sub, icd_label, model_prefix, mrr_file) {
  req <- c("COUNTY_FIPS","EQI_Period","Time_Period","Cancer_Type",
           "AAMR_Lower","AAMR_Upper","Smoking_Rate","RUCC",
           "EQI","EQI_Air","EQI_Water","EQI_Land","EQI_Built","EQI_Social","State_FIPS")
  miss <- setdiff(req, names(dt_sub))
  if (length(miss)) { message("[Skip] Missing cols: ", paste(miss, collapse=",")); return(invisible(NULL)) }

  lag5_ref     <- list()   # key: eqi_p — Lag5 Q1 mu draws (reference)
  lag_q5_store <- list()   # key: lag — EQI0005 Q5 beta draws for lag test

  for (sc in scenario_list) {
    scen_key <- sc$key; eqi_p <- sc$eqi; aamr_p <- sc$aamr; lagv <- sc$lag
    scen_dt  <- dt_sub[EQI_Period == eqi_p & Time_Period == aamr_p & Cancer_Type == icd_label]
    if (nrow(scen_dt) < opt$`min-n`) {
      message("[Skip] ", model_prefix, " ", scen_key, " n=", nrow(scen_dt)); next
    }
    eqi_out  <- gsub('-', '_', eqi_p)
    aamr_out <- gsub('-', '_', aamr_p)

    # ---- Overall EQI model ----
    des_o   <- build_design_overall(scen_dt)
    draws_o <- fit_stan(des_o, paste(model_prefix, "EQI", scen_key))

    if (!is.null(draws_o)) {
      sm_idx  <- match("Smoking_Rate", des_o$names)
      ms      <- mean(scen_dt$Smoking_Rate, na.rm=TRUE)
      mu_Q1   <- draws_o[["beta[1]"]] + draws_o[[paste0("beta[", sm_idx, "]")]] * ms
      if (lagv == 5) lag5_ref[[eqi_p]] <- mu_Q1

      if (!is.null(lag5_ref[[eqi_p]])) {
        mrr_df <- compute_mrr(draws_o, des_o$names, scen_dt,
                              icd_label, eqi_out, aamr_out, lagv,
                              ref_draws  = lag5_ref[[eqi_p]],
                              model_label = paste0(model_prefix, "_EQI"))
        if (nrow(mrr_df) > 0) append_rows(mrr_file, mrr_df)
        message("[OK] ", model_prefix, "_EQI ", scen_key)

        if (eqi_p == "2000-2005") {
          q5i <- match("EQI_factor5", des_o$names)
          if (!is.na(q5i)) lag_q5_store[[as.character(lagv)]] <- draws_o[[paste0("beta[", q5i, "]")]]
        }
      } else {
        message("[WARN] No lag5_ref for ", model_prefix, " eqi=", eqi_p, " lag=", lagv)
      }
    }

    # ---- Multi-domain models (Air/Water/Land/Built/Social) ----
    if (!is.null(lag5_ref[[eqi_p]])) {
      des_m   <- build_design_multi(scen_dt)
      draws_m <- fit_stan(des_m, paste(model_prefix, "Multi", scen_key))
      if (!is.null(draws_m)) {
        for (dm in domain_models) {
          mrr_df <- compute_mrr(draws_m, des_m$names, scen_dt,
                                icd_label, eqi_out, aamr_out, lagv,
                                ref_draws   = lag5_ref[[eqi_p]],
                                model_label = paste0(model_prefix, "_", dm$suffix),
                                eqi_prefix  = dm$prefix)
          if (nrow(mrr_df) > 0) append_rows(mrr_file, mrr_df)
        }
        message("[OK] ", model_prefix, " domains ", scen_key)
      }
    }
  }

  # Lag test: EQI0005 Q5 (lag5 vs lag10 vs lag15)
  pairs <- list(c("5","10"), c("10","15"), c("15","5"))
  lag_rows <- lapply(pairs, function(p) {
    la <- p[1]; lb <- p[2]
    if (!la %in% names(lag_q5_store) || !lb %in% names(lag_q5_store)) return(NULL)
    da <- lag_q5_store[[la]]; db <- lag_q5_store[[lb]]
    d  <- da - db
    tibble(ICD_Code=icd_label, Model=paste0(model_prefix,"_EQI"),
           comparison=paste0("lag",la,"_vs_lag",lb),
           diff_mean=round(mean(d,na.rm=TRUE),4),
           diff_lower=round(quantile(d,0.025,na.rm=TRUE),4),
           diff_upper=round(quantile(d,0.975,na.rm=TRUE),4),
           P_a_gt_b=format_p(mean(da>db,na.rm=TRUE)))
  })
  lag_result <- bind_rows(Filter(Negate(is.null), lag_rows))
  if (nrow(lag_result) > 0) {
    lag_file <- file.path(out_dir, paste0(icd_label, "_lag_test.csv"))
    append_rows(lag_file, lag_result)
  }
  invisible(NULL)
}

# ===========================================================================
# LOAD DATA (only what this section needs)
# ===========================================================================
sec <- opt$section
path_strat <- file.path(project_root, opt$`data-strat`)
path_main  <- file.path(project_root, opt$`data-main`)

sex_race_sections <- c("Male","Female","White","Black","Asian","Indian")
needs_strat <- sec %in% sex_race_sections
needs_main  <- sec %in% c("Typology","RUCC","Koppen","CensusRegion")

dt_strat <- NULL
dt_main  <- NULL

if (needs_strat) {
  if (!file.exists(path_strat)) stop("Stratified data not found: ", path_strat)
  dt_strat <- prep_data(fread(path_strat))
}
if (needs_main) {
  if (!file.exists(path_main)) stop("Main data not found: ", path_main)
  dt_main <- prep_data(fread(path_main))[Cancer_Type == base_cancer]
}

message("\n===== Section: ", sec, " | Disease: ", base_cancer, " =====")

# ===========================================================================
# Sex / Race  — one file per stratum
# ICD_Code = {base_cancer}_{stratum}   Model = Stratified_EQI / Stratified_EQI_Air …
# ===========================================================================
if (needs_strat) {
  # Stratified data uses short prefix "NDD" for the combined NDD ICD code
  strat_prefix_map <- c("G20_G30_G12.2_F01_F03" = "NDD")
  strat_prefix <- if (base_cancer %in% names(strat_prefix_map)) strat_prefix_map[[base_cancer]] else base_cancer
  cancer_full <- paste0(strat_prefix, "_", sec)
  if (!cancer_full %in% unique(dt_strat$Cancer_Type))
    stop(cancer_full, " not found in stratified data")
  mrr_out <- file.path(out_dir, paste0(base_cancer, "_", sec, "_MRR.csv"))
  message("  → ", cancer_full)
  dt_sub  <- dt_strat[Cancer_Type == cancer_full]
  run_strat_mrr(dt_sub, cancer_full, "Stratified", mrr_out)
}

# ===========================================================================
# Typology  (econdep 1-6)
# ICD_Code = base_cancer   Model = Typology_Farming_EQI / Typology_Farming_EQI_Air …
# ===========================================================================
if (sec == "Typology") {
  if (!"econdep" %in% names(dt_main)) stop("econdep column not found in main data")
  typo_values <- 1:6
  typo_labels <- c("Farming","Mining","Manufacturing","Government","Services","Nonspecialized")
  mrr_out <- file.path(out_dir, paste0(base_cancer, "_Typology_MRR.csv"))
  for (i in seq_along(typo_values)) {
    sv <- typo_values[i]; sl <- typo_labels[i]
    label  <- paste0("Typology_", sl)
    dt_sub <- dt_main[econdep == sv]
    message("  → ", label, " (n=", nrow(dt_sub), ")")
    if (nrow(dt_sub) < opt$`min-n`) { message("  [Skip] too few rows"); next }
    run_strat_mrr(dt_sub, base_cancer, label, mrr_out)
  }
}

# ===========================================================================
# RUCC  (1-4)
# ICD_Code = base_cancer   Model = RUCC1_EQI / RUCC1_EQI_Air …
# ===========================================================================
if (sec == "RUCC") {
  if (!"RUCC" %in% names(dt_main)) stop("RUCC column not found in main data")
  mrr_out <- file.path(out_dir, paste0(base_cancer, "_RUCC_MRR.csv"))
  for (rv in 1:4) {
    label  <- paste0("RUCC", rv)
    dt_sub <- dt_main[RUCC == rv]
    message("  → ", label, " (n=", nrow(dt_sub), ")")
    if (nrow(dt_sub) < opt$`min-n`) { message("  [Skip] too few rows"); next }
    run_strat_mrr(dt_sub, base_cancer, label, mrr_out)
  }
}

# ===========================================================================
# Köppen-Geiger  (koppen_major: B/C/D)
# ICD_Code = base_cancer   Model = koppen_major_B_EQI / koppen_major_B_EQI_Air …
# ===========================================================================
if (sec == "Koppen") {
  if (!"koppen_major" %in% names(dt_main)) stop("koppen_major column not found in main data")
  mrr_out <- file.path(out_dir, paste0(base_cancer, "_Koppen_MRR.csv"))
  for (kv in c("B","C","D")) {
    label  <- paste0("koppen_major_", kv)
    dt_sub <- dt_main[koppen_major == kv]
    message("  → ", label, " (n=", nrow(dt_sub), ")")
    if (nrow(dt_sub) < opt$`min-n`) { message("  [Skip] too few rows"); next }
    run_strat_mrr(dt_sub, base_cancer, label, mrr_out)
  }
}

# ===========================================================================
# Census Region  (census_region: 1-4)
# ICD_Code = base_cancer   Model = census_region_1_EQI / census_region_1_EQI_Air …
# ===========================================================================
if (sec == "CensusRegion") {
  if (!"census_region" %in% names(dt_main)) stop("census_region column not found in main data")
  mrr_out <- file.path(out_dir, paste0(base_cancer, "_CensusRegion_MRR.csv"))
  for (rv in 1:4) {
    label  <- paste0("census_region_", rv)
    dt_sub <- dt_main[census_region == rv]
    message("  → ", label, " (n=", nrow(dt_sub), ")")
    if (nrow(dt_sub) < opt$`min-n`) { message("  [Skip] too few rows"); next }
    run_strat_mrr(dt_sub, base_cancer, label, mrr_out)
  }
}

message("\n===== Done: ", sec, " | ", base_cancer, " =====")
message("Output: ", out_dir)
