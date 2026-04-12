#!/usr/bin/env Rscript
# Sensitivity analysis: EQI × healthcare / environment covariates
# Input:  Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Sensativity.csv
# Each invocation runs ONE model group (--model-group) for ONE cancer type (--cancer-types).
# Writes split files: {prefix}_sensitivity_EQI_{GroupTag}.csv
#                     {prefix}_sensitivity_Covar_{GroupTag}.csv
# A separate merge step (merge_sensativity_covar.R) combines splits into final CSVs.
#
# Model groups (5 total):
#   Baseline_Full  -> Baseline + Full (both in one task)
#   Insurance      -> +Insurance only
#   Doctor         -> +Doctor only
#   Forest         -> +Forest only
#   Monitoring     -> +Monitoring only

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
utils::globalVariables(c("EQI","Smoking_Rate","State_FIPS",
                         "uninsured_rate_2005","physician_density_per100k",
                         "forest_coverage","site_number_mean"))

option_list <- list(
  make_option(c("--data"),          type="character", default="Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Sensativity.csv"),
  make_option(c("--output-dir"),    type="character", default="Result/brms_Sensativity"),
  make_option(c("--cancer-types"),  type="character", default=NA, help="Single cancer type code"),
  make_option(c("--model-group"),   type="character", default=NA,
              help="One of: Baseline_Full | Insurance | Doctor | Forest | Monitoring"),
  make_option(c("--chains"),        type="integer",   default=4),
  make_option(c("--iter"),          type="integer",   default=2000),
  make_option(c("--warmup"),        type="integer",   default=1000),
  make_option(c("--adapt-delta"),   type="double",    default=0.95),
  make_option(c("--max-treedepth"), type="integer",   default=12),
  make_option(c("--min-n"),         type="integer",   default=50),
  make_option(c("--seed"),          type="integer",   default=1234),
  make_option(c("--test"),          action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))
if(opt$test){ opt$iter <- min(opt$iter,800); opt$warmup <- min(opt$warmup,300); message("[TEST MODE] iter=",opt$iter," warmup=",opt$warmup) }
set.seed(opt$seed)

if(is.na(opt$`cancer-types`)) stop("--cancer-types is required")
if(is.na(opt$`model-group`))  stop("--model-group is required")

cancer    <- str_trim(opt$`cancer-types`)
mod_group <- str_trim(opt$`model-group`)

cores_avail <- parallel::detectCores(logical=TRUE); cores_used <- max(1,floor(cores_avail*0.8)); options(mc.cores=cores_used)
message("Detected cores: ", cores_avail, " | Using: ", cores_used)
message("Cancer: ", cancer, "  |  Model group: ", mod_group)

# ---- Stan model (identical to cmdstan_main.R) ----
stan_code <- "data {\n  int<lower=1> N;\n  int<lower=1> S;\n  array[N] int<lower=1,upper=S> state;\n  vector[N] y_lower;\n  vector[N] y_upper;\n  array[N] int<lower=0,upper=2> cens;\n  int<lower=1> K;\n  matrix[N,K] X;\n} \nparameters {\n  vector[K] beta;\n  vector[S] z_u;\n  real<lower=0> sigma;\n  real<lower=0> sigma_u;\n} \ntransformed parameters {\n  vector[S] u = sigma_u * z_u;\n} \nmodel {\n  beta ~ normal(0,5);\n  z_u ~ normal(0,1);\n  sigma ~ exponential(1);\n  sigma_u ~ exponential(1);\n  for (i in 1:N) {\n    real mu = X[i] * beta + u[state[i]];\n    if (cens[i]==0) {\n      target += normal_lpdf(y_lower[i] | mu, sigma);\n    } else {\n      real p_up = normal_cdf(y_upper[i] | mu, sigma);\n      real p_lo = normal_cdf(y_lower[i] | mu, sigma);\n      real diff = fmax(p_up - p_lo, 1e-12);\n      target += log(diff);\n    }\n  }\n}"
stan_file <- file.path(tempdir(), "interval_mixed_model.stan"); writeLines(stan_code, stan_file)
mod <- cmdstan_model(stan_file)

# ---- Load & validate data ----
project_root <- normalizePath(".")
path <- file.path(project_root, opt$data); if(!file.exists(path)) stop("Data not found: ", path)
dt <- fread(path)

req <- c("COUNTY_FIPS","EQI_Period","Time_Period","Lag_Years","Cancer_Type",
         "AAMR_Lower","AAMR_Upper","Smoking_Rate","RUCC","EQI",
         "uninsured_rate_2005","physician_density_per100k","forest_coverage","site_number_mean")
miss <- setdiff(req, names(dt)); if(length(miss)) stop("Missing columns: ", paste(miss, collapse=","))

if(!"State_FIPS" %in% names(dt)) dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]

dt <- dt[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
dt[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]
dt <- dt[EQI_Period == "2000-2005"]
dt <- dt[RUCC %in% 1:4 | is.na(RUCC)]

if(!cancer %in% dt$Cancer_Type) stop("Cancer type not found in data: ", cancer)

# ---- Model group definitions ----
# covars: z-scored column names added beyond EQI_factor (intercept always included)
ALL_SENSITIVITY_MODELS <- list(
  Baseline   = list(tag="Baseline",    covars=character(0)),
  Full       = list(tag="Full",        covars=c("smoking_z","uninsured_z","physician_z","forest_z","monitor_z")),
  Insurance  = list(tag="+Insurance",  covars="uninsured_z"),
  Doctor     = list(tag="+Doctor",     covars="physician_z"),
  Forest     = list(tag="+Forest",     covars="forest_z"),
  Monitoring = list(tag="+Monitoring", covars="monitor_z")
)
GROUP_MODELS <- list(
  Baseline_Full = c("Baseline","Full"),
  Insurance     = "Insurance",
  Doctor        = "Doctor",
  Forest        = "Forest",
  Monitoring    = "Monitoring"
)
if(!mod_group %in% names(GROUP_MODELS))
  stop("Unknown --model-group: '", mod_group, "'. Valid: ", paste(names(GROUP_MODELS), collapse=", "))

run_models <- lapply(GROUP_MODELS[[mod_group]], function(k) ALL_SENSITIVITY_MODELS[[k]])

# ---- Covariate mapping: z-score col name -> original col name ----
COVAR_ORIG <- c(
  smoking_z   = "Smoking_Rate",
  uninsured_z = "uninsured_rate_2005",
  physician_z = "physician_density_per100k",
  forest_z    = "forest_coverage",
  monitor_z   = "site_number_mean"
)

# ---- Output paths ----
get_prefix <- function(x){
  if(x == "G20_G30_G12.2_F01_F03") return("NDD")
  if(x == "C00_C97")               return("Cancer")
  gsub("[^A-Za-z0-9]","_",x)
}
prefix  <- get_prefix(cancer)
out_dir <- file.path(project_root, opt$`output-dir`)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

# Split files — one per task, no concurrent-write collision
eqi_file   <- file.path(out_dir, paste0(prefix, "_sensitivity_EQI_",   mod_group, ".csv"))
covar_file <- file.path(out_dir, paste0(prefix, "_sensitivity_Covar_", mod_group, ".csv"))
if(file.exists(eqi_file))   file.remove(eqi_file)
if(file.exists(covar_file)) file.remove(covar_file)

# ---- Scenarios ----
scenario_list <- list(
  list(key="EQI0005_AAMR2006_2010", aamr="2006-2010", lag=5),
  list(key="EQI0005_AAMR2011_2015", aamr="2011-2015", lag=10),
  list(key="EQI0005_AAMR2016_2020", aamr="2016-2020", lag=15)
)

# ---- Helpers (identical to cmdstan_main.R) ----
sig_mark <- function(p){ if(is.na(p)) return(""); if(p<0.001) return("***"); if(p<0.01) return("**"); if(p<0.05) return("*"); "" }
format_cell <- function(draws){ if(length(draws)==0) return(""); ci <- quantile(draws,c(0.025,0.975),na.rm=TRUE); p <- 2*min(mean(draws>0),mean(draws<0)); sprintf("%0.2f(%0.2f,%0.2f)%s",mean(draws),ci[1],ci[2],sig_mark(p)) }
compute_p <- function(draws){
  if(length(draws)==0) return(NA_character_)
  pos <- sum(draws>0,na.rm=TRUE); neg <- sum(draws<0,na.rm=TRUE); n <- pos+neg
  if(n==0) return(NA_character_)
  p_pos <- (pos+0.5)/(n+1); p_neg <- (neg+0.5)/(n+1)
  sprintf("%.4f", 2*min(p_pos,p_neg))
}
append_rows <- function(path, df){ if(!file.exists(path)) write_csv(df,path) else suppressWarnings(write.table(df,path,sep=",",col.names=FALSE,row.names=FALSE,append=TRUE)) }

# ---- Design matrix ----
build_design_sensitivity <- function(d, extra_covars=character(0)){
  d <- d %>% mutate(EQI_factor = factor(EQI, levels=1:5))
  keep_cols <- c("EQI_factor","AAMR_Lower","AAMR_Upper","cens","State_FIPS", extra_covars)
  d <- d[complete.cases(d[, keep_cols, drop=FALSE]), ]
  rhs <- paste(c("EQI_factor", extra_covars), collapse=" + ")
  mm  <- model.matrix(as.formula(paste("~", rhs)), data=d,
                      contrasts.arg=list(EQI_factor=contr.treatment(5)))
  colnames(mm) <- make.names(colnames(mm))
  list(X=mm, names=colnames(mm), df=d)
}

# ---- Quintile extractors (same signature as cmdstan_main.R) ----
extract_quintiles <- function(draw_df, names_vec, prefix){
  out <- list(Q1="0.00", Q2="", Q3="", Q4="", Q5="")
  if(any(grepl(paste0(prefix,"\\.L"), names_vec))){ return(out) }
  for(q in 2:5){ nm <- paste0(prefix,q); idx <- match(nm,names_vec); out[[paste0("Q",q)]] <- if(is.na(idx)) "" else format_cell(draw_df[[paste0("beta[",idx,"]")]]) }
  out
}
extract_quintile_metrics <- function(draw_df, names_vec, prefix, summ_df){
  out <- list(Q2_p=NA_real_,Q3_p=NA_real_,Q4_p=NA_real_,Q5_p=NA_real_,
              Q2_rhat=NA_real_,Q3_rhat=NA_real_,Q4_rhat=NA_real_,Q5_rhat=NA_real_,
              Q2_ess_bulk=NA_real_,Q3_ess_bulk=NA_real_,Q4_ess_bulk=NA_real_,Q5_ess_bulk=NA_real_,
              Q2_ess_tail=NA_real_,Q3_ess_tail=NA_real_,Q4_ess_tail=NA_real_,Q5_ess_tail=NA_real_)
  if(any(grepl(paste0(prefix,"\\.L"), names_vec))){ return(out) }
  for(q in 2:5){
    nm <- paste0(prefix,q); idx <- match(nm,names_vec)
    if(!is.na(idx)){
      col <- paste0("beta[",idx,"]")
      out[[paste0("Q",q,"_p")]] <- compute_p(draw_df[[col]])
      sr <- summ_df[summ_df$variable==col,,drop=FALSE]
      if(nrow(sr)){ out[[paste0("Q",q,"_rhat")]] <- sr$rhat; out[[paste0("Q",q,"_ess_bulk")]] <- sr$ess_bulk; out[[paste0("Q",q,"_ess_tail")]] <- sr$ess_tail }
    }
  }
  out
}

# ---- Covariate extractor ----
extract_covar_cells <- function(draw_df, names_vec, model_covars){
  # names of cells are the original column names (values of COVAR_ORIG)
  orig_names <- unname(COVAR_ORIG)
  cells <- setNames(rep("", length(orig_names)), orig_names)
  for(z_col in model_covars){
    orig_col <- COVAR_ORIG[[z_col]]   # unname via [[, gives plain character
    idx <- match(z_col, names_vec)
    if(!is.na(idx)) cells[[orig_col]] <- format_cell(draw_df[[paste0("beta[",idx,"]")]])
  }
  cells
}

# ---- Prepare cancer data + z-scores ----
cancer_dt <- dt[Cancer_Type == cancer]

# Z-score on full cancer dataset; use [[ for safe column access, no get()
for(z_col in names(COVAR_ORIG)){
  orig_col <- COVAR_ORIG[[z_col]]          # plain character, no name attribute
  vals     <- cancer_dt[[orig_col]]
  mu_      <- mean(vals, na.rm=TRUE)
  sd_      <- sd(vals,   na.rm=TRUE)
  sd_safe  <- ifelse(is.na(sd_) | sd_==0, 1, sd_)
  cancer_dt[, (z_col) := (cancer_dt[[orig_col]] - mu_) / sd_safe]
}

# ---- Main loop ----
message("\n===== ", prefix, " | Group: ", mod_group, " =====")

for(sc in scenario_list){
  aamr_p   <- sc$aamr; lagv <- sc$lag; aamr_out <- gsub("-","_",aamr_p)
  scen_dt  <- cancer_dt[Time_Period == aamr_p]
  if(nrow(scen_dt) < opt$`min-n`){ message("[Skip] ",sc$key," n=",nrow(scen_dt)); next }

  for(sm in run_models){
    model_tag    <- sm$tag
    model_covars <- sm$covars

    des <- build_design_sensitivity(scen_dt, model_covars)
    if(nrow(des$df) < opt$`min-n`){ message("[Skip] ",sc$key," ",model_tag," n=",nrow(des$df)); next }

    states_s  <- sort(unique(des$df$State_FIPS)); state_idx <- match(des$df$State_FIPS, states_s)
    data_list <- list(N=nrow(des$df), S=length(states_s), state=state_idx,
                      y_lower=des$df$AAMR_Lower, y_upper=des$df$AAMR_Upper, cens=des$df$cens,
                      K=ncol(des$X), X=des$X)

    # init identical to cmdstan_main.R — plain function referencing data_list in enclosing scope
    init_fun <- function() list(beta=rep(0,data_list$K), z_u=rep(0,data_list$S), sigma=50, sigma_u=10)

    fit <- try(mod$sample(data=data_list, chains=opt$chains,
                          iter_sampling=opt$iter-opt$warmup, iter_warmup=opt$warmup,
                          adapt_delta=opt$`adapt-delta`, max_treedepth=opt$`max-treedepth`,
                          parallel_chains=min(opt$chains,cores_used), refresh=0, seed=opt$seed,
                          init=rep(list(init_fun()), opt$chains)), silent=FALSE)

    if(inherits(fit,"try-error")){ message("[Fail] ",sc$key," ",model_tag,": ",conditionMessage(attr(fit,"condition"))); next }

    draws <- as_draws_df(fit$draws("beta"))
    colnames(draws) <- paste0("beta[",seq_len(ncol(draws)),"]")
    summ  <- posterior::summarize_draws(fit$draws("beta"))

    # EQI quintile row
    q   <- extract_quintiles(draws, des$names, "EQI_factor")
    met <- extract_quintile_metrics(draws, des$names, "EQI_factor", summ)
    eqi_row <- tibble(
      ICD_Code=cancer, AAMR_Period=aamr_out, Lag=lagv, Sensitivity_Model=model_tag,
      Q1=q$Q1, Q2=q$Q2, Q3=q$Q3, Q4=q$Q4, Q5=q$Q5,
      Q2_p=met$Q2_p, Q3_p=met$Q3_p, Q4_p=met$Q4_p, Q5_p=met$Q5_p,
      Q2_rhat=sprintf("%.4f",met$Q2_rhat), Q3_rhat=sprintf("%.4f",met$Q3_rhat),
      Q4_rhat=sprintf("%.4f",met$Q4_rhat), Q5_rhat=sprintf("%.4f",met$Q5_rhat),
      Q2_ess_bulk=as.integer(round(met$Q2_ess_bulk)), Q3_ess_bulk=as.integer(round(met$Q3_ess_bulk)),
      Q4_ess_bulk=as.integer(round(met$Q4_ess_bulk)), Q5_ess_bulk=as.integer(round(met$Q5_ess_bulk)),
      Q2_ess_tail=as.integer(round(met$Q2_ess_tail)), Q3_ess_tail=as.integer(round(met$Q3_ess_tail)),
      Q4_ess_tail=as.integer(round(met$Q4_ess_tail)), Q5_ess_tail=as.integer(round(met$Q5_ess_tail))
    )
    append_rows(eqi_file, eqi_row)

    # Covariate coefficient row
    cells <- extract_covar_cells(draws, des$names, model_covars)
    covar_row <- tibble(
      ICD_Code=cancer, AAMR_Period=aamr_out, Lag=lagv, Sensitivity_Model=model_tag,
      Smoking_Rate              = cells[["Smoking_Rate"]],
      uninsured_rate_2005       = cells[["uninsured_rate_2005"]],
      physician_density_per100k = cells[["physician_density_per100k"]],
      forest_coverage           = cells[["forest_coverage"]],
      site_number_mean          = cells[["site_number_mean"]]
    )
    append_rows(covar_file, covar_row)

    message("[OK] ", sc$key, " ", model_tag)
  }
}

message("===== Done: ", prefix, " | Group: ", mod_group, " =====")
message("  EQI file:   ", eqi_file)
message("  Covar file: ", covar_file)
