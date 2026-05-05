#!/usr/bin/env Rscript
# Sensitivity analysis: EQI × covariate groups
# Input:  Data/Processed/df/df_Sensitivity.csv
# Output: Result/brms_Sensitivity/{shortname}_Sensitivity.csv  (one file per overall outcome)
#
# 13 models per scenario:
#   EQI+A:   EQI + Smoking_rate + Heavy_Drinking_rate + Physical_Activities_rate + Obesity_rate
#   EQI+B:   EQI + Forest_Coverage + Uninsured_rate + Physician_Density_per100k +
#                  AQS_Number + Homeownership_rate + Diabetes_Prevalence_rate
#   EQI+A+B: EQI + A + B
#   +{each A covariate separately}: 4 models
#   +{each B covariate separately}: 6 models

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(readr)
  library(cmdstanr)
  library(posterior)
})
utils::globalVariables(c("EQI", "State_FIPS", "AAMR_Lower", "AAMR_Upper", "RUCC",
                         "EQI_Period", "Time_Period", "Outcome", "cens"))

option_list <- list(
  make_option(c("--data"),          type="character", default="Data/Processed/df/df_Sensitivity.csv"),
  make_option(c("--output-dir"),    type="character", default="Result/brms_Sensitivity"),
  make_option(c("--outcomes"),      type="character", default=NA,
              help="Comma-separated ICD codes (default: all overall outcomes)"),
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

cores_avail <- parallel::detectCores(logical=TRUE); cores_used <- max(1,floor(cores_avail*0.8)); options(mc.cores=cores_used)
message("Detected cores: ", cores_avail, " | Using: ", cores_used)

# ---- Overall outcome shortnames ----
OVERALL_SHORTNAMES <- c(
  "K70_K76_C22"                  = "CLD",
  "J40_J47_J60_J70_J84_D86_C34" = "CRD",
  "N00_N29_C64_C65"              = "CKD",
  "I00_I99"                      = "CVD",
  "X60_X84_Y87.0"                = "Suicide",
  "G20_G30_G12.2_F01_F03"        = "NDD",
  "C00_C97"                      = "Cancer"
)

# ---- Covariate groups (original column names in df_Sensitivity) ----
GROUP_A <- c("Smoking_rate", "Heavy_Drinking_rate", "Physical_Activities_rate", "Obesity_rate")
GROUP_B <- c("Forest_Coverage", "Uninsured_rate", "Physician_Density_per100k",
             "AQS_Number", "Homeownership_rate", "Diabetes_Prevalence_rate")
ALL_COVARS <- c(GROUP_A, GROUP_B)

ALL_SENSITIVITY_MODELS <- c(
  list(list(tag="EQI+A",   covars=GROUP_A)),
  list(list(tag="EQI+B",   covars=GROUP_B)),
  list(list(tag="EQI+A+B", covars=c(GROUP_A, GROUP_B))),
  lapply(GROUP_A, function(v) list(tag=paste0("+", v), covars=v)),
  lapply(GROUP_B, function(v) list(tag=paste0("+", v), covars=v))
)

# ---- Stan model ----
stan_code <- "data {\n  int<lower=1> N;\n  int<lower=1> S;\n  array[N] int<lower=1,upper=S> state;\n  vector[N] y_lower;\n  vector[N] y_upper;\n  array[N] int<lower=0,upper=2> cens;\n  int<lower=1> K;\n  matrix[N,K] X;\n} \nparameters {\n  vector[K] beta;\n  vector[S] z_u;\n  real<lower=0> sigma;\n  real<lower=0> sigma_u;\n} \ntransformed parameters {\n  vector[S] u = sigma_u * z_u;\n} \nmodel {\n  beta ~ normal(0,5);\n  z_u ~ normal(0,1);\n  sigma ~ exponential(1);\n  sigma_u ~ exponential(1);\n  for (i in 1:N) {\n    real mu = X[i] * beta + u[state[i]];\n    if (cens[i]==0) {\n      target += normal_lpdf(y_lower[i] | mu, sigma);\n    } else {\n      real p_up = normal_cdf(y_upper[i] | mu, sigma);\n      real p_lo = normal_cdf(y_lower[i] | mu, sigma);\n      real diff = fmax(p_up - p_lo, 1e-12);\n      target += log(diff);\n    }\n  }\n}"
stan_file <- file.path(tempdir(), "interval_mixed_model.stan"); writeLines(stan_code, stan_file)
mod <- cmdstan_model(stan_file)

# ---- Load data ----
project_root <- normalizePath(".")
path <- file.path(project_root, opt$data); if(!file.exists(path)) stop("Data not found: ", path)
dt <- fread(path)

req <- c("COUNTY_FIPS","EQI_Period","Time_Period","Lag_Years","Outcome","AAMR_Lower","AAMR_Upper","RUCC","EQI")
miss <- setdiff(req, names(dt)); if(length(miss)) stop("Missing columns: ", paste(miss, collapse=","))

if(!"State_FIPS" %in% names(dt)) dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]
dt <- dt[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
dt[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]
dt <- dt[RUCC %in% 1:4 | is.na(RUCC)]

# ---- Select outcomes ----
available <- intersect(names(OVERALL_SHORTNAMES), unique(dt$Outcome))
selected <- if(is.na(opt$`outcomes`)) {
  available
} else {
  reqc <- str_trim(as.vector(str_split(opt$`outcomes`, ",", simplify=TRUE)))
  inv <- setdiff(reqc, available)
  if(length(inv)) stop("Invalid/unavailable outcomes: ", paste(inv, collapse=","))
  reqc
}
message("Outcomes: ", paste(selected, collapse=", "))

# ---- Output dir ----
out_dir <- file.path(project_root, opt$`output-dir`)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

# ---- Scenarios (same as cmdstan_main.R) ----
scenario_list <- list(
  list(key="EQI0005_AAMR2006_2010", eqi="2000-2005", aamr="2006-2010", lag=5),
  list(key="EQI0005_AAMR2011_2015", eqi="2000-2005", aamr="2011-2015", lag=10),
  list(key="EQI0005_AAMR2016_2020", eqi="2000-2005", aamr="2016-2020", lag=15),
  list(key="EQI0005_AAMR2021_2024", eqi="2000-2005", aamr="2021-2024", lag=20),
  list(key="EQI0610_AAMR2011_2015", eqi="2006-2010", aamr="2011-2015", lag=5),
  list(key="EQI0610_AAMR2016_2020", eqi="2006-2010", aamr="2016-2020", lag=10),
  list(key="EQI0610_AAMR2021_2024", eqi="2006-2010", aamr="2021-2024", lag=15)
)

# ---- Helpers ----
format_cell <- function(draws){
  if(length(draws)==0) return("")
  ci <- quantile(draws, c(0.025,0.975), na.rm=TRUE)
  sprintf("%0.2f(%0.2f,%0.2f)", mean(draws), ci[1], ci[2])
}
compute_p <- function(draws){
  if(length(draws)==0) return(NA_character_)
  pos <- sum(draws>0,na.rm=TRUE); neg <- sum(draws<0,na.rm=TRUE); n <- pos+neg
  if(n==0) return(NA_character_)
  sprintf("%.4f", 2*min((pos+0.5)/(n+1), (neg+0.5)/(n+1)))
}
append_rows <- function(path, df){
  if(!file.exists(path)) write_csv(df,path) else
    suppressWarnings(write.table(df,path,sep=",",col.names=FALSE,row.names=FALSE,append=TRUE))
}

# ---- Design matrix ----
build_design_sensitivity <- function(d, z_cols=character(0)){
  d <- d %>% mutate(EQI_factor = factor(EQI, levels=1:5))
  keep_cols <- c("EQI_factor","AAMR_Lower","AAMR_Upper","cens","State_FIPS", z_cols)
  d <- d[complete.cases(d[, ..keep_cols]), ]
  rhs <- paste(c("EQI_factor", z_cols), collapse=" + ")
  mm  <- model.matrix(as.formula(paste("~", rhs)), data=d,
                      contrasts.arg=list(EQI_factor=contr.treatment(5)))
  colnames(mm) <- make.names(colnames(mm))
  list(X=mm, names=colnames(mm), df=d)
}

# ---- Quintile extractors ----
extract_quintiles <- function(draw_df, names_vec, prefix){
  out <- list(Q1="0.00", Q2="", Q3="", Q4="", Q5="")
  if(any(grepl(paste0(prefix,"\\.L"), names_vec))) return(out)
  for(q in 2:5){
    nm <- paste0(prefix,q); idx <- match(nm,names_vec)
    out[[paste0("Q",q)]] <- if(is.na(idx)) "" else format_cell(draw_df[[paste0("beta[",idx,"]")]])
  }
  out
}
extract_quintile_metrics <- function(draw_df, names_vec, prefix, summ_df){
  out <- list(Q2_p=NA_character_,Q3_p=NA_character_,Q4_p=NA_character_,Q5_p=NA_character_,
              Q2_rhat=NA_real_,Q3_rhat=NA_real_,Q4_rhat=NA_real_,Q5_rhat=NA_real_,
              Q2_ess_bulk=NA_real_,Q3_ess_bulk=NA_real_,Q4_ess_bulk=NA_real_,Q5_ess_bulk=NA_real_,
              Q2_ess_tail=NA_real_,Q3_ess_tail=NA_real_,Q4_ess_tail=NA_real_,Q5_ess_tail=NA_real_)
  if(any(grepl(paste0(prefix,"\\.L"), names_vec))) return(out)
  for(q in 2:5){
    nm <- paste0(prefix,q); idx <- match(nm,names_vec)
    if(!is.na(idx)){
      col <- paste0("beta[",idx,"]")
      out[[paste0("Q",q,"_p")]] <- compute_p(draw_df[[col]])
      sr <- summ_df[summ_df$variable==col,,drop=FALSE]
      if(nrow(sr)){
        out[[paste0("Q",q,"_rhat")]]     <- sr$rhat
        out[[paste0("Q",q,"_ess_bulk")]] <- sr$ess_bulk
        out[[paste0("Q",q,"_ess_tail")]] <- sr$ess_tail
      }
    }
  }
  out
}

extract_covar_cell <- function(draw_df, names_vec, z_name){
  idx <- match(z_name, names_vec)
  if(is.na(idx)) return("")
  format_cell(draw_df[[paste0("beta[",idx,"]")]])
}

# ---- Main loop ----
for(outcome in selected){
  shortname <- OVERALL_SHORTNAMES[[outcome]]
  outfile   <- file.path(out_dir, paste0(shortname, "_Sensitivity.csv"))
  if(file.exists(outfile)) file.remove(outfile)

  message("\n===== ", shortname, " (", outcome, ") =====")
  out_dt <- dt[Outcome == outcome]

  # Z-score all covariates on the full outcome dataset
  z_name_map <- setNames(paste0(ALL_COVARS, "_z"), ALL_COVARS)
  for(orig_col in ALL_COVARS){
    z_col <- z_name_map[[orig_col]]
    if(!orig_col %in% names(out_dt)){
      out_dt[, (z_col) := NA_real_]; next
    }
    vals <- out_dt[[orig_col]]
    mu_ <- mean(vals, na.rm=TRUE); sd_ <- sd(vals, na.rm=TRUE)
    sd_safe <- ifelse(is.na(sd_)|sd_==0, 1, sd_)
    out_dt[, (z_col) := (out_dt[[orig_col]] - mu_) / sd_safe]
  }

  for(sc in scenario_list){
    eqi_p <- sc$eqi; aamr_p <- sc$aamr; lagv <- sc$lag
    eqi_out <- gsub("-","_",eqi_p); aamr_out <- gsub("-","_",aamr_p)
    scen_dt <- out_dt[EQI_Period == eqi_p & Time_Period == aamr_p]
    if(nrow(scen_dt) < opt$`min-n`){ message("[Skip] ",sc$key," n=",nrow(scen_dt)); next }

    for(sm in ALL_SENSITIVITY_MODELS){
      model_tag <- sm$tag
      z_cols    <- unname(z_name_map[sm$covars])

      des <- build_design_sensitivity(scen_dt, z_cols)
      if(nrow(des$df) < opt$`min-n`){ message("[Skip] ",sc$key," ",model_tag," n=",nrow(des$df)); next }

      states_s  <- sort(unique(des$df$State_FIPS)); state_idx <- match(des$df$State_FIPS, states_s)
      data_list <- list(N=nrow(des$df), S=length(states_s), state=state_idx,
                        y_lower=des$df$AAMR_Lower, y_upper=des$df$AAMR_Upper, cens=des$df$cens,
                        K=ncol(des$X), X=des$X)
      init_fun <- function() list(beta=rep(0,data_list$K), z_u=rep(0,data_list$S), sigma=50, sigma_u=10)

      fit <- try(mod$sample(data=data_list, chains=opt$chains,
                            iter_sampling=opt$iter-opt$warmup, iter_warmup=opt$warmup,
                            adapt_delta=opt$`adapt-delta`, max_treedepth=opt$`max-treedepth`,
                            parallel_chains=min(opt$chains,cores_used), refresh=0, seed=opt$seed,
                            init=rep(list(init_fun()), opt$chains)), silent=FALSE)

      if(inherits(fit,"try-error")){ message("[Fail] ",sc$key," ",model_tag); next }

      draws <- as_draws_df(fit$draws("beta"))
      colnames(draws) <- paste0("beta[",seq_len(ncol(draws)),"]")
      summ  <- posterior::summarize_draws(fit$draws("beta"))

      q   <- extract_quintiles(draws, des$names, "EQI_factor")
      met <- extract_quintile_metrics(draws, des$names, "EQI_factor", summ)

      # Covariate coefficients — one column per covariate, empty if not in this model
      all_z_names <- unname(z_name_map)
      covar_cells <- setNames(
        lapply(all_z_names, function(z) extract_covar_cell(draws, des$names, z)),
        ALL_COVARS
      )

      row <- tibble(
        ICD_Code         = outcome,
        EQI_Period       = eqi_out,
        AAMR_Period      = aamr_out,
        Lag              = lagv,
        Sensitivity_Model = model_tag,
        Q1=q$Q1, Q2=q$Q2, Q3=q$Q3, Q4=q$Q4, Q5=q$Q5,
        Q2_p=met$Q2_p, Q3_p=met$Q3_p, Q4_p=met$Q4_p, Q5_p=met$Q5_p,
        Q2_rhat    =sprintf("%.4f",met$Q2_rhat),
        Q3_rhat    =sprintf("%.4f",met$Q3_rhat),
        Q4_rhat    =sprintf("%.4f",met$Q4_rhat),
        Q5_rhat    =sprintf("%.4f",met$Q5_rhat),
        Q2_ess_bulk=as.integer(round(met$Q2_ess_bulk)),
        Q3_ess_bulk=as.integer(round(met$Q3_ess_bulk)),
        Q4_ess_bulk=as.integer(round(met$Q4_ess_bulk)),
        Q5_ess_bulk=as.integer(round(met$Q5_ess_bulk)),
        Q2_ess_tail=as.integer(round(met$Q2_ess_tail)),
        Q3_ess_tail=as.integer(round(met$Q3_ess_tail)),
        Q4_ess_tail=as.integer(round(met$Q4_ess_tail)),
        Q5_ess_tail=as.integer(round(met$Q5_ess_tail)),
        # Group A — behaviour and lifestyle
        Smoking_rate             = covar_cells[["Smoking_rate"]],
        Heavy_Drinking_rate      = covar_cells[["Heavy_Drinking_rate"]],
        Physical_Activities_rate = covar_cells[["Physical_Activities_rate"]],
        Obesity_rate             = covar_cells[["Obesity_rate"]],
        # Group B — area-level contextual
        Forest_Coverage           = covar_cells[["Forest_Coverage"]],
        Uninsured_rate            = covar_cells[["Uninsured_rate"]],
        Physician_Density_per100k = covar_cells[["Physician_Density_per100k"]],
        AQS_Number                = covar_cells[["AQS_Number"]],
        Homeownership_rate        = covar_cells[["Homeownership_rate"]],
        Diabetes_Prevalence_rate  = covar_cells[["Diabetes_Prevalence_rate"]]
      )
      append_rows(outfile, row)
      message("[OK] ", sc$key, " ", model_tag)
    }
  }
  message("===== Done: ", shortname, " → ", outfile, " =====")
}
message("All done. Output: ", out_dir)
