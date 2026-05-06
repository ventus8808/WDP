#!/usr/bin/env Rscript
# Sensitivity combination analysis: all 2^6 subsets of 6 control covariates added to EQI exposure.
# Models: +0 covariate (EQI only) = 1; +1 = C(6,1)=6; +2 = C(6,2)=15; +3 = C(6,3)=20;
#         +4 = C(6,4)=15; +5 = C(6,5)=6; +6 = 1 → total 64 models per scenario/layer.
# EQI 2000-2005 scenarios only (2006-2010 EQI removed).
# Output format: ICD_Code, EQI_Period, AAMR_Period, Lag, Model, Q1..Q5 + diagnostics.

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
utils::globalVariables(c('EQI','Smoking_rate','Physical_Activities_rate','Obesity_rate',
                         'Uninsured_rate','Physician_Density_per100k','Diabetes_Prevalence_rate','State_FIPS'))

option_list <- list(
  make_option(c("--data"), type="character", default="Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Cluster_Climate.csv", help="Input interval data"),
  make_option(c("--output-dir"), type="character", default="Result/brms_Sensitivity_Combination", help="Output directory"),
  make_option(c("--cancer-types"), type="character", default=NA, help="Comma separated ICD codes"),
  make_option(c("--chains"), type="integer", default=4),
  make_option(c("--iter"), type="integer", default=2000),
  make_option(c("--warmup"), type="integer", default=1000),
  make_option(c("--adapt-delta"), type="double", default=0.95),
  make_option(c("--max-treedepth"), type="integer", default=12),
  make_option(c("--min-n"), type="integer", default=50),
  make_option(c("--min-n-rucc"), type="integer", default=30),
  make_option(c("--seed"), type="integer", default=1234),
  make_option(c("--test"), action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))
if (opt$test) { opt$iter <- min(opt$iter, 800); opt$warmup <- min(opt$warmup, 300); message("[TEST MODE] iter=",opt$iter," warmup=",opt$warmup) }
set.seed(opt$seed)

cores_avail <- parallel::detectCores(logical=TRUE); cores_used <- max(1,floor(cores_avail*0.8)); options(mc.cores=cores_used)
message("Detected cores: ", cores_avail, " | Using: ", cores_used)

# Stan model (generic design matrix X, group random intercept u)
stan_code <- "data {\n  int<lower=1> N;\n  int<lower=1> S;\n  array[N] int<lower=1,upper=S> state;\n  vector[N] y_lower;\n  vector[N] y_upper;\n  array[N] int<lower=0,upper=2> cens;\n  int<lower=1> K;\n  matrix[N,K] X;\n} \nparameters {\n  vector[K] beta;\n  vector[S] z_u;\n  real<lower=0> sigma;\n  real<lower=0> sigma_u;\n} \ntransformed parameters {\n  vector[S] u = sigma_u * z_u;\n} \nmodel {\n  beta ~ normal(0,5);\n  z_u ~ normal(0,1);\n  sigma ~ exponential(1);\n  sigma_u ~ exponential(1);\n  for (i in 1:N) {\n    real mu = X[i] * beta + u[state[i]];\n    if (cens[i]==0) {\n      target += normal_lpdf(y_lower[i] | mu, sigma);\n    } else {\n      real p_up = normal_cdf(y_upper[i] | mu, sigma);\n      real p_lo = normal_cdf(y_lower[i] | mu, sigma);\n      real diff = fmax(p_up - p_lo, 1e-12);\n      target += log(diff);\n    }\n  }\n}";
stan_file <- file.path(tempdir(), "interval_mixed_model.stan"); writeLines(stan_code, stan_file)
mod <- cmdstan_model(stan_file)

# Load data
project_root <- normalizePath(".")
path <- file.path(project_root, opt$data); if(!file.exists(path)) stop("Data not found: ", path)
dt <- fread(path)

req <- c("COUNTY_FIPS","EQI_Period","Time_Period","Lag_Years","Cancer_Type","AAMR_Lower","AAMR_Upper","RUCC","EQI",
         "Smoking_rate","Physical_Activities_rate","Obesity_rate","Uninsured_rate","Physician_Density_per100k","Diabetes_Prevalence_rate")
miss <- setdiff(req, names(dt)); if(length(miss)) stop("Missing cols: ", paste(miss,collapse=","))

if(!"State_FIPS" %in% names(dt)) dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS),1,2)]

dt <- dt[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
dt[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]
dt <- dt[RUCC %in% 1:4 | is.na(RUCC)]

# EQI 2000-2005 scenarios only (2006-2010 EQI removed)
scenario_list <- list(
  list(key="EQI0005_AAMR2006_2010", eqi="2000-2005", aamr="2006-2010", lag=5),
  list(key="EQI0005_AAMR2011_2015", eqi="2000-2005", aamr="2011-2015", lag=10),
  list(key="EQI0005_AAMR2016_2020", eqi="2000-2005", aamr="2016-2020", lag=15)
)

# 6 control covariates and their abbreviations for model labels
all_covariates <- c("Smoking_rate","Physical_Activities_rate","Obesity_rate",
                    "Uninsured_rate","Physician_Density_per100k","Diabetes_Prevalence_rate")
covar_abbrev <- c(Smoking_rate="SM", Physical_Activities_rate="PA", Obesity_rate="OB",
                  Uninsured_rate="UN", Physician_Density_per100k="PD", Diabetes_Prevalence_rate="DB")

# All 2^6 = 64 covariate subsets (k = 0 to 6)
combo_list <- list(character(0))  # k=0: EQI only
for (k in seq_along(all_covariates)) {
  combo_list <- c(combo_list, combn(all_covariates, k, simplify=FALSE))
}
message("Total covariate combinations: ", length(combo_list),
        " (k=0:1, k=1:6, k=2:15, k=3:20, k=4:15, k=5:6, k=6:1)")

all_cancers <- sort(unique(dt$Cancer_Type))
selected <- if (is.na(opt$`cancer-types`)) all_cancers else {
  reqc <- str_split(opt$`cancer-types`,",",simplify=TRUE) |> as.vector() |> str_trim()
  inv <- setdiff(reqc, all_cancers); if(length(inv)) stop("Invalid cancer types: ", paste(inv,collapse=","))
  reqc
}
message("Cancer types to analyze: ", paste(selected,collapse=","))

out_dir <- file.path(project_root, opt$`output-dir`); if(!dir.exists(out_dir)) dir.create(out_dir,recursive=TRUE)

sig_mark <- function(p){ if(is.na(p)) return(""); if(p<0.001) return("***"); if(p<0.01) return("**"); if(p<0.05) return("*"); "" }
format_cell <- function(draws){ if(length(draws)==0) return(""); ci <- quantile(draws,c(0.025,0.975),na.rm=TRUE); p <- 2*min(mean(draws>0), mean(draws<0)); sprintf("%0.2f(%0.2f,%0.2f)%s", mean(draws), ci[1], ci[2], sig_mark(p)) }
append_rows <- function(path, df){ if(!file.exists(path)) write_csv(df,path) else suppressWarnings(write.table(df,path,sep=",",col.names=FALSE,row.names=FALSE,append=TRUE)) }

# p_posterior: two-sided posterior tail-area probability with Jeffreys correction.
compute_p <- function(draws){
  if(length(draws)==0) return(NA_character_)
  pos <- sum(draws > 0, na.rm=TRUE); neg <- sum(draws < 0, na.rm=TRUE); n <- pos + neg
  if(n == 0) return(NA_character_)
  p_pos <- (pos + 0.5) / (n + 1); p_neg <- (neg + 0.5) / (n + 1)
  sprintf("%.4f", 2 * min(p_pos, p_neg))
}

extract_quintile_metrics <- function(draw_df, names_vec, prefix, summ_df){
  out <- list(Q2_p=NA_real_, Q3_p=NA_real_, Q4_p=NA_real_, Q5_p=NA_real_,
              Q2_rhat=NA_real_, Q3_rhat=NA_real_, Q4_rhat=NA_real_, Q5_rhat=NA_real_,
              Q2_ess_bulk=NA_real_, Q3_ess_bulk=NA_real_, Q4_ess_bulk=NA_real_, Q5_ess_bulk=NA_real_,
              Q2_ess_tail=NA_real_, Q3_ess_tail=NA_real_, Q4_ess_tail=NA_real_, Q5_ess_tail=NA_real_)
  if(any(grepl(paste0(prefix,"\\.L"), names_vec))){ return(out) }
  for(q in 2:5){
    nm <- paste0(prefix, q); idx <- match(nm, names_vec)
    if(!is.na(idx)){
      col <- paste0("beta[", idx, "]"); draws_col <- draw_df[[col]]
      out[[paste0("Q",q,"_p")]] <- compute_p(draws_col)
      sr <- summ_df[summ_df$variable == col, , drop=FALSE]
      if(nrow(sr)){ out[[paste0("Q",q,"_rhat")]] <- sr$rhat; out[[paste0("Q",q,"_ess_bulk")]] <- sr$ess_bulk; out[[paste0("Q",q,"_ess_tail")]] <- sr$ess_tail }
    }
  }
  out
}

extract_quintiles <- function(draw_df, names_vec, prefix){
  out <- list(Q1="0.00", Q2="", Q3="", Q4="", Q5="")
  if(any(grepl(paste0(prefix,"\\.L"), names_vec))){ return(out) }
  for(q in 2:5){ nm <- paste0(prefix, q); idx <- match(nm, names_vec); out[[paste0("Q",q)]] <- if(is.na(idx)) "" else format_cell(draw_df[[paste0("beta[",idx,"]")]]) }
  out
}

# Build design matrix: intercept + optional covariates + EQI_factor (treatment contrasts, Q1 as ref).
build_design_combo <- function(d, covariates){
  d <- d %>% mutate(EQI_factor = factor(EQI, levels=1:5))
  needed <- c("EQI_factor","AAMR_Lower","AAMR_Upper","cens","State_FIPS", covariates)
  d <- d[complete.cases(d[, needed, drop=FALSE]), ]
  form <- if (length(covariates) == 0) ~ EQI_factor else
    as.formula(paste("~", paste(c(covariates, "EQI_factor"), collapse=" + ")))
  mm <- model.matrix(form, d, contrasts.arg=list(EQI_factor=contr.treatment(5)))
  colnames(mm) <- make.names(colnames(mm))
  list(X=mm, names=colnames(mm), df=d)
}

for(cancer in selected){
  message("===== Disease: ", cancer, " =====")
  outfile <- file.path(out_dir, paste0(cancer, "_main.csv"))
  for(sc in scenario_list){
    scen_key <- sc$key; eqi_p <- sc$eqi; aamr_p <- sc$aamr; lagv <- sc$lag
    scen_dt <- dt[EQI_Period==eqi_p & Time_Period==aamr_p & Cancer_Type==cancer]
    if(nrow(scen_dt) < opt$`min-n`){ message("[Skip] Scenario ", scen_key, " overall n=", nrow(scen_dt)); next }
    eqi_out <- gsub('-','_',eqi_p); aamr_out <- gsub('-','_',aamr_p)
    layers <- c("Overall", paste0("RUCC",1:4))
    for(lay in layers){
      if(lay=="Overall"){ layer_dt <- scen_dt; min_req <- opt$`min-n`; layer_tag <- "" } else { rv <- as.integer(sub("RUCC","",lay)); layer_dt <- scen_dt[RUCC==rv]; min_req <- opt$`min-n-rucc`; layer_tag <- paste0("RUCC",rv,"_") }
      if(nrow(layer_dt) < min_req){ message("[Skip] ", scen_key, " ", lay, " n=", nrow(layer_dt)); next }

      for(combo in combo_list){
        model_label <- if(length(combo)==0) "EQI" else paste(c("EQI", covar_abbrev[combo]), collapse="+")
        full_model_tag <- paste0(layer_tag, model_label)

        des <- build_design_combo(layer_dt, combo)
        if(nrow(des$df) < min_req){ message("[Skip] ", scen_key, " ", lay, " ", model_label, " n=", nrow(des$df)); next }

        states_c <- sort(unique(des$df$State_FIPS)); state_index_c <- match(des$df$State_FIPS, states_c)
        data_list <- list(N=nrow(des$df), S=length(states_c), state=state_index_c,
                          y_lower=des$df$AAMR_Lower, y_upper=des$df$AAMR_Upper, cens=des$df$cens,
                          K=ncol(des$X), X=des$X)
        init_fn <- function() list(beta=rep(0, data_list$K), z_u=rep(0, data_list$S), sigma=50, sigma_u=10)

        fit <- try(mod$sample(data=data_list, chains=opt$chains, iter_sampling=opt$iter-opt$warmup, iter_warmup=opt$warmup,
                              adapt_delta=opt$`adapt-delta`, max_treedepth=opt$`max-treedepth`,
                              parallel_chains=min(opt$chains, cores_used), refresh=0, seed=opt$seed,
                              init=rep(list(init_fn()), opt$chains)), silent=TRUE)

        if(inherits(fit,"try-error")){ message("[Fail] ", scen_key, " ", lay, " ", model_label); next }

        draws <- as_draws_df(fit$draws("beta")); colnames(draws) <- paste0("beta[",seq_len(ncol(draws)),"]")
        q_vals <- extract_quintiles(draws, des$names, "EQI_factor")
        summ <- posterior::summarize_draws(fit$draws("beta"))
        met <- extract_quintile_metrics(draws, des$names, "EQI_factor", summ)

        row <- tibble(
          ICD_Code=cancer, EQI_Period=eqi_out, AAMR_Period=aamr_out, Lag=lagv, Model=full_model_tag,
          Q1=q_vals$Q1, Q2=q_vals$Q2, Q3=q_vals$Q3, Q4=q_vals$Q4, Q5=q_vals$Q5,
          Q2_p=met$Q2_p, Q3_p=met$Q3_p, Q4_p=met$Q4_p, Q5_p=met$Q5_p,
          Q2_rhat=sprintf("%.4f", met$Q2_rhat), Q3_rhat=sprintf("%.4f", met$Q3_rhat),
          Q4_rhat=sprintf("%.4f", met$Q4_rhat), Q5_rhat=sprintf("%.4f", met$Q5_rhat),
          Q2_ess_bulk=as.integer(round(met$Q2_ess_bulk)), Q3_ess_bulk=as.integer(round(met$Q3_ess_bulk)),
          Q4_ess_bulk=as.integer(round(met$Q4_ess_bulk)), Q5_ess_bulk=as.integer(round(met$Q5_ess_bulk)),
          Q2_ess_tail=as.integer(round(met$Q2_ess_tail)), Q3_ess_tail=as.integer(round(met$Q3_ess_tail)),
          Q4_ess_tail=as.integer(round(met$Q4_ess_tail)), Q5_ess_tail=as.integer(round(met$Q5_ess_tail))
        )
        append_rows(outfile, row)
        message("[OK] ", scen_key, " ", lay, " ", model_label)
      }
    }
  }
  message("===== Completed: ", cancer, " =====")
}
message("All requested analyses complete. Output directory: ", out_dir)
