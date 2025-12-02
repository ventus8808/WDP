#!/usr/bin/env Rscript
# Pure cmdstanr interval-censored mixed model pipeline stratified by Typology and LandUse
# Output format mimics LMM results (ICD_Code, EQI_Period, AAMR_Period, Lag, Model, Q1..Q5)
# This script runs BOTH typology and landuse stratifications for a given disease

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
  make_option(c("--data"), type="character", default="Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Cluster_Climate_Typology_LandUse.csv", help="Input interval data with typology and landuse variables"),
  make_option(c("--output-dir"), type="character", default="Result/brms_Typology_LandUse", help="Output directory"),
  make_option(c("--cancer-types"), type="character", default="C00_C97", help="Comma separated ICD codes"),
  make_option(c("--chains"), type="integer", default=4),
  make_option(c("--iter"), type="integer", default=2000),
  make_option(c("--warmup"), type="integer", default=1000),
  make_option(c("--adapt-delta"), type="double", default=0.95),
  make_option(c("--max-treedepth"), type="integer", default=12),
  make_option(c("--min-n"), type="integer", default=50),
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
path <- file.path(project_root,opt$data); if(!file.exists(path)) stop("Data not found: ", path)
dt <- fread(path)

# Check required columns
req <- c("COUNTY_FIPS","EQI_Period","Time_Period","Lag_Years","Cancer_Type","AAMR_Lower","AAMR_Upper","Smoking_Rate","RUCC","EQI","EQI_Air","EQI_Water","EQI_Land","EQI_Built","EQI_Social","econdep","landuse_cluster")
miss <- setdiff(req, names(dt)); if(length(miss)) stop("Missing cols: ", paste(miss,collapse=","))

if(!"State_FIPS" %in% names(dt)) dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS),1,2)]

# interval censoring code
dt <- dt[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
dt[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]

# Define stratification configurations
# Typology: 1=Farming, 2=Mining, 3=Manufacturing, 4=Government, 5=Services, 6=Nonspecialized
# LandUse: 0=Natural, 1=Water-Sensitive, 2=Agricultural, 3=Urban
strat_configs <- list(
  list(name = "Typology", var = "econdep", values = 1:6,
       labels = c("Farming", "Mining", "Manufacturing", "Government", "Services", "Nonspecialized")),
  list(name = "LandUse", var = "landuse_cluster", values = 0:3,
       labels = c("Natural", "Water_Sensitive", "Agricultural", "Urban"))
)

scenario_list <- list(
  list(key="EQI0005_AAMR2006_2010", eqi="2000-2005", aamr="2006-2010", lag=5),
  list(key="EQI0005_AAMR2011_2015", eqi="2000-2005", aamr="2011-2015", lag=10),
  list(key="EQI0005_AAMR2016_2020", eqi="2000-2005", aamr="2016-2020", lag=15),
  list(key="EQI0610_AAMR2011_2015", eqi="2006-2010", aamr="2011-2015", lag=5),
  list(key="EQI0610_AAMR2016_2020", eqi="2006-2010", aamr="2016-2020", lag=10)
)

all_cancers <- sort(unique(dt$Cancer_Type))
selected <- if (is.na(opt$`cancer-types`)) all_cancers else { reqc <- str_split(opt$`cancer-types`,",",simplify=TRUE) |> as.vector() |> str_trim(); inv <- setdiff(reqc,all_cancers); if(length(inv)) stop("Invalid cancer types: ", paste(inv,collapse=",")); reqc }
# Force to single cancer for stratification analysis
if (length(selected) != 1) stop("For stratification analysis, exactly one cancer type must be specified")
message("Cancer type to analyze: ", paste(selected,collapse=","))
out_dir <- file.path(project_root,opt$`output-dir`); if(!dir.exists(out_dir)) dir.create(out_dir,recursive=TRUE)

sig_mark <- function(p){ if(is.na(p)) return(""); if(p<0.001) return("***"); if(p<0.01) return("**"); if(p<0.05) return("*"); "" }
format_cell <- function(draws){ if(length(draws)==0) return(""); ci <- quantile(draws,c(0.025,0.975),na.rm=TRUE); p <- 2*min(mean(draws>0), mean(draws<0)); sprintf("%0.2f(%0.2f,%0.2f)%s", mean(draws), ci[1], ci[2], sig_mark(p)) }
append_rows <- function(path, df){ if(!file.exists(path)) write_csv(df,path) else suppressWarnings(write.table(df,path,sep=",",col.names=FALSE,row.names=FALSE,append=TRUE)) }

# Posterior diagnostics helpers
compute_p <- function(draws){
  if(length(draws)==0) return(NA_character_)
  pos <- sum(draws > 0, na.rm = TRUE)
  neg <- sum(draws < 0, na.rm = TRUE)
  n <- pos + neg
  if(n == 0) return(NA_character_)
  p_pos <- (pos + 0.5) / (n + 1)
  p_neg <- (neg + 0.5) / (n + 1)
  p <- 2 * min(p_pos, p_neg)
  sprintf("%.4f", p)
}

extract_quintile_metrics <- function(draw_df, names_vec, prefix, summ_df){
  out <- list(
    Q2_p=NA_character_, Q3_p=NA_character_, Q4_p=NA_character_, Q5_p=NA_character_,
    Q2_rhat=NA_real_, Q3_rhat=NA_real_, Q4_rhat=NA_real_, Q5_rhat=NA_real_,
    Q2_ess_bulk=NA_real_, Q3_ess_bulk=NA_real_, Q4_ess_bulk=NA_real_, Q5_ess_bulk=NA_real_,
    Q2_ess_tail=NA_real_, Q3_ess_tail=NA_real_, Q4_ess_tail=NA_real_, Q5_ess_tail=NA_real_
  )
  if(any(grepl(paste0(prefix,"\\.L"), names_vec))){ return(out) }
  for(q in 2:5){
    nm <- paste0(prefix, q)
    idx <- match(nm, names_vec)
    if(!is.na(idx)){
      col <- paste0("beta[", idx, "]")
      draws_col <- draw_df[[col]]
      out[[paste0("Q",q,"_p")]] <- compute_p(draws_col)
      sr <- summ_df[summ_df$variable == col, , drop=FALSE]
      if(nrow(sr)){
        out[[paste0("Q",q,"_rhat")]] <- sr$rhat
        out[[paste0("Q",q,"_ess_bulk")]] <- sr$ess_bulk
        out[[paste0("Q",q,"_ess_tail")]] <- sr$ess_tail
      }
    }
  }
  out
}

build_design_overall <- function(d){ # Intercept + Smoking + EQI Q2..Q5
  d <- d %>% mutate(EQI_factor = factor(EQI, levels=1:5))
  d <- d[complete.cases(d[,c("Smoking_Rate","EQI_factor","AAMR_Lower","AAMR_Upper","cens","State_FIPS")]),]
  mm <- model.matrix(~ Smoking_Rate + EQI_factor, d, contrasts.arg = list(EQI_factor=contr.treatment(5)))
  colnames(mm) <- make.names(colnames(mm))
  list(X = mm, names = colnames(mm), df = d)
}

build_design_multi <- function(d){
  d <- d %>% mutate(
    EQI_Air_factor = factor(EQI_Air, levels=1:5),
    EQI_Water_factor = factor(EQI_Water, levels=1:5),
    EQI_Land_factor = factor(EQI_Land, levels=1:5),
    EQI_Built_factor = factor(EQI_Built, levels=1:5),
    EQI_Social_factor = factor(EQI_Social, levels=1:5)
  )
  d <- d[complete.cases(d[,c("Smoking_Rate","EQI_Air_factor","EQI_Water_factor","EQI_Land_factor","EQI_Built_factor","EQI_Social_factor","AAMR_Lower","AAMR_Upper","cens","State_FIPS")]),]
  form <- as.formula("~ Smoking_Rate + EQI_Air_factor + EQI_Water_factor + EQI_Land_factor + EQI_Built_factor + EQI_Social_factor")
  mm <- model.matrix(form, d,
                     contrasts.arg = list(
                       EQI_Air_factor=contr.treatment(5), EQI_Water_factor=contr.treatment(5),
                       EQI_Land_factor=contr.treatment(5), EQI_Built_factor=contr.treatment(5),
                       EQI_Social_factor=contr.treatment(5)))
  colnames(mm) <- make.names(colnames(mm))
  list(X=mm, names=colnames(mm), df = d)
}

extract_quintiles <- function(draw_df, names_vec, prefix){
  out <- list(Q1="0.00", Q2="", Q3="", Q4="", Q5="")
  if(any(grepl(paste0(prefix,"\\.L"), names_vec))){ return(out) }
  for(q in 2:5){ nm <- paste0(prefix, q); idx <- match(nm, names_vec); out[[paste0("Q",q)]] <- if(is.na(idx)) "" else format_cell(draw_df[[paste0("beta[",idx,"]")]]) }
  out
}

# Main analysis loop
for(cancer in selected){
  message("===== Disease: ", cancer, " =====")

  # Loop through both stratification types
  for(strat_config in strat_configs){
    strat_name <- strat_config$name
    strat_var <- strat_config$var
    strat_values <- strat_config$values
    strat_labels <- strat_config$labels

    message("===== Stratification: ", strat_name, " (", strat_var, ") =====")
    outfile <- file.path(out_dir, paste0(cancer, "_", strat_name, ".csv"))

    # Filter data to include only valid stratification values
    strat_dt <- dt[get(strat_var) %in% strat_values]
    message("Total observations with valid ", strat_name, ": ", nrow(strat_dt))

    # Loop through each stratum value
    for(i in seq_along(strat_values)){
      strat_value <- strat_values[i]
      strat_label <- strat_labels[i]
      strat_title <- paste0(strat_name, "_", strat_label)

      subset_dt <- strat_dt[get(strat_var) == strat_value]
      message("Analyzing ", strat_title, " (n=", nrow(subset_dt), ")")

      # Loop through scenarios
      for(sc in scenario_list){
        scen_key <- sc$key; eqi_p <- sc$eqi; aamr_p <- sc$aamr; lagv <- sc$lag
        scen_dt <- subset_dt[EQI_Period==eqi_p & Time_Period==aamr_p & Cancer_Type==cancer]
        if(nrow(scen_dt) < opt$`min-n`){ message("[Skip] Scenario ", scen_key, " n=", nrow(scen_dt)); next }
        eqi_out <- gsub('-', '_', eqi_p); aamr_out <- gsub('-', '_', aamr_p)

        # Overall model design
        des_overall <- build_design_overall(scen_dt)
        states_o <- sort(unique(des_overall$df$State_FIPS)); state_index_o <- match(des_overall$df$State_FIPS, states_o)
        data_list <- list(
          N = nrow(des_overall$df), S = length(states_o), state = state_index_o,
          y_lower = des_overall$df$AAMR_Lower, y_upper = des_overall$df$AAMR_Upper, cens = des_overall$df$cens,
          K = ncol(des_overall$X), X = des_overall$X
        )
        init_fun <- function() list(beta=rep(0, data_list$K), z_u=rep(0, data_list$S), sigma=50, sigma_u=10)
        fit_overall <- try(mod$sample(data=data_list, chains=opt$chains, iter_sampling=opt$iter-opt$warmup, iter_warmup=opt$warmup,
                  adapt_delta=opt$`adapt-delta`, max_treedepth=opt$`max-treedepth`, parallel_chains=min(opt$chains, cores_used), refresh=0, seed=opt$seed,
                  init=rep(list(init_fun()), opt$chains)), silent=TRUE)
        if(inherits(fit_overall,"try-error")){ message("[Fail] Overall EQI model ", scen_key); } else {
          draws <- as_draws_df(fit_overall$draws("beta"))
          colnames(draws) <- paste0("beta[",seq_len(ncol(draws)),"]")
          q_over <- extract_quintiles(draws, des_overall$names, "EQI_factor")
          summ_over <- posterior::summarize_draws(fit_overall$draws("beta"))
          met_over <- extract_quintile_metrics(draws, des_overall$names, "EQI_factor", summ_over)
          row_over <- tibble(
            ICD_Code=cancer, EQI_Period=eqi_out, AAMR_Period=aamr_out, Lag=lagv, Model=paste0(strat_title,"_EQI"),
            Q1=q_over$Q1, Q2=q_over$Q2, Q3=q_over$Q3, Q4=q_over$Q4, Q5=q_over$Q5,
            Q2_p=met_over$Q2_p, Q3_p=met_over$Q3_p, Q4_p=met_over$Q4_p, Q5_p=met_over$Q5_p,
            Q2_rhat=sprintf("%.4f", met_over$Q2_rhat), Q3_rhat=sprintf("%.4f", met_over$Q3_rhat), Q4_rhat=sprintf("%.4f", met_over$Q4_rhat), Q5_rhat=sprintf("%.4f", met_over$Q5_rhat),
            Q2_ess_bulk=as.integer(round(met_over$Q2_ess_bulk)), Q3_ess_bulk=as.integer(round(met_over$Q3_ess_bulk)), Q4_ess_bulk=as.integer(round(met_over$Q4_ess_bulk)), Q5_ess_bulk=as.integer(round(met_over$Q5_ess_bulk)),
            Q2_ess_tail=as.integer(round(met_over$Q2_ess_tail)), Q3_ess_tail=as.integer(round(met_over$Q3_ess_tail)), Q4_ess_tail=as.integer(round(met_over$Q4_ess_tail)), Q5_ess_tail=as.integer(round(met_over$Q5_ess_tail))
          )
          append_rows(outfile,row_over); message("[OK] ", scen_key, " Overall")
        }

        # Multi-domain design
        des_multi <- build_design_multi(scen_dt)
        states_m <- sort(unique(des_multi$df$State_FIPS)); state_index_m <- match(des_multi$df$State_FIPS, states_m)
        data_list2 <- list(N=nrow(des_multi$df), S=length(states_m), state=state_index_m, y_lower=des_multi$df$AAMR_Lower, y_upper=des_multi$df$AAMR_Upper, cens=des_multi$df$cens, K=ncol(des_multi$X), X=des_multi$X)
        init_fun2 <- function() list(beta=rep(0, data_list2$K), z_u=rep(0, data_list2$S), sigma=50, sigma_u=10)
        fit_multi <- try(mod$sample(data=data_list2, chains=opt$chains, iter_sampling=opt$iter-opt$warmup, iter_warmup=opt$warmup,
                    adapt_delta=opt$`adapt-delta`, max_treedepth=opt$`max-treedepth`, parallel_chains=min(opt$chains, cores_used), refresh=0, seed=opt$seed,
                    init=rep(list(init_fun2()), opt$chains)), silent=TRUE)
        if(inherits(fit_multi,"try-error")){ message("[Fail] Multi-domain model ", scen_key) } else {
          draws_m <- as_draws_df(fit_multi$draws("beta")); colnames(draws_m) <- paste0("beta[",seq_len(ncol(draws_m)),"]")
          domain_prefix <- c("EQI_Air_factor","EQI_Water_factor","EQI_Land_factor","EQI_Built_factor","EQI_Social_factor")
          summ_m <- posterior::summarize_draws(fit_multi$draws("beta"))
          for(dom in domain_prefix){
            qd <- extract_quintiles(draws_m, des_multi$names, dom)
            md <- extract_quintile_metrics(draws_m, des_multi$names, dom, summ_m)
            row_dom <- tibble(
              ICD_Code=cancer, EQI_Period=eqi_out, AAMR_Period=aamr_out, Lag=lagv, Model=paste0(strat_title, "_", sub("_factor","",dom)),
              Q1=qd$Q1, Q2=qd$Q2, Q3=qd$Q3, Q4=qd$Q4, Q5=qd$Q5,
              Q2_p=md$Q2_p, Q3_p=md$Q3_p, Q4_p=md$Q4_p, Q5_p=md$Q5_p,
              Q2_rhat=sprintf("%.4f", md$Q2_rhat), Q3_rhat=sprintf("%.4f", md$Q3_rhat), Q4_rhat=sprintf("%.4f", md$Q4_rhat), Q5_rhat=sprintf("%.4f", md$Q5_rhat),
              Q2_ess_bulk=as.integer(round(md$Q2_ess_bulk)), Q3_ess_bulk=as.integer(round(md$Q3_ess_bulk)), Q4_ess_bulk=as.integer(round(md$Q4_ess_bulk)), Q5_ess_bulk=as.integer(round(md$Q5_ess_bulk)),
              Q2_ess_tail=as.integer(round(md$Q2_ess_tail)), Q3_ess_tail=as.integer(round(md$Q3_ess_tail)), Q4_ess_tail=as.integer(round(md$Q4_ess_tail)), Q5_ess_tail=as.integer(round(md$Q5_ess_tail))
            )
            append_rows(outfile,row_dom)
          }
          message("[OK] ", scen_key, " Multi-domain")
        }
      }
    }
    message("===== Completed: ", strat_name, " for ", cancer, " =====")
  }
  message("===== Completed all stratifications for: ", cancer, " =====")
}
message("All requested analyses complete. Output directory: ", out_dir)
