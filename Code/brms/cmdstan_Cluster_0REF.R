#!/usr/bin/env Rscript
# Modified cmdstanr interval-censored mixed model pipeline using Cluster 0 as reference
# Output format: ICD_Code, Lag, K_Value, Model_Type, Intercept_Cluster0_Baseline, MRD_ClusterX_vs_Cluster0 (dynamic), Control_Smoking_Rate, N_Counties, Rhat_max
# Pools data across clusters, treats Cluster 0 as reference, calculates MRD as cluster betas

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
utils::globalVariables(c('Smoking_Rate','State_FIPS','Cluster'))

option_list <- list(
  make_option(c("--data"), type="character", default="Data/Processed/df_EQI_AAMR/EQI_AAMR_Interval_Clustered.csv", help="Input interval data with cluster IDs"),
  make_option(c("--output-dir"), type="character", default="Result/brms_cluster_0asREF", help="Output directory"),
  make_option(c("--cancer-types"), type="character", default="C00_C97", help="Comma separated ICD codes"),
  make_option(c("--k"), type="character", default="3,4,5,6", help="Comma separated k values to analyze (default: 3,4,5,6)"),
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

# Parse cluster IDs
k_values <- str_split(opt$k,",",simplify=TRUE) |> as.vector() |> str_trim() |> as.integer()
message("K values to analyze: ", paste(k_values,collapse=","))

req <- c("COUNTY_FIPS","EQI_Period","Time_Period","Lag_Years","Cancer_Type","AAMR_lower","AAMR_upper","Smoking_Rate")
miss <- setdiff(req, names(dt)); if(length(miss)) stop("Missing cols: ", paste(miss,collapse=","))

if(!"State_FIPS" %in% names(dt)) dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS),1,2)]

# interval censoring code
dt <- dt[!is.na(AAMR_lower) & !is.na(AAMR_upper)]
dt[, cens := ifelse(AAMR_lower == AAMR_upper, 0, 2)]

scenario_list <- list(
  list(key="EQI0005_AAMR2006_2010", eqi="2000-2005", aamr="2006-2010", lag=5),
  list(key="EQI0005_AAMR2011_2015", eqi="2000-2005", aamr="2011-2015", lag=10),
  list(key="EQI0610_AAMR2011_2015", eqi="2006-2010", aamr="2011-2015", lag=5),
  list(key="EQI0610_AAMR2016_2020", eqi="2006-2010", aamr="2016-2020", lag=10)
)

all_cancers <- sort(unique(dt$Cancer_Type))
selected <- if (is.na(opt$`cancer-types`)) all_cancers else { reqc <- str_split(opt$`cancer-types`,",",simplify=TRUE) |> as.vector() |> str_trim(); inv <- setdiff(reqc,all_cancers); if(length(inv)) stop("Invalid cancer types: ", paste(inv,collapse=",")); reqc }
message("Cancer types to analyze: ", paste(selected,collapse=","))
out_dir <- file.path(project_root,opt$`output-dir`); if(!dir.exists(out_dir)) dir.create(out_dir,recursive=TRUE)

sig_mark <- function(p){ if(is.na(p)) return(""); if(p<0.001) return("***"); if(p<0.01) return("**"); if(p<0.05) return("*"); "" }
format_cell <- function(draws){ if(length(draws)==0) return(""); ci <- quantile(draws,c(0.025,0.975),na.rm=TRUE); p <- 2*min(mean(draws>0), mean(draws<0)); sprintf("%0.2f(%0.2f,%0.2f)%s", mean(draws), ci[1], ci[2], sig_mark(p)) }
append_rows <- function(path, df){ if(!file.exists(path)) write_csv(df,path) else suppressWarnings(write.table(df,path,sep=",",col.names=FALSE,row.names=FALSE,append=TRUE)) }

build_design_cluster_main_effect <- function(d, k_val){ # Intercept + Smoking + Cluster_factor (0 as ref)
  cluster_col <- paste0("cluster_", k_val)
  d <- d %>% mutate(Cluster_factor = factor(get(cluster_col), levels=0:(k_val-1)))
  d <- d[complete.cases(d[,c("Smoking_Rate","Cluster_factor","AAMR_lower","AAMR_upper","cens","State_FIPS")]),]
  mm <- model.matrix(~ Smoking_Rate + Cluster_factor, d, contrasts.arg = list(Cluster_factor=contr.treatment(levels(d$Cluster_factor), base=1)))
  colnames(mm) <- make.names(colnames(mm))
  list(X = mm, names = colnames(mm), df = d)
}

extract_cluster_main_effect_results <- function(draws, names_vec, cancer, eqi_out, aamr_out, lagv, k_val, model_type, design_df, fit){
  # Extract intercept (Cluster 0 baseline)
  intercept_idx <- 1  # beta[1] is intercept
  intercept <- format_cell(draws[[paste0("beta[",intercept_idx,"]")]])

  # Extract MRD for each Cluster 1 to (k-1)
  mrd_list <- setNames(as.list(rep("", k_val-1)), paste0("MRD_Cluster",1:(k_val-1)))
  for(i in 1:(k_val-1)){
    cluster_nm <- paste0("Cluster_factor", levels(des_cluster$df$Cluster_factor)[i+1])
    cluster_idx <- match(cluster_nm, names_vec)
    if(!is.na(cluster_idx)){
      mrd_list[[paste0("MRD_Cluster",i)]] <- format_cell(draws[[paste0("beta[",cluster_idx,"]")]])
    }
  }

  # Diagnostics
  n_counties <- nrow(design_df)
  rhat_max <- max(fit$summary()$rhat, na.rm=TRUE)

  # Build tibble row
  row <- tibble(
    ICD_Code = cancer,
    EQI_Period = eqi_out,
    AAMR_Period = aamr_out,
    Lag = lagv,
    K_Value = k_val,
    Model_Type = model_type,
    Intercept = intercept
  )

  # Dynamically add MRD columns
  row <- bind_cols(row, as_tibble(mrd_list))

  # Add remaining fixed columns
  row <- row %>%
    mutate(
      N_Counties = n_counties,
      Rhat_max = sprintf("%.3f", rhat_max)
    )

  row
}

for(cancer in selected){
  message("===== Disease: ", cancer, " =====")

  for(k_val in k_values){
    cluster_ids <- 0:(k_val - 1)
    message("Processing k=", k_val, " with clusters: ", paste(cluster_ids,collapse=","))

    # Select the appropriate cluster column based on k
    cluster_col <- paste0("cluster_", k_val)
    if(!cluster_col %in% names(dt)) stop("Cluster column '", cluster_col, "' not found in data. Available columns: ", paste(names(dt), collapse=", "))
    dt_k <- dt[, Cluster := get(cluster_col)]

    # Single output file per cancer type and k value
    outfile <- file.path(out_dir, paste0(cancer, "_k", k_val, "_0asRef.csv"))
    message("Output file: ", outfile)

    for(sc in scenario_list){
      scen_key <- sc$key; eqi_p <- sc$eqi; aamr_p <- sc$aamr; lagv <- sc$lag
      scen_dt <- dt_k[EQI_Period==eqi_p & Time_Period==aamr_p & Cancer_Type==cancer]
      message("Pre-filter for ", scen_key, ": nrow=", nrow(scen_dt), ", unique FIPS=", length(unique(scen_dt$COUNTY_FIPS)), ", mean AAMR_lower=", mean(scen_dt$AAMR_lower, na.rm=TRUE))
      scen_dt <- scen_dt[!is.na(Cluster)]
      initial_counts <- table(scen_dt$Cluster, useNA="ifany")
      message("Post-filter for ", scen_key, ": nrow=", nrow(scen_dt), ", unique FIPS=", length(unique(scen_dt$COUNTY_FIPS)), ", mean AAMR_lower=", mean(scen_dt$AAMR_lower, na.rm=TRUE))
      message("Initial cluster counts for ", scen_key, ": ", paste(names(initial_counts), initial_counts, sep="=", collapse=", "))
      if(nrow(scen_dt) < opt$`min-n`){ message("[Skip] Scenario ", scen_key, " n=", nrow(scen_dt)); next }
      if(sum(scen_dt$Cluster == 0, na.rm=TRUE) == 0){ message("[Skip] No data in Cluster 0 for ", scen_key); next }
      eqi_out <- gsub('-', '_', eqi_p); aamr_out <- gsub('-', '_', aamr_p)

      # Cluster main effect model design (pooled)
      des_cluster <- build_design_cluster_main_effect(scen_dt, k_val)
      cluster_counts <- table(des_cluster$df$Cluster_factor, useNA="ifany")
      message("Cluster counts for ", scen_key, ": ", paste(names(cluster_counts), cluster_counts, sep="=", collapse=", "))
      message("Factor levels: ", paste(levels(des_cluster$df$Cluster_factor), collapse=", "))
      message("Unique clusters: ", paste(sort(unique(des_cluster$df$Cluster_factor)), collapse=", "))
      message("Model matrix columns: ", paste(des_cluster$names, collapse=", "))
      states_c <- sort(unique(des_cluster$df$State_FIPS)); state_index_c <- match(des_cluster$df$State_FIPS, states_c)
      data_list <- list(
        N = nrow(des_cluster$df), S = length(states_c), state = state_index_c,
        y_lower = des_cluster$df$AAMR_lower, y_upper = des_cluster$df$AAMR_upper, cens = des_cluster$df$cens,
        K = ncol(des_cluster$X), X = des_cluster$X
      )
      init_fun <- function() list(beta=rep(0, data_list$K), z_u=rep(0, data_list$S), sigma=50, sigma_u=10)
      fit_cluster <- try(mod$sample(data=data_list, chains=opt$chains, iter_sampling=opt$iter-opt$warmup, iter_warmup=opt$warmup,
                adapt_delta=opt$`adapt-delta`, max_treedepth=opt$`max-treedepth`, parallel_chains=min(opt$chains, cores_used), refresh=0, seed=opt$seed,
                init=rep(list(init_fun()), opt$chains)), silent=TRUE)
      if(inherits(fit_cluster,"try-error")){ message("[Fail] Cluster main effect model ", scen_key); next }

      draws <- as_draws_df(fit_cluster$draws("beta"))
      colnames(draws) <- paste0("beta[",seq_len(ncol(draws)),"]")
      row_cluster <- extract_cluster_main_effect_results(draws, des_cluster$names, cancer, eqi_out, aamr_out, lagv, k_val, "Cluster_Main_Effect", des_cluster$df, fit_cluster)
      append_rows(outfile, row_cluster)
      message("[OK] ", scen_key, " Cluster Main Effect")
    }
  }
  message("===== Completed: ", cancer, " =====")
}
message("All requested analyses complete. Output directory: ", out_dir)
