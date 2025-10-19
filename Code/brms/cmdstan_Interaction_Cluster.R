#!/usr/bin/env Rscript
# Pure cmdstanr interval-censored mixed model pipeline for Quintile * Cluster interactions
# Output format mimics LMM results (ICD_Code, EQI_Period, AAMR_Period, Lag, Model, Q1..Q5)
# Uses cluster IDs from clustering analysis, with Cluster2 (Advantageous) as baseline

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
utils::globalVariables(c('EQI','EQI_Air','EQI_Water','EQI_Land','EQI_Built','EQI_Social','Smoking_Rate','State_FIPS','Cluster'))

option_list <- list(
  make_option(c("--data"), type="character", default="Data/Processed/df_EQI_AAMR/EQI_AAMR_Interval_Clustered.csv", help="Input interval data with cluster IDs"),
  make_option(c("--output-dir"), type="character", default="Result/brms_interaction_cluster", help="Output directory"),
  make_option(c("--cancer-types"), type="character", default="C00_C97", help="Comma separated ICD codes (default C00_C97)"),
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

# Parse cancer types
all_cancers <- sort(unique(dt$Cancer_Type))
selected <- if (is.na(opt$`cancer-types`)) all_cancers else { reqc <- str_split(opt$`cancer-types`,",",simplify=TRUE) |> as.vector() |> str_trim(); inv <- setdiff(reqc,all_cancers); if(length(inv)) stop("Invalid cancer types: ", paste(inv,collapse=",")); reqc }
message("Cancer types to analyze: ", paste(selected,collapse=","))

req <- c("COUNTY_FIPS","EQI_Period","Time_Period","Lag_Years","Cancer_Type","AAMR_lower","AAMR_upper","Smoking_Rate","EQI","EQI_Air","EQI_Water","EQI_Land","EQI_Built","EQI_Social","Cluster")
miss <- setdiff(req, names(dt)); if(length(miss)) stop("Missing cols: ", paste(miss,collapse=","))

if(!"State_FIPS" %in% names(dt)) dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS),1,2)]

# Cluster mapping (Cluster2 as baseline)
dt[, Cluster := factor(Cluster, levels=c(2,1,0), labels=c("Advantageous","Unbalanced","Disadvantaged"))]

# interval censoring code
dt <- dt[!is.na(AAMR_lower) & !is.na(AAMR_upper)]
dt[, cens := ifelse(AAMR_lower == AAMR_upper, 0, 2)]

scenario_list <- list(
  list(key="EQI0005_AAMR2006_2010", eqi="2000-2005", aamr="2006-2010", lag=5),
  list(key="EQI0005_AAMR2011_2015", eqi="2000-2005", aamr="2011-2015", lag=10),
  list(key="EQI0610_AAMR2011_2015", eqi="2006-2010", aamr="2011-2015", lag=5),
  list(key="EQI0610_AAMR2016_2020", eqi="2006-2010", aamr="2016-2020", lag=10)
)

out_dir <- file.path(project_root,opt$`output-dir`); if(!dir.exists(out_dir)) dir.create(out_dir,recursive=TRUE)

sig_mark <- function(p){ if(is.na(p)) return(""); if(p<0.001) return("***"); if(p<0.01) return("**"); if(p<0.05) return("*"); "" }
format_cell <- function(draws){ if(length(draws)==0) return(""); ci <- quantile(draws,c(0.025,0.975),na.rm=TRUE); p <- 2*min(mean(draws>0), mean(draws<0)); sprintf("%0.2f(%0.2f,%0.2f)%s", mean(draws), ci[1], ci[2], sig_mark(p)) }
append_rows <- function(path, df){ if(!file.exists(path)) write_csv(df,path) else suppressWarnings(write.table(df,path,sep=",",col.names=FALSE,row.names=FALSE,append=TRUE)) }

# Function to build design matrix for interaction: Quintile * Cluster
build_design_interaction <- function(d, quintile_var){
  d <- d %>% mutate(Quintile = factor(.data[[quintile_var]], levels=1:5))
  d <- d[complete.cases(d[,c("Smoking_Rate","Quintile","Cluster","AAMR_lower","AAMR_upper","cens","State_FIPS")]),]
  mm <- model.matrix(~ Smoking_Rate + Quintile * Cluster, d, contrasts.arg = list(Quintile=contr.treatment(5), Cluster=contr.treatment(3)))
  colnames(mm) <- make.names(colnames(mm))
  list(X = mm, names = colnames(mm), df = d)
}

# Extract quintiles and interaction p-values (relative to Cluster2 baseline)
extract_interaction_results <- function(draws_df, names_vec, quintile_prefix, cluster_levels){
  message("Design matrix columns: ", paste(names_vec, collapse = ", "))
  
  # Cluster numeric mapping: Cluster1=baseline(Advantageous), Cluster2=Unbalanced, Cluster3=Disadvantaged
  # Note: After factor(), levels are c("Advantageous","Unbalanced","Disadvantaged")
  # But model.matrix uses numeric codes, so we need to map:
  # Cluster2 in colname = Unbalanced (2nd level)
  # Cluster3 in colname = Disadvantaged (3rd level)
  cluster_map <- list(
    "Unbalanced" = "Cluster2",
    "Disadvantaged" = "Cluster3"
  )
  
  # Quintile main effects (relative to Q1)
  q_main <- list()
  for(q in 2:5){
    nm <- paste0(quintile_prefix, q)
    idx <- match(nm, names_vec)
    if(!is.na(idx)) {
      q_main[[paste0("Q",q)]] <- draws_df[[paste0("beta[",idx,"]")]]
      message("Found main effect: ", nm, " at index ", idx)
    } else {
      message("WARNING: Main effect not found: ", nm)
    }
  }
  
  # Interaction effects: QuintileX.ClusterY (after make.names, colon becomes dot)
  int_effects <- list()
  for(q in 2:5){
    for(cl in names(cluster_map)){  # Loop through Unbalanced and Disadvantaged
      cluster_num <- cluster_map[[cl]]
      # After make.names: "Quintile2:Cluster2" becomes "Quintile2.Cluster2"
      int_nm <- paste0(quintile_prefix, q, ".", cluster_num)
      idx <- match(int_nm, names_vec)
      if(!is.na(idx)) {
        int_effects[[paste0("Q",q,".",cl)]] <- draws_df[[paste0("beta[",idx,"]")]]
        message("Found interaction: ", int_nm, " (", cl, ") at index ", idx)
      } else {
        message("WARNING: Interaction term not found: ", int_nm)
      }
    }
  }
  
  # 对于每个quintile和cluster，计算完整效应
  # 效应 = 主效应(相对Q1) + 交互效应(相对baseline cluster)
  mrd <- list()
  for(q in 2:5){
    mrd[[paste0("Q",q)]] <- list()
    baseline_q <- if(length(q_main[[paste0("Q",q)]])>0) q_main[[paste0("Q",q)]] else rep(0, nrow(draws_df))
    
    # 对于baseline cluster (Advantageous)，只有主效应
    mrd[[paste0("Q",q)]][["Advantageous"]] <- baseline_q
    
    # 对于其他clusters，主效应 + 交互效应
    for(cl in names(cluster_map)){
      int_draws <- int_effects[[paste0("Q",q,".",cl)]]
      if(length(int_draws)>0) {
        mrd[[paste0("Q",q)]][[cl]] <- baseline_q + int_draws
        message("Q", q, " ", cl, ": mean main=", round(mean(baseline_q),2), 
                ", mean interaction=", round(mean(int_draws),2),
                ", mean total=", round(mean(baseline_q + int_draws),2))
      } else {
        mrd[[paste0("Q",q)]][[cl]] <- baseline_q
        message("Q", q, " ", cl, ": no interaction found, using main effect only")
      }
    }
  }
  
  # Interaction p-values for Q5 (testing if interaction is significant)
  int_p <- list()
  int_effect_q5 <- list()
  for(cl in names(cluster_map)){
    int_draws <- int_effects[[paste0("Q5.",cl)]]
    if(length(int_draws)>0){
      p <- 2*min(mean(int_draws>0), mean(int_draws<0))
      int_p[[cl]] <- sprintf("%.3f", p)
      int_effect_q5[[cl]] <- format_cell(int_draws)  # Format the pure interaction effect
      message("Interaction p-value for Q5 x ", cl, ": ", p, ", effect: ", mean(int_draws))
    } else {
      int_p[[cl]] <- "NA"
      int_effect_q5[[cl]] <- ""
      message("No interaction p-value for ", cl)
    }
  }
  
  list(mrd_q2 = mrd$Q2, mrd_q3 = mrd$Q3, mrd_q4 = mrd$Q4, mrd_q5 = mrd$Q5, 
       int_p = int_p, int_effect_q5 = int_effect_q5)
}

# Domains
domains <- c("EQI", "EQI_Air", "EQI_Built", "EQI_Water", "EQI_Land", "EQI_Social")
domain_labels <- c("EQI", "Air", "Built", "Water", "Land", "Social")

for(cancer in selected){
  message("===== Disease: ", cancer, " =====")
  outfile <- file.path(out_dir, paste0(cancer, "_Interaction_Cluster.csv"))
  
  for(sc_idx in seq_along(scenario_list)){
    sc <- scenario_list[[sc_idx]]
    scen_key <- sc$key; eqi_p <- sc$eqi; aamr_p <- sc$aamr; lagv <- sc$lag
    scen_dt <- dt[EQI_Period==eqi_p & Time_Period==aamr_p & Cancer_Type==cancer]
    if(nrow(scen_dt) < opt$`min-n`){ message("[Skip] Scenario ", scen_key, " overall n=", nrow(scen_dt)); next }
    eqi_out <- gsub('-', '_', eqi_p); aamr_out <- gsub('-', '_', aamr_p)
    
    for(dom_idx in seq_along(domains)){
      dom <- domains[dom_idx]; dom_lab <- domain_labels[dom_idx]
      dom_dt <- scen_dt
      if(nrow(dom_dt) < opt$`min-n`){ message("[Skip] ", scen_key, " ", dom_lab, " n=", nrow(dom_dt)); next }
      
      # Build design
      des <- build_design_interaction(dom_dt, dom)
      states <- sort(unique(des$df$State_FIPS)); state_index <- match(des$df$State_FIPS, states)
      data_list <- list(
        N = nrow(des$df), S = length(states), state = state_index,
        y_lower = des$df$AAMR_lower, y_upper = des$df$AAMR_upper, cens = des$df$cens,
        K = ncol(des$X), X = des$X
      )
      
      # Initial values
      init_fun <- function() list(beta=rep(0, data_list$K), z_u=rep(0, data_list$S), sigma=50, sigma_u=10)
      fit <- try(mod$sample(data=data_list, chains=opt$chains, iter_sampling=opt$iter-opt$warmup, iter_warmup=opt$warmup,
                adapt_delta=opt$`adapt-delta`, max_treedepth=opt$`max-treedepth`, parallel_chains=min(opt$chains, cores_used), refresh=0, seed=opt$seed,
                init=rep(list(init_fun()), opt$chains)), silent=TRUE)
      
      # Check convergence
      if(inherits(fit,"try-error") || any(fit$summary()$rhat > 1.01)){
        message("[Retry] ", scen_key, " ", dom_lab, " with iter=4000")
        opt$iter <- 4000; opt$warmup <- 2000
        fit <- try(mod$sample(data=data_list, chains=opt$chains, iter_sampling=opt$iter-opt$warmup, iter_warmup=opt$warmup,
                  adapt_delta=opt$`adapt-delta`, max_treedepth=opt$`max-treedepth`, parallel_chains=min(opt$chains, cores_used), refresh=0, seed=opt$seed,
                  init=rep(list(init_fun()), opt$chains)), silent=TRUE)
        if(inherits(fit,"try-error")){ message("[Fail] ", scen_key, " ", dom_lab); next }
      }
      
      draws <- as_draws_df(fit$draws("beta"))
      colnames(draws) <- paste0("beta[",seq_len(ncol(draws)),"]")
      
      res <- extract_interaction_results(draws, des$names, "Quintile", c("Advantageous","Unbalanced","Disadvantaged"))
      
      # Output rows for each cluster (only Unbalanced and Disadvantaged, relative to Advantageous)
      clusters <- c("Unbalanced","Disadvantaged")
      for(cl_idx in seq_along(clusters)){
        cl <- clusters[cl_idx]
        int_p_val <- res$int_p[[cl]]
        int_effect <- res$int_effect_q5[[cl]]  # Pure interaction effect for Q5
        row <- tibble(
          ICD_Code = cancer,
          EQI_Period = eqi_out,
          AAMR_Period = aamr_out,
          Lag = lagv,
          Model = paste0(dom_lab, "_", cl),
          Interaction_Q5 = int_effect,
          Interaction_P_Value = int_p_val,
          Q1 = "0.00",
          Q2 = format_cell(res$mrd_q2[[cl]]),
          Q3 = format_cell(res$mrd_q3[[cl]]),
          Q4 = format_cell(res$mrd_q4[[cl]]),
          Q5 = format_cell(res$mrd_q5[[cl]])
        )
        append_rows(outfile, row)
      }
      message("[OK] ", scen_key, " ", dom_lab)
    }
  }
  message("===== Completed: ", cancer, " =====")
}
message("All requested analyses complete. Output directory: ", out_dir)