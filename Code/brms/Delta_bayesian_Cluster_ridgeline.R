
#!/usr/bin/env Rscript
# Delta cluster ridgeline runner
# Fits interval-censored mixed models for delta AAMR with EQI change categories,
# and saves posterior draws as RDS for ridgeline plots.
# Supports national and cluster-stratified runs (k=3) for C00_C97 and Lag in {5,10}.

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

option_list <- list(
  make_option(c("--data"), type="character",
              default="Data/Processed/df_EQI_AAMR_Triangulation/Delta_EQI_AAMR.csv",
              help="Input delta data CSV (must contain delta_AAMR_* and *_Change_Category columns)"),
  make_option(c("--output-dir"), type="character", default="Result/Delta_Ridgeline",
              help="Output directory for RDS files"),
  make_option(c("--cancer"), type="character", default="C00_C97",
              help="ICD code to analyze (default C00_C97)"),
  make_option(c("--lag"), type="integer", default=NA,
              help="Lag years to analyze (5 or 10). Default: run both 5 and 10"),
  make_option(c("--k"), type="integer", default=3,
              help="Number of clusters (default 3)"),
  make_option(c("--cluster"), type="character", default="national",
              help="Cluster to analyze: 'national' or one of {0,1,2} for k=3"),
  make_option(c("--chains"), type="integer", default=4),
  make_option(c("--iter"), type="integer", default=2000),
  make_option(c("--warmup"), type="integer", default=1000),
  make_option(c("--adapt-delta"), type="double", default=0.95),
  make_option(c("--max-treedepth"), type="integer", default=12),
  make_option(c("--seed"), type="integer", default=1234),
  make_option(c("--test"), action="store_true", default=FALSE,
              help="Reduce iterations for quick checks")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (opt$test) { opt$iter <- min(opt$iter, 800); opt$warmup <- min(opt$warmup, 300); message("[TEST MODE] iter=",opt$iter," warmup=",opt$warmup) }
set.seed(opt$seed)

# Threading
cores_avail <- parallel::detectCores(logical=TRUE); cores_used <- max(1, floor(cores_avail * 0.8)); options(mc.cores = cores_used)
message("Detected cores: ", cores_avail, " | Using: ", cores_used)

# Stan model (interval-censored mixed model with random intercept per state)
stan_code <- "data {\n  int<lower=1> N;\n  int<lower=1> S;\n  array[N] int<lower=1,upper=S> state;\n  vector[N] y_lower;\n  vector[N] y_upper;\n  array[N] int<lower=0,upper=2> cens;\n  int<lower=1> K;\n  matrix[N,K] X;\n}\nparameters {\n  vector[K] beta;\n  vector[S] z_u;\n  real<lower=0> sigma;\n  real<lower=0> sigma_u;\n}\ntransformed parameters {\n  vector[S] u = sigma_u * z_u;\n}\nmodel {\n  beta ~ normal(0,5);\n  z_u ~ normal(0,1);\n  sigma ~ exponential(1);\n  sigma_u ~ exponential(1);\n  for (i in 1:N) {\n    real mu = X[i] * beta + u[state[i]];\n    if (cens[i]==0) {\n      target += normal_lpdf(y_lower[i] | mu, sigma);\n    } else {\n      real p_up = normal_cdf(y_upper[i] | mu, sigma);\n      real p_lo = normal_cdf(y_lower[i] | mu, sigma);\n      real diff = fmax(p_up - p_lo, 1e-12);\n      target += log(diff);\n    }\n  }\n}\n"
stan_file <- file.path(tempdir(), "delta_interval_mixed_model.stan"); writeLines(stan_code, stan_file)
mod <- cmdstan_model(stan_file)

# Helpers
sig_mark <- function(p){ if(is.na(p)) return(""); if(p<0.001) return("***"); if(p<0.01) return("**"); if(p<0.05) return("*"); "" }
format_cell <- function(draws){ if(length(draws)==0) return(""); ci <- quantile(draws, c(0.025, 0.975), na.rm=TRUE); p <- 2 * min(mean(draws>0), mean(draws<0)); sprintf("%0.2f(%0.2f,%0.2f)%s", mean(draws), ci[1], ci[2], sig_mark(p)) }

compute_p <- function(draws){
  if(length(draws)==0) return(NA_real_)
  pos <- sum(draws > 0, na.rm = TRUE)
  neg <- sum(draws < 0, na.rm = TRUE)
  n <- pos + neg
  if(n == 0) return(NA_real_)
  p_pos <- (pos + 0.5) / (n + 1)
  p_neg <- (neg + 0.5) / (n + 1)
  2 * min(p_pos, p_neg)
}

extract_quintile_draws <- function(draw_df, names_vec, prefix){
  # Treatment contrasts expected: prefix2..prefix5 relative to Q1 baseline
  if(any(grepl(paste0(prefix,"\\.L"), names_vec))) {
    return(list(Q1=numeric(0), Q2=numeric(0), Q3=numeric(0), Q4=numeric(0), Q5=numeric(0)))
  }
  out <- list(Q1=numeric(0), Q2=numeric(0), Q3=numeric(0), Q4=numeric(0), Q5=numeric(0))
  for(q in 2:5){
    nm <- paste0(prefix, q)
    idx <- match(nm, names_vec)
    if(!is.na(idx)){
      col <- paste0("beta[", idx, "]")
      out[[paste0("Q", q)]] <- draw_df[[col]]
    }
  }
  out
}

build_design_overall <- function(d){ # Smoking + EQI change category Q2..Q5 vs Q1 baseline
  d <- d %>% mutate(EQI_change = factor(EQI_Change_Category, levels=1:5))
  d <- d[complete.cases(d[,c("delta_Smoking_Rate","EQI_change","delta_AAMR_lower","delta_AAMR_upper","cens","State_FIPS")]),]
  mm <- model.matrix(~ delta_Smoking_Rate + EQI_change, d, contrasts.arg=list(EQI_change=contr.treatment(5)))
  colnames(mm) <- make.names(colnames(mm))
  list(X=mm, names=colnames(mm), df=d)
}

build_design_single_domain <- function(d, domain){ # Smoking + domain change category Q2..Q5
  dom_col <- paste0(domain, "_Change_Category")
  fac <- factor(d[[dom_col]], levels=1:5)
  d2 <- d
  d2[[paste0(domain, "_change")]] <- fac
  d2 <- d2[complete.cases(d2[,c("delta_Smoking_Rate", paste0(domain,"_change"), "delta_AAMR_lower","delta_AAMR_upper","cens","State_FIPS")]),]
  form <- as.formula(paste0("~ delta_Smoking_Rate + ", domain, "_change"))
  mm <- model.matrix(form, d2, contrasts.arg=setNames(list(contr.treatment(5)), paste0(domain,"_change")))
  colnames(mm) <- make.names(colnames(mm))
  list(X=mm, names=colnames(mm), df=d2, prefix=paste0(domain, "_change"))
}

fit_and_extract <- function(design_info, chains, iter, warmup, adapt_delta, max_treedepth, seed){
  states <- sort(unique(design_info$df$State_FIPS)); state_index <- match(design_info$df$State_FIPS, states)
  stan_data <- list(
    N = nrow(design_info$df), S = length(states), state = state_index,
    y_lower = design_info$df$delta_AAMR_lower, y_upper = design_info$df$delta_AAMR_upper, cens = design_info$df$cens,
    K = ncol(design_info$X), X = design_info$X
  )
  init_fun <- function() list(beta=rep(0, stan_data$K), z_u=rep(0, stan_data$S), sigma=50, sigma_u=10)
  fit <- mod$sample(
    data=stan_data, chains=chains, iter_sampling=iter-warmup, iter_warmup=warmup,
    adapt_delta=adapt_delta, max_treedepth=max_treedepth,
    parallel_chains=min(chains, cores_used), refresh=0, seed=seed,
    init=rep(list(init_fun()), chains)
  )
  draws_df <- as_draws_df(fit$draws("beta"))
  colnames(draws_df) <- paste0("beta[", seq_len(ncol(draws_df)), "]")
  summ <- posterior::summarize_draws(fit$draws("beta"))
  list(fit=fit, draws=draws_df, summary=summ)
}

make_ridge_rds <- function(draws_df, names_vec, prefix, meta){
  # Collect draws
  q2 <- draws_df[[paste0("beta[", match(paste0(prefix, 2), names_vec), "]")]]
  q3 <- draws_df[[paste0("beta[", match(paste0(prefix, 3), names_vec), "]")]]
  q4 <- draws_df[[paste0("beta[", match(paste0(prefix, 4), names_vec), "]")]]
  q5 <- draws_df[[paste0("beta[", match(paste0(prefix, 5), names_vec), "]")]]
  n_draws <- max(length(q2), length(q3), length(q4), length(q5))
  draws_wide <- tibble::tibble(draw=q2) %>%
    mutate(Q2=q2, Q3=q3, Q4=q4, Q5=q5) %>%
    select(Q2,Q3,Q4,Q5)
  draws_long <- draws_wide %>%
    pivot_longer(cols=everything(), names_to="Quintile", values_to="Effect") %>%
    mutate(Model=meta$model, Domain=meta$domain, ICD_Code=meta$icd, Lag=meta$lag, Cluster=meta$cluster)
  # Summaries
  summary_df <- draws_long %>%
    group_by(Quintile, Model, Domain, ICD_Code, Lag, Cluster) %>%
    summarize(mean=mean(Effect, na.rm=TRUE),
              l95=quantile(Effect, 0.025, na.rm=TRUE),
              u95=quantile(Effect, 0.975, na.rm=TRUE),
              p=compute_p(Effect), .groups="drop")
  list(draws_long=draws_long, summary=summary_df, meta=meta)
}

# Load data
project_root <- normalizePath(".")
data_path <- file.path(project_root, opt$data)
if(!file.exists(data_path)) stop("Data not found: ", data_path)
dt <- fread(data_path)

# Validate and prepare columns
req <- c("COUNTY_FIPS","State_FIPS","Cancer_Type","Lag",
         "delta_AAMR_lower","delta_AAMR_upper","delta_Smoking_Rate",
         "EQI_Change_Category","Air_Change_Category","Water_Change_Category",
         "Land_Change_Category","Built_Change_Category","Social_Change_Category")
miss <- setdiff(req, names(dt)); if(length(miss)) stop("Missing columns: ", paste(miss,collapse=","))

# Interval codes
dt <- dt[!is.na(delta_AAMR_lower) & !is.na(delta_AAMR_upper)]
dt[, cens := ifelse(delta_AAMR_lower == delta_AAMR_upper, 0, 2)]

# Filter cancer and lag
dt <- dt[Cancer_Type == opt$cancer]
lags_to_run <- if (is.na(opt$lag)) c(5,10) else {
  if (!opt$lag %in% c(5,10)) stop("--lag must be 5 or 10"); c(opt$lag)
}

# Merge cluster info (k=opt$k)
cluster_path <- file.path(project_root, "Result/Cluster_Visualization/EQI_Clusters_All_K.csv")
if(!file.exists(cluster_path)) stop("Cluster data not found: ", cluster_path)
cluster_dt <- fread(cluster_path)
cluster_col <- paste0("cluster_", opt$k)
if(!cluster_col %in% names(cluster_dt)) stop("Cluster column missing in cluster file: ", cluster_col)
dt <- merge(dt, cluster_dt[, c("COUNTY_FIPS", cluster_col), with=FALSE], by="COUNTY_FIPS", all.x=TRUE)

# Subset by cluster selection
cluster_sel <- opt$cluster
if (cluster_sel == "national") {
  dt_sel <- dt
} else {
  if (!cluster_sel %in% as.character(0:(opt$k-1))) stop("--cluster must be 'national' or one of: ", paste(0:(opt$k-1), collapse=","))
  dt_sel <- dt[get(cluster_col) == as.integer(cluster_sel)]
}
if (nrow(dt_sel) == 0) stop("No rows for cluster selection: ", cluster_sel)

# Output directory
out_dir <- file.path(project_root, opt$`output-dir`)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

domains <- c("Air","Water","Land","Built","Social")

for (lag_run in lags_to_run) {
  message("Running models for Lag=", lag_run)
  dt_lag <- dt_sel[Lag == lag_run]
  if (nrow(dt_lag) == 0) { message("No rows for Lag=", lag_run, " after cluster filtering; skipping"); next }

  # Overall EQI change
  message("Running Overall EQI change model...")
  design_overall <- build_design_overall(dt_lag)
  res_overall <- fit_and_extract(design_overall, opt$chains, opt$iter, opt$warmup, opt$`adapt-delta`, opt$`max-treedepth`, opt$seed)
  meta_overall <- list(model="EQI", domain="Overall", icd=opt$cancer, lag=lag_run, cluster=cluster_sel)
  ridge_overall <- make_ridge_rds(res_overall$draws, design_overall$names, "EQI_change", meta_overall)

  # Single-domain change models
  ridge_domains <- list()
  for (dom in domains) {
    message("Running domain change model: ", dom)
    design_dom <- build_design_single_domain(dt_lag, dom)
    res_dom <- fit_and_extract(design_dom, opt$chains, opt$iter, opt$warmup, opt$`adapt-delta`, opt$`max-treedepth`, opt$seed)
    meta_dom <- list(model=dom, domain=dom, icd=opt$cancer, lag=lag_run, cluster=cluster_sel)
    ridge_domains[[dom]] <- make_ridge_rds(res_dom$draws, design_dom$names, paste0(dom, "_change"), meta_dom)
  }

  # Combine into a single RDS with 6 entries (EQI + 5 domains)
  combined <- list(
    EQI = ridge_overall,
    Air = ridge_domains[["Air"]],
    Water = ridge_domains[["Water"]],
    Land = ridge_domains[["Land"]],
    Built = ridge_domains[["Built"]],
    Social = ridge_domains[["Social"]]
  )
  out_file <- file.path(out_dir, paste0(opt$cancer, "_Lag", lag_run, "_", if (cluster_sel=="national") "National" else paste0("Cluster", cluster_sel), ".rds"))
  saveRDS(combined, out_file)
  message("✓ Saved combined ridgeline RDS: ", out_file)
}

message("✅ Completed ridgeline RDS generation for Cancer=", opt$cancer, " Lag=", opt$lag, " Cluster=", cluster_sel)
