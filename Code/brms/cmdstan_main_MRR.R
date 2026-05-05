#!/usr/bin/env Rscript
# cmdstanr MRR pipeline — Overall EQI only, EQI0005 lags only (5, 10, 15)
# Outputs (all in Result/brms_MRR_lag/):
#   {outcome}_main.csv     — MRD (mean difference, unchanged format)
#   {outcome}_MRR.csv      — Mortality Rate Ratio Q2-Q5 vs Q1
#   {outcome}_lag_test.csv — Pairwise lag comparison (Q5 draws)

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
  make_option(c("--data"), type="character",
              default="Data/Processed/df_EQI_AAMR_Triangulation/df_Main.csv",
              help="Input interval data"),
  make_option(c("--output-dir"), type="character", default="Result/brms_MRR_lag", help="Output directory for all results"),
  make_option(c("--outcomes"),        type="character", default=NA,                     help="Comma separated ICD codes"),
  make_option(c("--chains"),         type="integer",   default=4),
  make_option(c("--iter"),           type="integer",   default=2000),
  make_option(c("--warmup"),         type="integer",   default=1000),
  make_option(c("--adapt-delta"),    type="double",    default=0.95),
  make_option(c("--max-treedepth"),  type="integer",   default=12),
  make_option(c("--min-n"),          type="integer",   default=50),
  make_option(c("--seed"),           type="integer",   default=1234),
  make_option(c("--test"),           action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))
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
# Stan model
# ---------------------------------------------------------------------------
stan_code <- "data {\n  int<lower=1> N;\n  int<lower=1> S;\n  array[N] int<lower=1,upper=S> state;\n  vector[N] y_lower;\n  vector[N] y_upper;\n  array[N] int<lower=0,upper=2> cens;\n  int<lower=1> K;\n  matrix[N,K] X;\n} \nparameters {\n  vector[K] beta;\n  vector[S] z_u;\n  real<lower=0> sigma;\n  real<lower=0> sigma_u;\n} \ntransformed parameters {\n  vector[S] u = sigma_u * z_u;\n} \nmodel {\n  beta ~ normal(0,5);\n  z_u ~ normal(0,1);\n  sigma ~ exponential(1);\n  sigma_u ~ exponential(1);\n  for (i in 1:N) {\n    real mu = X[i] * beta + u[state[i]];\n    if (cens[i]==0) {\n      target += normal_lpdf(y_lower[i] | mu, sigma);\n    } else {\n      real p_up = normal_cdf(y_upper[i] | mu, sigma);\n      real p_lo = normal_cdf(y_lower[i] | mu, sigma);\n      real diff = fmax(p_up - p_lo, 1e-12);\n      target += log(diff);\n    }\n  }\n}"
stan_file <- file.path(tempdir(), "interval_mixed_model.stan")
writeLines(stan_code, stan_file)
mod <- cmdstan_model(stan_file)

# ---------------------------------------------------------------------------
# Load and prep data
# ---------------------------------------------------------------------------
project_root <- normalizePath(".")
path <- file.path(project_root, opt$data)
if (!file.exists(path)) stop("Data not found: ", path)
dt <- fread(path)

req  <- c("COUNTY_FIPS","EQI_Period","Time_Period","Lag_Years","Outcome",
          "AAMR_Lower","AAMR_Upper","Smoking_Rate","RUCC",
          "EQI","EQI_Air","EQI_Water","EQI_Land","EQI_Built","EQI_Social")
miss <- setdiff(req, names(dt))
if (length(miss)) stop("Missing cols: ", paste(miss, collapse=","))

if (!"State_FIPS" %in% names(dt))
  dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]

dt <- dt[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
dt[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]
dt <- dt[RUCC %in% 1:4 | is.na(RUCC)]

# ---------------------------------------------------------------------------
# Scenarios: EQI0005 only — lags 5, 10, 15
# ---------------------------------------------------------------------------
scenario_list <- list(
  list(key="EQI0005_AAMR2006_2010", eqi="2000-2005", aamr="2006-2010", lag=5),
  list(key="EQI0005_AAMR2011_2015", eqi="2000-2005", aamr="2011-2015", lag=10),
  list(key="EQI0005_AAMR2016_2020", eqi="2000-2005", aamr="2016-2020", lag=15)
)

all_outcomes <- sort(unique(dt$Outcome))
selected <- if (is.na(opt$`outcomes`)) {
  all_outcomes
} else {
  reqc <- str_split(opt$`outcomes`, ",", simplify=TRUE) |> as.vector() |> str_trim()
  inv  <- setdiff(reqc, all_outcomes)
  if (length(inv)) stop("Invalid outcomes: ", paste(inv, collapse=","))
  reqc
}
message("Outcomes to analyze: ", paste(selected, collapse=","))

# ---------------------------------------------------------------------------
# Output directory (single dir for all 3 file types)
# ---------------------------------------------------------------------------
out_dir <- file.path(project_root, opt$`output-dir`)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
sig_mark <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***"); if (p < 0.01) return("**"); if (p < 0.05) return("*"); ""
}

format_cell <- function(draws) {
  if (length(draws) == 0) return("")
  ci <- quantile(draws, c(0.025, 0.975), na.rm=TRUE)
  p  <- 2 * min(mean(draws > 0), mean(draws < 0))
  sprintf("%0.2f(%0.2f,%0.2f)%s", mean(draws), ci[1], ci[2], sig_mark(p))
}

append_rows <- function(path, df) {
  if (!file.exists(path)) write_csv(df, path) else
    suppressWarnings(write.table(df, path, sep=",", col.names=FALSE, row.names=FALSE, append=TRUE))
}

compute_p <- function(draws) {
  if (length(draws) == 0) return(NA_character_)
  pos <- sum(draws > 0, na.rm=TRUE); neg <- sum(draws < 0, na.rm=TRUE); n <- pos + neg
  if (n == 0) return(NA_character_)
  p_pos <- (pos + 0.5) / (n + 1); p_neg <- (neg + 0.5) / (n + 1)
  sprintf("%.4f", 2 * min(p_pos, p_neg))
}

# Format a probability value: "p<0.0001" if tiny, otherwise 4 decimal places.
format_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.0001) return("p<0.0001")
  sprintf("%.4f", p)
}

extract_quintile_metrics <- function(draw_df, names_vec, prefix, summ_df) {
  out <- list(Q2_p=NA_real_, Q3_p=NA_real_, Q4_p=NA_real_, Q5_p=NA_real_,
              Q2_rhat=NA_real_, Q3_rhat=NA_real_, Q4_rhat=NA_real_, Q5_rhat=NA_real_,
              Q2_ess_bulk=NA_real_, Q3_ess_bulk=NA_real_, Q4_ess_bulk=NA_real_, Q5_ess_bulk=NA_real_,
              Q2_ess_tail=NA_real_, Q3_ess_tail=NA_real_, Q4_ess_tail=NA_real_, Q5_ess_tail=NA_real_)
  if (any(grepl(paste0(prefix, "\\.L"), names_vec))) return(out)
  for (q in 2:5) {
    nm  <- paste0(prefix, q); idx <- match(nm, names_vec)
    if (!is.na(idx)) {
      col <- paste0("beta[", idx, "]")
      out[[paste0("Q", q, "_p")]] <- compute_p(draw_df[[col]])
      sr  <- summ_df[summ_df$variable == col, , drop=FALSE]
      if (nrow(sr)) {
        out[[paste0("Q", q, "_rhat")]]     <- sr$rhat
        out[[paste0("Q", q, "_ess_bulk")]] <- sr$ess_bulk
        out[[paste0("Q", q, "_ess_tail")]] <- sr$ess_tail
      }
    }
  }
  out
}

build_design_overall <- function(d) {
  d  <- d %>% mutate(EQI_factor = factor(EQI, levels=1:5))
  d  <- d[complete.cases(d[, c("Smoking_Rate","EQI_factor","AAMR_Lower","AAMR_Upper","cens","State_FIPS")]), ]
  mm <- model.matrix(~ Smoking_Rate + EQI_factor, d,
                     contrasts.arg = list(EQI_factor=contr.treatment(5)))
  colnames(mm) <- make.names(colnames(mm))
  list(X=mm, names=colnames(mm), df=d)
}

extract_quintiles <- function(draw_df, names_vec, prefix) {
  out <- list(Q1="0.00", Q2="", Q3="", Q4="", Q5="")
  if (any(grepl(paste0(prefix, "\\.L"), names_vec))) return(out)
  for (q in 2:5) {
    nm  <- paste0(prefix, q); idx <- match(nm, names_vec)
    out[[paste0("Q", q)]] <- if (is.na(idx)) "" else format_cell(draw_df[[paste0("beta[", idx, "]")]])
  }
  out
}

# Compute MRR (Q1-Q5 vs Lag5 Q1) from posterior draws.
# ref_draws: posterior draws of mu_Q1 from the Lag5 model (the universal reference).
#            If NULL (should not happen after lag5 is processed), falls back to current Q1.
compute_mrr <- function(draws, names_vec, layer_dt, outcome, eqi_out, aamr_out, lagv, ref_draws) {
  smoking_idx  <- match("Smoking_Rate", names_vec)
  mean_smoking <- mean(layer_dt$Smoking_Rate, na.rm=TRUE)
  mu_Q1_draws  <- draws[["beta[1]"]] + draws[[paste0("beta[", smoking_idx, "]")]] * mean_smoking

  mrr_rows <- lapply(1:5, function(q) {
    if (q == 1) {
      mu_Qq <- mu_Q1_draws
    } else {
      q_idx <- match(paste0("EQI_factor", q), names_vec)
      if (is.na(q_idx)) return(NULL)
      mu_Qq <- mu_Q1_draws + draws[[paste0("beta[", q_idx, "]")]]
    }
    MRR_q <- mu_Qq / ref_draws
    # Two-sided posterior p: P(MRR != 1) via Jeffreys-corrected proportion
    mrr_centered <- MRR_q - 1
    pos <- sum(mrr_centered > 0, na.rm=TRUE); neg <- sum(mrr_centered < 0, na.rm=TRUE); nn <- pos + neg
    p_raw <- if (nn == 0) NA_real_ else 2 * min((pos + 0.5) / (nn + 1), (neg + 0.5) / (nn + 1))
    tibble(
      ICD_Code    = outcome,
      EQI_Period  = eqi_out,
      AAMR_Period = aamr_out,
      Lag         = lagv,
      Quintile    = paste0("Q", q),
      MRR_mean    = round(mean(MRR_q,            na.rm=TRUE), 4),
      MRR_lower   = round(quantile(MRR_q, 0.025, na.rm=TRUE), 4),
      MRR_upper   = round(quantile(MRR_q, 0.975, na.rm=TRUE), 4),
      pct_diff    = round((mean(MRR_q,            na.rm=TRUE) - 1) * 100, 4),
      p           = format_p(p_raw)
    )
  })
  bind_rows(Filter(Negate(is.null), mrr_rows))
}

# Pairwise lag test on Q5 beta draws: lag5 vs lag10, lag10 vs lag15, lag15 vs lag5.
run_lag_test <- function(lag_q5_store, outcome, out_dir) {
  pairs <- list(c("5","10"), c("10","15"), c("15","5"))
  rows  <- lapply(pairs, function(p) {
    la <- p[1]; lb <- p[2]
    if (!la %in% names(lag_q5_store) || !lb %in% names(lag_q5_store)) return(NULL)
    da <- lag_q5_store[[la]]; db <- lag_q5_store[[lb]]
    diff_draws <- da - db
    tibble(
      ICD_Code   = outcome,
      comparison = paste0("lag", la, "_vs_lag", lb),
      diff_mean  = round(mean(diff_draws,            na.rm=TRUE), 4),
      diff_lower = round(quantile(diff_draws, 0.025, na.rm=TRUE), 4),
      diff_upper = round(quantile(diff_draws, 0.975, na.rm=TRUE), 4),
      P_a_gt_b   = format_p(mean(da > db,            na.rm=TRUE))
    )
  })
  result <- bind_rows(Filter(Negate(is.null), rows))
  if (nrow(result) > 0) {
    lag_file <- file.path(out_dir, paste0(outcome, "_lag_test.csv"))
    append_rows(lag_file, result)
    message("[LAG TEST] ", nrow(result), " comparisons written for ", outcome)
  }
}

# ===========================================================================
# Main loop
# ===========================================================================
for (outcome in selected) {
  message("===== Outcome: ", outcome, " =====")
  outfile        <- file.path(out_dir, paste0(outcome, "_main.csv"))
  mrr_file       <- file.path(out_dir, paste0(outcome, "_MRR.csv"))
  lag_q5_store   <- list()   # keyed by lag value: "5", "10", "15"
  lag5_ref_draws <- NULL     # Lag5 Q1 draws; used as universal MRR reference for all lags

  for (sc in scenario_list) {
    scen_key <- sc$key; eqi_p <- sc$eqi; aamr_p <- sc$aamr; lagv <- sc$lag
    scen_dt  <- dt[EQI_Period==eqi_p & Time_Period==aamr_p & Outcome==outcome]
    if (nrow(scen_dt) < opt$`min-n`) {
      message("[Skip] ", scen_key, " n=", nrow(scen_dt)); next
    }
    eqi_out  <- gsub('-', '_', eqi_p)
    aamr_out <- gsub('-', '_', aamr_p)

    # Overall layer only (no RUCC stratification)
    layer_dt <- scen_dt

    des_overall   <- build_design_overall(layer_dt)
    states_o      <- sort(unique(des_overall$df$State_FIPS))
    state_index_o <- match(des_overall$df$State_FIPS, states_o)
    data_list <- list(
      N       = nrow(des_overall$df),
      S       = length(states_o),
      state   = state_index_o,
      y_lower = des_overall$df$AAMR_Lower,
      y_upper = des_overall$df$AAMR_Upper,
      cens    = des_overall$df$cens,
      K       = ncol(des_overall$X),
      X       = des_overall$X
    )
    init_fun <- function() list(beta=rep(0, data_list$K), z_u=rep(0, data_list$S), sigma=50, sigma_u=10)

    fit_overall <- try(mod$sample(
      data            = data_list,
      chains          = opt$chains,
      iter_sampling   = opt$iter - opt$warmup,
      iter_warmup     = opt$warmup,
      adapt_delta     = opt$`adapt-delta`,
      max_treedepth   = opt$`max-treedepth`,
      parallel_chains = min(opt$chains, cores_used),
      refresh         = 0,
      seed            = opt$seed,
      init            = rep(list(init_fun()), opt$chains)
    ), silent=TRUE)

    if (inherits(fit_overall, "try-error")) {
      message("[Fail] ", scen_key); next
    }

    draws <- as_draws_df(fit_overall$draws("beta"))
    colnames(draws) <- paste0("beta[", seq_len(ncol(draws)), "]")

    # ---- Compute and cache Lag5 Q1 as universal MRR reference ----
    {
      smoking_idx_ref  <- match("Smoking_Rate", des_overall$names)
      mean_smoking_ref <- mean(layer_dt$Smoking_Rate, na.rm=TRUE)
      mu_Q1_current    <- draws[["beta[1]"]] + draws[[paste0("beta[", smoking_idx_ref, "]")]] * mean_smoking_ref
      if (lagv == 5) lag5_ref_draws <- mu_Q1_current
    }
    if (is.null(lag5_ref_draws)) {
      message("[WARN] lag5_ref_draws not yet available for lag", lagv, " — skipping MRR"); next
    }

    # ---- MRD output (original format, unchanged) ----
    q_over    <- extract_quintiles(draws, des_overall$names, "EQI_factor")
    summ_over <- posterior::summarize_draws(fit_overall$draws("beta"))
    met_over  <- extract_quintile_metrics(draws, des_overall$names, "EQI_factor", summ_over)
    row_over  <- tibble(
      ICD_Code=outcome, EQI_Period=eqi_out, AAMR_Period=aamr_out, Lag=lagv, Model="EQI",
      Q1=q_over$Q1, Q2=q_over$Q2, Q3=q_over$Q3, Q4=q_over$Q4, Q5=q_over$Q5,
      Q2_p=met_over$Q2_p, Q3_p=met_over$Q3_p, Q4_p=met_over$Q4_p, Q5_p=met_over$Q5_p,
      Q2_rhat=sprintf("%.4f", met_over$Q2_rhat),
      Q3_rhat=sprintf("%.4f", met_over$Q3_rhat),
      Q4_rhat=sprintf("%.4f", met_over$Q4_rhat),
      Q5_rhat=sprintf("%.4f", met_over$Q5_rhat),
      Q2_ess_bulk=as.integer(round(met_over$Q2_ess_bulk)),
      Q3_ess_bulk=as.integer(round(met_over$Q3_ess_bulk)),
      Q4_ess_bulk=as.integer(round(met_over$Q4_ess_bulk)),
      Q5_ess_bulk=as.integer(round(met_over$Q5_ess_bulk)),
      Q2_ess_tail=as.integer(round(met_over$Q2_ess_tail)),
      Q3_ess_tail=as.integer(round(met_over$Q3_ess_tail)),
      Q4_ess_tail=as.integer(round(met_over$Q4_ess_tail)),
      Q5_ess_tail=as.integer(round(met_over$Q5_ess_tail))
    )
    append_rows(outfile, row_over)
    message("[OK] ", scen_key, " MRD")

    # ---- MRR: Q1-Q5 vs Lag5 Q1 (universal reference for all lags) ----
    mrr_df <- compute_mrr(draws, des_overall$names, layer_dt, outcome, eqi_out, aamr_out, lagv, ref_draws=lag5_ref_draws)
    if (nrow(mrr_df) > 0) append_rows(mrr_file, mrr_df)
    message("[OK] ", scen_key, " MRR")

    # ---- Store Q5 beta draws for lag test ----
    q5_idx <- match("EQI_factor5", des_overall$names)
    if (!is.na(q5_idx))
      lag_q5_store[[as.character(lagv)]] <- draws[[paste0("beta[", q5_idx, "]")]]

  }  # end scenario loop

  # ---- Pairwise lag test (Q5): lag5 vs lag10, lag10 vs lag15, lag15 vs lag5 ----
  run_lag_test(lag_q5_store, outcome, out_dir)

  message("===== Completed: ", outcome, " =====")
}

message("All analyses complete. Output directory: ", out_dir)
