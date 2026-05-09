#!/usr/bin/env Rscript
# Ridgeline posterior extraction pipeline (new data structure)
# Scenarios: 2000-2005 EQI only, lags 5/10/15/20
# Models: EQI+SM+PA+OB (Overall) and domain-specific (Air/Water/Land/Built/Social) + SM+PA+OB
# Output: Result/brms_RL/{shortname}_Lag{N}_{Overall/Air/Water/Land/Built/Social}.rds

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(cmdstanr)
  library(posterior)
})

utils::globalVariables(c(
  "EQI", "EQI_Air", "EQI_Water", "EQI_Land", "EQI_Built", "EQI_Social",
  "State_FIPS", "Smoking_rate", "Physical_Activities_rate", "Obesity_rate"
))

option_list <- list(
  make_option(c("--data"), type = "character", default = "Data/Processed/df.csv"),
  make_option(c("--output-dir"), type = "character", default = "Result/brms_RL"),
  make_option(c("--outcomes"), type = "character", default = NA, help = "Comma-separated ICD codes"),
  make_option(c("--chains"), type = "integer", default = 16),
  make_option(c("--iter"), type = "integer", default = 1500),
  make_option(c("--warmup"), type = "integer", default = 1000),
  make_option(c("--adapt-delta"), type = "double", default = 0.95),
  make_option(c("--max-treedepth"), type = "integer", default = 12),
  make_option(c("--min-n"), type = "integer", default = 50),
  make_option(c("--seed"), type = "integer", default = 1234),
  make_option(c("--test"), action = "store_true", default = FALSE)
)
opt <- parse_args(OptionParser(option_list = option_list))


# ── Stan model ─────────────────────────────────────────────────────────────────
stan_code <- "data {
  int<lower=1> N;
  int<lower=1> S;
  array[N] int<lower=1,upper=S> state;
  vector[N] y_lower;
  vector[N] y_upper;
  array[N] int<lower=0,upper=2> cens;
  int<lower=1> K;
  matrix[N,K] X;
}
parameters {
  vector[K] beta;
  vector[S] z_u;
  real<lower=0> sigma;
  real<lower=0> sigma_u;
}
transformed parameters {
  vector[S] u = sigma_u * z_u;
}
model {
  beta ~ normal(0,5);
  z_u ~ normal(0,1);
  sigma ~ exponential(1);
  sigma_u ~ exponential(1);
  for (i in 1:N) {
    real mu = X[i] * beta + u[state[i]];
    if (cens[i]==0) {
      target += normal_lpdf(y_lower[i] | mu, sigma);
    } else {
      real p_up = normal_cdf(y_upper[i] | mu, sigma);
      real p_lo = normal_cdf(y_lower[i] | mu, sigma);
      real diff = fmax(p_up - p_lo, 1e-12);
      target += log(diff);
    }
  }
}"
stan_file <- file.path(tempdir(), "interval_mixed_rl.stan")
writeLines(stan_code, stan_file)
message("Compiling Stan model...")
mod <- cmdstan_model(stan_file, quiet = TRUE)
message("Model compiled.")

# ── Main diseases only (same as brms_Stratified_Demo.R) ───────────────────────
icd_to_name <- c(
  "I00_I99" = "CVD",
  "J40_J47_J60_J70_J84_D86_C34" = "CRD",
  "K70_K76_C22" = "CLD",
  "N00_N29_C64_C65" = "CKD",
  "X60_X84_Y87.0" = "Suicide",
  "G20_G30_G12.2_F01_F03" = "NDD",
  "C00_C97" = "Cancer"
)

get_shortname <- function(outcome) {
  nm <- icd_to_name[outcome]
  if (!is.na(nm)) unname(nm) else outcome
}

# ── Load data ──────────────────────────────────────────────────────────────────
project_root <- normalizePath(".")
path <- file.path(project_root, opt$data)
if (!file.exists(path)) stop("Data not found: ", path)
dt <- fread(path)

req <- c(
  "COUNTY_FIPS", "EQI_Period", "Time_Period", "Lag_Years", "Outcome",
  "AAMR_Lower", "AAMR_Upper",
  "EQI", "EQI_Air", "EQI_Water", "EQI_Land", "EQI_Built", "EQI_Social",
  "Smoking_rate", "Physical_Activities_rate", "Obesity_rate"
)
miss <- setdiff(req, names(dt))
if (length(miss)) stop("Missing cols: ", paste(miss, collapse = ","))

if (!"State_FIPS" %in% names(dt)) dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]
dt <- dt[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
dt[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]

# ── Scenarios: 2000-2005 EQI only, lags 5/10/15/20 ───────────────────────────
scenario_list <- list(
  list(key = "EQI0005_AAMR2006_2010", eqi = "2000-2005", aamr = "2006-2010", lag = 5),
  list(key = "EQI0005_AAMR2011_2015", eqi = "2000-2005", aamr = "2011-2015", lag = 10),
  list(key = "EQI0005_AAMR2016_2020", eqi = "2000-2005", aamr = "2016-2020", lag = 15),
  list(key = "EQI0005_AAMR2021_2024", eqi = "2000-2005", aamr = "2021-2024", lag = 20)
)

all_outcomes <- sort(unique(dt$Outcome))
selected <- if (is.na(opt$`outcomes`)) {
  intersect(names(icd_to_name), all_outcomes)
} else {
  reqc <- str_split(opt$`outcomes`, ",", simplify = TRUE) |>
    as.vector() |>
    str_trim()
  inv <- setdiff(reqc, all_outcomes)
  if (length(inv)) stop("Invalid outcomes: ", paste(inv, collapse = ","))
  reqc
}
message("Outcomes to analyze: ", paste(selected, collapse = ","))

out_dir <- file.path(project_root, opt$`output-dir`)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── Design matrix builders (identical to brms_Main.R) ─────────────────────────
build_design_overall <- function(d) {
  d <- d %>% mutate(EQI_factor = factor(EQI, levels = 1:5))
  d <- d[complete.cases(d[, c(
    "EQI_factor", "AAMR_Lower", "AAMR_Upper", "cens", "State_FIPS",
    "Smoking_rate", "Physical_Activities_rate", "Obesity_rate"
  )]), ]
  mm <- model.matrix(
    ~ Smoking_rate + Physical_Activities_rate + Obesity_rate + EQI_factor, d,
    contrasts.arg = list(EQI_factor = contr.treatment(5))
  )
  colnames(mm) <- make.names(colnames(mm))
  list(X = mm, names = colnames(mm), df = d)
}

build_design_multi <- function(d) {
  d <- d %>% mutate(
    EQI_Air_factor = factor(EQI_Air, levels = 1:5),
    EQI_Water_factor = factor(EQI_Water, levels = 1:5),
    EQI_Land_factor = factor(EQI_Land, levels = 1:5),
    EQI_Built_factor = factor(EQI_Built, levels = 1:5),
    EQI_Social_factor = factor(EQI_Social, levels = 1:5)
  )
  d <- d[complete.cases(d[, c(
    "EQI_Air_factor", "EQI_Water_factor", "EQI_Land_factor",
    "EQI_Built_factor", "EQI_Social_factor",
    "AAMR_Lower", "AAMR_Upper", "cens", "State_FIPS",
    "Smoking_rate", "Physical_Activities_rate", "Obesity_rate"
  )]), ]
  form <- as.formula("~ Smoking_rate + Physical_Activities_rate + Obesity_rate +
    EQI_Air_factor + EQI_Water_factor + EQI_Land_factor + EQI_Built_factor + EQI_Social_factor")
  mm <- model.matrix(form, d,
    contrasts.arg = list(
      EQI_Air_factor = contr.treatment(5), EQI_Water_factor = contr.treatment(5),
      EQI_Land_factor = contr.treatment(5), EQI_Built_factor = contr.treatment(5),
      EQI_Social_factor = contr.treatment(5)
    )
  )
  colnames(mm) <- make.names(colnames(mm))
  list(X = mm, names = colnames(mm), df = d)
}

# ── Build RDS object for ridgeline plotting ────────────────────────────────────
# Extracts Q2-Q5 posterior draws for a given design prefix and packages them
# into the standard ridgeline format consumed by check scripts and plot code.
build_ridge_rds <- function(draws_df, design_names, summ, prefix, model_label,
                            outcome, eqi_p, aamr_p, lagv, n_obs, n_states) {
  q_list <- lapply(2:5, function(q) {
    nm <- paste0(prefix, q)
    idx <- which(design_names == nm)
    if (length(idx) == 0) {
      return(NULL)
    }
    list(draws = draws_df[[paste0("beta[", idx, "]")]], idx = idx)
  })

  if (any(sapply(q_list, is.null))) {
    warning("Could not find all quintile parameters for prefix: ", prefix)
    return(NULL)
  }

  n_draws <- length(q_list[[1]]$draws)

  draws_wide <- data.frame(
    draw_id = seq_len(n_draws),
    Q2 = q_list[[1]]$draws, Q3 = q_list[[2]]$draws,
    Q4 = q_list[[3]]$draws, Q5 = q_list[[4]]$draws
  )

  draws_long <- data.frame(
    draw_id = rep(seq_len(n_draws), 4),
    quintile = factor(
      rep(c("Q2", "Q3", "Q4", "Q5"), each = n_draws),
      levels = c("Q2", "Q3", "Q4", "Q5")
    ),
    effect = c(
      q_list[[1]]$draws, q_list[[2]]$draws,
      q_list[[3]]$draws, q_list[[4]]$draws
    )
  )

  summary_df <- do.call(rbind, lapply(seq_along(q_list), function(i) {
    q <- i + 1
    d <- q_list[[i]]$draws
    col <- paste0("beta[", q_list[[i]]$idx, "]")
    sr <- summ[summ$variable == col, , drop = FALSE]
    data.frame(
      quintile = paste0("Q", q),
      mean = mean(d, na.rm = TRUE),
      sd = sd(d, na.rm = TRUE),
      q025 = quantile(d, 0.025, na.rm = TRUE),
      q975 = quantile(d, 0.975, na.rm = TRUE),
      rhat = if (nrow(sr)) sr$rhat else NA_real_,
      ess_bulk = if (nrow(sr)) sr$ess_bulk else NA_real_,
      ess_tail = if (nrow(sr)) sr$ess_tail else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  max_rhat <- max(summary_df$rhat, na.rm = TRUE)
  min_ess <- min(summary_df$ess_bulk, na.rm = TRUE)

  list(
    metadata = list(
      cancer_type = outcome,
      eqi_period  = eqi_p,
      aamr_period = aamr_p,
      lag         = lagv,
      model_type  = model_label,
      n_obs       = n_obs,
      n_states    = n_states,
      n_draws     = n_draws,
      timestamp   = Sys.time(),
      convergence = list(max_rhat = max_rhat, min_ess_bulk = min_ess)
    ),
    draws_wide = draws_wide,
    draws_long = draws_long,
    summary = summary_df
  )
}

# ── Helper: run Stan and return draws + summarize_draws ───────────────────────
run_stan <- function(des, scen_key, label) {
  states <- sort(unique(des$df$State_FIPS))
  si <- match(des$df$State_FIPS, states)
  dl <- list(
    N = nrow(des$df), S = length(states), state = si,
    y_lower = des$df$AAMR_Lower, y_upper = des$df$AAMR_Upper, cens = des$df$cens,
    K = ncol(des$X), X = des$X
  )
  init_fn <- function() list(beta = rep(0, dl$K), z_u = rep(0, dl$S), sigma = 50, sigma_u = 10)

  fit <- try(mod$sample(
    data            = dl,
    chains          = opt$chains,
    iter_sampling   = opt$iter - opt$warmup,
    iter_warmup     = opt$warmup,
    adapt_delta     = opt$`adapt-delta`,
    max_treedepth   = opt$`max-treedepth`,
    parallel_chains = min(opt$chains, cores_used),
    refresh         = 0,
    seed            = opt$seed,
    init            = rep(list(init_fn()), opt$chains)
  ), silent = TRUE)

  if (inherits(fit, "try-error")) {
    message("[Fail] ", scen_key, " ", label)
    return(NULL)
  }

  draws <- as_draws_df(fit$draws("beta"))
  colnames(draws) <- paste0("beta[", seq_len(ncol(draws)), "]")
  summ <- posterior::summarize_draws(fit$draws("beta"))
  list(draws = draws, summ = summ, n_states = length(states))
}

# ── Main loop ──────────────────────────────────────────────────────────────────
for (outcome in selected) {
  message("===== Outcome: ", outcome, " =====")
  sn <- get_shortname(outcome)

  for (sc in scenario_list) {
    scen_key <- sc$key
    eqi_p <- sc$eqi
    aamr_p <- sc$aamr
    lagv <- sc$lag

    scen_dt <- dt[EQI_Period == eqi_p & Time_Period == aamr_p & Outcome == outcome]
    if (nrow(scen_dt) < opt$`min-n`) {
      message("[Skip] ", scen_key, " n=", nrow(scen_dt))
      next
    }

    # ── Overall EQI + SM + PA + OB ──────────────────────────────────────────
    des_o <- build_design_overall(scen_dt)
    if (nrow(des_o$df) < opt$`min-n`) {
      message("[Skip] ", scen_key, " Overall n=", nrow(des_o$df))
    } else {
      res_o <- run_stan(des_o, scen_key, "Overall")
      if (!is.null(res_o)) {
        rds_o <- build_ridge_rds(
          res_o$draws, des_o$names, res_o$summ, "EQI_factor", "overall",
          outcome, eqi_p, aamr_p, lagv, nrow(des_o$df), res_o$n_states
        )
        if (!is.null(rds_o)) {
          out_file <- file.path(out_dir, sprintf("%s_Lag%d_Overall.rds", sn, lagv))
          saveRDS(rds_o, out_file)
          message("[OK] ", scen_key, " Overall → ", basename(out_file))
        }
      }
    }

    # ── Multi-domain (Air/Water/Land/Built/Social) + SM + PA + OB ───────────
    des_m <- build_design_multi(scen_dt)
    if (nrow(des_m$df) < opt$`min-n`) {
      message("[Skip] ", scen_key, " Multi-domain n=", nrow(des_m$df))
    } else {
      res_m <- run_stan(des_m, scen_key, "Multi-domain")
      if (!is.null(res_m)) {
        domain_map <- list(
          Air    = "EQI_Air_factor",
          Water  = "EQI_Water_factor",
          Land   = "EQI_Land_factor",
          Built  = "EQI_Built_factor",
          Social = "EQI_Social_factor"
        )
        for (dom_name in names(domain_map)) {
          rds_d <- build_ridge_rds(
            res_m$draws, des_m$names, res_m$summ, domain_map[[dom_name]],
            paste0("domain_", tolower(dom_name)),
            outcome, eqi_p, aamr_p, lagv, nrow(des_m$df), res_m$n_states
          )
          if (!is.null(rds_d)) {
            out_file <- file.path(out_dir, sprintf("%s_Lag%d_%s.rds", sn, lagv, dom_name))
            saveRDS(rds_d, out_file)
            message("[OK] ", scen_key, " ", dom_name, " → ", basename(out_file))
          }
        }
      }
    }
  }
  message("===== Completed: ", outcome, " =====")
}
message("All analyses complete. Output directory: ", out_dir)
