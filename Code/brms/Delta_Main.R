#!/usr/bin/env Rscript
# Interval-censored mixed model: EQI change → disease mortality change (national level).
# Outcome: delta AAMR (change between consecutive EQI periods, same lag).
# Predictor: EQI_Change_Category (Stable=reference, Improved, Worsened) per domain.
# Lags: 5, 10, 15 years (EQI0005 vs EQI0610, matched by lag).

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
utils::globalVariables(c(
  "State_FIPS", "delta_AAMR_lower", "delta_AAMR_upper", "delta_Smoking_Rate",
  "EQI_Change_Category", "Air_Change_Category", "Water_Change_Category",
  "Land_Change_Category", "Built_Change_Category", "Social_Change_Category"
))

option_list <- list(
  make_option(c("--data"), type = "character", default = "Data/Processed/df_Delta.csv"),
  make_option(c("--output-dir"), type = "character", default = "Result/Delta"),
  make_option(c("--ridgeline-dir"), type = "character", default = "Result/Delta_RL"),
  make_option(c("--outcomes"), type = "character", default = NA),
  make_option(c("--chains"), type = "integer", default = 6),
  make_option(c("--iter"), type = "integer", default = 1800),
  make_option(c("--warmup"), type = "integer", default = 1000),
  make_option(c("--adapt-delta"), type = "double", default = 0.95),
  make_option(c("--max-treedepth"), type = "integer", default = 12),
  make_option(c("--min-n"), type = "integer", default = 50),
  make_option(c("--seed"), type = "integer", default = 1234),
  make_option(c("--test"), action = "store_true", default = FALSE)
)
opt <- parse_args(OptionParser(option_list = option_list))

cores_avail <- parallel::detectCores(logical = TRUE)
slurm_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA)))
cores_used <- opt$chains
options(mc.cores = cores_used)

message("--- CPU Resource Report ---")
message("Environment: ", if (!is.na(slurm_cpus)) "Slurm (HPC)" else "Local Machine")
message("Total Cores Available: ", cores_avail)
message("Setting mc.cores to:   ", cores_used)
message("---------------------------")

set.seed(opt$seed)

# ── Stan model ───────────────────────────────────────────────────────────────
stan_code <- "data {\n  int<lower=1> N;\n  int<lower=1> S;\n  array[N] int<lower=1,upper=S> state;\n  vector[N] y_lower;\n  vector[N] y_upper;\n  array[N] int<lower=0,upper=2> cens;\n  int<lower=1> K;\n  matrix[N,K] X;\n} \nparameters {\n  vector[K] beta;\n  vector[S] z_u;\n  real<lower=0> sigma;\n  real<lower=0> sigma_u;\n} \ntransformed parameters {\n  vector[S] u = sigma_u * z_u;\n} \nmodel {\n  beta ~ normal(0,5);\n  z_u ~ normal(0,1);\n  sigma ~ exponential(1);\n  sigma_u ~ exponential(1);\n  for (i in 1:N) {\n    real mu = X[i] * beta + u[state[i]];\n    if (cens[i]==0) {\n      target += normal_lpdf(y_lower[i] | mu, sigma);\n    } else {\n      real p_up = normal_cdf(y_upper[i] | mu, sigma);\n      real p_lo = normal_cdf(y_lower[i] | mu, sigma);\n      real diff = fmax(p_up - p_lo, 1e-12);\n      target += log(diff);\n    }\n  }\n}"
stan_file <- file.path(tempdir(), "delta_main.stan")
writeLines(stan_code, stan_file)
mod <- cmdstan_model(stan_file)

# ── Data ─────────────────────────────────────────────────────────────────────
project_root <- normalizePath(".")
path <- file.path(project_root, opt$data)
if (!file.exists(path)) stop("Data not found: ", path)
dt <- fread(path)

req <- c(
  "COUNTY_FIPS", "State_FIPS", "Outcome", "Lag",
  "delta_AAMR_lower", "delta_AAMR_upper", "delta_Smoking_Rate",
  "EQI_Change_Category", "Air_Change_Category", "Water_Change_Category",
  "Land_Change_Category", "Built_Change_Category", "Social_Change_Category"
)
miss <- setdiff(req, names(dt))
if (length(miss)) stop("Missing cols: ", paste(miss, collapse = ","))

dt <- dt[!is.na(delta_AAMR_lower) & !is.na(delta_AAMR_upper)]
dt[, cens := ifelse(delta_AAMR_lower == delta_AAMR_upper, 0L, 2L)]

all_outcomes <- sort(unique(dt$Outcome))
selected <- if (is.na(opt$outcomes)) {
  all_outcomes
} else {
  reqc <- str_split(opt$outcomes, ",", simplify = TRUE) |>
    as.vector() |>
    str_trim()
  inv <- setdiff(reqc, all_outcomes)
  if (length(inv)) stop("Invalid outcomes: ", paste(inv, collapse = ","))
  reqc
}
message("Outcomes to analyze: ", paste(selected, collapse = ","))

out_dir <- file.path(project_root, opt$`output-dir`)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

rl_dir <- file.path(project_root, opt$`ridgeline-dir`)
if (!dir.exists(rl_dir)) dir.create(rl_dir, recursive = TRUE)

icd_to_name <- c(
  "I00_I99" = "CVD", "J40_J47_J60_J70_J84_D86_C34" = "CRD",
  "K70_K76_C22" = "CLD", "N00_N29_C64_C65" = "CKD",
  "X60_X84_Y87.0" = "Suicide", "G20_G30_G12.2_F01_F03" = "NDD",
  "C00_C97" = "Cancer"
)
disease_label <- function(icd) {
  nm <- icd_to_name[icd]
  if (is.na(nm)) icd else unname(nm)
}

# ── Helpers ───────────────────────────────────────────────────────────────────
sig_mark <- function(p) {
  if (is.na(p) || !is.numeric(p)) {
    return("")
  }
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
}

format_cell <- function(draws) {
  if (length(draws) == 0 || all(is.na(draws))) {
    return("")
  }
  ci <- quantile(draws, c(0.025, 0.975), na.rm = TRUE)
  p <- 2 * min(mean(draws > 0, na.rm = TRUE), mean(draws < 0, na.rm = TRUE))
  sprintf("%.2f(%.2f,%.2f)%s", mean(draws, na.rm = TRUE), ci[1], ci[2], sig_mark(p))
}

compute_p <- function(draws) {
  if (length(draws) == 0) {
    return(NA_character_)
  }
  pos <- sum(draws > 0, na.rm = TRUE)
  neg <- sum(draws < 0, na.rm = TRUE)
  n <- pos + neg
  if (n == 0) {
    return(NA_character_)
  }
  sprintf("%.4f", 2 * min((pos + 0.5) / (n + 1), (neg + 0.5) / (n + 1)))
}

append_rows <- function(path, df) {
  if (!file.exists(path)) {
    write_csv(df, path)
  } else {
    suppressWarnings(write.table(df, path,
      sep = ",", col.names = FALSE,
      row.names = FALSE, append = TRUE
    ))
  }
}

# ── Build RDS object for ridgeline plotting ────────────────────────────────────
build_ridge_rds <- function(draws_df, design_names, summ, pattern, model_label,
                            outcome, lag, n_obs, n_states) {
  # Extract I and W factor levels
  i_idx <- grep(paste0(pattern, "I$"), design_names)
  w_idx <- grep(paste0(pattern, "W$"), design_names)

  if (length(i_idx) == 0 || length(w_idx) == 0) {
    return(NULL)
  }

  i_draws <- draws_df[[paste0("beta[", i_idx[1], "]")]]
  w_draws <- draws_df[[paste0("beta[", w_idx[1], "]")]]

  if (length(i_draws) == 0 || length(w_draws) == 0) {
    return(NULL)
  }

  n_draws <- length(i_draws)

  draws_long <- data.frame(
    draw_id = rep(seq_len(n_draws), 2),
    category = rep(c("I", "W"), each = n_draws),
    effect = c(i_draws, w_draws)
  )

  summary_df <- do.call(rbind, lapply(list(
    list(idx = i_idx[1], cat = "I"),
    list(idx = w_idx[1], cat = "W")
  ), function(x) {
    sr <- summ[summ$variable == paste0("beta[", x$idx, "]"), , drop = FALSE]
    draws <- if (x$cat == "I") i_draws else w_draws
    data.frame(
      category = x$cat,
      mean = mean(draws, na.rm = TRUE),
      sd = sd(draws, na.rm = TRUE),
      q025 = quantile(draws, 0.025, na.rm = TRUE),
      q975 = quantile(draws, 0.975, na.rm = TRUE),
      rhat = if (nrow(sr)) sr$rhat else NA_real_,
      ess_bulk = if (nrow(sr)) sr$ess_bulk else NA_real_,
      ess_tail = if (nrow(sr)) sr$ess_tail else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  list(
    metadata = list(
      outcome = outcome,
      lag = lag,
      model_type = model_label,
      n_obs = n_obs,
      n_states = n_states,
      n_draws = n_draws,
      timestamp = Sys.time()
    ),
    draws_wide = data.frame(draw_id = seq_len(n_draws), I = i_draws, W = w_draws),
    draws_long = draws_long,
    summary = summary_df
  )
}


# ── Design & fitting ──────────────────────────────────────────────────────────
build_design <- function(d, cat_col) {
  original_n <- nrow(d)
  d[[cat_col]] <- factor(d[[cat_col]], levels = c("S", "I", "W"))
  keep <- c(
    "delta_AAMR_lower", "delta_AAMR_upper", "cens", "State_FIPS",
    "delta_Smoking_Rate", cat_col
  )
  d <- d[complete.cases(d[, keep, with = FALSE])]
  
  if (nrow(d) == 0) {
    warning("build_design: No complete cases after filtering. Original n=", original_n,
            ". Check if ", cat_col, " has unexpected values.")
    return(NULL)
  }
  
  mm <- model.matrix(as.formula(paste0("~ delta_Smoking_Rate + `", cat_col, "`")), d)
  colnames(mm) <- make.names(colnames(mm))
  list(X = mm, names = colnames(mm), df = d)
}

extract_effects <- function(fit, names_vec) {
  draws <- as_draws_df(fit$draws("beta"))
  colnames(draws) <- paste0("beta[", seq_len(ncol(draws)), "]")
  summ <- posterior::summarize_draws(fit$draws("beta"))
  get_draws <- function(pattern) {
    idx <- grep(pattern, names_vec)
    if (!length(idx)) {
      return(numeric(0))
    }
    draws[[paste0("beta[", idx[1], "]")]]
  }
  get_diag <- function(pattern) {
    idx <- grep(pattern, names_vec)
    if (!length(idx)) {
      return(list(rhat = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_))
    }
    sr <- summ[summ$variable == paste0("beta[", idx[1], "]"), , drop = FALSE]
    if (!nrow(sr)) {
      return(list(rhat = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_))
    }
    list(rhat = sr$rhat, ess_bulk = sr$ess_bulk, ess_tail = sr$ess_tail)
  }
  list(
    imp = get_draws("Improved"), wor = get_draws("Worsened"),
    inter = get_draws("Intercept"), sm = get_draws("delta_Smoking_Rate"),
    imp_d = get_diag("Improved"), wor_d = get_diag("Worsened"),
    rhat_max = max(summ$rhat, na.rm = TRUE)
  )
}

make_row <- function(outcome, lag, model_name, eff, n) {
  tibble(
    ICD_Code = outcome, Lag = lag, Model = model_name,
    MRD_Improved = format_cell(eff$imp),
    MRD_Worsened = format_cell(eff$wor),
    Intercept = format_cell(eff$inter),
    SM = format_cell(eff$sm),
    Improved_p = compute_p(eff$imp),
    Worsened_p = compute_p(eff$wor),
    Improved_rhat = sprintf("%.4f", eff$imp_d$rhat),
    Worsened_rhat = sprintf("%.4f", eff$wor_d$rhat),
    Improved_ess_bulk = as.integer(round(eff$imp_d$ess_bulk)),
    Worsened_ess_bulk = as.integer(round(eff$wor_d$ess_bulk)),
    Improved_ess_tail = as.integer(round(eff$imp_d$ess_tail)),
    Worsened_ess_tail = as.integer(round(eff$wor_d$ess_tail)),
    Rhat_max = round(eff$rhat_max, 4), N = n
  )
}

fit_one <- function(d, cat_col) {
  des <- build_design(d, cat_col)
  if (is.null(des) || nrow(des$df) < opt$`min-n`) {
    return(NULL)
  }
  states <- sort(unique(des$df$State_FIPS))
  state_idx <- match(des$df$State_FIPS, states)
  data_list <- list(
    N = nrow(des$df), S = length(states), state = state_idx,
    y_lower = des$df$delta_AAMR_lower, y_upper = des$df$delta_AAMR_upper,
    cens = des$df$cens, K = ncol(des$X), X = des$X
  )
  init_fn <- function() {
    list(
      beta = rep(0, data_list$K), z_u = rep(0, data_list$S),
      sigma = 5, sigma_u = 2
    )
  }
  fit <- try(mod$sample(
    data = data_list, chains = opt$chains,
    iter_sampling = opt$iter - opt$warmup, iter_warmup = opt$warmup,
    adapt_delta = opt$`adapt-delta`, max_treedepth = opt$`max-treedepth`,
    parallel_chains = opt$chains, refresh = 0, seed = opt$seed,
    init = rep(list(init_fn()), opt$chains)
  ), silent = FALSE)
  if (inherits(fit, "try-error")) {
    message("    [Stan Error] ", attr(fit, "condition")$message)
    return(NULL)
  }
  list(fit = fit, des = des)
}

domains <- list(
  list(col = "EQI_Change_Category", name = "EQI"),
  list(col = "Air_Change_Category", name = "Air"),
  list(col = "Water_Change_Category", name = "Water"),
  list(col = "Land_Change_Category", name = "Land"),
  list(col = "Built_Change_Category", name = "Built"),
  list(col = "Social_Change_Category", name = "Social")
)

# ── Main loop ─────────────────────────────────────────────────────────────────
for (outcome in selected) {
  dlabel <- disease_label(outcome)
  outfile <- file.path(out_dir, paste0(dlabel, "_Main.csv"))
  message("===== Outcome: ", dlabel, " (", outcome, ") =====")

  out_dt <- dt[Outcome == outcome]
  lag_vals <- sort(unique(out_dt$Lag))

  for (lag_v in lag_vals) {
    lag_dt <- out_dt[Lag == lag_v]
    message("  Lag=", lag_v, "  n=", nrow(lag_dt))
    if (nrow(lag_dt) < opt$`min-n`) {
      message("  [Skip]")
      next
    }

    for (dom in domains) {
      res <- fit_one(lag_dt, dom$col)
      if (is.null(res)) {
        message("  [Skip/Fail] ", dom$name, " Lag=", lag_v)
        next
      }
      eff <- extract_effects(res$fit, res$des$names)
      append_rows(outfile, make_row(outcome, lag_v, dom$name, eff, nrow(res$des$df)))

      # Extract and save ridgeline data
      rds <- build_ridge_rds(
        as_draws_df(res$fit$draws("beta")),
        res$des$names,
        posterior::summarize_draws(res$fit$draws("beta")),
        dom$col,
        dom$name,
        outcome, lag_v, nrow(res$des$df),
        length(sort(unique(res$des$df$State_FIPS)))
      )
      if (!is.null(rds)) {
        rl_file <- file.path(rl_dir, sprintf("%s_Lag%d_%s_Main.rds", dlabel, lag_v, dom$name))
        saveRDS(rds, rl_file)
        message("       → RL: ", basename(rl_file))
      }

      message("  [OK] ", dom$name, " Lag=", lag_v, "  Rhat_max=", round(eff$rhat_max, 3))
    }
  }
  message("===== Completed: ", outcome, " =====")
}
message("All analyses complete. Output: ", out_dir, " | Ridgeline: ", rl_dir)
