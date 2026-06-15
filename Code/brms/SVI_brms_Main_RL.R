#!/usr/bin/env Rscript
# SVI ridgeline posterior extraction pipeline.
#   - Exposure: SVI category A/B/C/D (A = reference); extracts B/C/D posterior draws.
#   - Disease-specific covariate adjustment (see DISEASE_COVSET); category resolved
#     from a hardcoded ICD->category map mirroring config.yaml (no yaml dependency).
#     Cancer takes precedence for cancer (sub)types.
#   - SVI is static -> four AAMR periods treated as four lags (5/10/15/20 yr).
#   - Single SVI model (no multi-domain).
# Output: Result/brms_SVI_Main_RL/{shortname}_Lag{N}_SVI.rds

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(cmdstanr)
  library(posterior)
})

utils::globalVariables(c(
  "SVI", "State_FIPS",
  "Smoking_rate", "Physical_Activities_rate", "Obesity_rate",
  "Uninsured_rate", "Physician_Density_per100k", "Diabetes_Prevalence_rate"
))

# ── Disease-specific covariate adjustment sets (config.yaml category keys) ───────
DISEASE_COVSET <- list(
  cancer      = c("Smoking_rate", "Uninsured_rate"),
  liver       = c("Smoking_rate", "Uninsured_rate"),
  cvd         = c("Smoking_rate", "Uninsured_rate"),
  kidney      = c("Smoking_rate", "Uninsured_rate"),
  ndd         = c("Physical_Activities_rate", "Uninsured_rate"),
  suicide     = c("Uninsured_rate", "Physician_Density_per100k", "Diabetes_Prevalence_rate"),
  respiratory = c("Uninsured_rate")
)
DISEASE_CODES <- list(
  liver       = c("K70_K76_C22", "K70_K76", "K70", "K71", "K73", "K74",
                  "K71_K73_K74", "K76", "K76.7", "C22"),
  respiratory = c("J40_J47_J60_J70_J84_D86_C34", "J40_J47_J60_J70_J84_D86",
                  "J43_J44", "J45", "J84_D86", "J60_J66", "C34"),
  kidney      = c("N00_N29_C64_C65", "N00_N29", "N18_N19", "N00_N15", "C64_C65"),
  cvd         = c("I00_I99", "I20_I25", "I60_I69", "I10_I15", "I50"),
  suicide     = c("X60_X84_Y87.0", "X60_X69", "X70_X84", "Y87.0"),
  ndd         = c("G20_G30_G12.2_F01_F03", "G30_F01_F03", "G20", "G10", "G12.2"),
  cancer      = c("C00_C97", "C18_C21", "C22", "C25", "C34", "C50", "C56",
                  "C61", "C64_C65", "C82_C85", "C91_C95")
)
icd_to_cats <- list()
for (cat in names(DISEASE_CODES)) {
  for (code in DISEASE_CODES[[cat]]) {
    icd_to_cats[[code]] <- unique(c(icd_to_cats[[code]], cat))
  }
}
resolve_covset <- function(outcome) {
  cats <- icd_to_cats[[outcome]]
  if (is.null(cats)) return(list(ok = FALSE, reason = paste0("no category for '", outcome, "'")))
  if ("cancer" %in% cats && "cancer" %in% names(DISEASE_COVSET)) {
    return(list(ok = TRUE, covset = DISEASE_COVSET[["cancer"]], cat = "cancer"))
  }
  defined <- cats[cats %in% names(DISEASE_COVSET)]
  if (length(defined) == 0) return(list(ok = FALSE, reason = paste0("no covariate set for: ", paste(cats, collapse = "/"))))
  uniq <- unique(lapply(defined, function(k) sort(DISEASE_COVSET[[k]])))
  if (length(uniq) > 1) return(list(ok = FALSE, reason = paste0("ambiguous covariate sets: ", paste(defined, collapse = "/"))))
  list(ok = TRUE, covset = DISEASE_COVSET[[defined[1]]], cat = defined[1])
}

option_list <- list(
  make_option(c("--data"), type = "character", default = "Data/Processed/df_SVI.csv"),
  make_option(c("--output-dir"), type = "character", default = "Result/brms_SVI_Main_RL"),
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
stan_file <- file.path(tempdir(), "interval_mixed_svi_rl.stan")
writeLines(stan_code, stan_file)
message("Compiling Stan model...")
mod <- cmdstan_model(stan_file, quiet = TRUE)
message("Model compiled.")

# ── Main diseases (overall categories) ─────────────────────────────────────────
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
  "COUNTY_FIPS", "Time_Period", "Outcome", "AAMR_Lower", "AAMR_Upper", "SVI",
  "Smoking_rate", "Physical_Activities_rate", "Obesity_rate",
  "Uninsured_rate", "Physician_Density_per100k", "Diabetes_Prevalence_rate"
)
miss <- setdiff(req, names(dt))
if (length(miss)) stop("Missing cols: ", paste(miss, collapse = ","))

if (!"State_FIPS" %in% names(dt)) dt[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]
dt <- dt[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
dt[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]

# ── SVI static -> 4 lag scenarios (filter by AAMR Time_Period only) ───────────
scenario_list <- list(
  list(key = "SVI_AAMR2006_2010", aamr = "2006-2010", lag = 5),
  list(key = "SVI_AAMR2011_2015", aamr = "2011-2015", lag = 10),
  list(key = "SVI_AAMR2016_2020", aamr = "2016-2020", lag = 15),
  list(key = "SVI_AAMR2021_2024", aamr = "2021-2024", lag = 20)
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

# ── Design builder: intercept + disease covariates + SVI_factor (A ref) ─────────
build_design <- function(d, covariates) {
  d <- as.data.frame(d)
  d$SVI_factor <- factor(d$SVI, levels = c("A", "B", "C", "D"))
  needed <- c("SVI_factor", "AAMR_Lower", "AAMR_Upper", "cens", "State_FIPS", covariates)
  d <- d[stats::complete.cases(d[, needed, drop = FALSE]), , drop = FALSE]
  form <- if (length(covariates) == 0) {
    ~SVI_factor
  } else {
    as.formula(paste("~", paste(c(covariates, "SVI_factor"), collapse = " + ")))
  }
  mm <- model.matrix(form, d, contrasts.arg = list(SVI_factor = contr.treatment(c("A", "B", "C", "D"))))
  colnames(mm) <- make.names(colnames(mm))
  list(X = mm, names = colnames(mm), df = d)
}

# ── Build RDS object for ridgeline plotting (B/C/D vs A posterior draws) ────────
build_ridge_rds <- function(draws_df, design_names, summ, model_label,
                            outcome, aamr_p, lagv, n_obs, n_states) {
  levels_vec <- c("B", "C", "D")
  lv_list <- lapply(levels_vec, function(lv) {
    nm <- paste0("SVI_factor", lv)
    idx <- which(design_names == nm)
    if (length(idx) == 0) return(NULL)
    list(draws = draws_df[[paste0("beta[", idx, "]")]], idx = idx)
  })
  if (any(sapply(lv_list, is.null))) {
    warning("Could not find all SVI category parameters for ", outcome)
    return(NULL)
  }
  n_draws <- length(lv_list[[1]]$draws)

  draws_wide <- data.frame(
    draw_id = seq_len(n_draws),
    B = lv_list[[1]]$draws, C = lv_list[[2]]$draws, D = lv_list[[3]]$draws
  )
  draws_long <- data.frame(
    draw_id  = rep(seq_len(n_draws), 3),
    category = factor(rep(levels_vec, each = n_draws), levels = levels_vec),
    effect   = c(lv_list[[1]]$draws, lv_list[[2]]$draws, lv_list[[3]]$draws)
  )
  summary_df <- do.call(rbind, lapply(seq_along(lv_list), function(i) {
    d <- lv_list[[i]]$draws
    col <- paste0("beta[", lv_list[[i]]$idx, "]")
    sr <- summ[summ$variable == col, , drop = FALSE]
    data.frame(
      category = levels_vec[i],
      mean = mean(d, na.rm = TRUE), sd = sd(d, na.rm = TRUE),
      q025 = quantile(d, 0.025, na.rm = TRUE), q975 = quantile(d, 0.975, na.rm = TRUE),
      rhat = if (nrow(sr)) sr$rhat else NA_real_,
      ess_bulk = if (nrow(sr)) sr$ess_bulk else NA_real_,
      ess_tail = if (nrow(sr)) sr$ess_tail else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  list(
    metadata = list(
      cancer_type = outcome,
      aamr_period = aamr_p,
      lag         = lagv,
      model_type  = model_label,
      n_obs       = n_obs,
      n_states    = n_states,
      n_draws     = n_draws,
      timestamp   = Sys.time(),
      convergence = list(max_rhat = max(summary_df$rhat, na.rm = TRUE),
                         min_ess_bulk = min(summary_df$ess_bulk, na.rm = TRUE))
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
    data = dl, chains = opt$chains,
    iter_sampling = opt$iter - opt$warmup, iter_warmup = opt$warmup,
    adapt_delta = opt$`adapt-delta`, max_treedepth = opt$`max-treedepth`,
    parallel_chains = opt$chains, refresh = 0, seed = opt$seed,
    init = rep(list(init_fn()), opt$chains)
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
  rc <- resolve_covset(outcome)
  if (!isTRUE(rc$ok)) {
    message("[Skip] Outcome ", outcome, " -> ", rc$reason)
    next
  }
  covset <- rc$covset
  model_label <- paste(c("SVI", names(covset)), collapse = "+")
  sn <- get_shortname(outcome)
  message("===== Outcome: ", outcome, "  [", rc$cat, "]  covset: ",
          paste(covset, collapse = "+"), " =====")

  for (sc in scenario_list) {
    scen_key <- sc$key
    aamr_p <- sc$aamr
    lagv <- sc$lag

    scen_dt <- dt[Time_Period == aamr_p & Outcome == outcome]
    if (nrow(scen_dt) < opt$`min-n`) {
      message("[Skip] ", scen_key, " n=", nrow(scen_dt))
      next
    }

    des <- build_design(scen_dt, covset)
    if (nrow(des$df) < opt$`min-n`) {
      message("[Skip] ", scen_key, " after design n=", nrow(des$df))
      next
    }
    res <- run_stan(des, scen_key, "SVI")
    if (is.null(res)) next

    rds <- build_ridge_rds(
      res$draws, des$names, res$summ, "SVI",
      outcome, aamr_p, lagv, nrow(des$df), res$n_states
    )
    if (!is.null(rds)) {
      out_file <- file.path(out_dir, sprintf("%s_Lag%d_SVI.rds", sn, lagv))
      saveRDS(rds, out_file)
      message("[OK] ", scen_key, " SVI -> ", basename(out_file))
    }
  }
  message("===== Completed: ", outcome, " =====")
}
message("All analyses complete. Output directory: ", out_dir)
