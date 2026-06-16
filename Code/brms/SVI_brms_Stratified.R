#!/usr/bin/env Rscript
# Comprehensive SVI stratified pipeline — runs all four stratification schemes
# in order (Cluster -> Demo -> Geo -> RUCC) in one script.
#   - Exposure: SVI category A/B/C/D (A = reference); output B/C/D MRD vs A and
#     within-stratum MRR (each scenario's adjusted category-A mean = reference).
#   - Disease-specific covariate adjustment (see DISEASE_COVSET); category resolved
#     from a hardcoded ICD->category map mirroring config.yaml (no yaml dependency).
#     Cancer takes precedence for cancer (sub)types.
#   - SVI is static -> four AAMR periods treated as four lags (5/10/15/20 yr).
#   - Single SVI index (no multi-domain).
#
# Data per phase (a phase is skipped with a warning if its file / strata columns
# are absent):
#   Cluster : Data/Processed/df_SVI.csv          strata Cluster_EQI, Cluster_NLCD
#   Demo    : Data/Processed/df_SVI_Stratified.csv  strata Stratum (sex/race)
#   Geo     : Data/Processed/df_SVI.csv          strata Census_Region, Climate_Zone,
#                                                       Economic_type, Homeownership_tertile
#   RUCC    : Data/Processed/df_SVI.csv          strata RUCC
#
# Output (same style as the EQI stratified scripts):
#   Result/brms_SVI_Stratified_{Cluster,Demo,Geo,RUCC}/{dlabel}_Stratified_{Tag}_{MRD,MRR}.csv

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
  "SVI", "State_FIPS", "Outcome", "Time_Period", "AAMR_Lower", "AAMR_Upper",
  "Smoking_rate", "Physical_Activities_rate", "Obesity_rate",
  "Uninsured_rate", "Physician_Density_per100k", "Diabetes_Prevalence_rate"
))

# ── Disease-specific covariate adjustment sets ──────────────────────────────────
DISEASE_COVSET <- list(
  cancer      = c("Smoking_rate", "Uninsured_rate"),
  liver       = c("Smoking_rate", "Uninsured_rate"),
  cvd         = c("Smoking_rate", "Uninsured_rate"),
  kidney      = c("Smoking_rate", "Uninsured_rate"),
  ndd         = c("Physical_Activities_rate", "Uninsured_rate"),
  suicide     = c("Uninsured_rate", "Physician_Density_per100k", "Diabetes_Prevalence_rate"),
  respiratory = c("Uninsured_rate")
)
covar_abbrev <- c(
  Smoking_rate = "SM", Physical_Activities_rate = "PA", Obesity_rate = "OB",
  Uninsured_rate = "UN", Physician_Density_per100k = "PD", Diabetes_Prevalence_rate = "DB"
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
  for (code in DISEASE_CODES[[cat]]) icd_to_cats[[code]] <- unique(c(icd_to_cats[[code]], cat))
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
  make_option(c("--output-base"), type = "character", default = "Result"),
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
stan_code <- "data {\n  int<lower=1> N;\n  int<lower=1> S;\n  array[N] int<lower=1,upper=S> state;\n  vector[N] y_lower;\n  vector[N] y_upper;\n  array[N] int<lower=0,upper=2> cens;\n  int<lower=1> K;\n  matrix[N,K] X;\n} \nparameters {\n  vector[K] beta;\n  vector[S] z_u;\n  real<lower=0> sigma;\n  real<lower=0> sigma_u;\n} \ntransformed parameters {\n  vector[S] u = sigma_u * z_u;\n} \nmodel {\n  beta ~ normal(0,5);\n  z_u ~ normal(0,1);\n  sigma ~ exponential(1);\n  sigma_u ~ exponential(1);\n  for (i in 1:N) {\n    real mu = X[i] * beta + u[state[i]];\n    if (cens[i]==0) {\n      target += normal_lpdf(y_lower[i] | mu, sigma);\n    } else {\n      real p_up = normal_cdf(y_upper[i] | mu, sigma);\n      real p_lo = normal_cdf(y_lower[i] | mu, sigma);\n      real diff = fmax(p_up - p_lo, 1e-12);\n      target += log(diff);\n    }\n  }\n}"
stan_file <- file.path(tempdir(), "interval_mixed_svi_strat.stan")
writeLines(stan_code, stan_file)
mod <- cmdstan_model(stan_file)

project_root <- normalizePath(".")

# ── Scenarios / disease labels ─────────────────────────────────────────────────
scenario_list <- list(
  list(key = "SVI_AAMR2006_2010", aamr = "2006-2010", lag = 5),
  list(key = "SVI_AAMR2011_2015", aamr = "2011-2015", lag = 10),
  list(key = "SVI_AAMR2016_2020", aamr = "2016-2020", lag = 15),
  list(key = "SVI_AAMR2021_2024", aamr = "2021-2024", lag = 20)
)
icd_to_name <- c(
  "I00_I99" = "CVD", "J40_J47_J60_J70_J84_D86_C34" = "CRD", "K70_K76_C22" = "CLD",
  "N00_N29_C64_C65" = "CKD", "X60_X84_Y87.0" = "Suicide",
  "G20_G30_G12.2_F01_F03" = "NDD", "C00_C97" = "Cancer"
)
disease_label <- function(icd) {
  nm <- icd_to_name[icd]
  if (is.na(nm)) icd else unname(nm)
}

# ── Stratification phases (run in this order) ──────────────────────────────────
PHASES <- list(
  list(tag = "Cluster", data = "Data/Processed/df_SVI.csv",
       out = "brms_SVI_Stratified_Cluster",
       vars = c("Cluster_EQI", "Cluster_NLCD"),
       tagfn = function(v, x) paste0(v, "_", x)),
  list(tag = "Demo", data = "Data/Processed/df_SVI_Stratified.csv",
       out = "brms_SVI_Stratified_Demo",
       vars = c("Stratum"),
       tagfn = function(v, x) as.character(x)),
  list(tag = "Geo", data = "Data/Processed/df_SVI.csv",
       out = "brms_SVI_Stratified_Geo",
       vars = c("Census_Region", "Climate_Zone", "Economic_type", "Homeownership_tertile"),
       tagfn = function(v, x) paste0(v, "_", x)),
  list(tag = "RUCC", data = "Data/Processed/df_SVI.csv",
       out = "brms_SVI_Stratified_RUCC",
       vars = c("RUCC"),
       tagfn = function(v, x) paste0("RUCC", x))
)

# ── Helpers ────────────────────────────────────────────────────────────────────
format_cell <- function(draws) {
  if (length(draws) == 0) return("")
  ci <- quantile(draws, c(0.025, 0.975), na.rm = TRUE)
  sprintf("%0.2f(%0.2f,%0.2f)", mean(draws), ci[1], ci[2])
}
append_rows <- function(path, df) {
  if (!file.exists(path)) write_csv(df, path) else suppressWarnings(write.table(df, path, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE))
}
compute_p <- function(draws) {
  if (length(draws) == 0) return(NA_character_)
  pos <- sum(draws > 0, na.rm = TRUE); neg <- sum(draws < 0, na.rm = TRUE)
  n <- pos + neg
  if (n == 0) return(NA_character_)
  sprintf("%.4f", 2 * min((pos + 0.5) / (n + 1), (neg + 0.5) / (n + 1)))
}
svi_est <- function(draws, names_vec) {
  out <- list(A = "0.00", B = "", C = "", D = "")
  for (lv in c("B", "C", "D")) {
    idx <- match(paste0("SVI_factor", lv), names_vec)
    out[[lv]] <- if (is.na(idx)) "" else format_cell(draws[[paste0("beta[", idx, "]")]])
  }
  out
}
svi_diag <- function(draws, names_vec, summ_df) {
  out <- list()
  for (lv in c("B", "C", "D")) {
    p <- NA_character_; rhat <- NA_real_; eb <- NA_real_; et <- NA_real_
    idx <- match(paste0("SVI_factor", lv), names_vec)
    if (!is.na(idx)) {
      col <- paste0("beta[", idx, "]")
      p <- compute_p(draws[[col]])
      sr <- summ_df[summ_df$variable == col, , drop = FALSE]
      if (nrow(sr)) { rhat <- sr$rhat; eb <- sr$ess_bulk; et <- sr$ess_tail }
    }
    out[[paste0(lv, "_p")]] <- p
    out[[paste0(lv, "_rhat")]] <- sprintf("%.4f", rhat)
    out[[paste0(lv, "_ess_bulk")]] <- as.integer(round(eb))
    out[[paste0(lv, "_ess_tail")]] <- as.integer(round(et))
  }
  out
}
compute_mu_A <- function(draws, names_vec, d, covariates) {
  mu <- draws[["beta[1]"]]
  for (cn in covariates) {
    idx <- match(cn, names_vec)
    if (!is.na(idx)) mu <- mu + draws[[paste0("beta[", idx, "]")]] * mean(d[[cn]], na.rm = TRUE)
  }
  mu
}
# Within-stratum (LagRef) MRR: divide each category mu by this scenario's mu_A.
mrr_block <- function(draws, names_vec, mu_A) {
  est <- list(A = "", B = "", C = "", D = "")
  pvl <- list(A_p = NA_character_, B_p = NA_character_, C_p = NA_character_, D_p = NA_character_)
  est$A <- format_cell(mu_A / mu_A)            # = 1.0000
  for (lv in c("B", "C", "D")) {
    idx <- match(paste0("SVI_factor", lv), names_vec)
    if (is.na(idx)) next
    rr <- (mu_A + draws[[paste0("beta[", idx, "]")]]) / mu_A
    est[[lv]] <- format_cell(rr)
    pvl[[paste0(lv, "_p")]] <- compute_p(rr - 1)
  }
  list(est = est, pvals = pvl)
}
build_design <- function(d, covariates) {
  d <- as.data.frame(d)
  d$SVI_factor <- factor(d$SVI, levels = c("A", "B", "C", "D"))
  needed <- c("SVI_factor", "AAMR_Lower", "AAMR_Upper", "cens", "State_FIPS", covariates)
  d <- d[stats::complete.cases(d[, needed, drop = FALSE]), , drop = FALSE]
  form <- if (length(covariates) == 0) ~SVI_factor else as.formula(paste("~", paste(c(covariates, "SVI_factor"), collapse = " + ")))
  mm <- model.matrix(form, d, contrasts.arg = list(SVI_factor = contr.treatment(c("A", "B", "C", "D"))))
  colnames(mm) <- make.names(colnames(mm))
  list(X = mm, names = colnames(mm), df = d)
}

# memoized data loader (adds State_FIPS + cens)
.data_cache <- new.env(parent = emptyenv())
load_data <- function(rel_path) {
  if (exists(rel_path, envir = .data_cache)) return(get(rel_path, envir = .data_cache))
  p <- file.path(project_root, rel_path)
  if (!file.exists(p)) { assign(rel_path, NULL, envir = .data_cache); return(NULL) }
  d <- fread(p)
  if (!"State_FIPS" %in% names(d)) d[, State_FIPS := substr(sprintf("%05s", COUNTY_FIPS), 1, 2)]
  d <- d[!is.na(AAMR_Lower) & !is.na(AAMR_Upper)]
  d[, cens := ifelse(AAMR_Lower == AAMR_Upper, 0, 2)]
  assign(rel_path, d, envir = .data_cache)
  d
}

# ── One stratification phase ───────────────────────────────────────────────────
run_phase <- function(phase, selected) {
  message("\n########## Phase: ", phase$tag, " ##########")
  dt <- load_data(phase$data)
  if (is.null(dt)) { message("[Skip phase] data not found: ", phase$data); return(invisible()) }

  base_req <- c("COUNTY_FIPS", "Time_Period", "Outcome", "AAMR_Lower", "AAMR_Upper", "SVI")
  if (length(setdiff(base_req, names(dt)))) {
    message("[Skip phase] ", phase$tag, " missing base cols: ",
            paste(setdiff(base_req, names(dt)), collapse = ",")); return(invisible())
  }
  present_vars <- phase$vars[phase$vars %in% names(dt)]
  if (length(present_vars) == 0) {
    message("[Skip phase] ", phase$tag, " — none of its strata columns present (",
            paste(phase$vars, collapse = ","), ")"); return(invisible())
  }
  if (length(present_vars) < length(phase$vars)) {
    message("[Note] ", phase$tag, " missing strata cols: ",
            paste(setdiff(phase$vars, present_vars), collapse = ","))
  }

  out_dir <- file.path(project_root, opt$`output-base`, phase$out)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  for (outcome in selected) {
    if (!(outcome %in% dt$Outcome)) next
    rc <- resolve_covset(outcome)
    if (!isTRUE(rc$ok)) { message("[Skip] ", outcome, " -> ", rc$reason); next }
    covset <- intersect(rc$covset, names(dt))     # only covariates present in this data
    if (length(covset) < length(rc$covset)) {
      message("[Note] ", outcome, " missing covariate(s): ",
              paste(setdiff(rc$covset, covset), collapse = ","))
    }
    model_core <- paste(c("SVI", covar_abbrev[covset]), collapse = "+")
    dlabel <- disease_label(outcome)
    outfile_mrd <- file.path(out_dir, paste0(dlabel, "_Stratified_", phase$tag, "_MRD.csv"))
    outfile_mrr <- file.path(out_dir, paste0(dlabel, "_Stratified_", phase$tag, "_MRR.csv"))
    message("== Outcome: ", dlabel, " (", outcome, ")  [", rc$cat, "]  Model core: ", model_core)

    for (sv in present_vars) {
      strat_vals <- sort(unique(na.omit(dt[Outcome == outcome][[sv]])))
      for (val in strat_vals) {
        strat_tag <- phase$tagfn(sv, val)
        strat_dt <- dt[Outcome == outcome & !is.na(get(sv)) & get(sv) == val]
        if (nrow(strat_dt) < opt$`min-n`) next

        for (sc in scenario_list) {
          scen_dt <- strat_dt[Time_Period == sc$aamr]
          if (nrow(scen_dt) < opt$`min-n`) next
          aamr_out <- gsub("-", "_", sc$aamr)

          des <- build_design(scen_dt, covset)
          if (nrow(des$df) < opt$`min-n`) next
          states <- sort(unique(des$df$State_FIPS))
          si <- match(des$df$State_FIPS, states)
          dl <- list(
            N = nrow(des$df), S = length(states), state = si,
            y_lower = des$df$AAMR_Lower, y_upper = des$df$AAMR_Upper, cens = des$df$cens,
            K = ncol(des$X), X = des$X
          )
          init_fn <- function() list(beta = rep(0, dl$K), z_u = rep(0, dl$S), sigma = 50, sigma_u = 10)
          fit <- try(mod$sample(
            data = dl, chains = opt$chains, iter_sampling = opt$iter - opt$warmup, iter_warmup = opt$warmup,
            adapt_delta = opt$`adapt-delta`, max_treedepth = opt$`max-treedepth`,
            parallel_chains = opt$chains, refresh = 0, seed = opt$seed,
            init = rep(list(init_fn()), opt$chains)
          ), silent = TRUE)
          if (inherits(fit, "try-error")) { message("[Fail] ", strat_tag, " ", sc$key); next }

          draws <- as_draws_df(fit$draws("beta"))
          colnames(draws) <- paste0("beta[", seq_len(ncol(draws)), "]")
          summ <- posterior::summarize_draws(fit$draws("beta"))
          model_label <- paste0(strat_tag, "_", model_core)
          est <- svi_est(draws, des$names)
          diag <- svi_diag(draws, des$names, summ)
          mu_A <- compute_mu_A(draws, des$names, des$df, covset)
          mrr <- mrr_block(draws, des$names, mu_A)

          row_mrd <- do.call(tibble, c(
            list(ICD_Code = outcome, AAMR_Period = aamr_out, Lag = sc$lag, Model = model_label),
            est, diag
          ))
          row_mrr <- do.call(tibble, c(
            list(ICD_Code = outcome, AAMR_Period = aamr_out, Lag = sc$lag, Model = model_label),
            list(A = mrr$est$A, B = mrr$est$B, C = mrr$est$C, D = mrr$est$D),
            list(A_p = mrr$pvals$A_p, B_p = mrr$pvals$B_p, C_p = mrr$pvals$C_p, D_p = mrr$pvals$D_p),
            list(B_rhat = diag$B_rhat, C_rhat = diag$C_rhat, D_rhat = diag$D_rhat,
                 B_ess_bulk = diag$B_ess_bulk, C_ess_bulk = diag$C_ess_bulk, D_ess_bulk = diag$D_ess_bulk,
                 B_ess_tail = diag$B_ess_tail, C_ess_tail = diag$C_ess_tail, D_ess_tail = diag$D_ess_tail)
          ))
          append_rows(outfile_mrd, row_mrd)
          append_rows(outfile_mrr, row_mrr)
          message("[OK] ", strat_tag, " ", sc$key)
        }
      }
    }
    message("== Completed: ", dlabel)
  }
}

# ── Determine selected outcomes (default: 7 overall) ───────────────────────────
ref_dt <- load_data("Data/Processed/df_SVI.csv")
all_outcomes <- if (is.null(ref_dt)) names(icd_to_name) else sort(unique(ref_dt$Outcome))
selected <- if (is.na(opt$`outcomes`)) {
  intersect(names(icd_to_name), all_outcomes)
} else {
  reqc <- str_split(opt$`outcomes`, ",", simplify = TRUE) |> as.vector() |> str_trim()
  reqc
}
message("Outcomes to analyze: ", paste(selected, collapse = ","))

# ── Run all four phases in order ───────────────────────────────────────────────
for (phase in PHASES) run_phase(phase, selected)
message("\nAll stratified analyses complete. Output base: ", file.path(opt$`output-base`))
