#!/usr/bin/env Rscript
# 区间AAMR使用 brms 拟合单批次模型 (一个场景 + 一个癌症类型 + 一个 RUCC 分层)
# 输入: --data-file 预过滤 CSV (含 AAMR_lower, AAMR_upper, EQI分位, Smoking_Rate, State)
# 输出: --json-out JSON (每行一个模型, 包含 Q1..Q5 估计/区间/p 值)
# 模型: 总 EQI + 五域 EQI (如果存在) 以及 RUCC 分层变量 (已在上游切分)

suppressPackageStartupMessages({
  library(brms)
  library(jsonlite)
  library(dplyr)
  ## cmdstanr 后端用于替代 rstan，解决符号加载问题
  suppressWarnings(suppressMessages(require(cmdstanr)))
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
getArg <- function(flag) {
  i <- which(args == flag)
  if (length(i) == 0 || i == length(args)) return(NA)
  args[i+1]
}

DATA_FILE <- getArg('--data-file')
JSON_OUT  <- getArg('--json-out')
SCENARIO  <- getArg('--scenario')
CANCER    <- getArg('--cancer-type')
RUCC_NAME <- getArg('--rucc-name')
CHAINS    <- as.integer(getArg('--chains'))
ITER      <- as.integer(getArg('--iter'))
WARMUP    <- as.integer(getArg('--warmup'))
CORES     <- as.integer(getArg('--cores'))
SEED      <- as.integer(getArg('--seed'))

if (is.na(DATA_FILE) || !file.exists(DATA_FILE)) {
  stop('数据文件不存在: ', DATA_FILE)
}

raw <- read.csv(DATA_FILE, stringsAsFactors = FALSE)
# 基本区间清洗: 去除上下界缺失与反向区间
raw <- raw %>% filter(!is.na(AAMR_lower), !is.na(AAMR_upper)) %>%
  filter(AAMR_lower <= AAMR_upper)
if (nrow(raw) < 50) {
  write(jsonlite::toJSON(list(), auto_unbox = TRUE), JSON_OUT)
  quit(save='no')
}

## ================= 数据预处理 =================
## cens_type: 0 精确点, 2 区间 (brms 需要 cbind(lower, upper) | cens(cens_type))
raw$cens_type <- ifelse(raw$AAMR_lower == raw$AAMR_upper, 0, 2)

## 可选: 对响应做 log 变换 (仅当全部>0) —— 防止尺度过大导致采样失败
use_log <- all(raw$AAMR_lower > 0)
if (use_log) {
  raw$AAMR_lower <- log(raw$AAMR_lower)
  raw$AAMR_upper <- log(raw$AAMR_upper)
}

raw$EQI <- factor(raw$EQI, levels = 1:5)
raw$Smoking_Rate_std <- scale(raw$Smoking_Rate, scale = FALSE)
raw$State <- factor(raw$State)

# 尝试域列存在性
domains <- c('EQI_Air','EQI_Water','EQI_Land','EQI_Built','EQI_Social')
for (d in domains) {
  if (d %in% names(raw)) {
    raw[[d]] <- factor(raw[[d]], levels = 1:5)
  }
}

## ================ 模型公式 =====================
## 正式区间删失语法: cbind(lower, upper) | cens(cens_type) ~ predictors
form_total <- bf(cbind(AAMR_lower, AAMR_upper) | cens(cens_type) ~ EQI + Smoking_Rate_std + (1|State))

has_all_domains <- all(domains %in% names(raw))
form_domains <- if (has_all_domains) bf(cbind(AAMR_lower, AAMR_upper) | cens(cens_type) ~ 
  EQI_Air + EQI_Water + EQI_Land + EQI_Built + EQI_Social + Smoking_Rate_std + (1|State)) else NULL

## 使用 cmdstanr 后端（需已安装 cmdstan）
options(brms.backend = 'cmdstanr')

## 先验: 提供稳定性，防止区间尺度过大导致发散
prior_list <- c(
  set_prior('normal(0, 5)', class='b'),
  set_prior('student_t(3, 0, 5)', class='sd'),
  set_prior('student_t(3, 0, 5)', class='sigma'),
  set_prior(if (use_log) 'normal(4.5, 1)' else 'normal(200,50)', class='Intercept')
)

fits <- list()

safe_fit <- function(formula, data, model_key) {
  if (is.null(formula)) return(NULL)
  tryCatch({
    fit <- brm(
      formula = formula,
      data = data,
      family = gaussian(),
      prior = prior_list,
      chains = CHAINS, iter = ITER, warmup = WARMUP,
      cores = CORES, seed = SEED,
      backend = 'cmdstanr',
      control = list(adapt_delta = 0.99, max_treedepth = 12)
    )
    fit$model_key <- model_key
    fit
  }, error = function(e) {
    message('模型失败: ', model_key, ' | ', e$message)
    NULL
  })
}

fits[['EQI']] <- safe_fit(form_total, raw, 'EQI')
if (!is.null(form_domains)) {
  fits[['EQI_domains_multi']] <- safe_fit(form_domains, raw, 'EQI_domains_multi')
}

# 如果 RUCC 分层, 需要将模型 key 改为 RUCC前缀 (由上游 Python 决定 RUCC_NAME)
prefix <- ifelse(is.na(RUCC_NAME) || RUCC_NAME == 'ALL', '', RUCC_NAME)

extract_quintiles <- function(fit, var_prefix, model_key_out) {
  if (is.null(fit)) return(NULL)
  summ <- suppressWarnings(posterior_summary(fit))
  # Quintile levels: 1 baseline, 2..5 have系数
  # brms 命名: EQI2, EQI3, EQI4, EQI5 (因为因子自动处理)
  mkrow <- function(q) {
    if (q == 1) {
      return(list(q = 'Q1', est = 0, l = 0, u = 0, p = NA))
    }
    term <- paste0(var_prefix, q)  # e.g. EQI2
    if (!(term %in% rownames(summ))) {
      return(list(q = paste0('Q', q), est = NA, l = NA, u = NA, p = NA))
    }
    est <- summ[term,'Estimate']; l <- summ[term,'Q2.5']; u <- summ[term,'Q97.5']
    # 后验近似 p 值: 双尾概率 (保守) -> 2 * min(P(beta>0), P(beta<0))
    draws <- as.numeric(as_draws_df(fit)[[term]])
    prob_pos <- mean(draws > 0); prob_neg <- mean(draws < 0)
    p_two <- 2 * min(prob_pos, prob_neg)
    list(q = paste0('Q', q), est = est, l = l, u = u, p = p_two)
  }
  quint_rows <- lapply(1:5, mkrow)
  out <- list(model_key = model_key_out)
  for (qr in quint_rows) {
    out[[paste0(qr$q, '_estimate')]] <- qr$est
    out[[paste0(qr$q, '_lower')]]   <- qr$l
    out[[paste0(qr$q, '_upper')]]   <- qr$u
    out[[paste0(qr$q, '_p')]]       <- qr$p
  }
  out
}

results_list <- list()

# 总 EQI
if (!is.null(fits[['EQI']])) {
  mk <- if (prefix == '') 'EQI' else paste0(prefix, '_EQI')
  results_list[[length(results_list)+1]] <- extract_quintiles(fits[['EQI']], 'EQI', mk)
}

# 域模型: 如果拟合成功, 按域拆成多行 (与 LMM 域结构对齐)
if (!is.null(fits[['EQI_domains_multi']])) {
  # brms 中每个域的因子系数命名: EQI_Air2, EQI_Air3, ...
  domains_map <- c(Air='EQI_Air', Water='EQI_Water', Land='EQI_Land', Built='EQI_Built', Social='EQI_Social')
  for (lab in names(domains_map)) {
    varpref <- domains_map[[lab]]
    mk <- if (prefix == '') paste0('EQI_', lab) else paste0(prefix, '_EQI_', lab)
    results_list[[length(results_list)+1]] <- extract_quintiles(fits[['EQI_domains_multi']], paste0(varpref), mk)
  }
}

jsonlite::write_json(results_list, JSON_OUT, auto_unbox = TRUE, pretty = TRUE)
message('写出 JSON: ', JSON_OUT)
if (use_log) message('注意: 使用 log(AAMR) 尺度拟合, 结果系数在 log 量纲上')
