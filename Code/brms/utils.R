#!/usr/bin/env Rscript
#'
#' BRMS Utilities - Centralized Configuration and Common Functions
#'
#' 提供所有brms脚本共享的功能：
#' - 项目配置加载
#' - 项目根目录查找
#' - 通用数据处理函数
#'

suppressPackageStartupMessages({
  library(yaml)
})

# Utility operator
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Find project root by looking for config.yaml
find_project_root <- function(start = getwd()) {
  cur <- normalizePath(start)
  for (i in 1:6) {
    cand <- file.path(cur, "config.yaml")
    if (file.exists(cand)) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  stop("config.yaml not found from starting directory: ", start)
}

#' Load project configuration - centralized function
load_project_config <- function() {
  # Try to find root from script location first, then from working directory
  script_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- script_args[grep("^--file=", script_args)]
  
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    root <- tryCatch(find_project_root(script_dir), 
                     error = function(e) find_project_root(getwd()))
  } else {
    root <- find_project_root(getwd())
  }
  
  config_path <- file.path(root, "config.yaml")
  message(sprintf("[brms] Reading config from: %s", config_path))
  config <- yaml::read_yaml(config_path)
  
  # Attach root path to config for convenience
  config$project_root <- root
  
  return(config)
}