#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# WDP BYM INLA Logging Utilities
# Advanced logging and error tracking system
# Author: WDP Analysis Team
# Date: 2024

# Global logger settings
.logger_settings <- list(
  log_level = "INFO",
  log_file = NULL,
  timestamp_format = "%Y-%m-%d %H:%M:%S",
  console_output = TRUE
)

#' Set up logger configuration
#' @param log_level Minimum logging level ("DEBUG", "INFO", "WARNING", "ERROR")
#' @param log_file Optional file path for log output
#' @param console_output Whether to also print to console
setup_logger <- function(log_level = "INFO", log_file = NULL, console_output = TRUE) {
  .logger_settings$log_level <<- toupper(log_level)
  .logger_settings$log_file <<- log_file
  .logger_settings$console_output <<- console_output
  
  if (!is.null(log_file)) {
    # Create log directory if it doesn't exist
    log_dir <- dirname(log_file)
    if (!dir.exists(log_dir)) {
      dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Initialize log file with header
    log_header <- sprintf("=== WDP INLA Analysis Log Started at %s ===\n", 
                         format(Sys.time(), .logger_settings$timestamp_format))
    cat(log_header, file = log_file, append = FALSE)
  }
  
  log_info(sprintf("Logger initialized - Level: %s, File: %s, Console: %s", 
                   log_level, ifelse(is.null(log_file), "None", log_file), console_output))
}

#' Get numeric log level for comparison
#' @param level Character log level
#' @return Numeric level
get_log_level_numeric <- function(level) {
  levels <- c("DEBUG" = 1, "INFO" = 2, "WARNING" = 3, "ERROR" = 4)
  return(levels[toupper(level)])
}

#' Generic logging function
#' @param level Log level
#' @param message Log message
#' @param context Optional context information
#' @param ... Additional parameters for sprintf formatting
log_message <- function(level, message, context = NULL, ...) {
  level <- toupper(level)
  
  # Check if this message should be logged based on current log level
  if (get_log_level_numeric(level) < get_log_level_numeric(.logger_settings$log_level)) {
    return(invisible())
  }
  
  # Format message
  if (length(list(...)) > 0) {
    message <- sprintf(message, ...)
  }
  
  # Create timestamp
  timestamp <- format(Sys.time(), .logger_settings$timestamp_format)
  
  # Create context string
  context_str <- if (!is.null(context)) sprintf(" [%s]", context) else ""
  
  # Format complete log line
  log_line <- sprintf("[%s] %s%s: %s\n", timestamp, level, context_str, message)
  
  # Output to console if enabled
  if (.logger_settings$console_output) {
    # Use different formatting for console
    level_emoji <- switch(level,
                         "DEBUG" = "🔍",
                         "INFO" = "ℹ️",
                         "WARNING" = "⚠️",
                         "ERROR" = "❌")
    
    console_line <- sprintf("%s %s%s: %s\n", level_emoji, level, context_str, message)
    cat(console_line)
  }
  
  # Write to log file if specified
  if (!is.null(.logger_settings$log_file)) {
    cat(log_line, file = .logger_settings$log_file, append = TRUE)
  }
}

#' Log debug message
#' @param message Debug message
#' @param context Optional context
#' @param ... Additional sprintf parameters
log_debug <- function(message, context = NULL, ...) {
  log_message("DEBUG", message, context, ...)
}

#' Log info message
#' @param message Info message
#' @param context Optional context
#' @param ... Additional sprintf parameters
log_info <- function(message, context = NULL, ...) {
  log_message("INFO", message, context, ...)
}

#' Log warning message
#' @param message Warning message
#' @param context Optional context
#' @param ... Additional sprintf parameters
log_warning <- function(message, context = NULL, ...) {
  log_message("WARNING", message, context, ...)
}

#' Log error message
#' @param message Error message
#' @param context Optional context
#' @param ... Additional sprintf parameters
log_error <- function(message, context = NULL, ...) {
  log_message("ERROR", message, context, ...)
}

#' Enhanced error logging with stack trace
#' @param error Error object from tryCatch
#' @param context Context where error occurred
#' @param additional_info Additional debugging information
log_error_with_trace <- function(error, context = NULL, additional_info = NULL) {
  error_msg <- conditionMessage(error)
  
  # Get call stack
  call_stack <- sys.calls()
  stack_trace <- paste(capture.output(traceback(error)), collapse = "\n")
  
  # Log main error
  log_error("Error occurred: %s", context, error_msg)
  
  # Log stack trace
  if (nchar(stack_trace) > 0) {
    log_debug("Stack trace:\n%s", context, stack_trace)
  }
  
  # Log additional information
  if (!is.null(additional_info)) {
    log_debug("Additional info: %s", context, additional_info)
  }
  
  # Return structured error info
  return(list(
    message = error_msg,
    stack_trace = stack_trace,
    context = context,
    timestamp = Sys.time(),
    additional_info = additional_info
  ))
}

#' Log function entry (for debugging)
#' @param function_name Name of the function being entered
#' @param parameters List of parameters
log_function_entry <- function(function_name, parameters = NULL) {
  param_str <- if (!is.null(parameters)) {
    paste(names(parameters), parameters, sep = "=", collapse = ", ")
  } else {
    "no parameters"
  }
  
  log_debug("Entering function: %s(%s)", function_name, param_str)
}

#' Log function exit (for debugging)
#' @param function_name Name of the function being exited
#' @param result Optional result summary
log_function_exit <- function(function_name, result = NULL) {
  result_str <- if (!is.null(result)) {
    sprintf(" with result: %s", result)
  } else {
    ""
  }
  
  log_debug("Exiting function: %s%s", function_name, result_str)
}

#' Log analysis progress
#' @param step Current step description
#' @param current Current progress number
#' @param total Total number of steps
#' @param additional_info Additional progress information
log_progress <- function(step, current = NULL, total = NULL, additional_info = NULL) {
  progress_str <- if (!is.null(current) && !is.null(total)) {
    sprintf(" [%d/%d]", current, total)
  } else {
    ""
  }
  
  info_str <- if (!is.null(additional_info)) {
    sprintf(" - %s", additional_info)
  } else {
    ""
  }
  
  log_info("Progress%s: %s%s", "", progress_str, step, info_str)
}

#' Create a performance timer
#' @param name Timer name
#' @return Timer object
create_timer <- function(name) {
  timer <- list(
    name = name,
    start_time = Sys.time(),
    checkpoints = list()
  )
  
  log_debug("Timer '%s' started", name)
  return(timer)
}

#' Add checkpoint to timer
#' @param timer Timer object
#' @param checkpoint_name Checkpoint description
add_checkpoint <- function(timer, checkpoint_name) {
  current_time <- Sys.time()
  elapsed <- as.numeric(difftime(current_time, timer$start_time, units = "secs"))
  
  timer$checkpoints[[checkpoint_name]] <- list(
    time = current_time,
    elapsed = elapsed
  )
  
  log_debug("Timer '%s' checkpoint '%s': %.2f seconds", timer$name, checkpoint_name, elapsed)
  return(timer)
}

#' Finalize timer and log total time
#' @param timer Timer object
finalize_timer <- function(timer) {
  end_time <- Sys.time()
  total_elapsed <- as.numeric(difftime(end_time, timer$start_time, units = "secs"))
  
  log_info("Timer '%s' completed: %.2f seconds total", timer$name, total_elapsed)
  
  return(list(
    name = timer$name,
    start_time = timer$start_time,
    end_time = end_time,
    total_elapsed = total_elapsed,
    checkpoints = timer$checkpoints
  ))
}