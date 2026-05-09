#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hit <- grep(file_arg, args)
  if (length(hit) > 0) {
    raw_path <- sub(file_arg, "", args[hit[1]])
    decoded_path <- gsub("~\\+~", " ", raw_path)
    candidate <- if (file.exists(raw_path)) raw_path else decoded_path
    if (file.exists(candidate)) {
      return(dirname(normalizePath(candidate)))
    }
  }
  getwd()
}

extract_domain <- function(path) {
  nm <- basename(path)
  nm <- sub("\\.rds$", "", nm)
  sub("^.*_", "", nm)
}

cli_args <- commandArgs(trailingOnly = TRUE)
target_dir <- if (length(cli_args) >= 1) cli_args[1] else get_script_dir()
target_dir <- normalizePath(target_dir, mustWork = TRUE)
rds_files <- list.files(target_dir, pattern = "\\.rds$", full.names = TRUE)

if (length(rds_files) == 0) {
  stop("No .rds files found in: ", target_dir)
}

message("Found ", length(rds_files), " RDS files")

draws_list <- list()
summary_list <- list()

for (f in rds_files) {
  obj <- readRDS(f)
  domain <- extract_domain(f)

  if (is.null(obj$draws_long) || !all(c("quintile", "effect") %in% names(obj$draws_long))) {
    warning("Skipping (missing draws_long/quintile/effect): ", basename(f))
    next
  }

  d <- obj$draws_long
  d$domain <- domain
  d$source_file <- basename(f)
  draws_list[[length(draws_list) + 1]] <- d

  if (!is.null(obj$summary)) {
    s <- obj$summary
    s$domain <- domain
    s$source_file <- basename(f)
    summary_list[[length(summary_list) + 1]] <- s
  }
}

if (length(draws_list) == 0) {
  stop("No usable draws_long data found in RDS files.")
}

draws <- do.call(rbind, draws_list)
draws$quintile <- factor(draws$quintile, levels = c("Q2", "Q3", "Q4", "Q5"))

# Keep plotting responsive if draws are very large
if (nrow(draws) > 300000) {
  set.seed(1234)
  draws <- draws[sample.int(nrow(draws), 300000), , drop = FALSE]
}

p <- ggplot(draws, aes(x = effect, color = quintile, fill = quintile)) +
  geom_density(alpha = 0.18, linewidth = 0.5, adjust = 1.1) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "black") +
  facet_wrap(~domain, scales = "free_y", ncol = 2) +
  labs(
    title = "RDS Check: Posterior Effect Densities by EQI Domain",
    x = "Effect",
    y = "Density",
    color = "Quintile",
    fill = "Quintile"
  ) +
  theme_bw(base_size = 12)

plot_file <- file.path(target_dir, "rds_preview_density.png")
ggsave(plot_file, p, width = 12, height = 8, dpi = 300)

if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, summary_list)
  write.csv(summary_df, file.path(target_dir, "rds_preview_summary.csv"), row.names = FALSE)
}

message("Saved plot: ", plot_file)
if (length(summary_list) > 0) {
  message("Saved summary: ", file.path(target_dir, "rds_preview_summary.csv"))
}
