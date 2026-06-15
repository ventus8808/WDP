#!/usr/bin/env Rscript
# ============================================================================
# RidgeLine Visualization for SVI brms_SVI_Main_RL Results (7 outcomes, 4 lags)
# ============================================================================
# SVI exposure: posterior distributions of the B/C/D mortality rate differences
# (MRD vs the A reference). One SVI model per outcome (no EQI domains).
#
# Input : Result/brms_SVI_Main_RL/{Outcome}_Lag{N}_SVI.rds
#         (draws_long: category in {B,C,D}, effect; metadata$lag)
# Output: Result/brms_SVI_RL_Visualization/{Outcome}/
#           SVI_Lag5_10_15_{Palette}.png        (3-lag, gradient fill)
#           SVI_Lag5_10_15_20_{Palette}.png     (4-lag, gradient fill)
#           SVI_Lag5_10_15_20_Categorical.png   (fixed B/C/D colours)
# Usage : Rscript Code/Visualization/RidgeLine_SVI_RL.R
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggridges)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
})

# ── Configuration ──────────────────────────────────────────────────────────────
OVERALL_GRADIENTS <- list(
  "ViridisStrong"   = c("#FDE725", "#7AD151", "#22A884", "#2A788E", "#414487", "#440154"),
  "Spectral_Reverse" = c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598",
                          "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"),
  "RedBlue"         = c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                          "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
  "GreenGold"       = c("#00441B", "#1B7837", "#5AAE61", "#A6DBA0", "#D9F0D3",
                          "#FEE08B", "#E6AB02", "#A6761D")
)
# Fixed categorical colours (B low-mid -> D high), muted academic green->orange.
CAT_COLORS <- c(B = "#9CC196", C = "#E6C58A", D = "#D58A4E")

FONT_FAMILY <- "Helvetica"
CAT_ORDER <- c("D", "C", "B")               # bottom -> top (B on top, D at bottom)
MIN_SUBPLOT_WIDTH_PX <- 1600
MIN_SUBPLOT_HEIGHT_PX <- 700
SUBPLOT_GAP_Y_IN <- MIN_SUBPLOT_HEIGHT_PX / 10 / 300

OUTCOMES <- c("Cancer", "CKD", "CLD", "CRD", "CVD", "NDD", "Suicide")

# ── Utilities ──────────────────────────────────────────────────────────────────
theme_publication <- function(base_size = 14, base_family = FONT_FAMILY) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.x = element_line(color = "#E0E0E0", linewidth = 0.3),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      text = element_text(family = base_family, color = "#000000"),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 8)),
      axis.text = element_text(size = 14, color = "#000000"),
      axis.line.x = element_line(color = "#000000", linewidth = 0.5),
      axis.ticks.x = element_line(color = "#000000", linewidth = 0.3),
      axis.ticks.length = unit(0.15, "cm"),
      strip.text = element_text(size = 14, face = "bold", color = "#000000"),
      strip.background = element_rect(fill = "#F5F5F5", color = "#CCCCCC", linewidth = 0.3),
      legend.position = "none",
      panel.spacing.y = unit(SUBPLOT_GAP_Y_IN, "in"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.margin = margin(15, 15, 10, 10)
    )
}

save_plot_px <- function(plot_obj, output_file, panel_rows) {
  ggsave(output_file, plot = plot_obj,
         width = MIN_SUBPLOT_WIDTH_PX, height = panel_rows * MIN_SUBPLOT_HEIGHT_PX,
         units = "px", dpi = 300, bg = "white", limitsize = FALSE)
}

# ── Data loading ───────────────────────────────────────────────────────────────
load_ridgeline_data <- function(files) {
  all_data <- list()
  for (f in files) {
    data <- readRDS(f)
    draws <- data$draws_long
    draws$lag <- data$metadata$lag
    all_data[[length(all_data) + 1]] <- draws
  }
  bind_rows(all_data)
}

# ── Plot: B/C/D ridgelines faceted by lag ──────────────────────────────────────
plot_svi_lags <- function(data, lags, output_file, gradient_colors = NULL,
                          title = NULL) {
  combined <- data %>%
    filter(lag %in% lags) %>%
    mutate(
      cat_label = factor(category, levels = CAT_ORDER),
      lag_label = factor(sprintf("Lag %d", lag), levels = sprintf("Lag %d", lags))
    )

  if (is.null(gradient_colors)) {
    # categorical fill by SVI category
    p <- ggplot(combined, aes(x = effect, y = cat_label, fill = cat_label)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
      geom_density_ridges(alpha = 0.85, scale = 1.3, rel_min_height = 0.005,
                          color = "black", linewidth = 0.45, quantile_lines = FALSE) +
      scale_fill_manual(values = CAT_COLORS, guide = "none")
  } else {
    # gradient fill by effect magnitude
    p <- ggplot(combined, aes(x = effect, y = cat_label, fill = after_stat(x))) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
      geom_density_ridges_gradient(alpha = 0.85, scale = 1.3, rel_min_height = 0.005,
                                   color = "black", linewidth = 0.45, quantile_lines = FALSE) +
      scale_fill_gradientn(colors = gradient_colors, name = "MRD", guide = "none")
  }

  p <- p +
    scale_x_continuous(name = "MRD vs A (deaths/100k), posterior distribution",
                       breaks = pretty_breaks(n = 6),
                       expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_discrete(name = NULL, expand = expansion(add = c(0.1, 1.6))) +
    facet_wrap(~lag_label, ncol = 1) +
    theme_publication()
  if (!is.null(title)) p <- p + ggtitle(title)

  save_plot_px(p, output_file, panel_rows = length(lags))
}

# ── Per-outcome processing ─────────────────────────────────────────────────────
process_outcome <- function(outcome_name, project_root) {
  cat(sprintf("Processing: %s\n", outcome_name))
  output_dir <- file.path(project_root, sprintf("Result/brms_SVI_RL_Visualization/%s", outcome_name))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  brms_dir <- file.path(project_root, "Result/brms_SVI_Main_RL")
  pattern <- sprintf("^%s_Lag[0-9]+_SVI\\.rds$", outcome_name)
  files <- list.files(brms_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    cat(sprintf("  ERROR: no RDS files for '%s'\n", outcome_name))
    return(FALSE)
  }
  cat(sprintf("  found %d RDS files\n", length(files)))

  all_data <- load_ridgeline_data(files)
  avail_lags <- sort(unique(all_data$lag))
  lags4 <- intersect(c(5, 10, 15, 20), avail_lags)
  lags3 <- intersect(c(5, 10, 15), avail_lags)

  # Gradient versions per palette (3-lag and 4-lag)
  for (palette_name in names(OVERALL_GRADIENTS)) {
    grad <- OVERALL_GRADIENTS[[palette_name]]
    plot_svi_lags(all_data, lags3,
                  file.path(output_dir, sprintf("SVI_Lag5_10_15_%s.png", palette_name)),
                  gradient_colors = grad, title = outcome_name)
    plot_svi_lags(all_data, lags4,
                  file.path(output_dir, sprintf("SVI_Lag5_10_15_20_%s.png", palette_name)),
                  gradient_colors = grad, title = outcome_name)
  }
  # Categorical version (4-lag)
  plot_svi_lags(all_data, lags4,
                file.path(output_dir, "SVI_Lag5_10_15_20_Categorical.png"),
                gradient_colors = NULL, title = outcome_name)

  cat(sprintf("  done (%d plots)\n", 2 * length(OVERALL_GRADIENTS) + 1))
  TRUE
}

main <- function() {
  project_root <- normalizePath(".")
  cat("RidgeLine Visualization — brms_SVI_Main_RL\n")
  results <- vapply(OUTCOMES, function(o) isTRUE(process_outcome(o, project_root)),
                    logical(1))
  cat(sprintf("\nCompleted: %d/%d outcomes -> Result/brms_SVI_RL_Visualization/\n",
              sum(results), length(OUTCOMES)))
}

if (!interactive()) main()
