#!/usr/bin/env Rscript
# ============================================================================
# RidgeLine Visualization for brms_RL Results (7 Outcomes, 4 Lags)
# ============================================================================
# Updated to handle new 2021-2024 data (Lag 20)
# Generates 6 plots per outcome: 4 lag-specific + 2 overall plots
#
# Structure:
#   - 4 lag-specific multi-domain plots (Lag 5, 10, 15, 20)
#   - 2 overall EQI plots (3 lags: 5,10,15 AND 4 lags: 5,10,15,20)
#
# Usage:
#   Rscript Code/Visualization/RidgeLine_brms_RL.R
#
# Output:
#   Result/brms_RL_Visualization/{Outcome}/ (7 outcome folders)
#   Each with 6 PNG files (42 total for all outcomes)
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggridges)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
})

# ============================================================================
# Configuration
# ============================================================================

# ============================================================================
# Configuration
# ============================================================================

DOMAIN_PALETTES <- list(
  "Overall EQI" = c("#E8E8E8", "#BDBDBD", "#8C8C8C", "#525252", "#1A1A1A"),
  "Air" = c("#FFF3E0", "#FFCC80", "#FFA726", "#FB8C00", "#E65100"),
  "Water" = c("#E3F2FD", "#90CAF9", "#42A5F5", "#1E88E5", "#0D47A1"),
  "Land" = c("#E8F5E9", "#A5D6A7", "#66BB6A", "#43A047", "#2E7D32"),
  "Built" = c("#F5F5F5", "#D9D9D9", "#BDBDBD", "#969696", "#525252"),
  "Social" = c("#D4E3F0", "#A8C7E0", "#7DAAD0", "#5281A8", "#2E4E74")
)

# Overall plot gradient colors for testing
OVERALL_GRADIENTS <- list(


  ## Viridis-like（高对比、Nature风格常见）
  "ViridisStrong" = c(
    "#FDE725", "#7AD151", "#22A884", "#2A788E", "#414487", "#440154"
  ),

  ## Plasma（亮黄 -> 紫红，视觉冲击强）
  "Plasma" = c(
    "#F0F921", "#FDB42F", "#ED7953", "#CC4778", "#9C179E", "#5C01A6", "#0D0887"
  ),


  ## Spectral（reverse）
  "Spectral_Reverse" = c(
    "#5E4FA2", "#3288BD", "#66C2A5",
    "#ABDDA4", "#E6F598", "#FEE08B",
    "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"
  ),

  ## RedBlue（适合风险增加/减少）
  "RedBlue" = c(
    "#053061", "#2166AC", "#4393C3", "#92C5DE",
    "#D1E5F0",
    "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"
  ),

  ## Ocean（深海蓝绿）
  "Ocean" = c(
    "#081D58", "#253494", "#225EA8",
    "#1D91C0", "#41B6C4", "#7FCDBB", "#C7E9B4"
  ),

  ## Sunset（适合高端论文图）
  "Sunset" = c(
    "#355C7D", "#6C5B7B", "#C06C84",
    "#F67280", "#F8B195"
  ),

  ## Neon（超高对比，适合 poster）
  "Neon" = c(
    "#00F5D4", "#00BBF9", "#9B5DE5",
    "#F15BB5", "#FEE440"
  ),

  ## CyanMagenta（现代感强）
  "CyanMagenta" = c(
    "#00B4D8", "#48CAE4", "#90E0EF",
    "#CAF0F8", "#FFCAD4", "#F4ACB7",
    "#C9184A"
  ),

  ## GreenGold（environmental epi 很适合）
  "GreenGold" = c(
    "#00441B", "#1B7837", "#5AAE61",
    "#A6DBA0", "#D9F0D3",
    "#FEE08B", "#E6AB02", "#A6761D"
  )
)


FONT_FAMILY <- "Helvetica"

# ============================================================================
# Utility Functions
# ============================================================================

theme_publication <- function(base_size = 14, base_family = FONT_FAMILY) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "#E0E0E0", linewidth = 0.3),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      text = element_text(family = base_family, color = "#000000"),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 8)),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 8)),
      axis.text = element_text(size = 14, color = "#000000"),
      axis.line.x = element_line(color = "#000000", linewidth = 0.5),
      axis.ticks.x = element_line(color = "#000000", linewidth = 0.3),
      axis.ticks.length = unit(0.15, "cm"),
      strip.text = element_text(size = 14, face = "bold", color = "#000000"),
      strip.background = element_rect(fill = "#F5F5F5", color = NA),
      legend.position = "none",
      plot.margin = margin(15, 15, 10, 10)
    )
}

# ============================================================================
# Data Loading
# ============================================================================

load_ridgeline_data <- function(files) {
  all_data <- list()

  for (f in files) {
    data <- readRDS(f)
    lag <- data$metadata$lag

    # Extract domain from filename (e.g., "Cancer_Lag5_Air.rds" → "Air")
    filename <- basename(f)
    domain <- sub(".*_Lag[0-9]+_([^.]+)\\.rds$", "\\1", filename)

    draws <- data$draws_long %>%
      mutate(
        domain = domain,
        lag = lag
      )

    all_data[[length(all_data) + 1]] <- draws
  }

  bind_rows(all_data)
}

# ============================================================================
# Plot Functions
# ============================================================================

plot_lag_all_domains <- function(data, lag, output_file) {
  cat(sprintf("  Processing Lag %d with all domains...\n", lag))

  combined <- data %>%
    filter(lag == !!lag, domain != "Overall") %>%
    mutate(
      quintile_label = factor(quintile, levels = c("Q5", "Q4", "Q3", "Q2")),
      quintile_num = as.numeric(factor(quintile, levels = c("Q2", "Q3", "Q4", "Q5"))),
      domain = factor(domain, levels = c("Air", "Water", "Land", "Built", "Social"))
    )

  # Add colors using vectorized approach
  get_fill_color <- function(domain_name, q_num) {
    palette_colors <- DOMAIN_PALETTES[[domain_name]]
    return(palette_colors[q_num])
  }

  combined <- combined %>%
    mutate(
      fill_color = mapply(get_fill_color, as.character(domain), quintile_num)
    )

  p <- ggplot(combined, aes(x = effect, y = quintile_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_density_ridges(
      aes(fill = interaction(domain, quintile_num)),
      alpha = 0.85, scale = 1.5, rel_min_height = 0.005,
      color = "black", linewidth = 0.5, quantile_lines = FALSE
    ) +
    scale_fill_manual(
      values = setNames(
        combined %>% distinct(domain, quintile_num, fill_color) %>% pull(fill_color),
        combined %>% distinct(domain, quintile_num) %>%
          mutate(key = paste(domain, quintile_num, sep = ".")) %>% pull(key)
      ),
      guide = "none"
    ) +
    scale_x_continuous(
      name = "MRD with Posterior Probability Distribution",
      breaks = pretty_breaks(n = 6),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    scale_y_discrete(
      name = NULL, labels = c("Q5", "Q4", "Q3", "Q2"),
      expand = expansion(add = c(0, 1.6))
    ) +
    facet_grid(domain ~ ., scales = "free") +
    theme_publication() +
    theme(
      strip.text.y = element_text(size = 14, face = "bold", margin = margin(t = 5, r = 5, b = 5, l = 5)),
      strip.background = element_rect(fill = "#F5F5F5", color = "#CCCCCC", linewidth = 0.3),
      panel.spacing.y = unit(1, "lines"),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
    )

  cat(sprintf("    Saving to: %s\n", output_file))
  ggsave(output_file, plot = p, width = 7, height = 12, dpi = 300, bg = "white")
}

plot_overall_lags <- function(data, lags, output_file, gradient_colors = NULL) {
  cat(sprintf("  Processing Overall EQI for lags: %s\n", paste(lags, collapse = ", ")))

  combined <- data %>%
    filter(lag %in% lags, domain == "Overall") %>%
    mutate(
      quintile_label = factor(quintile, levels = c("Q5", "Q4", "Q3", "Q2")),
      lag_label = factor(sprintf("Lag %d", lag), levels = sprintf("Lag %d", lags))
    )

  # Use provided gradient colors or default
  if (is.null(gradient_colors)) {
    gradient_colors <- c("#FDE725", "#FCA636", "#E16462", "#B12A90", "#6A00A8", "#320A5E")
  }

  p <- ggplot(combined, aes(x = effect, y = quintile_label, fill = after_stat(x))) +
    geom_density_ridges_gradient(
      alpha = 0.85, scale = 1.4, rel_min_height = 0.005,
      color = "black", linewidth = 0.5, quantile_lines = FALSE
    ) +
    scale_fill_gradientn(
      colors = gradient_colors,
      name = "MRD",
      guide = "none"
    ) +
    scale_x_continuous(
      name = "MRD with Posterior Probability Distribution",
      breaks = pretty_breaks(n = 6),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    scale_y_discrete(
      name = NULL,
      labels = c("Q5", "Q4", "Q3", "Q2"),
      expand = expansion(add = c(0.1, 1.8))
    ) +
    facet_wrap(~lag_label, ncol = 1) +
    theme_publication() +
    theme(
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      panel.spacing = unit(4, "lines"),
      strip.background = element_rect(fill = "#F5F5F5", color = "#CCCCCC", linewidth = 0.3),
      strip.text = element_text(size = 14, face = "bold", margin = margin(t = 5, b = 5))
    )

  cat(sprintf("    Saving to: %s\n", output_file))
  ggsave(output_file, plot = p, width = 7, height = 12, dpi = 300, bg = "white")
}

# ============================================================================
# Main Processing
# ============================================================================

process_outcome <- function(outcome_name, project_root) {
  cat(sprintf("Processing: %s\n", outcome_name))

  output_dir <- file.path(project_root, sprintf("Result/brms_RL_Visualization/%s", outcome_name))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  brms_dir <- file.path(project_root, "Result/brms_RL")
  pattern <- sprintf("^%s_Lag[0-9]+_(Air|Water|Land|Built|Social|Overall)\\.rds$", outcome_name)
  files <- list.files(brms_dir, pattern = pattern, full.names = TRUE)

  if (length(files) == 0) {
    cat(sprintf("ERROR: No RDS files found for '%s'\n", outcome_name))
    return(FALSE)
  }

  cat(sprintf("Found %d RDS files\n\n", length(files)))

  # Load data
  all_data <- load_ridgeline_data(files)

  # Generate 4 lag-specific plots
  cat("Generating lag-specific plots:\n")
  for (current_lag in c(5, 10, 15, 20)) {
    output_file <- file.path(output_dir, sprintf("Lag%d_AllDomains.png", current_lag))
    plot_lag_all_domains(all_data, current_lag, output_file)
  }

  # Generate overall plots with different gradients
  cat("\nGenerating overall EQI plots with multiple color palettes:\n")
  for (palette_name in names(OVERALL_GRADIENTS)) {
    gradient <- OVERALL_GRADIENTS[[palette_name]]

    # 3-lag version
    output_file_3 <- file.path(output_dir, sprintf("Overall_Lag5_10_15_%s.png", palette_name))
    plot_overall_lags(all_data, c(5, 10, 15), output_file_3, gradient_colors = gradient)

    # 4-lag version
    output_file_4 <- file.path(output_dir, sprintf("Overall_Lag5_10_15_20_%s.png", palette_name))
    plot_overall_lags(all_data, c(5, 10, 15, 20), output_file_4, gradient_colors = gradient)
  }

  cat(sprintf("\n✓ Completed %s (4 lag + 8 overall plots = 12 plots)\n", outcome_name))
  return(TRUE)
}

main <- function() {
  project_root <- normalizePath(".")

  cat("\n")
  cat("╔════════════════════════════════════════════════════════════╗\n")
  cat("║   RidgeLine Visualization - brms_RL Results                ║\n")
  cat("║   7 Outcomes × 12 Plots = 84 PNG Files                    ║\n")
  cat("║   (4 lag plots + 8 overall color variants)                ║\n")
  cat("╚════════════════════════════════════════════════════════════╝\n")

  outcomes <- c("Cancer", "CKD", "CLD", "CRD", "CVD", "NDD", "Suicide")

  results <- list()
  for (outcome in outcomes) {
    results[[outcome]] <- process_outcome(outcome, project_root)
  }

  cat("\n")
  cat("════════════════════════════════════════════════════════════\n")
  cat("✅ Generation Complete\n")
  cat("════════════════════════════════════════════════════════════\n")
  cat("Output structure:\n")
  cat("  Result/brms_RL_Visualization/\n")
  cat("    ├── Cancer/\n")
  cat("    │   ├── Lag5_AllDomains.png\n")
  cat("    │   ├── Lag10_AllDomains.png\n")
  cat("    │   ├── Lag15_AllDomains.png\n")
  cat("    │   ├── Lag20_AllDomains.png\n")
  cat("    │   ├── Overall_Lag5_10_15_{Blue,Green,Purple,Blue2}.png\n")
  cat("    │   └── Overall_Lag5_10_15_20_{Blue,Green,Purple,Blue2}.png\n")
  cat("    └── ... (6 more outcomes)\n\n")

  total_success <- sum(unlist(results))
  cat(sprintf("Completed: %d/%d outcomes\n", total_success, length(outcomes)))
  cat("Expected: 84 PNG files total\n\n")
}

if (!interactive()) {
  main()
}
