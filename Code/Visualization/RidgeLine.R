#!/usr/bin/env Rscript
# ============================================================================
# Publication-Quality Gradient Ridgeline Plots
# ============================================================================
# Design for Nature / The Lancet level manuscripts
#
# Design Philosophy:
#   - Aesthetic minimalism: Single-color gradients, clean typography
#   - Information maximalism: Bayesian probabilities, dose-response curves
#   - Visual metaphor: "Risk waves" moving rightward with worsening EQI
#
# Key Features:
#   1. Monochromatic gradient (light → dark) for Q2 → Q5
#   2. Vertical dashed line at MRD = 0 (null effect)
#   3. Posterior probability annotations: P(effect > 0)
#   4. Faceting by lag period or cluster type
#   5. Y-axis: Environmental domains + Overall EQI
#   6. Minimalist theme with Nature/Lancet aesthetics
#
# Usage:
#   Rscript RidgeLine_Publication.R --mode <mode> [options]
#
# Modes:
#   domains_by_lag    All 5 domains + Overall, faceted by lag
#   domains_single    Single lag, all domains stacked
#   overall_evolution Overall EQI across lags (dose-response over time)
#
# Options:
#   --lag       Specific lag (5/10/15) for single mode
#   --output    Output PNG file
#   --width     Plot width in inches (default: 10)
#   --height    Plot height in inches (default: 8)
#   --dpi       Resolution (default: 600 for publication)
#   --palette   Color palette: "warm" (orange-red), "cool" (blue), "neutral" (gray)
#
# Examples:
#   # All domains, faceted by lag (recommended for main figure)
#   Rscript RidgeLine_Publication.R --mode domains_by_lag \
#     --output Result/Ridgeline/Fig_Domains_ByLag.png
#
#   # Single lag, all domains (for supplementary)
#   Rscript RidgeLine_Publication.R --mode domains_single --lag 10 \
#     --output Result/Ridgeline/Fig_Lag10_Domains.png
#
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggridges)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(grid)
  library(gridExtra)
  library(gtable)
  library(patchwork)
})

# ============================================================================
# Configuration
# ============================================================================

# Nature/Lancet color palettes
# Domain-specific color palettes (light to dark gradients)
DOMAIN_PALETTES <- list(
  "Overall EQI" = c("#E8E8E8", "#BDBDBD", "#8C8C8C", "#525252", "#1A1A1A"), # Gray scale
  "Air" = c("#FFF3E0", "#FFCC80", "#FFA726", "#FB8C00", "#E65100"), # Orange
  "Water" = c("#E3F2FD", "#90CAF9", "#42A5F5", "#1E88E5", "#0D47A1"), # Blue
  "Land" = c("#E8F5E9", "#A5D6A7", "#66BB6A", "#43A047", "#2E7D32"), # Green
  "Built" = c("#F5F5F5", "#D9D9D9", "#BDBDBD", "#969696", "#525252"), # Grey
  "Social" = c("#D4E3F0", "#A8C7E0", "#7DAAD0", "#5281A8", "#2E4E74") # Dark Blue
)

# Legacy palettes (kept for backward compatibility)
PALETTES <- list(
  warm = c("#FED9B7", "#F4A261", "#E76F51", "#D62828", "#9B2226"),
  cool = c("#BDE0FE", "#A2D2FF", "#6A9FD4", "#4A7BA7", "#2C5F8D"),
  neutral = c("#D9D9D9", "#BDBDBD", "#969696", "#636363", "#252525"),
  green = c("#D4EDDA", "#A8D5BA", "#6EBF8B", "#3D9970", "#1E6F5C"),
  blue = c("#E3F2FD", "#90CAF9", "#42A5F5", "#1E88E5", "#0D47A1")
)

# Typography (Nature/Lancet standards)
FONT_FAMILY <- "Helvetica"
BASE_SIZE <- 14
TITLE_SIZE <- 14
AXIS_SIZE <- 14
LABEL_SIZE <- 14

# ============================================================================
# Utility Functions
# ============================================================================

#' Calculate posterior probability of positive effect
#' @param draws Numeric vector of posterior draws
#' @return Probability that effect > 0
calc_prob_positive <- function(draws) {
  sum(draws > 0) / length(draws)
}

#' Format probability for annotation (vectorized)
#' @param p Probability value or vector
#' @return Formatted string or character vector
format_prob <- function(p) {
  sapply(p, function(x) {
    if (x >= 0.999) {
      return("P>0.999")
    }
    if (x >= 0.99) {
      return(sprintf("P=%.3f", x))
    }
    if (x >= 0.95) {
      return(sprintf("P=%.2f", x))
    }
    if (x <= 0.001) {
      return("P<0.001")
    }
    if (x <= 0.01) {
      return(sprintf("P=%.3f", x))
    }
    return(sprintf("P=%.2f", x))
  }, USE.NAMES = FALSE)
}

#' Load and process ridgeline data
#' @param files Vector of .rds file paths
#' @param include_overall Logical, include overall EQI model
#' @return Processed data frame
load_ridgeline_data <- function(files, include_overall = TRUE) {
  all_data <- list()

  for (f in files) {
    data <- readRDS(f)
    lag <- data$metadata$lag

    # Determine if this is Overall or MultiDomain
    is_overall <- is.null(data$draws_long$domain)

    if (is_overall) {
      if (!include_overall) next

      draws <- data$draws_long %>%
        mutate(
          domain = "Overall EQI",
          lag = lag
        )
    } else {
      draws <- data$draws_long %>%
        mutate(lag = lag)
    }

    all_data[[length(all_data) + 1]] <- draws
  }

  combined <- bind_rows(all_data)
  return(combined)
}

#' Calculate summary statistics with probabilities
#' @param data Data frame with effect draws
#' @return Summary data frame
calc_summary_with_prob <- function(data) {
  data %>%
    group_by(lag, domain, quintile) %>%
    summarise(
      mean = mean(effect),
      median = median(effect),
      sd = sd(effect),
      q025 = quantile(effect, 0.025),
      q975 = quantile(effect, 0.975),
      prob_positive = calc_prob_positive(effect),
      .groups = "drop"
    )
}

# ============================================================================
# Theme: Nature/Lancet Publication Style
# ============================================================================

theme_publication <- function(base_size = BASE_SIZE, base_family = FONT_FAMILY) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # Remove clutter
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "#E0E0E0", linewidth = 0.3),

      # Clean background
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),

      # Typography
      text = element_text(family = base_family, color = "#000000"),
      plot.title = element_text(
        size = TITLE_SIZE,
        face = "bold",
        hjust = 0,
        margin = margin(b = 10)
      ),
      plot.subtitle = element_text(
        size = BASE_SIZE,
        color = "#333333",
        hjust = 0,
        margin = margin(b = 15)
      ),
      plot.caption = element_text(
        size = LABEL_SIZE,
        color = "#666666",
        hjust = 0,
        margin = margin(t = 10)
      ),

      # Axes
      axis.title.x = element_text(
        size = AXIS_SIZE,
        face = "bold",
        margin = margin(t = 8)
      ),
      axis.title.y = element_text(
        size = AXIS_SIZE,
        face = "bold",
        margin = margin(r = 8)
      ),
      axis.text = element_text(size = LABEL_SIZE, color = "#000000"),
      axis.line.x = element_line(color = "#000000", linewidth = 0.5),
      axis.ticks.x = element_line(color = "#000000", linewidth = 0.3),
      axis.ticks.length = unit(0.15, "cm"),

      # Facets (clean, minimal)
      strip.text = element_text(
        size = BASE_SIZE,
        face = "bold",
        color = "#000000",
        margin = margin(b = 5)
      ),
      strip.background = element_rect(fill = "#F5F5F5", color = NA),

      # Legend
      legend.position = "right",
      legend.title = element_text(size = AXIS_SIZE, face = "bold"),
      legend.text = element_text(size = LABEL_SIZE),
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.8, "cm"),

      # Margins
      plot.margin = margin(15, 15, 10, 10)
    )
}

# ============================================================================
# Mode 1: All Domains by Lag (Main Figure)
# ============================================================================

plot_domains_by_lag <- function(project_root, output_file, palette = "warm") {
  cat("\n════════════════════════════════════════════════════════════\n")
  cat("Mode: Domains by Lag (Publication Main Figure)\n")
  cat("════════════════════════════════════════════════════════════\n\n")

  # Load all MultiDomain files
  ridge_dir <- file.path(project_root, "Result/Ridgeline")
  files_multi <- list.files(ridge_dir, pattern = "MultiDomain\\.rds$", full.names = TRUE)
  files_overall <- list.files(ridge_dir, pattern = "Overall\\.rds$", full.names = TRUE)

  if (length(files_multi) == 0) {
    stop("No MultiDomain files found")
  }

  cat("Loading MultiDomain models:", length(files_multi), "files\n")
  cat("Loading Overall models:", length(files_overall), "files\n\n")

  # Load and combine
  data_multi <- load_ridgeline_data(files_multi, include_overall = FALSE)
  data_overall <- load_ridgeline_data(files_overall, include_overall = TRUE)

  combined <- bind_rows(data_multi, data_overall) %>%
    mutate(
      # Quintile labels: Q2 at TOP (reverse order for Y-axis)
      quintile_label = factor(quintile, levels = c("Q2", "Q3", "Q4", "Q5")),
      quintile_num = as.numeric(factor(quintile, levels = c("Q2", "Q3", "Q4", "Q5"))),
      # Domain order for X-axis (columns)
      domain = factor(domain, levels = c("Overall EQI", "Air", "Water", "Land", "Built", "Social")),
      # Lag labels: 5 at TOP (reverse order for Y-axis)
      lag_label = factor(sprintf("Lag %d", lag),
        levels = sprintf("Lag %d", c(15, 10, 5)) # Reversed so 5 appears at top
      )
    )

  # Calculate summary with probabilities
  summary_stats <- calc_summary_with_prob(combined) %>%
    mutate(
      prob_label = format_prob(prob_positive),
      domain = factor(domain, levels = c("Overall EQI", "Air", "Water", "Land", "Built", "Social")),
      lag_label = factor(sprintf("Lag %d", lag),
        levels = sprintf("Lag %d", c(15, 10, 5))
      )
    )

  cat("Summary statistics calculated\n")
  cat("Domains:", paste(levels(combined$domain), collapse = ", "), "\n")
  cat("Lags:", paste(unique(combined$lag), collapse = ", "), "\n\n")

  # Create domain-specific color mapping
  cat("Creating publication-quality ridgeline plot with domain-specific colors...\n")

  # Add domain-specific colors to the data
  combined <- combined %>%
    rowwise() %>%
    mutate(
      domain_colors = list(DOMAIN_PALETTES[[as.character(domain)]]),
      fill_color = domain_colors[[quintile_num]]
    ) %>%
    ungroup()

  # Create separate plots for each lag period (rows)
  lag_plots <- list()

  for (lag_val in c(5, 10, 15)) {
    lag_data <- combined %>% filter(lag == lag_val)

    p_lag <- ggplot(lag_data, aes(x = effect, y = quintile_label)) +
      geom_density_ridges(
        aes(fill = interaction(domain, quintile_num)),
        alpha = 0.85,
        scale = 1.1,
        rel_min_height = 0.005,
        color = "white",
        linewidth = 0.4,
        quantile_lines = TRUE,
        quantiles = 0.5,
        vline_color = "white",
        vline_width = 0.3
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "#000000",
        linewidth = 0.6,
        alpha = 0.7
      ) +
      scale_fill_manual(
        values = setNames(
          lag_data %>% distinct(domain, quintile_num, fill_color) %>% pull(fill_color),
          lag_data %>% distinct(domain, quintile_num) %>%
            mutate(key = paste(domain, quintile_num, sep = ".")) %>% pull(key)
        ),
        guide = "none"
      ) +
      scale_x_continuous(
        name = if (lag_val == 15) "Effect on Age-Adjusted Mortality Rate (AAMR)" else NULL,
        breaks = pretty_breaks(n = 6),
        expand = expansion(mult = c(0.05, 0.05))
      ) +
      scale_y_discrete(
        name = NULL,
        expand = expansion(add = c(0.2, 0.5))
      ) +
      facet_grid(. ~ domain, scales = "free", switch = "x") +
      theme_publication() +
      theme(
        strip.text.x = element_text(size = 11, face = "bold", margin = margin(b = 5, t = 5)),
        strip.background.x = element_rect(fill = "#F5F5F5", color = "#CCCCCC", linewidth = 0.3),
        panel.spacing.x = unit(0.3, "lines"),
        axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        legend.position = "none",
        plot.margin = margin(5, 10, 5, 10)
      )

    # Add horizontal lag label strip
    lag_label_grob <- textGrob(
      sprintf("Lag %d years", lag_val),
      gp = gpar(fontsize = 14, fontface = "bold", col = "#000000")
    )

    lag_strip <- rectGrob(
      gp = gpar(fill = "gray", col = "black", lwd = 1)
    )

    # Create a frame with the lag strip
    lag_frame <- gTree(children = gList(lag_strip, lag_label_grob))

    lag_plots[[as.character(lag_val)]] <- list(
      label = lag_frame,
      plot = p_lag
    )
  }

  # Combine plots with patchwork
  p <- (
    wrap_elements(lag_plots[["5"]]$label) / lag_plots[["5"]]$plot /
      wrap_elements(lag_plots[["10"]]$label) / lag_plots[["10"]]$plot /
      wrap_elements(lag_plots[["15"]]$label) / lag_plots[["15"]]$plot
  ) +
    plot_layout(heights = c(0.6, 5, 0.6, 5, 0.6, 5)) +
    plot_annotation(
      title = "Environmental Quality Index and Cancer Mortality: Domain-Specific Dose-Response",
      subtitle = "Posterior distributions of quintile effects (Q2-Q5 vs. Q1 reference) across lag periods",
      caption = "Note: Dashed lines indicate null effect (MRD = 0). Color gradients represent exposure quintiles (lighter = lower exposure).\nWhite lines show posterior median. Gray = Overall EQI; Colors = domain-specific effects.",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        plot.caption = element_text(size = 8, hjust = 0)
      )
    )

  # Save
  cat("\nSaving to:", output_file, "\n")
  ggsave(
    output_file,
    plot = p,
    width = 20,
    height = 12,
    dpi = 600,
    bg = "white"
  )

  cat("✓ Publication figure saved\n")

  # Print summary table
  cat("\n════════════════════════════════════════════════════════════\n")
  cat("Summary Statistics\n")
  cat("════════════════════════════════════════════════════════════\n\n")

  summary_table <- summary_stats %>%
    arrange(lag, domain, quintile) %>%
    mutate(
      ci_95 = sprintf("[%.1f, %.1f]", q025, q975)
    ) %>%
    select(lag, domain, quintile, median, ci_95, prob_positive)

  print(summary_table, n = 50)

  cat("\n")
}

# ============================================================================
# Mode 2: Single Lag, All Domains (Supplementary Figure)
# ============================================================================

plot_domains_single_lag <- function(project_root, lag, output_file, palette = "warm") {
  cat("\n════════════════════════════════════════════════════════════\n")
  cat("Mode: Single Lag, All Domains (Supplementary Figure)\n")
  cat("════════════════════════════════════════════════════════════\n\n")

  ridge_dir <- file.path(project_root, "Result/Ridgeline")

  # Find files for this lag
  file_multi <- list.files(ridge_dir,
    pattern = sprintf("Lag%d_MultiDomain\\.rds$", lag),
    full.names = TRUE
  )
  file_overall <- list.files(ridge_dir,
    pattern = sprintf("Lag%d_Overall\\.rds$", lag),
    full.names = TRUE
  )

  if (length(file_multi) == 0) {
    stop(sprintf("No MultiDomain file found for Lag %d", lag))
  }

  cat("Loading Lag", lag, "models\n\n")

  # Load data
  data_multi <- load_ridgeline_data(file_multi, include_overall = FALSE)
  data_overall <- load_ridgeline_data(file_overall, include_overall = TRUE)

  combined <- bind_rows(data_multi, data_overall) %>%
    mutate(
      quintile_label = factor(quintile, levels = c("Q5", "Q4", "Q3", "Q2")),
      quintile_num = as.numeric(factor(quintile, levels = c("Q2", "Q3", "Q4", "Q5"))),
      domain = factor(domain, levels = c("Overall EQI", "Air", "Water", "Land", "Built", "Social"))
    )

  # Summary with probabilities
  summary_stats <- calc_summary_with_prob(combined) %>%
    mutate(
      prob_label = format_prob(prob_positive),
      domain = factor(domain, levels = c("Overall EQI", "Air", "Water", "Land", "Built", "Social"))
    )

  # Get palette
  colors <- PALETTES[[palette]]

  # Create plot
  cat("Creating publication-quality ridgeline plot...\n")

  p <- ggplot(combined, aes(x = effect, y = domain, fill = quintile_num)) +

    # Ridgelines
    geom_density_ridges(
      aes(height = after_stat(density)),
      alpha = 0.85,
      scale = 0.95,
      rel_min_height = 0.005,
      color = "white",
      linewidth = 0.4,
      quantile_lines = TRUE,
      quantiles = 0.5,
      vline_color = "white",
      vline_width = 0.3
    ) +

    # Zero line
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "#000000",
      linewidth = 0.6,
      alpha = 0.7
    ) +

    # Probability annotations
    geom_text(
      data = summary_stats,
      aes(
        x = median, y = as.numeric(domain) + 0.3 * (as.numeric(factor(quintile, levels = c("Q2", "Q3", "Q4", "Q5"))) - 2.5),
        label = prob_label
      ),
      size = 2.2,
      color = "#000000",
      fontface = "italic",
      family = FONT_FAMILY,
      hjust = 0.5,
      inherit.aes = FALSE
    ) +

    # Gradient fill
    scale_fill_gradientn(
      colors = colors,
      name = "Quintile",
      breaks = 1:4,
      labels = c("Q2", "Q3", "Q4", "Q5"),
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barheight = unit(3, "cm"),
        barwidth = unit(0.3, "cm")
      )
    ) +

    # Axes
    scale_x_continuous(
      name = "Effect on Age-Adjusted Mortality Rate (AAMR)",
      breaks = pretty_breaks(n = 10)
    ) +
    scale_y_discrete(
      name = "EQI Domain",
      expand = expansion(add = c(0.3, 1.2))
    ) +

    # Theme
    theme_publication() +
    theme(
      legend.position = "right"
    ) +

    # Labels
    labs(
      title = sprintf("Environmental Quality Domains and Cancer Mortality (Lag %d Years)", lag),
      subtitle = "Posterior distributions showing dose-response gradients across EQI quintiles",
      caption = "Note: Dashed line indicates null effect. P-values show posterior probability of positive effect.\nGradient shading represents exposure quintiles (Q2-Q5 vs. Q1 reference). White line shows median."
    )

  # Save
  cat("\nSaving to:", output_file, "\n")
  ggsave(
    output_file,
    plot = p,
    width = 10,
    height = 8,
    dpi = 600,
    bg = "white"
  )

  cat("✓ Supplementary figure saved\n\n")
}

# ============================================================================
# Mode 3: Overall Evolution (Dose-Response Over Time)
# ============================================================================

plot_overall_evolution <- function(project_root, output_file, palette = "blue") {
  cat("\n════════════════════════════════════════════════════════════\n")
  cat("Mode: Overall EQI Evolution Across Lags\n")
  cat("════════════════════════════════════════════════════════════\n\n")

  ridge_dir <- file.path(project_root, "Result/Ridgeline")
  files_overall <- list.files(ridge_dir, pattern = "Overall\\.rds$", full.names = TRUE)

  if (length(files_overall) == 0) {
    stop("No Overall model files found")
  }

  cat("Loading Overall models:", length(files_overall), "files\n\n")

  # Load data
  combined <- load_ridgeline_data(files_overall, include_overall = TRUE) %>%
    mutate(
      quintile_label = factor(quintile, levels = c("Q5", "Q4", "Q3", "Q2")),
      quintile_num = as.numeric(factor(quintile, levels = c("Q2", "Q3", "Q4", "Q5"))),
      lag_label = factor(sprintf("Lag %d", lag),
        levels = sprintf("Lag %d", c(5, 10, 15))
      )
    )

  # Summary
  summary_stats <- calc_summary_with_prob(combined) %>%
    mutate(
      prob_label = format_prob(prob_positive),
      lag_label = factor(sprintf("Lag %d", lag),
        levels = sprintf("Lag %d", c(5, 10, 15))
      )
    )

  # Palette
  colors <- PALETTES[[palette]]

  # Plot
  cat("Creating temporal evolution plot...\n")

  p <- ggplot(combined, aes(x = effect, y = quintile_label, fill = quintile_num)) +
    geom_density_ridges(
      alpha = 0.85,
      scale = 2.5,
      rel_min_height = 0.005,
      color = "white",
      linewidth = 0.5,
      quantile_lines = FALSE
    ) +
    scale_fill_gradientn(
      colors = colors,
      name = "Quintile",
      breaks = 4:1,
      labels = c("Q5", "Q4", "Q3", "Q2"),
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barheight = unit(2.5, "cm"),
        barwidth = unit(0.3, "cm")
      )
    ) +
    scale_x_continuous(
      name = "Mortality Rate Difference",
      limits = c(0, 16),
      breaks = seq(0, 16, 2),
      expand = c(0, 0)
    ) +
    scale_y_discrete(
      name = "Posterior Probability",
      labels = NULL,
      expand = expansion(add = c(0.2, 0.8))
    ) +
    facet_wrap(~lag_label, ncol = 1) +
    theme_publication() +
    theme(
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      axis.text.y = element_blank()
    ) +
    labs(
      title = NULL,
      subtitle = NULL,
      caption = NULL
    )

  # Save
  cat("\nSaving to:", output_file, "\n")
  ggsave(
    output_file,
    plot = p,
    width = 10,
    height = 10,
    dpi = 300,
    bg = "white"
  )

  cat("✓ Evolution figure saved\n\n")
}

# ============================================================================
# Mode 3b: Overall Evolution v2 (with Y-axis labels and MRD color mapping)
# ============================================================================

plot_overall_evolution_v2 <- function(project_root, output_file) {
  cat("\n════════════════════════════════════════════════════════════\n")
  cat("Mode: Overall EQI Evolution Across Lags (Version 2)\n")
  cat("════════════════════════════════════════════════════════════\n\n")

  ridge_dir <- file.path(project_root, "Result/Ridgeline")
  files_overall <- list.files(ridge_dir, pattern = "Overall\\.rds$", full.names = TRUE)

  if (length(files_overall) == 0) {
    stop("No Overall model files found")
  }

  cat("Loading Overall models:", length(files_overall), "files\n\n")

  # Load data
  combined <- load_ridgeline_data(files_overall, include_overall = TRUE) %>%
    mutate(
      quintile_label = factor(quintile, levels = c("Q5", "Q4", "Q3", "Q2")),
      quintile_num = as.numeric(factor(quintile, levels = c("Q2", "Q3", "Q4", "Q5"))),
      lag_label = factor(sprintf("Lag %d", lag),
        levels = sprintf("Lag %d", c(5, 10, 15))
      )
    )

  # Summary
  summary_stats <- calc_summary_with_prob(combined) %>%
    mutate(
      prob_label = format_prob(prob_positive),
      lag_label = factor(sprintf("Lag %d", lag),
        levels = sprintf("Lag %d", c(5, 10, 15))
      )
    )

  # Plot
  cat("Creating temporal evolution plot (v2)...\n")

  p <- ggplot(combined, aes(x = effect, y = quintile_label, fill = after_stat(x))) +
    geom_density_ridges_gradient(
      alpha = 0.85,
      scale = 1.4,
      rel_min_height = 0.005,
      color = "black",
      linewidth = 0.5,
      quantile_lines = FALSE
    ) +
    scale_fill_gradientn(
      colors = c("#FDE725", "#FCA636", "#E16462", "#B12A90", "#6A00A8", "#320A5E"),
      name = "MRD",
      limits = c(0, 16),
      breaks = seq(0, 16, 4),
      guide = "none"
    ) +
    scale_x_continuous(
      name = "MRD with Posterior Probability Distribution",
      limits = c(0, 16),
      breaks = seq(0, 16, 2),
      expand = c(0, 0)
    ) +
    scale_y_discrete(
      name = NULL,
      labels = c("Q5", "Q4", "Q3", "Q2"),
      expand = expansion(add = c(0, 1.8))
    ) +
    facet_wrap(~lag_label, ncol = 1) +
    theme_publication() +
    theme(
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      axis.text.y = element_text(size = 14, family = FONT_FAMILY),
      panel.spacing = unit(4, "lines"),
      plot.margin = margin(20, 20, 20, 20),
      strip.background = element_rect(fill = "#F5F5F5", color = "#CCCCCC", linewidth = 0.3),
      strip.text = element_text(size = 14, face = "bold", margin = margin(t = 5, b = 5)),
      plot.title = element_text(size = 14, face = "bold", hjust = 0)
    ) +
    labs(
      title = "(A)",
      subtitle = NULL,
      caption = NULL
    )

  # Save
  cat("\nSaving to:", output_file, "\n")
  ggsave(
    output_file,
    plot = p,
    width = 7,
    height = 12,
    dpi = 300,
    bg = "white"
  )

  cat("✓ Evolution figure v2 saved\n\n")
}

# ============================================================================
# Mode 3c: Overall Evolution v2 Low Saturation (with lower saturation colors)
# ============================================================================

plot_overall_evolution_v2_low_sat <- function(project_root, output_file) {
  cat("\n════════════════════════════════════════════════════════════\n")
  cat("Mode: Overall EQI Evolution Across Lags (Version 2 - Low Saturation)\n")
  cat("════════════════════════════════════════════════════════════\n\n")

  ridge_dir <- file.path(project_root, "Result/Ridgeline")
  files_overall <- list.files(ridge_dir, pattern = "Overall\\.rds$", full.names = TRUE)

  if (length(files_overall) == 0) {
    stop("No Overall model files found")
  }

  cat("Loading Overall models:", length(files_overall), "files\n\n")

  # Load data
  combined <- load_ridgeline_data(files_overall, include_overall = TRUE) %>%
    mutate(
      quintile_label = factor(quintile, levels = c("Q5", "Q4", "Q3", "Q2")),
      quintile_num = as.numeric(factor(quintile, levels = c("Q2", "Q3", "Q4", "Q5"))),
      lag_label = factor(sprintf("Lag %d", lag),
        levels = sprintf("Lag %d", c(5, 10, 15))
      )
    )

  # Summary
  summary_stats <- calc_summary_with_prob(combined) %>%
    mutate(
      prob_label = format_prob(prob_positive),
      lag_label = factor(sprintf("Lag %d", lag),
        levels = sprintf("Lag %d", c(5, 10, 15))
      )
    )

  # Plot with lower saturation colors
  cat("Creating temporal evolution plot (v2 - low saturation)...\n")

  p <- ggplot(combined, aes(x = effect, y = quintile_label, fill = after_stat(x))) +
    geom_density_ridges_gradient(
      alpha = 0.85,
      scale = 1.4,
      rel_min_height = 0.005,
      color = "black",
      linewidth = 0.5,
      quantile_lines = FALSE
    ) +
    scale_fill_gradientn(
      colors = c("#FFF9E0", "#FFD9A8", "#F0B8B0", "#D8A8C8", "#B8A0D0", "#A090B8"),
      name = "MRD",
      limits = c(0, 16),
      breaks = seq(0, 16, 4),
      guide = "none"
    ) +
    scale_x_continuous(
      name = "MRD with Posterior Probability Distribution",
      limits = c(0, 16),
      breaks = seq(0, 16, 2),
      expand = c(0, 0)
    ) +
    scale_y_discrete(
      name = NULL,
      labels = c("Q5", "Q4", "Q3", "Q2"),
      expand = expansion(add = c(0.2, 1.2))
    ) +
    facet_wrap(~lag_label, ncol = 1) +
    theme_publication() +
    theme(
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      axis.text.y = element_text(size = 14, family = FONT_FAMILY),
      panel.spacing = unit(5, "lines"),
      plot.margin = margin(20, 20, 20, 20),
      strip.background = element_rect(fill = "#F5F5F5", color = "#CCCCCC", linewidth = 0.3),
      strip.text = element_text(size = 14, face = "bold", margin = margin(t = 10, b = 10)),
      plot.title = element_text(size = 14, face = "bold", hjust = 0)
    ) +
    labs(
      title = "(A)",
      subtitle = NULL,
      caption = NULL
    )

  # Save
  cat("\nSaving to:", output_file, "\n")
  ggsave(
    output_file,
    plot = p,
    width = 7,
    height = 12,
    dpi = 300,
    bg = "white"
  )

  cat("✓ Evolution figure v2 (low saturation) saved\n\n")
}

# ============================================================================
# Mode 4: Lags by Domains (Transposed Figure)
# ============================================================================

plot_lags_by_domains <- function(project_root, output_file, palette = "green") {
  cat("\n════════════════════════════════════════════════════════════\n")
  cat("Mode: Lags by Domains (Split by Lag)\n")
  cat("════════════════════════════════════════════════════════════\n\n")

  # Load all MultiDomain files
  ridge_dir <- file.path(project_root, "Result/Ridgeline")
  files_multi <- list.files(ridge_dir, pattern = "MultiDomain\\.rds$", full.names = TRUE)
  files_overall <- list.files(ridge_dir, pattern = "Overall\\.rds$", full.names = TRUE)

  if (length(files_multi) == 0) {
    stop("No MultiDomain files found")
  }

  cat("Loading MultiDomain models:", length(files_multi), "files\n")
  cat("Loading Overall models:", length(files_overall), "files\n\n")

  # Load and combine
  data_multi <- load_ridgeline_data(files_multi, include_overall = FALSE)

  combined_all <- data_multi %>%
    filter(domain != "Overall EQI") %>%
    mutate(
      # Quintile labels: Q5 at TOP (reverse order for Y-axis)
      quintile_label = factor(quintile, levels = c("Q5", "Q4", "Q3", "Q2")),
      quintile_num = as.numeric(factor(quintile, levels = c("Q2", "Q3", "Q4", "Q5"))),
      # Domain order for Y-axis (rows)
      domain = factor(domain, levels = c("Air", "Water", "Land", "Built", "Social"))
    )

  # Calculate summary with probabilities
  summary_stats <- calc_summary_with_prob(combined_all) %>%
    mutate(
      prob_label = format_prob(prob_positive),
      domain = factor(domain, levels = c("Air", "Water", "Land", "Built", "Social"))
    )

  cat("Summary statistics calculated\n")

  # Add domain-specific colors
  combined_all <- combined_all %>%
    rowwise() %>%
    mutate(
      domain_colors = list(DOMAIN_PALETTES[[as.character(domain)]]),
      fill_color = domain_colors[[quintile_num]]
    ) %>%
    ungroup()

  # Loop over lags
  for (current_lag in c(5, 10, 15)) {
    cat(sprintf("\nProcessing Lag %d...\n", current_lag))

    combined <- combined_all %>% filter(lag == current_lag)

    # Set x-axis limits based on lag
    x_limits <- if (current_lag == 15) c(-10, 15) else c(-10, 20)
    x_breaks <- seq(x_limits[1], x_limits[2], 5)

    p <- ggplot(combined, aes(x = effect, y = quintile_label)) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "black",
        linewidth = 0.5
      ) +
      geom_density_ridges(
        aes(fill = interaction(domain, quintile_num)),
        alpha = 0.85,
        scale = 1.5,
        rel_min_height = 0.005,
        color = "black",
        linewidth = 0.5,
        quantile_lines = FALSE
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
        limits = x_limits,
        breaks = x_breaks,
        expand = c(0, 0)
      ) +
      scale_y_discrete(
        name = NULL,
        labels = c("Q5", "Q4", "Q3", "Q2"),
        expand = expansion(add = c(0, 1.6))
      ) +
      facet_grid(domain ~ ., scales = "free") +
      theme_publication() +
      theme(
        strip.text.y = element_text(size = 14, face = "bold", margin = margin(t = 5, r = 5, b = 5, l = 5)),
        strip.background = element_rect(fill = "#F5F5F5", color = "#CCCCCC", linewidth = 0.3),
        panel.spacing.y = unit(1, "lines"),
        panel.border = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, family = FONT_FAMILY),
        legend.position = "none",
        plot.margin = margin(20, 20, 20, 20),
        plot.title = element_text(size = 14, face = "bold", hjust = 0)
      ) +
      labs(
        title = "(B)",
        subtitle = NULL,
        caption = NULL
      )

    # Save
    current_output <- file.path(project_root, sprintf("Result/Ridgeline/Ridgeline_lag%d_Domain.png", current_lag))
    cat("Saving to:", current_output, "\n")
    ggsave(
      current_output,
      plot = p,
      width = 7,
      height = 12,
      dpi = 300,
      bg = "white"
    )
  }

  cat("✓ Split lag figures saved\n")

  # Print summary table
  cat("\n════════════════════════════════════════════════════════════\n")
  cat("Summary Statistics\n")
  cat("════════════════════════════════════════════════════════════\n\n")

  summary_table <- summary_stats %>%
    arrange(lag, domain, quintile) %>%
    mutate(
      ci_95 = sprintf("[%.1f, %.1f]", q025, q975)
    ) %>%
    select(lag, domain, quintile, median, ci_95, prob_positive)

  print(summary_table, n = 50)

  cat("\n")
}

# ============================================================================
# Command-Line Argument Parsing
# ============================================================================

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  parsed <- list(
    mode = "overall_evolution",
    lag = NULL,
    output = NULL,
    width = 10,
    height = 10,
    dpi = 600,
    palette = "blue"
  )

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--mode" && i < length(args)) {
      parsed$mode <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--lag" && i < length(args)) {
      parsed$lag <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--output" && i < length(args)) {
      parsed$output <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--width" && i < length(args)) {
      parsed$width <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--height" && i < length(args)) {
      parsed$height <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--dpi" && i < length(args)) {
      parsed$dpi <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--palette" && i < length(args)) {
      parsed$palette <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }

  return(parsed)
}

# ============================================================================
# Main Execution
# ============================================================================

main <- function() {
  project_root <- normalizePath(".")

  cat("\n")
  cat("╔════════════════════════════════════════════════════════════╗\n")
  cat("║   Publication-Quality Ridgeline Plots                     ║\n")
  cat("║   Design: Nature / The Lancet Standards                   ║\n")
  cat("║   Running all modes automatically                         ║\n")
  cat("╚════════════════════════════════════════════════════════════╝\n")

  # Run Ridgeline_lag2.png (overall evolution v2)
  plot_overall_evolution_v2(
    project_root,
    file.path(project_root, "Result/Ridgeline/Ridgeline_lag2.png")
  )

  # Run Ridgeline_lag2_low_saturation.png (overall evolution v2 low saturation)
  plot_overall_evolution_v2_low_sat(
    project_root,
    file.path(project_root, "Result/Ridgeline/Ridgeline_lag2_low_saturation.png")
  )

  # Run Ridgeline_lag{5/10/15}_Domain.png (lags by domains)
  plot_lags_by_domains(project_root, NULL, "blue")

  cat("\n")
  cat("════════════════════════════════════════════════════════════\n")
  cat("✅ All publication figures generated successfully\n")
  cat("════════════════════════════════════════════════════════════\n\n")
  cat("Generated files:\n")
  cat("  - Result/Ridgeline/Ridgeline_lag2.png\n")
  cat("  - Result/Ridgeline/Ridgeline_lag2_low_saturation.png\n")
  cat("  - Result/Ridgeline/Ridgeline_lag5_Domain.png\n")
  cat("  - Result/Ridgeline/Ridgeline_lag10_Domain.png\n")
  cat("  - Result/Ridgeline/Ridgeline_lag15_Domain.png\n\n")
}

# Run if not in interactive mode
if (!interactive()) {
  main()
}
