#!/usr/bin/env Rscript
# Ridgeline plot visualization for posterior distributions
# Visualizes Q2-Q5 quintile effects from Bayesian analysis
# Input:  Result/Ridgeline/C00_C97_Ridge_Test.rds
# Output: Result/Ridgeline/C00_C97_Ridge_Plot.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggridges)
  library(dplyr)
  library(scales)
})

cat("========================================\n")
cat("Ridgeline Plot Visualization\n")
cat("========================================\n\n")

# ============================================================================
# Configuration
# ============================================================================
project_root <- normalizePath(".")
input_file <- file.path(project_root, "Result/Ridgeline/C00_C97_Ridge_Test.rds")
output_file <- file.path(project_root, "Result/Ridgeline/C00_C97_Ridge_Plot.png")

# ============================================================================
# Load data
# ============================================================================
cat("Loading data...\n")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

ridge_data <- readRDS(input_file)

cat("Data loaded successfully\n")
cat("Metadata:\n")
cat("  Cancer type:", ridge_data$metadata$cancer_type, "\n")
cat("  EQI period: ", ridge_data$metadata$eqi_period, "\n")
cat("  AAMR period:", ridge_data$metadata$aamr_period, "\n")
cat("  Lag:        ", ridge_data$metadata$lag, "\n")
cat("  N draws:    ", ridge_data$metadata$n_draws, "\n")
cat("  N obs:      ", ridge_data$metadata$n_obs, "\n")
cat("\n")

# ============================================================================
# Check convergence
# ============================================================================
cat("Convergence diagnostics:\n")
cat("  Max R-hat:       ", sprintf("%.4f", ridge_data$metadata$convergence$max_rhat), "\n")
cat("  Min ESS (bulk):  ", sprintf("%.0f", ridge_data$metadata$convergence$min_ess_bulk), "\n")

if (ridge_data$metadata$convergence$max_rhat > 1.05) {
  cat("  ⚠️  Warning: Some chains may not have converged well\n")
} else {
  cat("  ✓ Convergence looks good\n")
}
cat("\n")

# ============================================================================
# Summary statistics
# ============================================================================
cat("Summary statistics:\n")
print(ridge_data$summary)
cat("\n")

# ============================================================================
# Prepare data for plotting
# ============================================================================
cat("Preparing plot data...\n")

plot_data <- ridge_data$draws_long %>%
  mutate(
    quintile_label = case_when(
      quintile == "Q2" ~ "Q2 vs Q1",
      quintile == "Q3" ~ "Q3 vs Q1",
      quintile == "Q4" ~ "Q4 vs Q1",
      quintile == "Q5" ~ "Q5 vs Q1",
      TRUE ~ as.character(quintile)
    ),
    quintile_label = factor(quintile_label,
      levels = c("Q5 vs Q1", "Q4 vs Q1", "Q3 vs Q1", "Q2 vs Q1")
    )
  )

# Get summary stats for annotations
summary_stats <- ridge_data$summary %>%
  mutate(
    quintile_label = case_when(
      quintile == "Q2" ~ "Q2 vs Q1",
      quintile == "Q3" ~ "Q3 vs Q1",
      quintile == "Q4" ~ "Q4 vs Q1",
      quintile == "Q5" ~ "Q5 vs Q1",
      TRUE ~ as.character(quintile)
    ),
    quintile_label = factor(quintile_label,
      levels = c("Q5 vs Q1", "Q4 vs Q1", "Q3 vs Q1", "Q2 vs Q1")
    ),
    label = sprintf("%.2f [%.2f, %.2f]", mean, q025, q975)
  )

# ============================================================================
# Create ridgeline plot with enhanced aesthetics
# ============================================================================
cat("Creating ridgeline plot...\n")

# Define custom color palette (gradient from light to dark, cool to warm)
custom_colors <- c(
  "Q2 vs Q1" = "#3B9AB2", # Teal blue
  "Q3 vs Q1" = "#78B7C5", # Sky blue
  "Q4 vs Q1" = "#EBCC2A", # Golden yellow
  "Q5 vs Q1" = "#E1AF00" # Dark gold/amber
)

p <- ggplot(plot_data, aes(x = effect, y = quintile_label, fill = quintile_label)) +

  # Main ridgeline densities with gradient shading
  geom_density_ridges(
    alpha = 0.85,
    scale = 1.3,
    rel_min_height = 0.01,
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.5, 0.975),
    color = "white",
    linewidth = 0.6,
    vline_width = 0.4,
    vline_color = "gray30"
  ) +

  # Add MRD = 0 reference line (no effect)
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "#D55E00", # Vermillion/orange-red
    linewidth = 1.0,
    alpha = 0.8
  ) +

  # Add subtle background shading for negative/positive regions
  annotate(
    "rect",
    xmin = -Inf, xmax = 0,
    ymin = -Inf, ymax = Inf,
    fill = "#E8F4F8", # Light blue tint
    alpha = 0.15
  ) +
  annotate(
    "rect",
    xmin = 0, xmax = Inf,
    ymin = -Inf, ymax = Inf,
    fill = "#FFF4E6", # Light warm tint
    alpha = 0.15
  ) +

  # Add "MRD = 0" label
  annotate(
    "text",
    x = 0,
    y = 4.6,
    label = "MRD = 0\n(No Effect)",
    color = "#D55E00",
    size = 3.5,
    fontface = "bold",
    hjust = -0.1,
    vjust = 0.5
  ) +

  # Add text annotations with mean and 95% CI
  geom_text(
    data = summary_stats,
    aes(x = Inf, y = quintile_label, label = label),
    hjust = 1.05,
    vjust = 1.8,
    size = 3.8,
    family = "mono",
    color = "#2C3E50", # Dark blue-gray
    fontface = "bold"
  ) +

  # Custom color scale
  scale_fill_manual(values = custom_colors) +

  # Enhanced x-axis
  scale_x_continuous(
    name = "Effect on AAMR (Age-Adjusted Mortality Rate)",
    breaks = pretty_breaks(n = 10),
    expand = expansion(mult = c(0.02, 0.02))
  ) +

  # Enhanced y-axis
  scale_y_discrete(
    name = "EQI Quintile",
    expand = expansion(add = c(0.2, 0.8))
  ) +

  # Enhanced theme
  theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
  theme(
    # Plot titles
    plot.title = element_text(
      size = 17,
      face = "bold",
      hjust = 0.5,
      color = "#1A1A1A",
      margin = margin(b = 8)
    ),
    plot.subtitle = element_text(
      size = 12.5,
      hjust = 0.5,
      color = "#4A4A4A",
      margin = margin(b = 15)
    ),
    plot.caption = element_text(
      size = 9.5,
      color = "#666666",
      hjust = 1,
      margin = margin(t = 10),
      lineheight = 1.2
    ),

    # Axes
    axis.title.x = element_text(
      size = 13,
      face = "bold",
      color = "#2C3E50",
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      size = 13,
      face = "bold",
      color = "#2C3E50",
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(size = 11, color = "#4A4A4A"),
    axis.text.y = element_text(size = 11, color = "#2C3E50", face = "bold"),

    # Grid
    panel.grid.major.x = element_line(color = "#D0D0D0", linewidth = 0.35),
    panel.grid.minor.x = element_line(color = "#E8E8E8", linewidth = 0.2),

    # Background
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),

    # Legend
    legend.position = "none",

    # Margins
    plot.margin = margin(20, 25, 15, 15)
  ) +

  # Title and caption
  labs(
    title = sprintf(
      "%s: Overall EQI Effect Distribution",
      ridge_data$metadata$cancer_type
    ),
    subtitle = sprintf(
      "Lag %d: %s EQI → %s AAMR | Bayesian Interval-Censored Mixed Model",
      ridge_data$metadata$lag,
      ridge_data$metadata$eqi_period,
      ridge_data$metadata$aamr_period
    ),
    caption = sprintf(
      "N = %d observations | %d MCMC draws per quintile\nWhite lines show median and 95%% credible intervals | Numbers show mean [95%% CI]\nDashed line (MRD = 0) indicates no mortality rate difference from reference quintile (Q1)",
      ridge_data$metadata$n_obs,
      ridge_data$metadata$n_draws
    )
  )

# ============================================================================
# Save plot
# ============================================================================
cat("Saving plot to:", output_file, "\n")

ggsave(
  output_file,
  plot = p,
  width = 14,
  height = 9,
  dpi = 350,
  bg = "white"
)

cat("✓ Plot saved successfully\n")
cat("\n")

# ============================================================================
# Additional summary output
# ============================================================================
cat("========================================\n")
cat("Posterior Summary\n")
cat("========================================\n\n")

for (i in 1:nrow(ridge_data$summary)) {
  row <- ridge_data$summary[i, ]
  cat(sprintf(
    "%-10s Mean: %7.2f  SD: %6.2f  95%% CI: [%7.2f, %7.2f]  R-hat: %.3f  ESS: %5.0f\n",
    row$quintile,
    row$mean,
    row$sd,
    row$q025,
    row$q975,
    row$rhat,
    row$ess_bulk
  ))
}

cat("\n")
cat("========================================\n")
cat("Interpretation Guide\n")
cat("========================================\n")
cat("- Distributions show full posterior uncertainty for each quintile effect\n")
cat("- Reference: Q1 (lowest EQI/best environmental quality) = 0\n")
cat("- Positive values indicate HIGHER mortality risk\n")
cat("- Negative values indicate LOWER mortality risk\n")
cat("- Wider distributions = more uncertainty\n")
cat("- Dashed orange line (MRD = 0) indicates no effect\n")
cat("- Light blue background (left) = protective effect region\n")
cat("- Light yellow background (right) = harmful effect region\n")
cat("========================================\n\n")

cat("✅ Visualization complete!\n")
cat("Output: ", output_file, "\n")
