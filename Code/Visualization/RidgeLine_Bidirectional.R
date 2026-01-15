#!/usr/bin/env Rscript


suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(ggridges)
  library(grid)
  library(ggtext)
  library(patchwork)
  library(cowplot)
})

# ==============================================================================
# CONFIGURATION & THEME
# ==============================================================================

# Colors from Visualization_Delta_Cluster.py
COLOR_IMPROVED <- "#7FB3D3"  # Ocean blue
COLOR_WORSENED <- "#F2AFAF"  # Cherry blossom pink

# Font settings
FONT_FAMILY <- "sans"
BASE_SIZE <- 15
TITLE_SIZE <- 17
AXIS_SIZE <- 15
LABEL_SIZE <- 13

# Publication Theme (Adapted from RidgeLine.R)
theme_publication <- function(base_size = BASE_SIZE, base_family = FONT_FAMILY) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # Text elements
      plot.title = element_text(size = TITLE_SIZE, face = "bold", hjust = 0, margin = margin(b = 10)),
      plot.subtitle = element_text(size = AXIS_SIZE, margin = margin(b = 10), color = "grey30"),
      plot.caption = element_text(size = LABEL_SIZE, color = "grey50", hjust = 1, margin = margin(t = 10)),

      # Axis elements
      axis.title = element_text(size = AXIS_SIZE, face = "bold"),
      axis.text = element_text(size = AXIS_SIZE, color = "black"),
      axis.text.y = element_text(margin = margin(r = 5)),

      # Grid elements
      panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", size = 0.2),
      panel.grid.minor.y = element_blank(),

      # Facet elements
      strip.text = element_text(size = AXIS_SIZE, face = "bold", hjust = 0),
      strip.background = element_rect(fill = "grey95", color = NA),
      panel.spacing = unit(3, "lines"),


      # Legend elements
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = AXIS_SIZE),
      legend.key.size = unit(0.8, "cm")
    )
}

# ==============================================================================
# DATA PROCESSING
# ==============================================================================

parse_key <- function(key) {
  # Expected formats:
  # Lag5_National_Overall
  # Lag5_National_Air
  # Lag5_K3_C0_Overall
  # Lag5_K3_C0_Air

  parts <- str_split(key, "_")[[1]]

  lag_str <- parts[1] # "Lag5"
  lag_val <- as.integer(gsub("Lag", "", lag_str))

  if (parts[2] == "National") {
    cluster_label <- "National"
    cluster_id <- -1
    domain <- parts[3]
  } else {
    # K3_C0
    k_val <- parts[2]
    c_val <- parts[3]
    cluster_id <- as.integer(gsub("C", "", c_val))
    cluster_label <- paste0("Cluster ", cluster_id)
    domain <- parts[4]
  }

  list(
    Lag = lag_val,
    Cluster_Label = cluster_label,
    Cluster_ID = cluster_id,
    Domain = domain
  )
}

load_and_process_data <- function(rds_path) {
  if (!file.exists(rds_path)) stop("File not found: ", rds_path)

  raw_list <- readRDS(rds_path)
  set.seed(42) # For reproducible downsampling

  all_draws <- list()

  for (key in names(raw_list)) {
    item <- raw_list[[key]]
    if (is.null(item) || is.null(item$draws)) next

    meta <- parse_key(key)

    # Extract draws
    d <- item$draws

    # Downsample for visualization performance (max 2000 draws)
    if (nrow(d) > 2000) {
      d <- d[sample(nrow(d), 2000), ]
    }

    d$Lag <- meta$Lag
    d$Cluster_Label <- meta$Cluster_Label
    d$Cluster_ID <- meta$Cluster_ID
    d$Domain <- meta$Domain

    all_draws[[key]] <- d
  }

  combined <- bind_rows(all_draws)
  rm(raw_list, all_draws)
  gc()

  # Shift distributions for visualization separation
  # Improved -> Left (-1), Worsened -> Right (+1)
  combined <- combined %>%
    mutate(
      effect = case_when(
        category == "Improved" ~ effect - 1,
        category == "Worsened" ~ effect + 1,
        TRUE ~ effect
      )
    )

  # Factor ordering

  # 1. Domain: Rename Overall -> EQI and set levels
  combined$Domain <- recode(combined$Domain, "Overall" = "EQI")
  domain_levels <- c("EQI", "Air", "Water", "Land", "Built", "Social")
  combined$Domain <- factor(combined$Domain, levels = rev(domain_levels)) # Rev for plotting bottom-to-top

  # 2. Cluster: National, then C0, C1...
  # Create a sorting index
  combined <- combined %>%
    mutate(
      Cluster_Order = ifelse(Cluster_Label == "National", -1, Cluster_ID),
      Cluster_Factor = factor(Cluster_Label, levels = unique(Cluster_Label[order(Cluster_Order)]))
    )

  # 3. Category: Improved, Deteriorated
  combined$category <- recode(combined$category, "Worsened" = "Deteriorated")
  combined$category <- factor(combined$category, levels = c("Improved", "Deteriorated"))

  return(combined)
}

# ==============================================================================
# PLOTTING
# ==============================================================================

plot_bidirectional_ridges <- function(data, output_path, cancer_code) {

  # Calculate x-axis limits (symmetric)
  max_val <- max(abs(quantile(data$effect, c(0.005, 0.995))), na.rm = TRUE)
  limit <- ceiling(max_val * 1.2)

  # Create Lag Labels and ensure factor ordering
  data$Lag_Label <- factor(paste0(data$Lag, "-Year Lag"), levels = c("5-Year Lag", "10-Year Lag"))

  # Split data by Lag
  d5 <- data %>% filter(Lag == 5)
  d10 <- data %>% filter(Lag == 10)

  # Common layers
  common_layers <- list(
    geom_vline(xintercept = 0, linetype = "solid", color = "grey40", size = 0.5),
    geom_density_ridges(
      alpha = 0.6,
      scale = 0.9,
      rel_min_height = 0.01,
      color = "white",
      size = 0.2
    ),
    scale_fill_manual(
      values = c("Improved" = COLOR_IMPROVED, "Deteriorated" = COLOR_WORSENED),
      labels = c("Improved Environment", "Deteriorated Environment")
    ),
    scale_x_continuous(
      limits = c(-limit, limit),
      breaks = scales::pretty_breaks(n = 5),
      expand = c(0, 0)
    ),
    theme_publication(),
    theme(
      panel.spacing.y = unit(1, "lines"),
      axis.text.y = element_text(vjust = 0),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.margin = margin(25, 25, 25, 25),
      legend.position = "none"
    )
  )

  # Plot Lag 5 (Left) - Axis on Right, Cluster Labels on Left
  p5 <- ggplot(d5, aes(x = effect, y = Domain, fill = category)) +
    common_layers +
    facet_grid(Cluster_Factor ~ Lag_Label, scales = "free_y", space = "free_y", switch = "y") +
    scale_y_discrete(position = "right") + # Put labels on right (middle of combined plot)
    labs(x = NULL, y = NULL) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 90, face = "bold", size = AXIS_SIZE, hjust = 0.5, vjust = 0.5),
      strip.text.x = element_text(size = TITLE_SIZE, face = "bold", hjust = 0.5),
      axis.text.y.right = element_text(hjust = 0.5, margin = margin(l = 10, r = 0)) # Center align labels
    )

  # Plot Lag 10 (Right) - No Axis, No Cluster Strip
  p10 <- ggplot(d10, aes(x = effect, y = Domain, fill = category)) +
    common_layers +
    facet_grid(Cluster_Factor ~ Lag_Label, scales = "free_y", space = "free_y") +
    labs(x = NULL, y = NULL) +
    theme(
      axis.text.y = element_blank(), # Hide axis text
      axis.ticks.y = element_blank(),
      strip.text.y = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = element_text(size = TITLE_SIZE, face = "bold", hjust = 0.5)
    )

  # Combine with Cowplot for precise alignment
  # align = "h" aligns rows, axis = "bt" aligns top and bottom axes

  p5 <- p5 + theme(plot.margin = margin(20, 5, 20, 25))   # 右边距从25减到15
  p10 <- p10 + theme(plot.margin = margin(20, 25, 20, 0)) # 左边距从15加到35

  p_combined <- plot_grid(p5, p10, align = "h", axis = "bt", ncol = 2, rel_widths = c(1.1, 0.9))

  # Add shared X axis label with legend info (using markdown for colors)
  x_label_text <- paste0(
    "MRD with Posterior Probability Distribution (",
    "<span style='color:", COLOR_IMPROVED, "'>Improved Environment</span>, ",
    "<span style='color:", COLOR_WORSENED, "'>Deteriorated Environment</span>)"
  )

  # Create a text grob for the label
  # We use ggdraw to place the label at the bottom
  label_grob <- gridtext::richtext_grob(
    x_label_text,
    x = 0.5, y = 1.2,
    gp = gpar(fontsize = AXIS_SIZE, fontface = "bold")
  )

  # Combine plot and label
  final_plot <- plot_grid(
    p_combined,
    label_grob,
    ncol = 1,
    rel_heights = c(1, 0.05) # Adjust height ratio for label
  )

  ggsave(output_path, plot = final_plot, width = 12, height = 12, dpi = 300, bg = "white")
  message("✅ Saved plot to: ", output_path)
}

# ==============================================================================
# MAIN
# ==============================================================================

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to input .rds file"),
  make_option(c("-o", "--output-dir"), type="character", default="Result/brms_delta_cluster_ridgeline",
              help="Output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input)) {
  stop("Please provide input file with --input")
}

# Extract cancer code from filename
filename <- basename(opt$input)
cancer_code <- gsub("_ridgeline.rds", "", filename)

# Setup output
if (!dir.exists(opt$`output-dir`)) dir.create(opt$`output-dir`, recursive = TRUE)
output_file <- file.path(opt$`output-dir`, paste0(cancer_code, "_Bidirectional_Ridge.png"))

message("📊 Processing ", cancer_code, "...")
message("   Input: ", opt$input)

# Load Data
dt <- load_and_process_data(opt$input)

# Plot
plot_bidirectional_ridges(dt, output_file, cancer_code)

message("Done.")
