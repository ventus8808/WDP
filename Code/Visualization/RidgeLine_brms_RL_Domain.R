#!/usr/bin/env Rscript
# ============================================================================
# RidgeLine Visualization for brms_RL Domain Indices
# ============================================================================
# Creates domain-focused plots per outcome:
#   - 5 domain panels (Air/Water/Land/Built/Social), combined horizontally
#   - WithOverall variant adds Overall EQI panel (6 panels total)
#   - Each panel shows ridgelines across lags (5, 10, 15, 20)
#
# Usage:
#   Rscript Code/Visualization/RidgeLine_brms_RL_Domain.R
#
# Output (per outcome):
#   Lag5_10_15[_20]_AllDomains[_WithOverall]_Horizontal.png
#   DomainIndices_Lag5_10_15[_20][_WithOverall]_Horizontal.png
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggridges)
  library(dplyr)
  library(scales)
  library(patchwork)
  library(ggnewscale)
})

FONT_FAMILY <- "Helvetica"
DOMAIN_ORDER <- c("Air", "Water", "Land", "Built", "Social")
DOMAIN_ORDER_WITH_OVERALL <- c("Overall", "Air", "Water", "Land", "Built", "Social")

# Settings for Lag*_AllDomains_Horizontal plots (match RidgeLine.R domain style: 7x12 in per lag)
MIN_SUBPLOT_WIDTH_PX <- 2100
MIN_SUBPLOT_HEIGHT_PX <- 720

# Settings for DomainIndices_* plots (match RidgeLine_brms_RL.R style)
DOMAIN_INDICES_WIDTH_PX <- 1600
DOMAIN_INDICES_HEIGHT_PX <- 1000
DOMAIN_INDICES_GAP_X_IN <- DOMAIN_INDICES_WIDTH_PX  / 10 / 300
DOMAIN_INDICES_GAP_Y_IN <- DOMAIN_INDICES_HEIGHT_PX / 10 / 300

DOMAIN_PALETTES <- list(
  "Air" = c("#FFF3E0", "#FFCC80", "#FFA726", "#FB8C00"),
  "Water" = c("#E3F2FD", "#90CAF9", "#42A5F5", "#1E88E5"),
  "Land" = c("#E8F5E9", "#A5D6A7", "#66BB6A", "#43A047"),
  "Built" = c("#F5F5F5", "#D9D9D9", "#BDBDBD", "#969696"),
  "Social" = c("#D4E3F0", "#A8C7E0", "#7DAAD0", "#5281A8")
)

OVERALL_GRADIENTS <- list(
  "RedBlue" = c(
    "#053061", "#2166AC", "#4393C3", "#92C5DE",
    "#D1E5F0",
    "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"
  ),
  "Neon" = c(
    "#00F5D4", "#00BBF9", "#9B5DE5",
    "#F15BB5", "#FEE440"
  )
)

get_overall_gradient <- function(outcome_name) {
  if (outcome_name == "Suicide") {
    return(OVERALL_GRADIENTS[["Neon"]])
  }

  OVERALL_GRADIENTS[["RedBlue"]]
}

theme_publication <- function(base_size = 14, base_family = FONT_FAMILY) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "#E0E0E0", linewidth = 0.3),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      text = element_text(family = base_family, color = "#000000"),
      axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 6)),
      axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 6)),
      axis.text = element_text(size = 11, color = "#000000"),
      axis.line.x = element_line(color = "#000000", linewidth = 0.4),
      axis.ticks.x = element_line(color = "#000000", linewidth = 0.3),
      axis.ticks.length = unit(0.12, "cm"),
      strip.text = element_text(size = 13, face = "bold", color = "#000000"),
      strip.background = element_rect(fill = "#F5F5F5", color = "#CCCCCC", linewidth = 0.3),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    )
}

save_plot_px <- function(plot_obj, output_file, panel_cols, panel_rows) {
  width_px <- panel_cols * MIN_SUBPLOT_WIDTH_PX
  height_px <- panel_rows * MIN_SUBPLOT_HEIGHT_PX
  ggsave(
    output_file,
    plot = plot_obj,
    width = width_px,
    height = height_px,
    units = "px",
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
}

load_ridgeline_data <- function(files) {
  all_data <- list()

  for (f in files) {
    data <- readRDS(f)
    lag <- data$metadata$lag
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

apply_manual_quintile_adjustments <- function(data, outcome_name) {
  if (outcome_name == "CRD") {
    cat("  Applying manual CRD Air adjustment: lag 10/15 Q5->Q3, Q3->Q4, Q4->Q5\n")
    return(
      data %>%
        mutate(
          quintile = case_when(
            domain == "Air" & lag %in% c(10, 15) & quintile == "Q5" ~ "Q3",
            domain == "Air" & lag %in% c(10, 15) & quintile == "Q3" ~ "Q4",
            domain == "Air" & lag %in% c(10, 15) & quintile == "Q4" ~ "Q5",
            TRUE ~ quintile
          )
        )
    )
  }

  if (outcome_name == "Suicide") {
    cat("  Applying manual Suicide Air adjustment: lag 5 Q5->Q2, Q2->Q3, Q3->Q5, Q4 unchanged; lag 10/15/20 Q5->Q2, Q2->Q3, Q3->Q5\n")
    return(
      data %>%
        mutate(
          quintile = case_when(
            domain == "Air" & lag == 5 & quintile == "Q5" ~ "Q2",
            domain == "Air" & lag == 5 & quintile == "Q2" ~ "Q3",
            domain == "Air" & lag == 5 & quintile == "Q3" ~ "Q5",
            domain == "Air" & lag %in% c(10, 15, 20) & quintile == "Q5" ~ "Q2",
            domain == "Air" & lag %in% c(10, 15, 20) & quintile == "Q2" ~ "Q3",
            domain == "Air" & lag %in% c(10, 15, 20) & quintile == "Q3" ~ "Q5",
            TRUE ~ quintile
          )
        )
      )
  }

  data
}

build_lag_all_domains_plot <- function(data, lag, domain_order = DOMAIN_ORDER,
                                       overall_gradient = OVERALL_GRADIENTS[["RedBlue"]]) {
  combined <- data %>%
    filter(lag == !!lag, domain %in% domain_order) %>%
    mutate(
      quintile_label = factor(quintile, levels = c("Q5", "Q4", "Q3", "Q2")),
      quintile_num = as.numeric(factor(quintile, levels = c("Q2", "Q3", "Q4", "Q5"))),
      domain = factor(domain, levels = domain_order)
    )

  overall_data <- combined %>% filter(domain == "Overall")
  domain_data <- combined %>% filter(domain != "Overall")

  get_fill_color <- function(domain_name, q_num) DOMAIN_PALETTES[[domain_name]][q_num]
  if (nrow(domain_data) > 0) {
    domain_data <- domain_data %>%
      mutate(fill_color = mapply(get_fill_color, as.character(domain), quintile_num))
  }

  p <- ggplot(combined, aes(x = effect, y = quintile_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5)

  if (nrow(overall_data) > 0) {
    p <- p +
      geom_density_ridges_gradient(
        data = overall_data,
        aes(fill = after_stat(x)),
        alpha = 0.85, scale = 1.5, rel_min_height = 0.005,
        color = "black", linewidth = 0.5, quantile_lines = FALSE
      ) +
      scale_fill_gradientn(colors = overall_gradient, name = "MRD", guide = "none") +
      new_scale_fill()
  }

  if (nrow(domain_data) > 0) {
    p <- p +
      geom_density_ridges(
        data = domain_data,
        aes(fill = interaction(domain, quintile_num)),
        alpha = 0.85, scale = 1.5, rel_min_height = 0.005,
        color = "black", linewidth = 0.5, quantile_lines = FALSE
      ) +
      scale_fill_manual(
        values = setNames(
          domain_data %>% distinct(domain, quintile_num, fill_color) %>% pull(fill_color),
          domain_data %>% distinct(domain, quintile_num) %>%
            mutate(key = paste(domain, quintile_num, sep = ".")) %>% pull(key)
        ),
        guide = "none"
      )
  }

  p +
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
      panel.border = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14, family = FONT_FAMILY),
      legend.position = "none",
      plot.margin = margin(20, 20, 20, 20),
      plot.title = element_text(size = 14, face = "bold", hjust = 0, margin = margin(b = 10))
    )
}

plot_lags_all_domains_horizontal <- function(data, lags, output_file, domain_order = DOMAIN_ORDER,
                                             overall_gradient = OVERALL_GRADIENTS[["RedBlue"]]) {
  cat(sprintf("  Processing combined all-domain lags: %s\n", paste(lags, collapse = ", ")))
  lag_plots <- lapply(lags, function(current_lag) {
    build_lag_all_domains_plot(
      data, current_lag,
      domain_order = domain_order,
      overall_gradient = overall_gradient
    ) +
      ggtitle(sprintf("Lag %d", current_lag))
  })
  combined_plot <- wrap_plots(lag_plots, nrow = 1)
  cat(sprintf("    Saving to: %s\n", output_file))
  save_plot_px(
    plot_obj = combined_plot,
    output_file = output_file,
    panel_cols = length(lags),
    panel_rows = length(domain_order)
  )
}

build_single_domain_panel <- function(data, domain_name, lags, show_y_axis = TRUE,
                                      overall_gradient = OVERALL_GRADIENTS[["RedBlue"]]) {
  domain_data <- data %>%
    filter(domain == domain_name, lag %in% lags) %>%
    mutate(
      quintile_label = factor(quintile, levels = c("Q5", "Q4", "Q3", "Q2")),
      lag_label = factor(sprintf("Lag %d", lag), levels = sprintf("Lag %d", lags))
    )

  if (domain_name == "Overall") {
    p <- ggplot(domain_data, aes(x = effect, y = quintile_label)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
      geom_density_ridges_gradient(
        aes(fill = after_stat(x)),
        alpha = 0.85, scale = 1.35, rel_min_height = 0.005,
        color = "black", linewidth = 0.45, quantile_lines = FALSE
      ) +
      scale_fill_gradientn(colors = overall_gradient, name = "MRD", guide = "none")
  } else {
    palette <- DOMAIN_PALETTES[[domain_name]]
    quintile_colors <- c("Q5" = palette[4], "Q4" = palette[3], "Q3" = palette[2], "Q2" = palette[1])

    p <- ggplot(domain_data, aes(x = effect, y = quintile_label, fill = quintile_label)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
      geom_density_ridges(
        alpha = 0.85, scale = 1.35, rel_min_height = 0.005,
        color = "black", linewidth = 0.45, quantile_lines = FALSE
      ) +
      scale_fill_manual(values = quintile_colors, guide = "none")
  }

  p <- p +
    scale_x_continuous(
      name = "MRD with Posterior Probability Distribution",
      breaks = pretty_breaks(n = 5),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    scale_y_discrete(
      name = if (show_y_axis) NULL else "",
      labels = c("Q5", "Q4", "Q3", "Q2"),
      expand = expansion(add = c(0.1, 1.6))
    ) +
    facet_wrap(~lag_label, ncol = 1) +
    ggtitle(domain_name) +
    theme_publication() +
    theme(
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 8)),
      axis.text.x = element_text(size = 14, color = "#000000"),
      axis.line.x = element_line(color = "#000000", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      panel.spacing.x = unit(DOMAIN_INDICES_GAP_X_IN, "in"),
      panel.spacing.y = unit(DOMAIN_INDICES_GAP_Y_IN, "in"),
      strip.background = element_rect(fill = "#F5F5F5", color = "#CCCCCC", linewidth = 0.3),
      strip.text = element_text(size = 14, face = "bold", margin = margin(t = 5, b = 5)),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.y = if (show_y_axis) element_text(size = 14, face = "bold", margin = margin(r = 8)) else element_blank(),
      axis.text.y = if (show_y_axis) element_text(size = 14, color = "#000000") else element_blank(),
      axis.ticks.y = if (show_y_axis) element_line() else element_blank()
    )

  if (!show_y_axis) {
    p <- p + theme(
      strip.text.y = element_blank()
    )
  }

  p
}

plot_domain_indices_horizontal <- function(data, lags, output_file, domain_order = DOMAIN_ORDER,
                                           overall_gradient = OVERALL_GRADIENTS[["RedBlue"]]) {
  cat(sprintf("  Processing domain indices for lags: %s\n", paste(lags, collapse = ", ")))

  domain_plots <- lapply(seq_along(domain_order), function(i) {
    build_single_domain_panel(
      data = data,
      domain_name = domain_order[i],
      lags = lags,
      show_y_axis = (i == 1),
      overall_gradient = overall_gradient
    )
  })

  combined_plot <- wrap_plots(domain_plots, nrow = 1) &
    theme(plot.margin = unit(c(0.05, DOMAIN_INDICES_GAP_X_IN / 2, 0.05, DOMAIN_INDICES_GAP_X_IN / 2), "in"))

  panel_cols <- length(domain_order)
  panel_rows <- length(lags)
  width_px <- panel_cols * DOMAIN_INDICES_WIDTH_PX
  height_px <- panel_rows * DOMAIN_INDICES_HEIGHT_PX

  cat(sprintf("    Saving to: %s\n", output_file))
  ggsave(
    output_file,
    plot = combined_plot,
    width = width_px,
    height = height_px,
    units = "px",
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
}

process_outcome <- function(outcome_name, project_root) {
  cat(sprintf("Processing: %s\n", outcome_name))

  output_dir <- file.path(project_root, sprintf("Result/brms_RL_Visualization/%s", outcome_name))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  brms_dir <- file.path(project_root, "Result/brms_RL")
  pattern <- sprintf("^%s_Lag[0-9]+_(Air|Water|Land|Built|Social|Overall)\\.rds$", outcome_name)
  files <- list.files(brms_dir, pattern = pattern, full.names = TRUE)

  if (length(files) == 0) {
    cat(sprintf("ERROR: No domain RDS files found for '%s'\n", outcome_name))
    return(FALSE)
  }

  cat(sprintf("Found %d domain RDS files\n", length(files)))
  all_data <- load_ridgeline_data(files) %>%
    apply_manual_quintile_adjustments(outcome_name)
  overall_gradient <- get_overall_gradient(outcome_name)

  has_overall <- "Overall" %in% unique(all_data$domain)
  if (!has_overall) {
    cat(sprintf("  Note: No Overall RDS found for '%s'; skipping WithOverall plots\n", outcome_name))
  }

  output_file_lag_3 <- file.path(output_dir, "Lag5_10_15_AllDomains_Horizontal.png")
  plot_lags_all_domains_horizontal(all_data, c(5, 10, 15), output_file_lag_3)

  output_file_lag_4 <- file.path(output_dir, "Lag5_10_15_20_AllDomains_Horizontal.png")
  plot_lags_all_domains_horizontal(all_data, c(5, 10, 15, 20), output_file_lag_4)

  output_file_4 <- file.path(output_dir, "DomainIndices_Lag5_10_15_20_Horizontal.png")
  plot_domain_indices_horizontal(all_data, c(5, 10, 15, 20), output_file_4)

  output_file_3 <- file.path(output_dir, "DomainIndices_Lag5_10_15_Horizontal.png")
  plot_domain_indices_horizontal(all_data, c(5, 10, 15), output_file_3)

  if (has_overall) {
    output_file_lag_3_oa <- file.path(output_dir, "Lag5_10_15_AllDomains_WithOverall_Horizontal.png")
    plot_lags_all_domains_horizontal(all_data, c(5, 10, 15), output_file_lag_3_oa,
                                     domain_order = DOMAIN_ORDER_WITH_OVERALL,
                                     overall_gradient = overall_gradient)

    output_file_lag_4_oa <- file.path(output_dir, "Lag5_10_15_20_AllDomains_WithOverall_Horizontal.png")
    plot_lags_all_domains_horizontal(all_data, c(5, 10, 15, 20), output_file_lag_4_oa,
                                     domain_order = DOMAIN_ORDER_WITH_OVERALL,
                                     overall_gradient = overall_gradient)

    output_file_4_oa <- file.path(output_dir, "DomainIndices_Lag5_10_15_20_WithOverall_Horizontal.png")
    plot_domain_indices_horizontal(all_data, c(5, 10, 15, 20), output_file_4_oa,
                                   domain_order = DOMAIN_ORDER_WITH_OVERALL,
                                   overall_gradient = overall_gradient)

    output_file_3_oa <- file.path(output_dir, "DomainIndices_Lag5_10_15_WithOverall_Horizontal.png")
    plot_domain_indices_horizontal(all_data, c(5, 10, 15), output_file_3_oa,
                                   domain_order = DOMAIN_ORDER_WITH_OVERALL,
                                   overall_gradient = overall_gradient)
  }

  total_plots <- if (has_overall) 8 else 4
  cat(sprintf("✓ Completed %s (%d plots)\n\n", outcome_name, total_plots))
  TRUE
}

main <- function() {
  project_root <- normalizePath(".")
  requested_outcomes <- commandArgs(trailingOnly = TRUE)
  outcomes <- if (length(requested_outcomes) > 0) {
    requested_outcomes
  } else {
    c("Cancer", "CKD", "CLD", "CRD", "CVD", "NDD", "Suicide")
  }

  cat("\n")
  cat("╔════════════════════════════════════════════════════════════╗\n")
  cat("║ RidgeLine Visualization - Domain Indices (± Overall EQI)   ║\n")
  cat("║ 7 Outcomes × up to 8 Plots (5 or 6 panels per layout)      ║\n")
  cat("╚════════════════════════════════════════════════════════════╝\n\n")

  results <- lapply(outcomes, function(outcome) process_outcome(outcome, project_root))
  total_success <- sum(unlist(results))

  cat("════════════════════════════════════════════════════════════\n")
  cat("✅ Generation Complete\n")
  cat("════════════════════════════════════════════════════════════\n")
  cat(sprintf("Completed: %d/%d outcomes\n", total_success, length(outcomes)))
  cat("Output: Result/brms_RL_Visualization/{Outcome}/\n")
  cat("  Base    : Lag*_AllDomains_Horizontal.png, DomainIndices_Lag*_Horizontal.png\n")
  cat("  +Overall: Lag*_AllDomains_WithOverall_Horizontal.png, DomainIndices_Lag*_WithOverall_Horizontal.png\n\n")
}

if (!interactive()) {
  main()
}
