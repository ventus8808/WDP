#!/usr/bin/env Rscript
# Ridgeline_Delta.R — Bidirectional ridgeline plots for Delta model results
#
# Input:  Result/Delta_RL/{Disease}_Lag{N}_{suffix}.rds
# Output: Result/Delta_Ridgeline/{Disease}_{strat_type}[_Simple].png
#
# Per disease, full PNG files are produced for each strat type:
#   Full    — 6 domains on Y, 3 lag columns, all strat rows
# Additionally, one combined simplified PNG is produced:
#   Combined_Simple — EQI domain only, combining RUCC + Cluster_EQI + Cluster_NLCD
#
# Strat level labels come from mapping.yaml:
#   Cluster_EQI : A=Low-burden, B=Mixed-burden, C=High-burden
#   Cluster_NLCD: A=Natural, B=Water-sensitive, C=Agricultural, D=Urban
#   RUCC        : 1=Metropolitan urban, 2=Nonmetropolitan urban,
#                 3=Less urban, 4=Thinly populated
#
# RDS structure (from Delta_Main.R / Delta_Cluster.R / Delta_RUCC.R):
#   $draws_long — data.frame: draw_id | category ("I"/"W") | effect

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggridges)
  library(cowplot)
  library(scales)
})

# ── Constants ─────────────────────────────────────────────────────────────────
COLOR_IMPROVED <- "#7494BC"
COLOR_WORSENED <- "#DE9F83"
BASE_SIZE <- 13
TITLE_SIZE <- 14
STRIP_SIZE <- 9

DOMAIN_LEVELS <- c("EQI", "Air", "Water", "Land", "Built", "Social")
DOMAIN_PAT <- "EQI|Air|Water|Land|Built|Social"

DISEASES_ALL <- c("CVD", "CRD", "CLD", "CKD", "Suicide", "NDD", "Cancer")
STRAT_TYPES <- c("Main", "Cluster_EQI", "Cluster_NLCD", "RUCC")
COMBINED_SIMPLE_TYPES <- c("RUCC", "Cluster_EQI", "Cluster_NLCD")

# Label maps derived from mapping.yaml
LABEL_EQI <- c(A = "Low-burden", B = "Mixed-burden", C = "High-burden")
LABEL_NLCD <- c(
  A = "Natural", B = "Water",
  C = "Agricultural", D = "Urban"
)
LABEL_RUCC <- c(
  `1` = "Metropolitan", `2` = "Nonmetropolitan",
  `3` = "Semi-rural", `4` = "Rural"
)

# ── Theme ─────────────────────────────────────────────────────────────────────
theme_pub <- function() {
  theme_minimal(base_size = BASE_SIZE) +
    theme(
      plot.title = element_text(
        size = TITLE_SIZE, face = "bold",
        hjust = 0.5, margin = margin(b = 6)
      ),
      axis.title = element_blank(),
      axis.text = element_text(size = BASE_SIZE, color = "black"),
      axis.text.y = element_text(margin = margin(r = 4)),
      axis.line.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.2),
      panel.grid.minor.y = element_blank(),
      strip.text = element_text(size = BASE_SIZE, face = "bold"),
      strip.background = element_rect(fill = "grey93", color = NA),
      panel.spacing = unit(0.8, "lines"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = BASE_SIZE),
      legend.key.size = unit(0.7, "cm")
    )
}

# ── CLI ───────────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--input-dir",
    type = "character", default = "Result/Delta_RL",
    help = "Directory containing .rds files"
  ),
  make_option("--output-dir",
    type = "character", default = "Result/Delta_Ridgeline_Visualization",
    help = "Output directory for PNG files"
  ),
  make_option("--diseases",
    type = "character", default = NA,
    help = "Comma-separated disease labels (default: all 7)"
  ),
  make_option("--stratify-by",
    type = "character", default = "All",
    help = "All | Main | Cluster_EQI | Cluster_NLCD | RUCC"
  ),
  make_option("--lags",
    type = "character", default = "5,10,15",
    help = "Comma-separated lag values"
  ),
  make_option("--max-draws",
    type = "integer", default = 2000,
    help = "Max draws per category per file (for performance)"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

# ── Filename parser ───────────────────────────────────────────────────────────
# Returns list(lag, strat_level, domain) or NULL if file doesn't match.
parse_filename <- function(fn, disease, strat_type) {
  dp <- DOMAIN_PAT
  d <- disease

  if (strat_type == "Main") {
    m <- regmatches(fn, regexec(
      sprintf("^%s_Lag(\\d+)_(%s)_Main\\.rds$", d, dp), fn
    ))[[1]]
    if (length(m) == 3) {
      return(list(lag = as.integer(m[2]), strat_level = "National", domain = m[3]))
    }
  }

  if (strat_type %in% c("Cluster_EQI", "Cluster_NLCD")) {
    # National row: Cluster_{Domain}.rds (Phase 1 of Delta_Cluster.R)
    m <- regmatches(fn, regexec(
      sprintf("^%s_Lag(\\d+)_Cluster_(%s)\\.rds$", d, dp), fn
    ))[[1]]
    if (length(m) == 3) {
      return(list(lag = as.integer(m[2]), strat_level = "National", domain = m[3]))
    }

    if (strat_type == "Cluster_EQI") {
      m <- regmatches(fn, regexec(
        sprintf("^%s_Lag(\\d+)_Cluster_EQI_([ABC])_(%s)\\.rds$", d, dp), fn
      ))[[1]]
      if (length(m) == 4) {
        return(list(
          lag = as.integer(m[2]),
          strat_level = unname(LABEL_EQI[m[3]]), domain = m[4]
        ))
      }
    }

    if (strat_type == "Cluster_NLCD") {
      m <- regmatches(fn, regexec(
        sprintf("^%s_Lag(\\d+)_Cluster_NLCD_([ABCD])_(%s)\\.rds$", d, dp), fn
      ))[[1]]
      if (length(m) == 4) {
        return(list(
          lag = as.integer(m[2]),
          strat_level = unname(LABEL_NLCD[m[3]]), domain = m[4]
        ))
      }
    }
  }

  if (strat_type == "RUCC") {
    # National row: RUCC_{Domain}.rds (Phase 1 of Delta_RUCC.R)
    m <- regmatches(fn, regexec(
      sprintf("^%s_Lag(\\d+)_RUCC_(%s)\\.rds$", d, dp), fn
    ))[[1]]
    if (length(m) == 3) {
      return(list(lag = as.integer(m[2]), strat_level = "National", domain = m[3]))
    }

    m <- regmatches(fn, regexec(
      sprintf("^%s_Lag(\\d+)_RUCC_([1-4])_(%s)\\.rds$", d, dp), fn
    ))[[1]]
    if (length(m) == 4) {
      return(list(
        lag = as.integer(m[2]),
        strat_level = unname(LABEL_RUCC[m[3]]), domain = m[4]
      ))
    }
  }

  NULL
}

# ── Data loader ───────────────────────────────────────────────────────────────
load_strat_data <- function(input_dir, disease, strat_type, lags, domains, max_draws) {
  files <- list.files(input_dir,
    pattern = sprintf("^%s_", disease),
    full.names = TRUE
  )
  rows <- list()

  for (fp in files) {
    fn <- basename(fp)
    meta <- parse_filename(fn, disease, strat_type)
    if (is.null(meta)) next
    if (!meta$lag %in% lags) next
    if (!meta$domain %in% domains) next

    rds <- tryCatch(readRDS(fp), error = function(e) NULL)
    if (is.null(rds) || is.null(rds$draws_long)) next

    dl <- rds$draws_long

    # Downsample each category independently (density doesn't need draw pairing)
    for (cat in c("I", "W")) {
      idx <- which(dl$category == cat)
      if (length(idx) > max_draws) dl <- dl[-sample(idx, length(idx) - max_draws), ]
    }

    dl$Lag <- meta$lag
    dl$Strat_Level <- meta$strat_level
    dl$Domain <- meta$domain
    rows[[fn]] <- dl
  }

  if (length(rows) == 0) {
    return(NULL)
  }
  bind_rows(rows)
}

# ── Factor ordering & distribution shift ─────────────────────────────────────
prepare_data <- function(df, strat_order, domain_levels) {
  df$category <- case_when(
    df$category == "I" ~ "Improved",
    df$category == "W" ~ "Deteriorated",
    TRUE ~ df$category
  )
  df$category <- factor(df$category, levels = c("Improved", "Deteriorated"))

  # Shift distributions apart for visual separation
  df <- mutate(df, effect = case_when(
    category == "Improved" ~ effect - 1,
    category == "Deteriorated" ~ effect + 1,
    TRUE ~ effect
  ))

  df$Domain <- factor(df$Domain, levels = rev(domain_levels))
  df$Strat_Level <- factor(df$Strat_Level, levels = strat_order)
  df
}

# ── Per-lag ggplot panel ──────────────────────────────────────────────────────
make_panel <- function(df_lag, lag_val, has_strat,
                       show_y_labels, show_strat_strip, xlim_val) {
  p <- ggplot(df_lag, aes(x = effect, y = Domain, fill = category)) +
    geom_vline(
      xintercept = 0, linetype = "solid",
      color = "grey35", linewidth = 0.5
    ) +
    geom_density_ridges(
      alpha = 0.65, scale = 0.9, rel_min_height = 0.01,
      color = "white", linewidth = 0.2
    ) +
    scale_fill_manual(
      values = c(
        "Improved" = COLOR_IMPROVED,
        "Deteriorated" = COLOR_WORSENED
      ),
      labels = c("Improved Environment", "Deteriorated Environment")
    ) +
    scale_x_continuous(
      limits = xlim_val,
      breaks = pretty(xlim_val, n = 7),
      labels = label_number(accuracy = 1),
      expand = c(0.02, 0)
    ) +
    scale_y_discrete(
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(title = paste0(lag_val, "-Year Lag"), x = NULL, y = NULL) +
    theme_pub() +
    theme(plot.margin = margin(12, 8, 12, 8), legend.position = "none")

  if (has_strat) {
    p <- p + facet_grid(
      Strat_Level ~ .,
      scales = "free_y", space = "free_y",
      switch = if (show_strat_strip) "y" else NULL
    )
    if (show_strat_strip) {
      p <- p + theme(
        strip.placement = "outside",
        strip.text.y.left = element_text(
          angle = 90, face = "bold", size = STRIP_SIZE, hjust = 0.5
        )
      )
    } else {
      p <- p + theme(
        strip.text.y       = element_blank(),
        strip.background.y = element_blank()
      )
    }
  }

  if (!show_y_labels) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }

  p
}

# ── Strat level ordering ──────────────────────────────────────────────────────
strat_order_for <- function(strat_type) {
  switch(strat_type,
    Main         = "National",
    Cluster_EQI  = c("National", unname(LABEL_EQI)),
    Cluster_NLCD = c("National", unname(LABEL_NLCD)),
    RUCC         = c("National", unname(LABEL_RUCC)),
    "National"
  )
}

auto_height <- function(n_strat, n_domains) {
  max(3, n_strat * n_domains * 0.375 + 1.5)
}

# ── Assemble & save ───────────────────────────────────────────────────────────
build_and_save <- function(df, strat_type, simplified,
                           lags_used, output_path, width, height,
                           show_all_y_labels = FALSE) {
  domain_levels <- if (simplified) "EQI" else DOMAIN_LEVELS
  strat_order <- strat_order_for(strat_type)
  has_strat <- strat_type != "Main"

  df <- prepare_data(df, strat_order, domain_levels)

  q <- quantile(df$effect, c(0.005, 0.995), na.rm = TRUE)
  lim <- ceiling(max(abs(q)) / 5) * 5
  xlim <- c(-lim, lim)

  avail_lags <- sort(intersect(lags_used, unique(df$Lag)))
  if (length(avail_lags) == 0) {
    message("  [Skip] no lags found in data")
    return(invisible(NULL))
  }

  panels <- lapply(seq_along(avail_lags), function(i) {
    lag_v <- avail_lags[i]
    sub <- df[df$Lag == lag_v, ]
    if (nrow(sub) == 0) {
      return(NULL)
    }
    make_panel(
      df_lag = sub,
      lag_val = lag_v,
      has_strat = has_strat,
      show_y_labels = show_all_y_labels || (i == 1),
      show_strat_strip = (i == 1) && has_strat,
      xlim_val = xlim
    )
  })
  panels <- Filter(Negate(is.null), panels)
  if (length(panels) == 0) {
    return(invisible(NULL))
  }

  # Set outer margins
  n_p <- length(panels)
  panels[[1]] <- panels[[1]] + theme(plot.margin = margin(12, 5, 12, 15))
  if (n_p >= 2) panels[[n_p]] <- panels[[n_p]] + theme(plot.margin = margin(12, 15, 12, 5))
  if (n_p == 3) panels[[2]] <- panels[[2]] + theme(plot.margin = margin(12, 5, 12, 5))

  rel_w <- switch(as.character(n_p),
    "3" = c(1.2, 1, 1),
    "2" = c(1.2, 1),
    1
  )

  top_row <- plot_grid(
    plotlist   = panels,
    align      = "h",
    axis       = "bt",
    ncol       = n_p,
    rel_widths = rel_w
  )

  # Shared legend — use a minimal rect plot to avoid ggridges grob issues
  legend <- get_legend(
    ggplot(
      data.frame(
        x = c(0, 1),
        y = c(0, 1),
        cat = factor(c("Improved", "Deteriorated"),
          levels = c("Improved", "Deteriorated")
        )
      ),
      aes(xmin = x, xmax = x + 0.4, ymin = y, ymax = y + 0.4, fill = cat)
    ) +
      geom_rect() +
      scale_fill_manual(
        values = c(
          "Improved" = COLOR_IMPROVED,
          "Deteriorated" = COLOR_WORSENED
        ),
        labels = c("Improved Environment", "Deteriorated Environment")
      ) +
      theme_pub() +
      theme(legend.position = "bottom")
  )

  # Shared x-axis label
  x_label <- draw_label(
    "MRD with Posterior Probability Distribution",
    fontface = "bold", size = BASE_SIZE
  )

  final <- plot_grid(
    top_row, legend, x_label,
    ncol = 1,
    rel_heights = c(1, 0.04, 0.04)
  )

  ggsave(output_path,
    plot = final,
    width = width, height = height, dpi = 300, bg = "white"
  )
  message("  ✓ ", basename(output_path))
}

load_combined_simple_data <- function(input_dir, disease, lags, max_draws, strat_types) {
  all_specs <- list(
    Cluster_NLCD = list(strat_type = "Cluster_NLCD", group = "Cluster NLCD"),
    Cluster_EQI  = list(strat_type = "Cluster_EQI", group = "Cluster EQI"),
    RUCC         = list(strat_type = "RUCC", group = "RUCC")
  )
  specs <- all_specs[intersect(names(all_specs), strat_types)]

  rows <- list()
  for (spec in specs) {
    df <- load_strat_data(
      input_dir   = input_dir,
      disease     = disease,
      strat_type  = spec$strat_type,
      lags        = lags,
      domains     = "EQI",
      max_draws   = max_draws
    )
    if (!is.null(df) && nrow(df) > 0) {
      df <- df[df$Strat_Level != "National", ]
      if (nrow(df) == 0) next

      df$Strat_Label <- df$Strat_Level

      df$Strat_Group <- spec$group
      rows[[spec$strat_type]] <- df
    }
  }

  if (length(rows) == 0) {
    return(NULL)
  }
  bind_rows(rows)
}

load_main_combined_data <- function(input_dir, disease, lags, max_draws) {
  rows <- list()

  df_main <- load_strat_data(
    input_dir = input_dir, disease = disease, strat_type = "Main",
    lags = lags, domains = DOMAIN_LEVELS, max_draws = max_draws
  )
  if (!is.null(df_main) && nrow(df_main) > 0) {
    df_main <- df_main[df_main$Strat_Level == "National", ]
    if (nrow(df_main) > 0) {
      df_main$Strat_Group <- "National"
      df_main$Strat_Label <- as.character(df_main$Domain)
      rows[["Main"]] <- df_main
    }
  }

  df_eqi <- load_strat_data(
    input_dir = input_dir, disease = disease, strat_type = "Cluster_EQI",
    lags = lags, domains = "EQI", max_draws = max_draws
  )
  if (!is.null(df_eqi) && nrow(df_eqi) > 0) {
    df_eqi <- df_eqi[df_eqi$Strat_Level != "National", ]
    if (nrow(df_eqi) > 0) {
      df_eqi$Strat_Group <- "Regime EQI"
      df_eqi$Strat_Label <- df_eqi$Strat_Level
      rows[["Cluster_EQI"]] <- df_eqi
    }
  }

  df_nlcd <- load_strat_data(
    input_dir = input_dir, disease = disease, strat_type = "Cluster_NLCD",
    lags = lags, domains = "EQI", max_draws = max_draws
  )
  if (!is.null(df_nlcd) && nrow(df_nlcd) > 0) {
    df_nlcd <- df_nlcd[df_nlcd$Strat_Level != "National", ]
    if (nrow(df_nlcd) > 0) {
      df_nlcd$Strat_Group <- "Land Use EQI"
      df_nlcd$Strat_Label <- df_nlcd$Strat_Level
      rows[["Cluster_NLCD"]] <- df_nlcd
    }
  }

  if (length(rows) == 0) {
    return(NULL)
  }
  bind_rows(rows)
}

prepare_combined_simple_data <- function(df, strat_group_levels, strat_label_levels) {
  df$category <- case_when(
    df$category == "I" ~ "Improved",
    df$category == "W" ~ "Deteriorated",
    TRUE ~ df$category
  )
  df$category <- factor(df$category, levels = c("Improved", "Deteriorated"))

  df <- mutate(df, effect = case_when(
    category == "Improved" ~ effect - 1,
    category == "Deteriorated" ~ effect + 1,
    TRUE ~ effect
  ))

  df$Strat_Group <- factor(df$Strat_Group, levels = strat_group_levels)
  df$Strat_Label <- factor(df$Strat_Label, levels = strat_label_levels)
  df
}

make_panel_combined_simple <- function(df_lag, lag_val, show_y_labels, show_group_strip, xlim_val) {
  p <- ggplot(df_lag, aes(x = effect, y = Strat_Label, fill = category)) +
    geom_vline(
      xintercept = 0, linetype = "solid",
      color = "grey35", linewidth = 0.5
    ) +
    geom_density_ridges(
      alpha = 0.65, scale = 0.9, rel_min_height = 0.01,
      color = "white", linewidth = 0.2
    ) +
    scale_fill_manual(
      values = c(
        "Improved" = COLOR_IMPROVED,
        "Deteriorated" = COLOR_WORSENED
      ),
      labels = c("Improved Environment", "Deteriorated Environment")
    ) +
    scale_x_continuous(
      limits = xlim_val,
      breaks = pretty(xlim_val, n = 7),
      labels = label_number(accuracy = 1),
      expand = c(0.02, 0)
    ) +
    scale_y_discrete(
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(title = paste0(lag_val, "-Year Lag"), x = NULL, y = NULL) +
    theme_pub() +
    theme(
      plot.margin = margin(12, 8, 12, 8),
      legend.position = "none",
      panel.spacing.y = unit(1.6, "lines")
    ) +
    facet_grid(
      Strat_Group ~ .,
      scales = "free_y", space = "free_y",
      switch = if (show_group_strip) "y" else NULL
    )

  if (show_group_strip) {
    p <- p + theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(
        angle = 90, face = "bold", size = STRIP_SIZE, hjust = 0.5
      ),
      axis.text.y = element_text(vjust = 0, margin = margin(r = 4))
    )
  } else {
    p <- p + theme(
      strip.text.y       = element_blank(),
      strip.background.y = element_blank()
    )
  }

  if (!show_y_labels) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }

  p
}

build_and_save_combined_simple <- function(df, lags_used, output_path, width, height,
                                           strat_group_levels, strat_label_levels,
                                           show_all_y_labels = FALSE) {
  df <- prepare_combined_simple_data(df, strat_group_levels, strat_label_levels)

  q <- quantile(df$effect, c(0.005, 0.995), na.rm = TRUE)
  lim <- ceiling(max(abs(q)) / 5) * 5
  xlim <- c(-lim, lim)

  avail_lags <- sort(intersect(lags_used, unique(df$Lag)))
  if (length(avail_lags) == 0) {
    message("  [Skip] no lags found in combined data")
    return(invisible(NULL))
  }

  panels <- lapply(seq_along(avail_lags), function(i) {
    lag_v <- avail_lags[i]
    sub <- df[df$Lag == lag_v, ]
    if (nrow(sub) == 0) {
      return(NULL)
    }
    make_panel_combined_simple(
      df_lag           = sub,
      lag_val          = lag_v,
      show_y_labels    = show_all_y_labels || (i == 1),
      show_group_strip = (i == 1),
      xlim_val         = xlim
    )
  })
  panels <- Filter(Negate(is.null), panels)
  if (length(panels) == 0) {
    return(invisible(NULL))
  }

  n_p <- length(panels)
  panels[[1]] <- panels[[1]] + theme(plot.margin = margin(12, 5, 12, 15))
  if (n_p >= 2) panels[[n_p]] <- panels[[n_p]] + theme(plot.margin = margin(12, 15, 12, 5))
  if (n_p == 3) panels[[2]] <- panels[[2]] + theme(plot.margin = margin(12, 5, 12, 5))

  rel_w <- switch(as.character(n_p),
    "3" = c(1.2, 1, 1),
    "2" = c(1.2, 1),
    1
  )
  top_row <- plot_grid(plotlist = panels, align = "h", axis = "bt", ncol = n_p, rel_widths = rel_w)

  legend <- get_legend(
    ggplot(
      data.frame(
        x = c(0, 1),
        y = c(0, 1),
        cat = factor(c("Improved", "Deteriorated"),
          levels = c("Improved", "Deteriorated")
        )
      ),
      aes(xmin = x, xmax = x + 0.4, ymin = y, ymax = y + 0.4, fill = cat)
    ) +
      geom_rect() +
      scale_fill_manual(
        values = c(
          "Improved" = COLOR_IMPROVED,
          "Deteriorated" = COLOR_WORSENED
        ),
        labels = c("Improved Environment", "Deteriorated Environment")
      ) +
      theme_pub() +
      theme(legend.position = "bottom")
  )

  x_label <- draw_label(
    "MRD with Posterior Probability Distribution",
    fontface = "bold", size = BASE_SIZE
  )

  final <- plot_grid(top_row, legend, x_label, ncol = 1, rel_heights = c(1, 0.04, 0.04))
  ggsave(output_path, plot = final, width = width, height = height, dpi = 300, bg = "white")
  message("  ✓ ", basename(output_path))
}

# ── Main ──────────────────────────────────────────────────────────────────────
project_root <- normalizePath(".")
input_dir <- file.path(project_root, opt$`input-dir`)
output_dir <- file.path(project_root, opt$`output-dir`)

if (!dir.exists(input_dir)) stop("Input directory not found: ", input_dir)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

lags_used <- as.integer(str_trim(str_split(opt$lags, ",")[[1]]))

diseases <- if (is.na(opt$diseases)) {
  DISEASES_ALL
} else {
  str_trim(str_split(opt$diseases, ",")[[1]])
}

strat_types <- if (opt$`stratify-by` == "All") {
  STRAT_TYPES
} else {
  str_trim(str_split(opt$`stratify-by`, ",")[[1]])
}

message("Diseases:    ", paste(diseases, collapse = ", "))
message("Strat types: ", paste(strat_types, collapse = ", "))
message("Lags:        ", paste(lags_used, collapse = ", "))

for (disease in diseases) {
  message("\n===== ", disease, " =====")

  for (st in strat_types) {
    strat_order <- strat_order_for(st)

    # ── Full plot (6 domains) ────────────────────────────────────────────────
    message("  [", st, "] loading...")
    df <- load_strat_data(
      input_dir, disease, st, lags_used,
      DOMAIN_LEVELS, opt$`max-draws`
    )

    if (is.null(df) || nrow(df) == 0) {
      message("  [Skip] ", st, " — no matching files")
    } else {
      h <- auto_height(length(strat_order), 6)
      build_and_save(
        df,
        strat_type = st, simplified = FALSE,
        lags_used = lags_used,
        output_path = file.path(output_dir, paste0(disease, "_", st, ".png")),
        width = 15, height = h,
        show_all_y_labels = TRUE
      )
    }
  }

  message("  [Combined] loading...")
  df_main_combined <- load_main_combined_data(
    input_dir = input_dir, disease = disease, lags = lags_used,
    max_draws = opt$`max-draws`
  )
  if (!is.null(df_main_combined) && nrow(df_main_combined) > 0) {
    h_mc <- auto_height(13, 1)
    build_and_save_combined_simple(
      df_main_combined,
      lags_used = lags_used,
      output_path = file.path(output_dir, paste0(disease, "_Combined.png")),
      width = 12, height = h_mc,
      strat_group_levels = c("National", "Regime EQI", "Land Use EQI"),
      strat_label_levels = unique(c(rev(DOMAIN_LEVELS), rev(unname(LABEL_EQI)), rev(unname(LABEL_NLCD)))),
      show_all_y_labels = TRUE
    )
  }
}

message("\nDone. Output: ", output_dir)
