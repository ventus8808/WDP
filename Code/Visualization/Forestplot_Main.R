library(grid)
library(forestploter)

pdf(NULL)

# ─── Disease switches ──────────────────────────────────────────────────────────
CVD <- 1
CLD <- 1
CRD <- 1
CKD <- 1
Suicide <- 1
NDD <- 1
Cancer <- 1

# ─── Lag switches ─────────────────────────────────────────────────────────────
Lag_5 <- 1
Lag_10 <- 1
Lag_15 <- 1
Lag_20 <- 1

# ─── MRR type switches (generates one plot per enabled type) ──────────────────
MRR_LagRef <- 0
MRR_SameRef <- 1

# ─── Overall-disease panel switch (1 = include top panel, 0 = others only) ────
Show_Overall <- 1

# ─── MRR panel switch (1 = show MRR column, 0 = MRD only) ─────────────────────
Show_MRR <- 1

# ─── Paths ────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- if (length(script_file) > 0) dirname(normalizePath(script_file)) else getwd()
base_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
extract_dir <- file.path(base_dir, "Result", "Extraction")
mrr_dir <- file.path(base_dir, "Result", "brms_Main_MRR")
out_dir <- file.path(base_dir, "Result", "Forest_Main")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ─── Derived constants ────────────────────────────────────────────────────────
QUINTILES <- c("Q2", "Q3", "Q4", "Q5")
DISEASE_NAMES <- c("CVD", "CLD", "CRD", "CKD", "Suicide", "NDD", "Cancer")
DISEASE_SW <- c(CVD, CLD, CRD, CKD, Suicide, NDD, Cancer)
SELECTED <- DISEASE_NAMES[DISEASE_SW == 1]
MRR_TYPES <- if (isTRUE(Show_MRR == 1)) c("LagRef", "SameRef")[c(MRR_LagRef, MRR_SameRef) == 1] else character(0)

ALL_LAGS <- c(5L, 10L, 15L, 20L)
LAG_ON <- c(Lag_5, Lag_10, Lag_15, Lag_20)
LAGS <- ALL_LAGS[LAG_ON == 1]
N_LAGS <- length(LAGS)
LAG_LABELS <- paste0(LAGS, "-year lag")
LAG_COLORS <- c("#2A9D8F", "#a98467", "#1685a9", "#1D3557")[seq_len(N_LAGS)]

# ─── Helpers ──────────────────────────────────────────────────────────────────
parse_cell <- function(x) {
  if (is.na(x) || trimws(x) %in% c("", "0.0", "0.00", "0")) {
    return(c(NA_real_, NA_real_, NA_real_))
  }
  nums <- as.numeric(regmatches(x, gregexpr("-?[0-9]+\\.?[0-9]*", x, perl = TRUE))[[1]])
  if (length(nums) != 3) {
    return(c(NA_real_, NA_real_, NA_real_))
  }
  nums
}

blank_row <- function(has_mrr) {
  r <- list(label = "")
  for (i in seq_len(N_LAGS)) {
    r[[paste0("mrd_est_", i)]] <- NA_real_
    r[[paste0("mrd_ll_", i)]] <- NA_real_
    r[[paste0("mrd_hl_", i)]] <- NA_real_
    if (has_mrr) {
      r[[paste0("mrr_est_", i)]] <- NA_real_
      r[[paste0("mrr_ll_", i)]] <- NA_real_
      r[[paste0("mrr_hl_", i)]] <- NA_real_
    }
  }
  r
}

# ─── Data loading ─────────────────────────────────────────────────────────────
load_mrd <- function(group) {
  path <- file.path(extract_dir, paste0(group, "_Main.csv"))
  if (!file.exists(path)) {
    message("MRD not found: ", path)
    return(NULL)
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  df <- df[df$Model == "EQI" & df$Lag %in% LAGS, ]
  if (nrow(df) == 0) NULL else df
}

load_mrr <- function(icd, group, mrr_type) {
  # Try per-ICD file first, fall back to group-level file
  for (stem in c(icd, group)) {
    path <- file.path(mrr_dir, paste0(stem, "_MRR_", mrr_type, ".csv"))
    if (file.exists(path)) {
      df <- read.csv(path, stringsAsFactors = FALSE)
      df <- df[df$ICD_Code == icd & df$EQI_Period == "2000_2005" & df$Lag %in% LAGS, ]
      if (nrow(df) > 0) {
        return(df)
      }
    }
  }
  NULL
}

# ─── Build combined table ─────────────────────────────────────────────────────
build_table <- function(group, mrd_df, mrr_type, icds = NULL) {
  has_mrr <- !is.null(mrr_type)
  rows <- list()
  icd_list <- if (is.null(icds)) unique(mrd_df$ICD_Code) else icds

  for (icd in icd_list) {
    mrd_sub <- mrd_df[mrd_df$ICD_Code == icd, ]
    mrr_sub <- if (has_mrr) load_mrr(icd, group, mrr_type) else NULL

    hdr <- list(label = mrd_sub$Outcome[1])
    for (i in seq_len(N_LAGS)) {
      hdr[[paste0("mrd_est_", i)]] <- NA_real_
      hdr[[paste0("mrd_ll_", i)]] <- NA_real_
      hdr[[paste0("mrd_hl_", i)]] <- NA_real_
      if (has_mrr) {
        hdr[[paste0("mrr_est_", i)]] <- NA_real_
        hdr[[paste0("mrr_ll_", i)]] <- NA_real_
        hdr[[paste0("mrr_hl_", i)]] <- NA_real_
      }
    }
    rows[[length(rows) + 1]] <- hdr

    for (q in QUINTILES) {
      row <- list(label = paste0("   ", q))
      for (i in seq_len(N_LAGS)) {
        mrd_lag <- mrd_sub[mrd_sub$Lag == LAGS[i], ]
        mrd_vals <- if (nrow(mrd_lag) == 0) c(NA_real_, NA_real_, NA_real_) else parse_cell(mrd_lag[[q]][1])
        row[[paste0("mrd_est_", i)]] <- mrd_vals[1]
        row[[paste0("mrd_ll_", i)]] <- mrd_vals[2]
        row[[paste0("mrd_hl_", i)]] <- mrd_vals[3]

        if (has_mrr) {
          mrr_lag <- if (!is.null(mrr_sub)) mrr_sub[mrr_sub$Lag == LAGS[i], ] else data.frame()
          mrr_vals <- if (nrow(mrr_lag) == 0) c(NA_real_, NA_real_, NA_real_) else parse_cell(mrr_lag[[q]][1])
          row[[paste0("mrr_est_", i)]] <- mrr_vals[1]
          row[[paste0("mrr_ll_", i)]] <- mrr_vals[2]
          row[[paste0("mrr_hl_", i)]] <- mrr_vals[3]
        }
      }
      rows[[length(rows) + 1]] <- row
    }
  }

  dt <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  num_cols <- grep("^(mrd|mrr)_(est|ll|hl)_", names(dt), value = TRUE)
  dt[num_cols] <- lapply(dt[num_cols], as.numeric)
  dt
}

# ─── Auto axes ────────────────────────────────────────────────────────────────
auto_axis_mrd <- function(dt, include_zero = FALSE) {
  vals <- unlist(dt[, grep("^mrd_(ll|hl)_", names(dt), value = TRUE)])
  vals <- vals[!is.na(vals) & is.finite(vals)]
  if (length(vals) == 0) {
    return(list(xlim = c(-10, 50), ticks = seq(-10, 50, 10)))
  }
  lo <- min(vals)
  hi <- max(vals)
  if (include_zero) {
    lo <- min(lo, 0)
    hi <- max(hi, 0)
  }
  tks <- pretty(c(lo, hi), n = 6)
  list(xlim = c(tks[1], tks[length(tks)]), ticks = tks)
}

auto_axis_mrr <- function(dt) {
  vals <- unlist(dt[, grep("^mrr_(ll|hl)_", names(dt), value = TRUE)])
  vals <- vals[!is.na(vals) & is.finite(vals)]
  if (length(vals) == 0) {
    return(list(xlim = c(0.5, 2.0), ticks = c(0.5, 1.0, 1.5, 2.0)))
  }
  lo <- min(c(vals, 1))
  hi <- max(c(vals, 1)) # ensure ref=1 is included
  rng <- hi - lo
  step <- if (rng > 3) 1 else if (rng > 1) 0.5 else if (rng > 0.5) 0.2 else 0.1
  xlim <- c(floor(lo / step) * step, ceiling(hi / step) * step)
  list(xlim = xlim, ticks = round(seq(xlim[1], xlim[2], by = step), 4))
}

# ─── Legend constants ─────────────────────────────────────────────────────────
LEGEND_ITEM_W <- 0.98 # 30% narrower than original 1.4
LEGEND_SEG_W <- 0.21 # 30% narrower than original 0.3
LEGEND_PAD <- 0.08
LEGEND_H <- 0.32
LEGEND_BOTTOM_MARGIN <- 0.18
LEGEND_AREA <- LEGEND_BOTTOM_MARGIN + LEGEND_H + 0.15
LEGEND_PTSIZE <- 0.32
LEGEND_LWD <- 1.5

draw_legend <- function(plot_w) {
  legend_w <- N_LAGS * LEGEND_ITEM_W + 2 * LEGEND_PAD
  legend_x <- (plot_w - legend_w) / 2
  legend_y <- LEGEND_BOTTOM_MARGIN
  item_y <- legend_y + LEGEND_H / 2

  grid.rect(
    x = unit(legend_x, "in"), y = unit(legend_y, "in"),
    width = unit(legend_w, "in"), height = unit(LEGEND_H, "in"),
    just = c("left", "bottom"),
    gp = gpar(fill = "white", col = "black", lwd = 0.8)
  )
  for (i in seq_along(LAG_LABELS)) {
    item_x <- legend_x + LEGEND_PAD + (i - 1) * LEGEND_ITEM_W
    seg_x0 <- item_x
    seg_x1 <- item_x + LEGEND_SEG_W
    pt_x <- (seg_x0 + seg_x1) / 2
    grid.segments(
      x0 = unit(seg_x0, "in"), y0 = unit(item_y, "in"),
      x1 = unit(seg_x1, "in"), y1 = unit(item_y, "in"),
      gp = gpar(col = LAG_COLORS[i], lwd = LEGEND_LWD)
    )
    grid.points(
      x = unit(pt_x, "in"), y = unit(item_y, "in"),
      pch = 16, size = unit(LEGEND_PTSIZE, "char"),
      gp = gpar(col = LAG_COLORS[i])
    )
    grid.text(
      label = LAG_LABELS[i],
      x = unit(seg_x1 + 0.06, "in"), y = unit(item_y, "in"),
      just = c("left", "center"),
      gp = gpar(fontsize = 9)
    )
  }
}

# ─── Forest theme (shared) ────────────────────────────────────────────────────
make_theme <- function() {
  forest_theme(
    base_size = 10,
    ci_pch = 16,
    ci_col = LAG_COLORS,
    ci_fill = LAG_COLORS,
    ci_alpha = 0.9,
    ci_lty = 1,
    ci_lwd = 1.5,
    ci_Theight = NA,
    refline_gp = gpar(lty = "dashed", col = "grey30", lwd = 1),
    legend_name = "Lag",
    legend_value = LAG_LABELS,
    legend_position = "none",
    core = list(bg = c("white", "#EBEBEB"), padding = unit(c(2, 3), "mm"))
  )
}

# ─── MRD-only plot ────────────────────────────────────────────────────────────
make_forest_mrd <- function(dt, mrd_xlim, mrd_ticks) {
  dt$plot_col <- paste(rep(" ", 30), collapse = " ")
  display <- dt[, c("label", "plot_col")]
  colnames(display) <- c("Outcome", "  MRD")

  forest(
    data = display,
    est = lapply(seq_len(N_LAGS), function(i) dt[[paste0("mrd_est_", i)]]),
    lower = lapply(seq_len(N_LAGS), function(i) dt[[paste0("mrd_ll_", i)]]),
    upper = lapply(seq_len(N_LAGS), function(i) dt[[paste0("mrd_hl_", i)]]),
    ci_column = 2, ref_line = 0,
    xlim = mrd_xlim, ticks_at = mrd_ticks,
    xlab = "MRD",
    nudge_y = 0.22, sizes = 0.32, theme = make_theme()
  )
}

# ─── MRD + MRR two-panel plot ─────────────────────────────────────────────────
make_forest_two <- function(dt, mrd_xlim, mrr_xlim, mrd_ticks, mrr_ticks) {
  dt$mrd_col <- paste(rep(" ", 20), collapse = " ")
  dt$gap <- paste(rep(" ", 6), collapse = " ")
  dt$mrr_col <- paste(rep(" ", 20), collapse = " ")
  display <- dt[, c("label", "mrd_col", "gap", "mrr_col")]
  colnames(display) <- c("Outcome", "               MRD               ", "          ", "               MRR               ")

  est_l <- lower_l <- upper_l <- vector("list", 2 * N_LAGS)
  for (i in seq_len(N_LAGS)) {
    est_l[[2 * i - 1]] <- dt[[paste0("mrd_est_", i)]]
    est_l[[2 * i]] <- dt[[paste0("mrr_est_", i)]]
    lower_l[[2 * i - 1]] <- dt[[paste0("mrd_ll_", i)]]
    lower_l[[2 * i]] <- dt[[paste0("mrr_ll_", i)]]
    upper_l[[2 * i - 1]] <- dt[[paste0("mrd_hl_", i)]]
    upper_l[[2 * i]] <- dt[[paste0("mrr_hl_", i)]]
  }

  forest(
    data = display,
    est = est_l, lower = lower_l, upper = upper_l,
    ci_column = c(2, 4), ref_line = c(0, 1),
    xlim = list(mrd_xlim, mrr_xlim),
    ticks_at = list(mrd_ticks, mrr_ticks),
    xlab = c("MRD", "MRR"),
    nudge_y = 0.22, sizes = 0.32, theme = make_theme()
  )
}

# ─── Save helper ──────────────────────────────────────────────────────────────
save_plot <- function(p, path) {
  wh <- get_wh(p, unit = "in")
  plot_w <- wh[1]
  plot_h <- wh[2]
  total_h <- plot_h + LEGEND_AREA

  png(path, width = plot_w, height = total_h, res = 300, units = "in")
  grid.newpage()
  pushViewport(viewport(
    x = 0, y = unit(total_h, "in"),
    width = unit(plot_w, "in"), height = unit(plot_h, "in"),
    just = c("left", "top")
  ))
  grid.draw(p)
  popViewport()
  draw_legend(plot_w)
  dev.off()
  cat("Saved ->", path, "\n")
}

# ─── Combined save helper (Overall stacked above Others) ──────────────────────
PANEL_GAP <- -0.3

# Make per-column widths identical between two forestploter gtables so the
# panels line up exactly, even when the bottom plot has longer outcome labels.
align_widths <- function(p_top, p_bot) {
  n <- min(length(p_top$widths), length(p_bot$widths))
  for (i in seq_len(n)) {
    w_t <- convertWidth(p_top$widths[i], "in", valueOnly = TRUE)
    w_b <- convertWidth(p_bot$widths[i], "in", valueOnly = TRUE)
    mw <- max(w_t, w_b)
    p_top$widths[i] <- unit(mw, "in")
    p_bot$widths[i] <- unit(mw, "in")
  }
  list(top = p_top, bot = p_bot)
}

save_plot_combined <- function(p_top, p_bot, path) {
  aligned <- align_widths(p_top, p_bot)
  p_top <- aligned$top
  p_bot <- aligned$bot

  wh_t <- get_wh(p_top, unit = "in")
  wh_b <- get_wh(p_bot, unit = "in")
  plot_w <- max(wh_t[1], wh_b[1])
  h_t <- wh_t[2]
  h_b <- wh_b[2]
  total_h <- h_t + PANEL_GAP + h_b + LEGEND_AREA

  png(path, width = plot_w, height = total_h, res = 300, units = "in")
  grid.newpage()
  # Top panel: overall — centered within max width so panels align visually
  pushViewport(viewport(
    x = unit(plot_w / 2, "in"), y = unit(total_h, "in"),
    width = unit(wh_t[1], "in"), height = unit(h_t, "in"),
    just = c("centre", "top")
  ))
  grid.draw(p_top)
  popViewport()
  # Bottom panel: others — centered within max width
  pushViewport(viewport(
    x = unit(plot_w / 2, "in"), y = unit(total_h - h_t - PANEL_GAP, "in"),
    width = unit(wh_b[1], "in"), height = unit(h_b, "in"),
    just = c("centre", "top")
  ))
  grid.draw(p_bot)
  popViewport()
  draw_legend(plot_w)
  dev.off()
  cat("Saved ->", path, "\n")
}

# ─── Filename helper ──────────────────────────────────────────────────────────
# Builds filenames like:
#   CVD_SameRef_Overall_MRR.png        (all 4 lags, MRR + Overall on, SameRef)
#   CVD_Overall.png                    (all 4 lags, MRR off, Overall on)
#   CVD_Lag5_15_SameRef_MRR.png        (lags 5+15 only, MRR on, Overall off)
make_filename <- function(group, mrr_type = NULL) {
  parts <- c(group)
  if (length(LAGS) < length(ALL_LAGS)) {
    parts <- c(parts, paste0("Lag", paste(LAGS, collapse = "_")))
  }
  if (isTRUE(Show_MRR == 1) && !is.null(mrr_type)) parts <- c(parts, mrr_type)
  if (isTRUE(Show_Overall == 1)) parts <- c(parts, "Overall")
  if (isTRUE(Show_MRR == 1)) parts <- c(parts, "MRR")
  paste0(paste(parts, collapse = "_"), ".png")
}

# ─── Run: one plot per disease × MRR type ─────────────────────────────────────
for (group in SELECTED) {
  mrd_df <- load_mrd(group)
  if (is.null(mrd_df)) next

  all_icds <- unique(mrd_df$ICD_Code)
  overall_icd <- all_icds[1] # first ICD is the overall disease (e.g., I00_I99 for CVD)
  other_icds <- all_icds[-1]
  has_others <- length(other_icds) > 0

  show_overall <- isTRUE(Show_Overall == 1)

  if (length(MRR_TYPES) == 0) {
    out_path <- file.path(out_dir, make_filename(group))
    if (show_overall) {
      dt_o <- build_table(group, mrd_df, NULL, icds = overall_icd)
      ax_o <- auto_axis_mrd(dt_o, include_zero = TRUE)
      p_o <- make_forest_mrd(dt_o, ax_o$xlim, ax_o$ticks)
    }
    if (!has_others) {
      if (show_overall) {
        save_plot(p_o, out_path)
      } else {
        message("Nothing to plot for ", group, " (no sub-outcomes and Show_Overall=0)")
        next
      }
    } else {
      dt_r <- build_table(group, mrd_df, NULL, icds = other_icds)
      ax_r <- auto_axis_mrd(dt_r)
      p_r <- make_forest_mrd(dt_r, ax_r$xlim, ax_r$ticks)
      if (show_overall) {
        save_plot_combined(p_o, p_r, out_path)
      } else {
        save_plot(p_r, out_path)
      }
    }
  } else {
    for (mrr_type in MRR_TYPES) {
      out_path <- file.path(out_dir, make_filename(group, mrr_type))
      if (show_overall) {
        dt_o <- build_table(group, mrd_df, mrr_type, icds = overall_icd)
        mrd_ax_o <- auto_axis_mrd(dt_o, include_zero = TRUE)
        mrr_ax_o <- auto_axis_mrr(dt_o)
        p_o <- make_forest_two(dt_o, mrd_ax_o$xlim, mrr_ax_o$xlim, mrd_ax_o$ticks, mrr_ax_o$ticks)
      }
      if (!has_others) {
        if (show_overall) {
          save_plot(p_o, out_path)
        } else {
          message("Nothing to plot for ", group, " (no sub-outcomes and Show_Overall=0)")
          next
        }
      } else {
        dt_r <- build_table(group, mrd_df, mrr_type, icds = other_icds)
        mrd_ax_r <- auto_axis_mrd(dt_r)
        mrr_ax_r <- auto_axis_mrr(dt_r)
        p_r <- make_forest_two(dt_r, mrd_ax_r$xlim, mrr_ax_r$xlim, mrd_ax_r$ticks, mrr_ax_r$ticks)
        if (show_overall) {
          save_plot_combined(p_o, p_r, out_path)
        } else {
          save_plot(p_r, out_path)
        }
      }
    }
  }
}
