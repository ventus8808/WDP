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

# ─── MRR type switches (one plot per enabled type) ────────────────────────────
MRR_LagRef <- 1
MRR_SameRef <- 0

# ─── Overall-disease panel switch (1 = include top panel, 0 = subtypes only) ──
Show_Overall <- 1

# ─── MRR panel switch (1 = show MRR column, 0 = MRD only) ─────────────────────
Show_MRR <- 1

# ─── Paths ────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- if (length(script_file) > 0) dirname(normalizePath(script_file)) else getwd()
base_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
mrr_dir <- file.path(base_dir, "Result", "brms_SVI_Main_MRR")
out_dir <- file.path(base_dir, "Result", "Forest_SVI_Main")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ─── Disease -> overall ICD + ordered subtypes (icd -> readable name) ─────────
# Mirrors config.yaml diseases. Each outcome's result file is named by its
# "stem": the disease shortname for the overall, the ICD code for subtypes.
DISEASE_CONFIG <- list(
  CVD = list(overall = "I00_I99", subtypes = c(
    I20_I25 = "IHD", I60_I69 = "Stroke", I10_I15 = "HD", I50 = "HF"
  )),
  CLD = list(overall = "K70_K76_C22", subtypes = c(
    K70_K76 = "Non-cancer CLD", K70 = "Alcoholic", K71 = "K71", K73 = "K73",
    K74 = "K74", K71_K73_K74 = "Toxic & chronic", K76 = "Non-alcoholic",
    "K76.7" = "HRS", C22 = "Liver cancer"
  )),
  CRD = list(overall = "J40_J47_J60_J70_J84_D86_C34", subtypes = c(
    "J40_J47_J60_J70_J84_D86" = "Non-cancer CRD", J43_J44 = "COPD", J45 = "Asthma",
    J84_D86 = "ILD", J60_J66 = "Pneumoconiosis", C34 = "Lung cancer"
  )),
  CKD = list(overall = "N00_N29_C64_C65", subtypes = c(
    N00_N29 = "Non-cancer CKD", N18_N19 = "RF", N00_N15 = "GRTI",
    C64_C65 = "Kidney cancer"
  )),
  Suicide = list(overall = "X60_X84_Y87.0", subtypes = c(
    X60_X69 = "NVS", X70_X84 = "VS", "Y87.0" = "SI"
  )),
  NDD = list(overall = "G20_G30_G12.2_F01_F03", subtypes = c(
    G30_F01_F03 = "Dementia", G20 = "PD", G10 = "HD", G12.2 = "ALS"
  )),
  Cancer = list(overall = "C00_C97", subtypes = c(
    C18_C21 = "Colorectal", C22 = "Liver", C25 = "Pancreatic", C34 = "Lung",
    C50 = "Breast", C56 = "Ovarian", C61 = "Prostate", C64_C65 = "Kidney",
    C82_C85 = "NHL", C91_C95 = "Leukemia"
  ))
)

# ─── Derived constants ────────────────────────────────────────────────────────
CATS <- c("B", "C", "D")
CAT_LABELS <- c(B = "   B (low-mid)", C = "   C (mid-high)", D = "   D (high)")
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

# ─── Data loading (one file per outcome, named by stem) ───────────────────────
load_mrd <- function(stem) {
  path <- file.path(mrr_dir, paste0(stem, "_MRD.csv"))
  if (!file.exists(path)) {
    return(NULL)
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  df <- df[df$Lag %in% LAGS, ]
  if (nrow(df) == 0) NULL else df
}
load_mrr <- function(stem, icd, mrr_type) {
  path <- file.path(mrr_dir, paste0(stem, "_MRR_", mrr_type, ".csv"))
  if (!file.exists(path)) {
    return(NULL)
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  df <- df[df$ICD_Code == icd & df$Lag %in% LAGS, ]
  if (nrow(df) == 0) NULL else df
}

# entries: list of list(icd, name, stem). Returns a table or NULL if all missing.
build_table <- function(entries, mrr_type) {
  has_mrr <- !is.null(mrr_type)
  rows <- list()
  for (e in entries) {
    mrd <- load_mrd(e$stem)
    if (is.null(mrd)) next
    mrr <- if (has_mrr) load_mrr(e$stem, e$icd, mrr_type) else NULL

    hdr <- list(label = e$name)
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

    for (cat in CATS) {
      row <- list(label = CAT_LABELS[[cat]])
      for (i in seq_len(N_LAGS)) {
        mrd_lag <- mrd[mrd$Lag == LAGS[i] & mrd$ICD_Code == e$icd, ]
        mrd_vals <- if (nrow(mrd_lag) == 0) c(NA_real_, NA_real_, NA_real_) else parse_cell(mrd_lag[[cat]][1])
        row[[paste0("mrd_est_", i)]] <- mrd_vals[1]
        row[[paste0("mrd_ll_", i)]] <- mrd_vals[2]
        row[[paste0("mrd_hl_", i)]] <- mrd_vals[3]
        if (has_mrr) {
          mrr_lag <- if (!is.null(mrr)) mrr[mrr$Lag == LAGS[i], ] else data.frame()
          mrr_vals <- if (nrow(mrr_lag) == 0) c(NA_real_, NA_real_, NA_real_) else parse_cell(mrr_lag[[cat]][1])
          row[[paste0("mrr_est_", i)]] <- mrr_vals[1]
          row[[paste0("mrr_ll_", i)]] <- mrr_vals[2]
          row[[paste0("mrr_hl_", i)]] <- mrr_vals[3]
        }
      }
      rows[[length(rows) + 1]] <- row
    }
  }
  if (length(rows) == 0) {
    return(NULL)
  }
  dt <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  num_cols <- grep("^(mrd|mrr)_(est|ll|hl)_", names(dt), value = TRUE)
  dt[num_cols] <- lapply(dt[num_cols], as.numeric)
  dt
}

# ─── Auto axes ────────────────────────────────────────────────────────────────
auto_axis_mrd <- function(dt, include_zero = TRUE) {
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
  est <- unlist(dt[, grep("^mrr_est_", names(dt), value = TRUE)])
  est <- est[!is.na(est) & is.finite(est)]
  if (length(est) == 0) {
    return(list(xlim = c(0.5, 2.0), ticks = c(0.5, 1.0, 1.5, 2.0)))
  }
  lo <- min(c(est, 1))
  hi <- max(c(est, 1))
  pad <- 0.15 * (hi - lo)
  if (pad == 0) pad <- 0.5
  lo <- max(0, lo - pad)
  hi <- hi + pad
  tks <- pretty(c(lo, hi), n = 5)
  list(xlim = c(tks[1], tks[length(tks)]), ticks = tks)
}

# ─── Legend ───────────────────────────────────────────────────────────────────
LEGEND_ITEM_W <- 0.98
LEGEND_SEG_W <- 0.21
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
    just = c("left", "bottom"), gp = gpar(fill = "white", col = "black", lwd = 0.8)
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
      x = unit(pt_x, "in"), y = unit(item_y, "in"), pch = 16,
      size = unit(LEGEND_PTSIZE, "char"), gp = gpar(col = LAG_COLORS[i])
    )
    grid.text(
      label = LAG_LABELS[i], x = unit(seg_x1 + 0.06, "in"), y = unit(item_y, "in"),
      just = c("left", "center"), gp = gpar(fontsize = 9)
    )
  }
}

# ─── Forest theme + plotters ──────────────────────────────────────────────────
make_theme <- function() {
  forest_theme(
    base_size = 10, ci_pch = 16,
    ci_col = LAG_COLORS, ci_fill = LAG_COLORS, ci_alpha = 0.9,
    ci_lty = 1, ci_lwd = 1.5, ci_Theight = NA,
    refline_gp = gpar(lty = "dashed", col = "grey30", lwd = 1),
    legend_name = "Lag", legend_value = LAG_LABELS, legend_position = "none",
    core = list(bg = c("white", "#EBEBEB"), padding = unit(c(2, 3), "mm"))
  )
}
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
    xlim = mrd_xlim, ticks_at = mrd_ticks, xlab = "MRD (deaths/100k vs A)",
    nudge_y = 0.22, sizes = 0.32, theme = make_theme()
  )
}
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
    data = display, est = est_l, lower = lower_l, upper = upper_l,
    ci_column = c(2, 4), ref_line = c(0, 1),
    xlim = list(mrd_xlim, mrr_xlim), ticks_at = list(mrd_ticks, mrr_ticks),
    xlab = c("MRD (deaths/100k vs A)", "MRR vs A"),
    nudge_y = 0.22, sizes = 0.32, theme = make_theme()
  )
}

# ─── Save helpers ─────────────────────────────────────────────────────────────
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

PANEL_GAP <- -0.3
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
  pushViewport(viewport(
    x = unit(plot_w / 2, "in"), y = unit(total_h, "in"),
    width = unit(wh_t[1], "in"), height = unit(h_t, "in"),
    just = c("centre", "top")
  ))
  grid.draw(p_top)
  popViewport()
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

make_filename <- function(group, mrr_type = NULL) {
  parts <- c(group)
  if (length(LAGS) < length(ALL_LAGS)) parts <- c(parts, paste0("Lag", paste(LAGS, collapse = "_")))
  if (isTRUE(Show_MRR == 1) && !is.null(mrr_type)) parts <- c(parts, mrr_type)
  if (isTRUE(Show_Overall == 1)) parts <- c(parts, "Overall")
  if (isTRUE(Show_MRR == 1)) parts <- c(parts, "MRR") else parts <- c(parts, "MRD")
  paste0("SVI_", paste(parts, collapse = "_"), ".png")
}

# entries (overall + subtypes that have result files) for a disease group
group_entries <- function(group) {
  cfg <- DISEASE_CONFIG[[group]]
  ents <- list(list(icd = cfg$overall, name = group, stem = group))
  subs <- cfg$subtypes
  for (i in seq_along(subs)) {
    icd <- names(subs)[i]
    ents[[length(ents) + 1]] <- list(icd = icd, name = unname(subs[i]), stem = icd)
  }
  ents
}

# ─── Run: one plot per disease (overall + subtypes), per MRR type ─────────────
mrr_loop <- if (length(MRR_TYPES) == 0) list(NULL) else as.list(MRR_TYPES)

for (group in SELECTED) {
  ents <- group_entries(group)
  overall_ent <- ents[1]
  sub_ents <- ents[-1]
  sub_ents <- Filter(function(e) file.exists(file.path(mrr_dir, paste0(e$stem, "_MRD.csv"))), sub_ents)
  show_overall <- isTRUE(Show_Overall == 1)

  for (mrr_type in mrr_loop) {
    out_path <- file.path(out_dir, make_filename(group, mrr_type))
    has_mrr <- !is.null(mrr_type)

    p_o <- NULL
    if (show_overall) {
      dt_o <- build_table(overall_ent, mrr_type)
      if (!is.null(dt_o)) {
        if (has_mrr) {
          p_o <- make_forest_two(
            dt_o, auto_axis_mrd(dt_o, TRUE)$xlim, auto_axis_mrr(dt_o)$xlim,
            auto_axis_mrd(dt_o, TRUE)$ticks, auto_axis_mrr(dt_o)$ticks
          )
        } else {
          ax <- auto_axis_mrd(dt_o, TRUE)
          p_o <- make_forest_mrd(dt_o, ax$xlim, ax$ticks)
        }
      }
    }

    p_r <- NULL
    if (length(sub_ents) > 0) {
      dt_r <- build_table(sub_ents, mrr_type)
      if (!is.null(dt_r)) {
        if (has_mrr) {
          p_r <- make_forest_two(
            dt_r, auto_axis_mrd(dt_r)$xlim, auto_axis_mrr(dt_r)$xlim,
            auto_axis_mrd(dt_r)$ticks, auto_axis_mrr(dt_r)$ticks
          )
        } else {
          ax <- auto_axis_mrd(dt_r)
          p_r <- make_forest_mrd(dt_r, ax$xlim, ax$ticks)
        }
      }
    }

    if (!is.null(p_o) && !is.null(p_r)) {
      save_plot_combined(p_o, p_r, out_path)
    } else if (!is.null(p_o)) {
      save_plot(p_o, out_path)
    } else if (!is.null(p_r)) {
      save_plot(p_r, out_path)
    } else {
      message("Nothing to plot for ", group)
    }
  }
}
