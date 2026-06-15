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
MRR_LagRef <- 0
MRR_SameRef <- 1

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

# ─── Derived constants ────────────────────────────────────────────────────────
# SVI trajectory categories vs the A reference (A = least vulnerable, omitted).
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

# ─── Data loading ─────────────────────────────────────────────────────────────
load_mrd <- function(group) {
  path <- file.path(mrr_dir, paste0(group, "_MRD.csv"))
  if (!file.exists(path)) {
    message("MRD not found: ", path)
    return(NULL)
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  df <- df[df$Lag %in% LAGS, ]
  if (nrow(df) == 0) NULL else df
}

load_mrr <- function(group, icd, mrr_type) {
  path <- file.path(mrr_dir, paste0(group, "_MRR_", mrr_type, ".csv"))
  if (!file.exists(path)) {
    return(NULL)
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  df <- df[df$ICD_Code == icd & df$Lag %in% LAGS, ]
  if (nrow(df) == 0) NULL else df
}

# ─── Build one combined table (all selected diseases) ─────────────────────────
build_table <- function(mrr_type) {
  has_mrr <- !is.null(mrr_type)
  rows <- list()

  for (group in SELECTED) {
    mrd <- load_mrd(group)
    if (is.null(mrd)) next
    icd <- mrd$ICD_Code[1]
    mrr <- if (has_mrr) load_mrr(group, icd, mrr_type) else NULL

    # disease header row (no estimates)
    hdr <- list(label = group)
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
        mrd_lag <- mrd[mrd$Lag == LAGS[i], ]
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

  dt <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  num_cols <- grep("^(mrd|mrr)_(est|ll|hl)_", names(dt), value = TRUE)
  dt[num_cols] <- lapply(dt[num_cols], as.numeric)
  dt
}

# ─── Auto axes ────────────────────────────────────────────────────────────────
auto_axis_mrd <- function(dt, include_zero = TRUE) {
  vals <- unlist(dt[, grep("^mrd_(ll|hl)_", names(dt), value = TRUE)])
  vals <- vals[!is.na(vals) & is.finite(vals)]
  if (length(vals) == 0) return(list(xlim = c(-10, 50), ticks = seq(-10, 50, 10)))
  lo <- min(vals); hi <- max(vals)
  if (include_zero) { lo <- min(lo, 0); hi <- max(hi, 0) }
  tks <- pretty(c(lo, hi), n = 6)
  list(xlim = c(tks[1], tks[length(tks)]), ticks = tks)
}

auto_axis_mrr <- function(dt) {
  # Use point estimates (means) for the axis range; wide ratio CIs (near-zero
  # denominators) would otherwise blow up the shared axis. Over-long CIs clip.
  est <- unlist(dt[, grep("^mrr_est_", names(dt), value = TRUE)])
  est <- est[!is.na(est) & is.finite(est)]
  if (length(est) == 0) return(list(xlim = c(0.5, 2.0), ticks = c(0.5, 1.0, 1.5, 2.0)))
  lo <- min(c(est, 1)); hi <- max(c(est, 1))
  pad <- 0.15 * (hi - lo); if (pad == 0) pad <- 0.5
  lo <- max(0, lo - pad); hi <- hi + pad
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
    seg_x0 <- item_x; seg_x1 <- item_x + LEGEND_SEG_W
    pt_x <- (seg_x0 + seg_x1) / 2
    grid.segments(x0 = unit(seg_x0, "in"), y0 = unit(item_y, "in"),
                  x1 = unit(seg_x1, "in"), y1 = unit(item_y, "in"),
                  gp = gpar(col = LAG_COLORS[i], lwd = LEGEND_LWD))
    grid.points(x = unit(pt_x, "in"), y = unit(item_y, "in"), pch = 16,
                size = unit(LEGEND_PTSIZE, "char"), gp = gpar(col = LAG_COLORS[i]))
    grid.text(label = LAG_LABELS[i], x = unit(seg_x1 + 0.06, "in"), y = unit(item_y, "in"),
              just = c("left", "center"), gp = gpar(fontsize = 9))
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

# ─── Save helper ──────────────────────────────────────────────────────────────
save_plot <- function(p, path) {
  wh <- get_wh(p, unit = "in")
  plot_w <- wh[1]; plot_h <- wh[2]
  total_h <- plot_h + LEGEND_AREA
  png(path, width = plot_w, height = total_h, res = 300, units = "in")
  grid.newpage()
  pushViewport(viewport(x = 0, y = unit(total_h, "in"),
                        width = unit(plot_w, "in"), height = unit(plot_h, "in"),
                        just = c("left", "top")))
  grid.draw(p)
  popViewport()
  draw_legend(plot_w)
  dev.off()
  cat("Saved ->", path, "\n")
}

make_filename <- function(mrr_type = NULL) {
  parts <- c("SVI_Main")
  if (length(LAGS) < length(ALL_LAGS)) parts <- c(parts, paste0("Lag", paste(LAGS, collapse = "_")))
  if (isTRUE(Show_MRR == 1) && !is.null(mrr_type)) parts <- c(parts, mrr_type)
  if (isTRUE(Show_MRR == 1)) parts <- c(parts, "MRR") else parts <- c(parts, "MRD")
  paste0(paste(parts, collapse = "_"), ".png")
}

# ─── Run ──────────────────────────────────────────────────────────────────────
if (length(MRR_TYPES) == 0) {
  dt <- build_table(NULL)
  ax <- auto_axis_mrd(dt, include_zero = TRUE)
  p <- make_forest_mrd(dt, ax$xlim, ax$ticks)
  save_plot(p, file.path(out_dir, make_filename()))
} else {
  for (mrr_type in MRR_TYPES) {
    dt <- build_table(mrr_type)
    mrd_ax <- auto_axis_mrd(dt, include_zero = TRUE)
    mrr_ax <- auto_axis_mrr(dt)
    p <- make_forest_two(dt, mrd_ax$xlim, mrr_ax$xlim, mrd_ax$ticks, mrr_ax$ticks)
    save_plot(p, file.path(out_dir, make_filename(mrr_type)))
  }
}
