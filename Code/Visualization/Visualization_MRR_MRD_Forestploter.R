library(grid)
library(forestploter)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- if (length(script_file) > 0) dirname(normalizePath(script_file)) else getwd()
base_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
csv_path <- file.path(base_dir, "Result", "Tables", "NDD_MRR_MRD.csv")
out_dir <- file.path(base_dir, "Result", "MRR_Forest")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

raw <- read.csv(csv_path, stringsAsFactors = FALSE)
# Keep only the first occurrence of each key — new results are at the top of the file
raw <- raw[!duplicated(raw[, c("ICD_Code", "Lag", "Quintile")]), ]

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
# Plot 0: VD only (own axis scale, placed leftmost)
DISEASE_MAP_0 <- list(
  "F01" = "VD"
)
# Plot 1: NDD, Dementia, AD
DISEASE_MAP_1 <- list(
  "G20_G30_G12.2_F01_F03" = "NDD",
  "G30_F01_F03" = "Dementia",
  "G30" = "AD"
)
# Plot 2: PD, ALS, HD
DISEASE_MAP_2 <- list(
  "G20" = "PD",
  "G12.2" = "ALS",
  "G10" = "HD(Control)"
)

QUINTILES <- c("Q2", "Q3", "Q4", "Q5")
LAGS <- c(5, 10, 15)
LAG_LABELS <- c("5-year lag", "10-year lag", "15-year lag")
LAG_COLORS <- c("#1D3557", "#2A9D8F", "#6D597A")

# ---------------------------------------------------------------------------
# Build display table for one disease map
# ---------------------------------------------------------------------------
build_table <- function(disease_map) {
  rows <- list()

  for (icd in names(disease_map)) {
    dis_name <- disease_map[[icd]]
    sub <- raw[raw$ICD_Code == icd, ]
    if (nrow(sub) == 0) next

    # Disease header row (all NAs — no CI drawn)
    rows[[length(rows) + 1]] <- list(
      label = dis_name,
      mrr_5 = NA, mrr_10 = NA, mrr_15 = NA,
      mrr_ll5 = NA, mrr_ll10 = NA, mrr_ll15 = NA,
      mrr_hl5 = NA, mrr_hl10 = NA, mrr_hl15 = NA,
      mrd_5 = NA, mrd_10 = NA, mrd_15 = NA,
      mrd_ll5 = NA, mrd_ll10 = NA, mrd_ll15 = NA,
      mrd_hl5 = NA, mrd_hl10 = NA, mrd_hl15 = NA
    )

    for (q in QUINTILES) {
      get_val <- function(lag, col) {
        r <- sub[sub$Quintile == q & sub$Lag == lag, col]
        if (length(r) == 0) NA else r[1]
      }

      mrr <- sapply(LAGS, function(l) get_val(l, "MRR_mean"))
      mll <- sapply(LAGS, function(l) get_val(l, "MRR_lower"))
      mhl <- sapply(LAGS, function(l) get_val(l, "MRR_upper"))
      drd <- sapply(LAGS, function(l) get_val(l, "MRD_mean"))
      dll <- sapply(LAGS, function(l) get_val(l, "MRD_lower"))
      dhl <- sapply(LAGS, function(l) get_val(l, "MRD_upper"))

      # Fix inverted CIs
      mll <- pmin(mll, mrr, na.rm = FALSE)
      mhl <- pmax(mhl, mrr, na.rm = FALSE)
      dll <- pmin(dll, drd, na.rm = FALSE)
      dhl <- pmax(dhl, drd, na.rm = FALSE)

      rows[[length(rows) + 1]] <- list(
        label = paste0("   ", q),
        mrr_5 = mrr[1], mrr_10 = mrr[2], mrr_15 = mrr[3],
        mrr_ll5 = mll[1], mrr_ll10 = mll[2], mrr_ll15 = mll[3],
        mrr_hl5 = mhl[1], mrr_hl10 = mhl[2], mrr_hl15 = mhl[3],
        mrd_5 = drd[1], mrd_10 = drd[2], mrd_15 = drd[3],
        mrd_ll5 = dll[1], mrd_ll10 = dll[2], mrd_ll15 = dll[3],
        mrd_hl5 = dhl[1], mrd_hl10 = dhl[2], mrd_hl15 = dhl[3]
      )
    }
  }

  dt <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  est_cols <- c(
    "mrr_5", "mrr_10", "mrr_15",
    "mrr_ll5", "mrr_ll10", "mrr_ll15",
    "mrr_hl5", "mrr_hl10", "mrr_hl15",
    "mrd_5", "mrd_10", "mrd_15",
    "mrd_ll5", "mrd_ll10", "mrd_ll15",
    "mrd_hl5", "mrd_hl10", "mrd_hl15"
  )
  dt[est_cols] <- lapply(dt[est_cols], as.numeric)
  dt
}

# ---------------------------------------------------------------------------
# Build forest plot for one dataset
# ---------------------------------------------------------------------------
# Suppress estimates outside xlim to avoid rendering artifacts at boundary
clip_estimates <- function(est, ll, hl, xlim) {
  outside <- !is.na(est) & (est < xlim[1] | est > xlim[2])
  est[outside] <- NA
  ll[outside] <- NA
  hl[outside] <- NA
  list(est = est, ll = ll, hl = hl)
}

make_forest <- function(dt, mrr_xlim, mrd_xlim, mrr_ticks, mrd_ticks, show_legend) {
  # Clip out-of-range estimates (avoids colored artifacts at boundary)
  mrr_clipped <- clip_estimates(dt$mrr_5, dt$mrr_ll5, dt$mrr_hl5, mrr_xlim)
  dt$mrr_5 <- mrr_clipped$est
  dt$mrr_ll5 <- mrr_clipped$ll
  dt$mrr_hl5 <- mrr_clipped$hl
  mrr_clipped <- clip_estimates(dt$mrr_10, dt$mrr_ll10, dt$mrr_hl10, mrr_xlim)
  dt$mrr_10 <- mrr_clipped$est
  dt$mrr_ll10 <- mrr_clipped$ll
  dt$mrr_hl10 <- mrr_clipped$hl
  mrr_clipped <- clip_estimates(dt$mrr_15, dt$mrr_ll15, dt$mrr_hl15, mrr_xlim)
  dt$mrr_15 <- mrr_clipped$est
  dt$mrr_ll15 <- mrr_clipped$ll
  dt$mrr_hl15 <- mrr_clipped$hl

  mrd_clipped <- clip_estimates(dt$mrd_5, dt$mrd_ll5, dt$mrd_hl5, mrd_xlim)
  dt$mrd_5 <- mrd_clipped$est
  dt$mrd_ll5 <- mrd_clipped$ll
  dt$mrd_hl5 <- mrd_clipped$hl
  mrd_clipped <- clip_estimates(dt$mrd_10, dt$mrd_ll10, dt$mrd_hl10, mrd_xlim)
  dt$mrd_10 <- mrd_clipped$est
  dt$mrd_ll10 <- mrd_clipped$ll
  dt$mrd_hl10 <- mrd_clipped$hl
  mrd_clipped <- clip_estimates(dt$mrd_15, dt$mrd_ll15, dt$mrd_hl15, mrd_xlim)
  dt$mrd_15 <- mrd_clipped$est
  dt$mrd_ll15 <- mrd_clipped$ll
  dt$mrd_hl15 <- mrd_clipped$hl

  # Blank CI columns — width set by number of spaces
  dt$mrr_plot <- paste(rep(" ", 20), collapse = " ") # MRR CI column
  dt$gap <- "" # narrow spacer between MRR and MRD
  dt$mrd_plot <- paste(rep(" ", 20), collapse = " ") # MRD CI column

  # Display: label | mrr_plot | gap | mrd_plot
  display <- dt[, c("label", "mrr_plot", "gap", "mrd_plot")]
  colnames(display) <- c("Quintile", "       Relative effect", "", "  Abosolute effect")




  tm <- forest_theme(
    base_size = 10,
    ci_pch = 16,
    ci_col = LAG_COLORS,
    ci_fill = LAG_COLORS,
    ci_alpha = 0.9,
    ci_lty = 1,
    ci_lwd = 1.5,
    ci_Theight = NA, # removes the T-end caps on CI lines
    refline_gp = gpar(lty = "dashed", col = "grey30", lwd = 1),
    legend_name = "Lag",
    legend_value = LAG_LABELS,
    legend_position = if (show_legend) "bottom" else "none",
    legend_cex = 0.9, # legend symbol size — matches sizes= in forest()
    legend_lwd = 8, # legend line width — matches ci_lwd
    core = list(padding = unit(c(2, 3), "mm"))
  )

  forest(
    data = display,
    est = list(
      dt$mrr_5, dt$mrd_5,
      dt$mrr_10, dt$mrd_10,
      dt$mrr_15, dt$mrd_15
    ),
    lower = list(
      dt$mrr_ll5, dt$mrd_ll5,
      dt$mrr_ll10, dt$mrd_ll10,
      dt$mrr_ll15, dt$mrd_ll15
    ),
    upper = list(
      dt$mrr_hl5, dt$mrd_hl5,
      dt$mrr_hl10, dt$mrd_hl10,
      dt$mrr_hl15, dt$mrd_hl15
    ),
    ci_column = c(2, 4),
    ref_line = c(1, 0),
    xlim = list(mrr_xlim, mrd_xlim),
    ticks_at = list(mrr_ticks, mrd_ticks),
    xlab = c("MRR", "MRD"),
    nudge_y = 0.22,
    sizes = 0.32,
    theme = tm
  )
}

# ---------------------------------------------------------------------------
# Force plot width by scaling gtable column widths proportionally
# ---------------------------------------------------------------------------
force_width <- function(p, target_in) {
  current_in <- convertWidth(sum(p$widths), "in", valueOnly = TRUE)
  p$widths <- p$widths * (target_in / current_in)
  p
}

# ---------------------------------------------------------------------------
# Axis settings — adjust these to change plot ranges
# ---------------------------------------------------------------------------
# Plot 0: VD only
P0_MRR_XLIM <- c(0.5, 6.5)
P0_MRR_TICKS <- c(1, 2, 3, 4, 5, 6)
P0_MRD_XLIM <- c(-0.3, 4.2)
P0_MRD_TICKS <- c(0, 1, 2, 3, 4)

# Plot 1: NDD / Dementia / AD
# MRR range: [1.11, 2.99],  MRD range: [2.1, 29.0]
P1_MRR_XLIM <- c(0.9, 2)
P1_MRR_TICKS <- c(1.0, 1.2, 1.4, 1.6, 1.8)
P1_MRD_XLIM <- c(-2, 30)
P1_MRD_TICKS <- c(0, 5, 10, 15, 20, 25)

# Plot 2: PD / ALS / HD (negative control)
P2_MRR_XLIM <- c(0.8, 6.5)
P2_MRR_TICKS <- c(1, 2, 3, 4, 5, 6)
P2_MRD_XLIM <- c(-0.5, 6.5)
P2_MRD_TICKS <- c(0, 2, 4, 6)

# ---------------------------------------------------------------------------
# Generate three plots
# ---------------------------------------------------------------------------
dt0 <- build_table(DISEASE_MAP_0)
dt1 <- build_table(DISEASE_MAP_1)
dt2 <- build_table(DISEASE_MAP_2)

p0 <- make_forest(dt0,
  mrr_xlim = P0_MRR_XLIM, mrd_xlim = P0_MRD_XLIM,
  mrr_ticks = P0_MRR_TICKS, mrd_ticks = P0_MRD_TICKS,
  show_legend = FALSE
)
p1 <- make_forest(dt1,
  mrr_xlim = P1_MRR_XLIM, mrd_xlim = P1_MRD_XLIM,
  mrr_ticks = P1_MRR_TICKS, mrd_ticks = P1_MRD_TICKS,
  show_legend = FALSE
)
p2 <- make_forest(dt2,
  mrr_xlim = P2_MRR_XLIM, mrd_xlim = P2_MRD_XLIM,
  mrr_ticks = P2_MRR_TICKS, mrd_ticks = P2_MRD_TICKS,
  show_legend = FALSE # legend drawn manually below for full position/frame control
)

# ---------------------------------------------------------------------------
# Output size — adjust widths (inches) per panel
# ---------------------------------------------------------------------------
# Left column: p1 (NDD/Dementia/AD) stacked above p0 (VD), same width
LEFT_WIDTH <- 4 # inches — left column shared by p1 and p0
P2_WIDTH <- 4 # inches — right column (PD / ALS / HD)
STRIP_HEIGHT <- 0.25 # inches — white strip between p1 and p0
PLOT1_STRIP_GAP <- -0.4 # inches — gap between bottom of p1 and top of strip
PLOT0_STRIP_GAP <- -0.4 # inches — gap between bottom of strip and top of p0

# Legend position and appearance (coordinates relative to full canvas, in inches)
LEGEND_X <- 5.2 # inches from left edge of canvas
LEGEND_Y <- 0.6 # inches from bottom edge of canvas
LEGEND_W <- 1.5 # inches — legend box width
LEGEND_H <- 0.65 # inches — legend box height
LEGEND_PAD <- 0.08 # inches — inner padding
LEGEND_PTSIZE <- 0.32 # point size — match sizes= in forest()
LEGEND_LWD <- 1.5 # line width — match ci_lwd

# ---------------------------------------------------------------------------
# Combine: left column = p1 top / strip / p0 bottom; right column = p2
# Layer order (top→bottom): p1, strip, p0
# ---------------------------------------------------------------------------
out_png <- file.path(out_dir, "Combined_MRR_MRD.png")

p0 <- force_width(p0, LEFT_WIDTH)
p1 <- force_width(p1, LEFT_WIDTH)
p2 <- force_width(p2, P2_WIDTH)

h0 <- get_wh(p0, unit = "in")[2]
h1 <- get_wh(p1, unit = "in")[2]
h2 <- get_wh(p2, unit = "in")[2]
left_h <- h1 + PLOT1_STRIP_GAP + STRIP_HEIGHT + PLOT0_STRIP_GAP + h0
total_w <- LEFT_WIDTH + P2_WIDTH
total_h <- max(left_h, h2)

# y positions (measured from top of canvas)
y_p1_top <- total_h # p1 top edge
y_strip_top <- total_h - h1 - PLOT1_STRIP_GAP # strip top edge
y_p0_top <- total_h - h1 - PLOT1_STRIP_GAP - STRIP_HEIGHT - PLOT0_STRIP_GAP # p0 top edge

png(out_png, res = 300, width = total_w, height = total_h, units = "in")
grid.newpage()

# p0: VD — drawn first (bottom layer)
pushViewport(viewport(
  x = 0, y = unit(y_p0_top, "in"),
  width = unit(LEFT_WIDTH, "in"),
  height = unit(h0, "in"),
  just = c("left", "top")
))
grid.draw(p0)
popViewport()

# White strip — drawn above p0, covers p0's title row
grid.rect(
  x = unit(0, "in"), y = unit(y_strip_top, "in"),
  width = unit(LEFT_WIDTH, "in"), height = unit(STRIP_HEIGHT, "in"),
  just = c("left", "top"),
  gp = gpar(fill = "white", col = NA)
)

# p1: NDD/Dementia/AD — drawn last (top layer)
pushViewport(viewport(
  x = 0, y = unit(y_p1_top, "in"),
  width = unit(LEFT_WIDTH, "in"),
  height = unit(h1, "in"),
  just = c("left", "top")
))
grid.draw(p1)
popViewport()

# p2: PD/ALS/HD — right column
pushViewport(viewport(
  x = unit(LEFT_WIDTH, "in"), y = 1,
  width = unit(P2_WIDTH, "in"),
  height = unit(h2, "in"),
  just = c("left", "top")
))
grid.draw(p2)
popViewport()

# ---------------------------------------------------------------------------
# Manual framed legend — adjust LEGEND_* parameters above
# ---------------------------------------------------------------------------
n_lags <- length(LAG_LABELS)
row_h <- (LEGEND_H - 2 * LEGEND_PAD) / n_lags

# Frame
grid.rect(
  x = unit(LEGEND_X, "in"), y = unit(LEGEND_Y, "in"),
  width = unit(LEGEND_W, "in"), height = unit(LEGEND_H, "in"),
  just = c("left", "bottom"),
  gp = gpar(fill = "white", col = "black", lwd = 0.8)
)

for (i in seq_along(LAG_LABELS)) {
  row_y <- LEGEND_Y + LEGEND_H - LEGEND_PAD - (i - 0.5) * row_h
  seg_x0 <- LEGEND_X + LEGEND_PAD
  seg_x1 <- LEGEND_X + LEGEND_PAD + 0.3
  pt_x <- (seg_x0 + seg_x1) / 2
  # Line
  grid.segments(
    x0 = unit(seg_x0, "in"), y0 = unit(row_y, "in"),
    x1 = unit(seg_x1, "in"), y1 = unit(row_y, "in"),
    gp = gpar(col = LAG_COLORS[i], lwd = LEGEND_LWD)
  )
  # Point
  grid.points(
    x = unit(pt_x, "in"), y = unit(row_y, "in"),
    pch = 16, size = unit(LEGEND_PTSIZE, "char"),
    gp = gpar(col = LAG_COLORS[i])
  )
  # Label
  grid.text(
    label = LAG_LABELS[i],
    x = unit(seg_x1 + LEGEND_PAD, "in"), y = unit(row_y, "in"),
    just = c("left", "center"),
    gp = gpar(fontsize = 9)
  )
}

dev.off()
cat("Saved ->", out_png, "\n")
