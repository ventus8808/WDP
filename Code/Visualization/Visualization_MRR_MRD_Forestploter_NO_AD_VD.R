library(grid)
library(forestploter)

pdf(NULL) # suppress automatic Rplots.pdf creation

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
# Plot 1: NDD, Dementia
DISEASE_MAP_1 <- list(
  "G20_G30_G12.2_F01_F03" = "NDD",
  "G30_F01_F03" = "Dementia"
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
# LAG_COLORS <- c("#1D3557", "#2A9D8F", "#FFA54F")


LAG_COLORS <- c("#1D3557", "#2A9D8F", "#a98467")

# ---------------------------------------------------------------------------
# Build display table for one disease map
# ---------------------------------------------------------------------------
build_table <- function(disease_map) {
  rows <- list()

  for (icd in names(disease_map)) {
    dis_name <- disease_map[[icd]]
    sub <- raw[raw$ICD_Code == icd, ]
    if (nrow(sub) == 0) next

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
# Suppress estimates outside xlim
# ---------------------------------------------------------------------------
clip_estimates <- function(est, ll, hl, xlim) {
  # hide only when the entire CI is off-scale on one side
  hide <- (!is.na(ll) & ll > xlim[2]) | (!is.na(hl) & hl < xlim[1])
  # cap point and CI endpoints at axis boundaries (partial CIs remain visible)
  est <- pmax(pmin(est, xlim[2]), xlim[1])
  ll  <- pmax(ll, xlim[1])
  hl  <- pmin(hl, xlim[2])
  est[hide] <- NA; ll[hide] <- NA; hl[hide] <- NA
  list(est = est, ll = ll, hl = hl)
}

# ---------------------------------------------------------------------------
# Build forest plot for one dataset
# ---------------------------------------------------------------------------
make_forest <- function(dt, mrr_xlim, mrd_xlim, mrr_ticks, mrd_ticks, show_legend) {
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

  dt$mrr_plot <- paste(rep(" ", 20), collapse = " ")
  dt$gap <- ""
  dt$mrd_plot <- paste(rep(" ", 20), collapse = " ")

  # ✅ SWAPPED ORDER: MRD first, MRR second
  display <- dt[, c("label", "mrd_plot", "gap", "mrr_plot")]
  colnames(display) <- c("Quintile", "     Absolute effect", "", "    Relative effect")

  tm <- forest_theme(
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
    legend_position = if (show_legend) "bottom" else "none",
    legend_cex = 0.9,
    legend_lwd = 8,
    core = list(bg = c("white", "#EBEBEB"), padding = unit(c(2, 3), "mm"))
  )

  forest(
    data = display,
    est = list(
      dt$mrd_5, dt$mrr_5,
      dt$mrd_10, dt$mrr_10,
      dt$mrd_15, dt$mrr_15
    ),
    lower = list(
      dt$mrd_ll5, dt$mrr_ll5,
      dt$mrd_ll10, dt$mrr_ll10,
      dt$mrd_ll15, dt$mrr_ll15
    ),
    upper = list(
      dt$mrd_hl5, dt$mrr_hl5,
      dt$mrd_hl10, dt$mrr_hl10,
      dt$mrd_hl15, dt$mrr_hl15
    ),
    ci_column = c(2, 4),
    ref_line = c(0, 1),
    xlim = list(mrd_xlim, mrr_xlim),
    ticks_at = list(mrd_ticks, mrr_ticks),
    xlab = c("MRD", "MRR"),
    nudge_y = 0.22,
    sizes = 0.32,
    theme = tm
  )
}

# ---------------------------------------------------------------------------
# Force plot width
# ---------------------------------------------------------------------------
force_width <- function(p, target_in) {
  current_in <- convertWidth(sum(p$widths), "in", valueOnly = TRUE)
  p$widths <- p$widths * (target_in / current_in)
  p
}

# ---------------------------------------------------------------------------
# Axis settings
# ---------------------------------------------------------------------------
P1_MRR_XLIM <- c(0.9, 2)
P1_MRR_TICKS <- c(1.0, 1.2, 1.4, 1.6, 1.8)
P1_MRD_XLIM <- c(-5, 32)
P1_MRD_TICKS <- c(0, 5, 10, 15, 20, 25, 30)

P2_MRR_XLIM <- c(0.7, 6.5)
P2_MRR_TICKS <- c(1, 2, 3, 4, 5, 6)
P2_MRD_XLIM <- c(-0.45, 6)
P2_MRD_TICKS <- c(0, 1, 2, 3, 4, 5)

# ---------------------------------------------------------------------------
# Generate plots
# ---------------------------------------------------------------------------
dt1 <- build_table(DISEASE_MAP_1)
dt2 <- build_table(DISEASE_MAP_2)

p1 <- make_forest(dt1,
  mrr_xlim = P1_MRR_XLIM, mrd_xlim = P1_MRD_XLIM,
  mrr_ticks = P1_MRR_TICKS, mrd_ticks = P1_MRD_TICKS,
  show_legend = FALSE
)

p2 <- make_forest(dt2,
  mrr_xlim = P2_MRR_XLIM, mrd_xlim = P2_MRD_XLIM,
  mrr_ticks = P2_MRR_TICKS, mrd_ticks = P2_MRD_TICKS,
  show_legend = FALSE
)

# ---------------------------------------------------------------------------
# Layout & output
# ---------------------------------------------------------------------------
PLOT_WIDTH <- 5
STRIP_HEIGHT <- 0.25
INTER_GAP <- -0.4 # negative = plots overlap slightly; white strip covers seam

LEGEND_ITEM_W <- 1.4 # width per item (seg + label)
LEGEND_SEG_W <- 0.3 # line segment width
LEGEND_PAD <- 0.12 # padding left of each item
LEGEND_H <- 0.35 # box height
LEGEND_W <- length(LAG_LABELS) * LEGEND_ITEM_W + 2 * LEGEND_PAD
LEGEND_BOTTOM_MARGIN <- 0.18
LEGEND_AREA <- LEGEND_BOTTOM_MARGIN + LEGEND_H + 0.15
LEGEND_PTSIZE <- 0.32
LEGEND_LWD <- 1.5

out_png <- file.path(out_dir, "Combined_MRR_MRD_No_AD_VD.png")

p1 <- force_width(p1, PLOT_WIDTH)
p2 <- force_width(p2, PLOT_WIDTH)

h1 <- get_wh(p1, unit = "in")[2]
h2 <- get_wh(p2, unit = "in")[2]

# Total height: stacked plots + legend area at bottom
plot_stack_h <- h1 + INTER_GAP + STRIP_HEIGHT + INTER_GAP + h2
total_w <- PLOT_WIDTH
total_h <- plot_stack_h + LEGEND_AREA

# Y positions from canvas bottom (just="top" means y is the top edge)
y_p1_top <- total_h
y_strip_top <- y_p1_top - h1 - INTER_GAP
y_p2_top <- y_strip_top - STRIP_HEIGHT - INTER_GAP

png(out_png, res = 300, width = total_w, height = total_h, units = "in")
grid.newpage()

# p2 (PD/ALS/HD) below p1
pushViewport(viewport(
  x = 0, y = unit(y_p2_top, "in"),
  width = unit(PLOT_WIDTH, "in"),
  height = unit(h2, "in"),
  just = c("left", "top")
))
grid.draw(p2)
popViewport()

# White strip to hide axis overlap between p1 and p2
grid.rect(
  x = unit(0, "in"), y = unit(y_strip_top, "in"),
  width = unit(PLOT_WIDTH, "in"), height = unit(STRIP_HEIGHT, "in"),
  just = c("left", "top"),
  gp = gpar(fill = "white", col = NA)
)

# p1 (NDD/Dementia) on top
pushViewport(viewport(
  x = 0, y = unit(y_p1_top, "in"),
  width = unit(PLOT_WIDTH, "in"),
  height = unit(h1, "in"),
  just = c("left", "top")
))
grid.draw(p1)
popViewport()

# Centered horizontal legend below both plots
legend_x <- (PLOT_WIDTH - LEGEND_W) / 2
legend_y <- LEGEND_BOTTOM_MARGIN
item_y <- legend_y + LEGEND_H / 2 # vertical center of box

grid.rect(
  x = unit(legend_x, "in"), y = unit(legend_y, "in"),
  width = unit(LEGEND_W, "in"), height = unit(LEGEND_H, "in"),
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

dev.off()
cat("Saved ->", out_png, "\n")
