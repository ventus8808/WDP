library(grid)
library(forestploter)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- if (length(script_file) > 0) dirname(normalizePath(script_file)) else getwd()
base_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
csv_path <- file.path(base_dir, "Result", "Tables", "main.csv")
out_dir <- file.path(base_dir, "Result", "MRR_Forest")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

raw <- read.csv(csv_path, stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
# Plot 1: NDD, Dementia, AD, PD
DISEASE_MAP_1 <- list(
  "G20_G30_G12.2_F01_F03" = "NDD",
  "G30_F01_F03" = "Dementia",
  "G30" = "AD",
  "G20" = "PD"
)
# Plot 2: VD, ALS, HD
DISEASE_MAP_2 <- list(
  "F01" = "VD",
  "G12.2" = "ALS",
  "G10" = "HD (Control)"
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
  colnames(display) <- c("Quintile", "Relative effect", "", "Abosolute effect")

  tm <- forest_theme(
    base_size = 9,
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
    nudge_y = 0.25,
    sizes = 0.35,
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
# Plot 1: NDD / Dementia / AD / PD
# MRR range: [1.11, 2.99],  MRD range: [2.1, 29.0]
P1_MRR_XLIM <- c(0.8, 3.2)
P1_MRR_TICKS <- c(1.0, 1.5, 2.0, 2.5, 3.0)
P1_MRD_XLIM <- c(-2, 32)
P1_MRD_TICKS <- c(0, 5, 10, 15, 20, 25, 30)

# Plot 2: VD / ALS / HD (negative control)
# MRR range: [1.21, 5.53] for VD+ALS; HD extreme values clipped
# MRD range: [0.0, 3.0]
P2_MRR_XLIM <- c(0.8, 6.0)
P2_MRR_TICKS <- c(1, 2, 3, 4, 5)
P2_MRD_XLIM <- c(-0.5, 3.5)
P2_MRD_TICKS <- c(0, 1, 2, 3)

# ---------------------------------------------------------------------------
# Generate both plots
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
  show_legend = TRUE
)

# ---------------------------------------------------------------------------
# Output size — adjust P1_WIDTH and P2_WIDTH (inches) to force panel widths
# ---------------------------------------------------------------------------
P1_WIDTH <- 4 # inches — left panel (NDD / Dementia / AD / PD)
P2_WIDTH <- 4 # inches — right panel (VD / ALS / HD)

# ---------------------------------------------------------------------------
# Combine side by side and save
# ---------------------------------------------------------------------------
out_png <- file.path(out_dir, "Combined_MRR_MRD.png")

p1 <- force_width(p1, P1_WIDTH)
p2 <- force_width(p2, P2_WIDTH)

h1 <- get_wh(p1, unit = "in")[2]
h2 <- get_wh(p2, unit = "in")[2]
total_w <- P1_WIDTH + P2_WIDTH
total_h <- max(h1, h2)

png(out_png, res = 300, width = total_w, height = total_h, units = "in")
grid.newpage()

pushViewport(viewport(
  x = 0, y = 1,
  width = unit(P1_WIDTH, "in"),
  height = unit(h1, "in"),
  just = c("left", "top")
))
grid.draw(p1)
popViewport()

pushViewport(viewport(
  x = unit(P1_WIDTH, "in"), y = 1,
  width = unit(P2_WIDTH, "in"),
  height = unit(h2, "in"),
  just = c("left", "top")
))
grid.draw(p2)
popViewport()

dev.off()
cat("Saved ->", out_png, "\n")
