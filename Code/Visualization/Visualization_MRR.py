"""
Forest plot — MRR (transposed)
Layout : 2 column groups via GridSpecFromSubplotSpec
  Left:  NDD / Dementia / HD  — height ratios 3 : 3 : 2
  Right: AD  / PD      / ALS  — height ratios 2 : 2 : 2
Y-axis: fixed ranges; tick marks on right; label contra-rotated (-90°) on right of Q5
Disease labels: rotation=90 (reads bottom → top) along the plot height
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from matplotlib.transforms import blended_transform_factory

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
plt.rcParams.update(
    {
        "font.family": "Helvetica",
        "font.size": 14,
        "axes.titlesize": 12,
        "axes.labelsize": 12,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
    }
)

DPI = 300
FRAME_COLOR = "#BBBBBB"
ROUND_SIZE = 0.06

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path(__file__).resolve().parents[2]
MRR_CSV = BASE / "Result" / "Tables" / "MRR.csv"
OUT_DIR = BASE / "Result" / "MRR_Forest"
OUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(MRR_CSV)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
LEFT_DISEASES = [
    ("G20_G30_G12.2_F01_F03", "(A) NDD"),
    ("G30", "(C) AD"),
    ("G12.2", "(E) ALS"),
]
RIGHT_DISEASES = [
    ("G30_F01_F03", "(B) Dementia"),
    ("G20", "(D) PD"),
    ("G10", "(F) HD (negative control)"),
]

QUINTILES = ["Q2", "Q3", "Q4", "Q5"]
N_COLS = len(QUINTILES)

LAGS = [5, 10, 15]

# Okabe-Ito — color-blind safe
LAG_COLORS = {
    5: "#2171B5",  # 深蓝
    10: "#6BAED6",  # 中蓝
    15: "#B2D8F0",  # 浅蓝
}
# LAG_COLORS = {5: "#1F77B4",10: "#FF7F0E",15: "#2CA02C",}
# LAG_COLORS = {5: "#E69F00", 10: "#0072B2", 15: "#009E73"}
LAG_LABELS = {5: "5-year lag", 10: "10-year lag", 15: "15-year lag"}

LAG_X = {5: 0.0, 10: 1.0, 15: 2.0}
X_LIM = (-0.6, 2.6)

# Fixed y-limits  (start, end of axis)
DISEASE_YLIM = {
    "G20_G30_G12.2_F01_F03": (0.8, 1.8),  # NDD
    "G30_F01_F03": (0.8, 1.8),  # Dementia
    "G10": (0.0, 6),  # HD  (0,1,10)
    "G30": (0.8, 1.8),  # AD
    "G20": (0.5, 3.0),  # PD
    "G12.2": (0, 6),  # ALS
}


DISEASE_YTICKS = {
    "G20_G30_G12.2_F01_F03": [1, 1.2, 1.4, 1.6],
    "G30_F01_F03": [1, 1.2, 1.4, 1.6],
    "G10": [1, 2, 3, 4, 5],
    "G30": [1, 1.2, 1.4, 1.6],
    "G20": [1, 1.5, 2.0, 2.5],
    "G12.2": [1, 2, 3, 4, 5],
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def style_ax(ax):
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_facecolor("white")
    # Compute aspect correction so the corner radius is equal in physical units.
    # rounding_size is in axes coords (0–1); without correction the corner is
    # 6% of width in x but 6% of height in y, which looks distorted for tall panels.
    fig = ax.get_figure()
    fw, fh = fig.get_size_inches()
    pos = ax.get_position()
    ax_w = pos.width * fw  # axes width  in inches
    ax_h = pos.height * fh  # axes height in inches
    mut_aspect = ax_w / ax_h if ax_h > 0 else 1.0
    frame = FancyBboxPatch(
        (0, 0),
        1,
        1,
        boxstyle=f"round,pad=0,rounding_size={ROUND_SIZE}",
        mutation_aspect=mut_aspect,
        facecolor="none",
        edgecolor=FRAME_COLOR,
        linewidth=1.1,
        transform=ax.transAxes,
        zorder=5,
        clip_on=False,
    )
    ax.add_patch(frame)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.yaxis.set_ticks_position("left")
    ax.tick_params(
        axis="y",
        which="both",
        direction="in",
        length=3.5,
        width=0.7,
        right=False,
        left=True,
    )
    ax.set_axisbelow(True)


def apply_yticks(ax, icd_code, col_idx, show_ylabel=True):
    """
    Tick marks on left of all panels.
    Numeric labels only on Q2 (col_idx == 0).
    Y-axis text label shown on Q2 only when show_ylabel is True.
    """
    ticks = DISEASE_YTICKS[icd_code]
    ax.set_yticks(ticks)

    if col_idx == 0:  # Q2 — show numeric labels on left
        labels = []
        for t in ticks:
            if abs(t - round(t)) < 1e-9:
                labels.append(str(int(round(t))))
            else:
                labels.append(f"{t:.1f}")
        ax.set_yticklabels(labels, fontsize=8.5)
        if show_ylabel:
            ax.text(
                -0.3,
                0.5,
                "MRR (95% CrI)",
                transform=ax.transAxes,
                rotation=90,
                va="center",
                ha="center",
                fontsize=9.5,
            )
    else:
        ax.set_yticklabels([])


def draw_ref_line(ax):
    """Horizontal dashed line at MRR = 1."""
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    ax.plot(
        [0, 1],
        [1.0, 1.0],
        color="#444444",
        linewidth=0.85,
        linestyle="--",
        transform=trans,
        clip_on=False,
        zorder=2,
        alpha=0.55,
    )


def draw_x_separators(ax):
    """Faint dotted vertical lines between lag positions."""
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    for sep_x in [0.5, 1.5]:
        ax.plot(
            [sep_x, sep_x],
            [0, 1],
            color=FRAME_COLOR,
            linewidth=0.5,
            linestyle=":",
            transform=trans,
            clip_on=False,
            zorder=1,
        )


def plot_point(ax, xpos, mean, lower, upper, color, sig, ymin, ymax):
    """
    Vertical errorbar at xpos.
    CI bounds clipped to [ymin, ymax]; out-of-range ends get an arrow.
    """
    mean_c = np.clip(mean, ymin, ymax)
    lower_c = np.clip(lower, ymin, ymax)
    upper_c = np.clip(upper, ymin, ymax)

    ax.errorbar(
        xpos,
        mean_c,
        yerr=[[max(mean_c - lower_c, 0)], [max(upper_c - mean_c, 0)]],
        fmt="o",
        color=color,
        markerfacecolor=color if sig else "white",
        markeredgecolor=color,
        markeredgewidth=1.5,
        markersize=5.0,
        capsize=3.0,
        capthick=1.2,
        linewidth=1.2,
        zorder=4,
        clip_on=True,
    )

    rng = ymax - ymin
    if lower < ymin:
        ax.annotate(
            "",
            xy=(xpos, ymin),
            xytext=(xpos, ymin + rng * 0.10),
            arrowprops=dict(arrowstyle="-|>", color=color, lw=1.1, mutation_scale=7),
            zorder=5,
        )
    if upper > ymax:
        ax.annotate(
            "",
            xy=(xpos, ymax),
            xytext=(xpos, ymax - rng * 0.10),
            arrowprops=dict(arrowstyle="-|>", color=color, lw=1.1, mutation_scale=7),
            zorder=5,
        )


# ---------------------------------------------------------------------------
# Figure
# Top-level GridSpec: [left_group | spacer | right_group]
# Left  group — height ratios 3:3:2 (NDD, Dementia, HD)
# Right group — height ratios 2:2:2 (AD, PD, ALS)
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(10, 10))

gs_top = GridSpec(
    1,
    3,
    figure=fig,
    left=0.03,
    right=0.96,
    top=0.96,
    bottom=0.07,
    width_ratios=[5.0, 0.6, 5.0],
    wspace=0.0,
)

gs_left = GridSpecFromSubplotSpec(
    3,
    N_COLS,
    subplot_spec=gs_top[0, 0],
    hspace=0.25,
    wspace=0.08,
    height_ratios=[1, 1, 1],
)

gs_right = GridSpecFromSubplotSpec(
    3,
    N_COLS,
    subplot_spec=gs_top[0, 2],
    hspace=0.25,
    wspace=0.08,
    height_ratios=[2, 2, 2],
)

left_axes = [[fig.add_subplot(gs_left[r, c]) for c in range(N_COLS)] for r in range(3)]
right_axes = [
    [fig.add_subplot(gs_right[r, c]) for c in range(N_COLS)] for r in range(3)
]


# ---------------------------------------------------------------------------
# Fill panels
# ---------------------------------------------------------------------------
N_ROWS = 3


def fill_group(axes_grid, diseases, show_ylabel=True):
    for row_idx, (icd_code, label) in enumerate(diseases):
        sub = df[df["ICD_Code"] == icd_code].copy()
        ymin, ymax = DISEASE_YLIM[icd_code]

        for col_idx, q in enumerate(QUINTILES):
            ax = axes_grid[row_idx][col_idx]
            style_ax(ax)
            draw_ref_line(ax)
            draw_x_separators(ax)
            ax.set_xlim(*X_LIM)
            ax.set_ylim(ymin, ymax)
            apply_yticks(ax, icd_code, col_idx, show_ylabel=show_ylabel)
            ax.set_xticks(list(LAG_X.values()))
            ax.set_xticklabels([])

            # Disease name — just above the top-left of Q2, outside the subplot
            if col_idx == 0:
                ax.text(
                    0.0,
                    1.02,
                    label,
                    transform=ax.transAxes,
                    ha="left",
                    va="bottom",
                    fontsize=10,
                    fontweight="bold",
                    clip_on=False,
                )

            # Quintile label — bottom of last row only
            if row_idx == N_ROWS - 1:
                ax.set_xlabel(q, fontsize=11, fontweight="bold", labelpad=4)

            q_sub = sub[sub["Quintile"] == q] if not sub.empty else pd.DataFrame()

            for lag in LAGS:
                xpos = LAG_X[lag]
                color = LAG_COLORS[lag]
                row_data = (
                    q_sub[q_sub["Lag"] == lag] if not q_sub.empty else pd.DataFrame()
                )

                if row_data.empty:
                    ax.plot(
                        xpos,
                        (ymin + ymax) / 2,
                        "x",
                        color="#CCCCCC",
                        markersize=4,
                        zorder=3,
                    )
                    continue

                r_ = row_data.iloc[0]
                mean = float(r_["MRR_mean"])
                lower = float(r_["MRR_lower"])
                upper = float(r_["MRR_upper"])

                p_raw = str(r_.get("p", "1"))
                try:
                    sig = float(p_raw) < 0.05 if p_raw != "p<0.0001" else True
                except ValueError:
                    sig = p_raw == "p<0.0001"

                plot_point(ax, xpos, mean, lower, upper, color, sig, ymin, ymax)


fill_group(left_axes, LEFT_DISEASES, show_ylabel=True)
fill_group(right_axes, RIGHT_DISEASES, show_ylabel=False)

# ---------------------------------------------------------------------------
# Legend
# ---------------------------------------------------------------------------
lag_handles = [
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="none",
        color=LAG_COLORS[lag],
        markerfacecolor=LAG_COLORS[lag],
        markersize=6.5,
        label=LAG_LABELS[lag],
    )
    for lag in LAGS
]
sig_handles = [
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="none",
        color="#555555",
        markerfacecolor="#555555",
        markersize=6.5,
        label="p < 0.05",
    ),
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="none",
        color="#555555",
        markerfacecolor="white",
        markeredgecolor="#555555",
        markeredgewidth=1.5,
        markersize=6.5,
        label="Non-significant",
    ),
]

fig.legend(
    handles=lag_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.0),
    ncol=6,
    fontsize=14,
    frameon=True,
    framealpha=0.97,
    edgecolor="#DDDDDD",
    handletextpad=0.5,
    columnspacing=1.3,
)

plt.savefig(OUT_DIR / "MRR_Forest.png", dpi=DPI, bbox_inches="tight", facecolor="white")
plt.savefig(OUT_DIR / "MRR_Forest.pdf", bbox_inches="tight", facecolor="white")
plt.close()
print("Saved → MRR_Forest.png / .pdf")
