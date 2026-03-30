"""
Forest plot — MRR with MRD background shading
Same layout as Visualization_MRR.py, plus:
  - For each lag position in every panel, a background colour block
    coloured by the MRD value using the same blue-gradient colormap
    as Visualization_Heatmap_NDD_Helvetica.py (vmin=0, vmax=30).
Source: Result/Tables/main.csv  (MRR + MRD merged)
Output: Result/MRR_Forest/MRR_MRD_Forest.png / .pdf
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
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
MRD_ALPHA = 0.50  # fixed alpha for MRD background blocks

# MRD colormap — same blue gradient as Visualization_Heatmap_NDD_Helvetica.py
_MRD_COLORS = [
    "#F5F8FA",
    "#BFDAE9",
    "#89BED9",
    "#408FBF",
    "#2165AA",
    "#073162",
    "#073162",
]
_MRD_POSITIONS = [0.0, 0.067, 0.133, 0.200, 0.3, 0.5, 1.0]
MRD_CMAP = LinearSegmentedColormap.from_list(
    "mrd_cmap", list(zip(_MRD_POSITIONS, _MRD_COLORS)), N=256
)
MRD_NORM = plt.Normalize(vmin=0, vmax=30)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path(__file__).resolve().parents[2]
MAIN_CSV = BASE / "Result" / "Tables" / "main.csv"
OUT_DIR = BASE / "Result" / "MRR_Forest"
OUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(MAIN_CSV)

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
N_ROWS = 3

LAGS = [5, 10, 15]

# LAG_COLORS = {5: "#3F007D", 10: "#756BB1", 15: "#BCBDDC"}
LAG_COLORS = {
    5: "#1D3557",  # deep muted blue
    10: "#2A9D8F",  # desaturated teal
    15: "#6D597A",  # muted purple-grey
}
LAG_LABELS = {5: "5-year lag", 10: "10-year lag", 15: "15-year lag"}

LAG_X = {5: 0.0, 10: 1.0, 15: 2.0}
X_LIM = (-0.6, 2.6)

DISEASE_YLIM = {
    "G20_G30_G12.2_F01_F03": (0.8, 1.8),
    "G30_F01_F03": (0.8, 1.8),
    "G10": (0.0, 6.0),
    "G30": (0.8, 1.8),
    "G20": (0.5, 3.0),
    "G12.2": (0.0, 6.0),
}

DISEASE_YTICKS = {
    "G20_G30_G12.2_F01_F03": [1, 1.2, 1.4, 1.6],
    "G30_F01_F03": [1, 1.2, 1.4, 1.6],
    "G10": [1, 2, 3, 4, 5],
    "G30": [1, 1.2, 1.4, 1.6],
    "G20": [1, 1.5, 2.0, 2.5],
    "G12.2": [1, 2, 3, 4, 5],
}

# Scalar mappable for the colorbar legend
_mrd_sm = plt.cm.ScalarMappable(cmap=MRD_CMAP, norm=MRD_NORM)
_mrd_sm.set_array([])


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def style_ax(ax):
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_facecolor("white")
    fig = ax.get_figure()
    fw, fh = fig.get_size_inches()
    pos = ax.get_position()
    ax_w = pos.width * fw
    ax_h = pos.height * fh
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
    ticks = DISEASE_YTICKS[icd_code]
    ax.set_yticks(ticks)
    if col_idx == 0:
        labels = []
        for t in ticks:
            labels.append(
                str(int(round(t))) if abs(t - round(t)) < 1e-9 else f"{t:.1f}"
            )
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


def draw_mrd_cap(ax, q_sub, ymin, ymax):
    """
    Draw three coloured bands at the bottom of the subplot, covering half
    of the y < 1 range.  Each band corresponds to one lag (left→right:
    5, 10, 15-year lag) and is filled with the MRD colour for that lag.
    The bands are clipped to the subplot's rounded corners and separated
    by gray dashed vertical lines.  MRR content is drawn on top (higher z).
    """
    from matplotlib.patches import Rectangle as _Rect

    y_span = ymax - ymin
    # half the <1 area expressed in axes fraction
    cap_h_ax = (1.0 - ymin) * 0.5 / y_span
    section_w = 1.0 / 3.0

    fig = ax.get_figure()
    fw, fh = fig.get_size_inches()
    pos = ax.get_position()
    ax_w = pos.width * fw
    ax_h = pos.height * fh
    mut_asp = ax_w / ax_h if ax_h > 0 else 1.0

    # Use the full subplot shape as clip path so fills respect the rounded
    # bottom corners of the subplot frame (same ROUND_SIZE / mutation_aspect).
    clip_frame = FancyBboxPatch(
        (0.0, 0.0),
        1.0,
        1.0,
        boxstyle=f"round,pad=0,rounding_size={ROUND_SIZE}",
        mutation_aspect=mut_asp,
        transform=ax.transAxes,
        facecolor="none",
        edgecolor="none",
        linewidth=0,
        zorder=0,
        clip_on=False,
    )
    ax.add_patch(clip_frame)

    # Coloured fill for each lag section
    for i, lag in enumerate(LAGS):
        row = q_sub[q_sub["Lag"] == lag]
        if row.empty or "MRD_mean" not in row.columns:
            continue
        mrd = row.iloc[0]["MRD_mean"]
        if pd.isna(mrd):
            continue
        rgba = MRD_CMAP(MRD_NORM(float(mrd)))
        fill = _Rect(
            (i * section_w, 0.0),
            section_w,
            cap_h_ax,
            transform=ax.transAxes,
            facecolor=rgba[:3],
            alpha=MRD_ALPHA,
            linewidth=0,
            zorder=1,
        )
        fill.set_clip_path(clip_frame)
        fill.set_clip_on(True)
        ax.add_patch(fill)

    # Dashed gray vertical separators between the three sections
    for sep_x in [section_w, 2 * section_w]:
        ax.plot(
            [sep_x, sep_x],
            [0.0, cap_h_ax],
            transform=ax.transAxes,
            color=FRAME_COLOR,
            linewidth=0.7,
            linestyle="--",
            zorder=2,
            clip_on=True,
        )


def plot_point(ax, xpos, mean, lower, upper, color, sig, ymin, ymax):
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
        markersize=3.5,
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
def fill_group(axes_grid, diseases, show_ylabel=True):
    for row_idx, (icd_code, label) in enumerate(diseases):
        sub = df[df["ICD_Code"] == icd_code].copy()
        ymin, ymax = DISEASE_YLIM[icd_code]

        for col_idx, q in enumerate(QUINTILES):
            ax = axes_grid[row_idx][col_idx]
            style_ax(ax)
            ax.set_xlim(*X_LIM)
            ax.set_ylim(ymin, ymax)

            q_sub = sub[sub["Quintile"] == q] if not sub.empty else pd.DataFrame()

            # MRD bands at the bottom (drawn before MRR content)
            if not q_sub.empty:
                draw_mrd_cap(ax, q_sub, ymin, ymax)

            draw_ref_line(ax)
            draw_x_separators(ax)
            apply_yticks(ax, icd_code, col_idx, show_ylabel=show_ylabel)
            ax.set_xticks(list(LAG_X.values()))
            ax.set_xticklabels([])

            # Disease name above subplot
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

            # Quintile label at bottom of last row
            if row_idx == N_ROWS - 1:
                ax.set_xlabel(q, fontsize=11, fontweight="bold", labelpad=4)

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

                p_col = "p_MRR" if "p_MRR" in r_.index else "p"
                p_raw = str(r_.get(p_col, "1"))
                try:
                    sig = float(p_raw) < 0.05 if p_raw != "p<0.0001" else True
                except ValueError:
                    sig = p_raw == "p<0.0001"

                plot_point(ax, xpos, mean, lower, upper, color, sig, ymin, ymax)


fill_group(left_axes, LEFT_DISEASES, show_ylabel=True)
fill_group(right_axes, RIGHT_DISEASES, show_ylabel=False)

# ---------------------------------------------------------------------------
# Legend — horizontal line + centre dot (parallel line style)
# ---------------------------------------------------------------------------
lag_handles = [
    Line2D(
        [0],
        [0],
        color=LAG_COLORS[lag],
        linewidth=1.2,
        marker="o",
        markerfacecolor=LAG_COLORS[lag],
        markeredgecolor=LAG_COLORS[lag],
        markeredgewidth=1.2,
        markersize=3.5,
        label=LAG_LABELS[lag],
    )
    for lag in LAGS
]

fig.legend(
    handles=lag_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.0),
    ncol=3,
    fontsize=11,
    frameon=True,
    framealpha=0.97,
    edgecolor="#DDDDDD",
    handletextpad=0.5,
    columnspacing=1.0,
)

# MRD colorbar below the lag legend
cbar_ax = fig.add_axes([0.15, -0.04, 0.70, 0.018])
cbar = fig.colorbar(_mrd_sm, cax=cbar_ax, orientation="horizontal")
cbar.set_label("Point Estimate of Mortality Rate Difference", fontsize=12)
cbar.ax.tick_params(labelsize=9)

plt.savefig(
    OUT_DIR / "MRR_MRD_Forest_V.png", dpi=DPI, bbox_inches="tight", facecolor="white"
)
plt.savefig(OUT_DIR / "MRR_MRD_Forest_V.pdf", bbox_inches="tight", facecolor="white")
plt.close()
print("Saved → MRR_MRD_Forest_V.png / .pdf")
