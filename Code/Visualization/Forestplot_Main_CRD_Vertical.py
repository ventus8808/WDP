#!/usr/bin/env python3
"""CRD LagRef forest plot with vertical credible-interval measurements.

The layout follows the panel style in Visualization.py: measurements are on the
vertical axis, categories run across the horizontal axis, and colours identify
the four EQI quintiles.  It produces one figure only: COPD, ILD, and lung
cancer form the columns; MRD and MRR form the rows.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Arial", "DejaVu Sans"],
        "axes.titleweight": "semibold",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXTRACTION_DIR = PROJECT_ROOT / "Result" / "Extraction"
MRR_DIR = PROJECT_ROOT / "Result" / "brms_Main_MRR"
OUTPUT_DIR = PROJECT_ROOT / "Result" / "Forest_Main"
OUTPUT_FILE = OUTPUT_DIR / "CRD_LagRef_SubtypeGrid_Vertical.png"

GROUP = "CRD"
MRR_TYPE = "LagRef"
SUBTYPES = {
    "C34": "Lung cancer",
    "J43_J44": "COPD",
    "J84_D86": "ILD",
}
LAGS = [5, 10, 15, 20]
QUINTILES = ["Q2", "Q3", "Q4", "Q5"]

# Okabe–Ito colour-blind-safe palette, commonly used for scientific figures.
QUINTILE_COLOURS = {
    "Q2": "#0072B2",  # blue
    "Q3": "#009E73",  # bluish green
    "Q4": "#D55E00",  # vermillion
    "Q5": "#CC79A7",  # reddish purple
}

EFFECT_RE = re.compile(
    r"^\s*([+-]?\d+(?:\.\d+)?)\s*\(\s*([+-]?\d+(?:\.\d+)?)\s*,\s*([+-]?\d+(?:\.\d+)?)\s*\)"
)


def parse_effect(value: object) -> tuple[float, float, float] | None:
    """Return estimate, lower CI, and upper CI from an extracted result cell."""
    if pd.isna(value):
        return None
    match = EFFECT_RE.match(str(value))
    if not match:
        return None
    return tuple(float(v) for v in match.groups())


def load_mrd() -> pd.DataFrame:
    path = EXTRACTION_DIR / f"{GROUP}_Main.csv"
    if not path.exists():
        raise FileNotFoundError(f"MRD input not found: {path}")
    df = pd.read_csv(path)
    return df[
        (df["Model"] == "EQI")
        & (df["Lag"].isin(LAGS))
        & (df["ICD_Code"].isin(SUBTYPES))
    ].copy()


def load_mrr(icd: str) -> pd.DataFrame:
    """Load the per-ICD LagRef result, with CRD-level fallback."""
    for stem in (icd, GROUP):
        path = MRR_DIR / f"{stem}_MRR_{MRR_TYPE}.csv"
        if path.exists():
            df = pd.read_csv(path)
            subset = df[
                (df["ICD_Code"] == icd)
                & (df["EQI_Period"] == "2000_2005")
                & (df["Lag"].isin(LAGS))
            ].copy()
            if not subset.empty:
                return subset
    raise FileNotFoundError(f"MRR LagRef input not found for {icd}")


def collect_measurements(df: pd.DataFrame) -> dict[tuple[int, str], tuple[float, float, float]]:
    """Shape a result frame into {(lag, quintile): (estimate, lower, upper)}."""
    values = {}
    for lag in LAGS:
        lag_rows = df[df["Lag"] == lag]
        if lag_rows.empty:
            continue
        row = lag_rows.iloc[0]
        for quintile in QUINTILES:
            effect = parse_effect(row[quintile])
            if effect is not None:
                values[(lag, quintile)] = effect
    return values


def shared_limits(
    datasets: list[dict[tuple[int, str], tuple[float, float, float]]], reference: float
) -> tuple[float, float]:
    values = [reference]
    for dataset in datasets:
        for _, lower, upper in dataset.values():
            values.extend((lower, upper))
    lo, hi = min(values), max(values)
    pad = max((hi - lo) * 0.10, 0.15 if reference == 1 else 0.75)
    return lo - pad, hi + pad


def draw_panel(
    ax: plt.Axes,
    measurements: dict[tuple[int, str], tuple[float, float, float]],
    *,
    measure: str,
    ylim: tuple[float, float],
) -> None:
    """Draw one subtype/measure panel with vertical credible intervals."""
    offsets = np.linspace(-0.24, 0.24, len(QUINTILES))
    x_positions = np.arange(len(LAGS))

    for q_index, quintile in enumerate(QUINTILES):
        for x, lag in zip(x_positions, LAGS):
            effect = measurements.get((lag, quintile))
            if effect is None:
                continue
            estimate, lower, upper = effect
            ax.errorbar(
                x + offsets[q_index],
                estimate,
                yerr=[[estimate - lower], [upper - estimate]],
                fmt="o",
                color=QUINTILE_COLOURS[quintile],
                ecolor=QUINTILE_COLOURS[quintile],
                markersize=3.8,
                capsize=3,
                capthick=1.25,
                elinewidth=1.35,
                markeredgecolor="white",
                markeredgewidth=0.5,
                zorder=3,
            )

    ax.set_ylim(ylim)
    ax.set_xlim(-0.55, len(LAGS) - 0.45)
    ax.set_xticks(x_positions, [f"{lag}-year" for lag in LAGS])
    if measure == "MRR":
        ax.set_yscale("log")
        ax.set_yticks([0.9, 1, 1.25, 1.5, 2, 2.5])
        ax.set_yticklabels(["0.9", "1.0", "1.25", "1.5", "2.0", "2.5"])
    else:
        ax.set_yticks([0, 3, 6, 9, 12, 15])
    ax.grid(axis="y", color="#E2E2E2", linewidth=0.65, zorder=0)
    ax.set_axisbelow(True)
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#555555")
        spine.set_linewidth(0.75)
    ax.tick_params(axis="both", colors="#333333", labelsize=9, length=3.5, width=0.75, pad=5)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    mrd = load_mrd()

    mrd_data = {}
    mrr_data = {}
    for icd in SUBTYPES:
        mrd_data[icd] = collect_measurements(mrd[mrd["ICD_Code"] == icd])
        mrr_data[icd] = collect_measurements(load_mrr(icd))

    # Current CRD values span −3.06 to 12.61 for MRD CrIs and 0.94 to 2.40 for
    # MRR CrIs. These bounds retain all positive values, provide headroom, and
    # make the common axes easy to read across the three subtypes.
    mrd_limits = (0, 15)
    mrr_limits = (0.9, 2.5)

    # A shared y-axis within each row makes estimates directly comparable across
    # all three CRD subtypes.
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.6), sharex=True, sharey="row")
    fig.patch.set_facecolor("white")

    for col, (icd, subtype) in enumerate(SUBTYPES.items()):
        axes[0, col].set_title(subtype, fontsize=12, fontweight="bold", color="#111111", pad=9)
        draw_panel(axes[0, col], mrd_data[icd], measure="MRD", ylim=mrd_limits)
        draw_panel(axes[1, col], mrr_data[icd], measure="MRR", ylim=mrr_limits)

    for ax in axes[0]:
        ax.tick_params(labelbottom=False)
    for row in range(2):
        for col in range(1, 3):
            axes[row, col].tick_params(labelleft=False, left=False)
    axes[0, 0].set_ylabel("MRD\n95% CrI", fontsize=11, fontweight="semibold", color="#3A3A3C", labelpad=16)
    axes[1, 0].set_ylabel("MRR\n95% CrI", fontsize=11, fontweight="semibold", color="#3A3A3C", labelpad=16)

    for letter, ax in zip("ABCDEF", axes.flat):
        ax.text(0.02, 0.95, letter, transform=ax.transAxes, fontsize=10,
                fontweight="bold", va="top", ha="left")

    handles = [
        Line2D([0], [0], marker="o", color=QUINTILE_COLOURS[q], label=q,
               markersize=4.6, linewidth=1.4, markeredgecolor="white")
        for q in QUINTILES
    ]
    fig.legend(
        handles=handles,
        title="EQI quintile",
        loc="lower center",
        ncol=4,
        frameon=True,
        fancybox=True,
        bbox_to_anchor=(0.5, 0.02),
        fontsize=10.5,
        title_fontsize=10.5,
    )
    legend = fig.legends[0]
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_edgecolor("#444444")
    legend.get_frame().set_linewidth(0.7)
    fig.tight_layout(rect=(0.06, 0.14, 0.99, 0.98), h_pad=2.2, w_pad=2.2)
    fig.savefig(OUTPUT_FILE, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved -> {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
