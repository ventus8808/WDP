"""
Bubble (dot-heatmap) plots of the EQI-Air × SVI joint sensitivity results.

Reads the +SM (smoking-adjusted) model row from each disease CSV in
    Result/brms_SVI_Air_22_Sensitivity_Combination/   (2x2 joint)
    Result/brms_SVI_Air_33_Sensitivity_Combination/   (3x3 joint)
and draws, per disease, a square matrix of joint cells a{air}s{svi}:
x = EQI-Air level, y = SVI level. Colour & size encode the mortality rate
difference (MRD, deaths/100k) vs the a1s1 reference; a solid ring marks
posterior p < 0.05.

Output: Result/SVI_Bubble_Visualization/{22,33}_{disease}.png
"""

import re
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

try:
    import seaborn as sns
    CMAP = sns.color_palette("vlag", as_cmap=True)
except Exception:                                    # pragma: no cover
    CMAP = plt.get_cmap("RdBu_r")

RESULT_DIRS = {
    "22": Path("Result/brms_SVI_Air_22_Sensitivity_Combination"),
    "33": Path("Result/brms_SVI_Air_33_Sensitivity_Combination"),
}
OUT_DIR = Path("Result/SVI_Bubble_Visualization")
MODEL = "AirSVI+SM"
CELL_RE = re.compile(r"^a(\d+)s(\d+)$")
EST_RE = re.compile(r"^\s*(-?[\d.]+)\s*\(\s*(-?[\d.]+)\s*,\s*(-?[\d.]+)\s*\)")

DISEASE_FULL = {
    "Cancer": "Cancer", "CVD": "Cardiovascular disease",
    "CRD": "Chronic respiratory disease", "CKD": "Chronic kidney disease",
    "CLD": "Chronic liver disease", "NDD": "Neurodegenerative disease",
    "Suicide": "Suicide",
}

AREA_MIN, AREA_SPAN = 260.0, 3200.0   # bubble area: ref -> AREA_MIN, max|MRD| -> +SPAN

plt.rcParams.update({
    "font.family": "Georgia",
    "font.size": 13,
    "axes.edgecolor": "#888888",
    "axes.linewidth": 0.8,
    "figure.dpi": 150,
})


def _parse_est(cell):
    if not isinstance(cell, str):
        return (float(cell), np.nan, np.nan) if pd.notna(cell) else (np.nan,) * 3
    m = EST_RE.match(cell)
    if m:
        return float(m.group(1)), float(m.group(2)), float(m.group(3))
    try:
        return float(cell), np.nan, np.nan
    except ValueError:
        return np.nan, np.nan, np.nan


def _stars(p):
    if pd.isna(p):
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def _text_color(rgba):
    lum = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
    return "white" if lum < 0.55 else "#1a1a1a"


def make_bubble(csv_path, grid, disease):
    df = pd.read_csv(csv_path)
    row = df[df["Model"] == MODEL]
    if row.empty:
        print(f"  [skip] {csv_path.name}: no {MODEL} row")
        return
    row = row.iloc[0]

    cells = []
    for col in df.columns:
        m = CELL_RE.match(col)
        if not m:
            continue
        a, s = int(m.group(1)), int(m.group(2))
        mean, lo, hi = _parse_est(row[col])
        p = row.get(f"{col}_p", np.nan)
        cells.append(dict(a=a, s=s, mean=mean, lo=lo, hi=hi,
                          p=float(p) if pd.notna(p) else np.nan,
                          is_ref=(a == 1 and s == 1)))
    cells = pd.DataFrame(cells)
    if cells.empty:
        print(f"  [skip] {csv_path.name}: no joint cells")
        return
    n_a, n_s = int(cells["a"].max()), int(cells["s"].max())

    maxabs = float(np.nanmax(np.abs(cells["mean"]))) or 1.0
    vlim = np.ceil(maxabs / 5) * 5                     # rounded symmetric limit
    norm = mcolors.Normalize(-vlim, vlim)
    cells["area"] = AREA_MIN + (np.abs(cells["mean"]) / maxabs) * AREA_SPAN

    fig, ax = plt.subplots(figsize=(3.4 + 1.55 * n_a, 2.9 + 1.55 * n_s),
                           constrained_layout=True)

    # faint cell separators for a matrix feel
    for x in np.arange(0.5, n_a + 1):
        ax.axvline(x, color="#e6e6e6", lw=0.8, zorder=0)
    for y in np.arange(0.5, n_s + 1):
        ax.axhline(y, color="#e6e6e6", lw=0.8, zorder=0)

    sig = (cells["p"] < 0.05).fillna(False)
    sc = ax.scatter(
        cells["a"], cells["s"], s=cells["area"], c=cells["mean"],
        cmap=CMAP, norm=norm,
        edgecolors=np.where(sig, "#222222", "#bdbdbd"),
        linewidths=np.where(sig, 1.6, 0.8),
        zorder=3,
    )

    for _, c in cells.iterrows():
        if c["is_ref"]:
            txt, tcol = "ref", "#444444"
        else:
            tcol = _text_color(CMAP(norm(c["mean"])))
            ci = "" if np.isnan(c["lo"]) else f"\n({c['lo']:.0f}, {c['hi']:.0f})"
            txt = f"{c['mean']:.1f}{_stars(c['p'])}{ci}"
        ax.annotate(txt, (c["a"], c["s"]), ha="center", va="center",
                    fontsize=9.5, color=tcol, zorder=4, linespacing=1.25)

    ax.set_xticks(range(1, n_a + 1))
    ax.set_yticks(range(1, n_s + 1))
    ax.set_xticklabels([f"Air\nQ{i}" for i in range(1, n_a + 1)])
    ax.set_yticklabels([f"SVI Q{i}" for i in range(1, n_s + 1)])
    ax.set_xlabel("EQI-Air  (low $\\rightarrow$ high air pollution)", labelpad=8)
    ax.set_ylabel("SVI  (low $\\rightarrow$ high vulnerability)", labelpad=8)
    ax.set_xlim(0.5, n_a + 0.5)
    ax.set_ylim(0.5, n_s + 0.5)
    ax.set_aspect("equal")
    ax.tick_params(length=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    fig.suptitle(f"{DISEASE_FULL.get(disease, disease)}",
                 fontsize=17, fontweight="bold", y=1.06)
    ax.set_title(f"Joint EQI-Air × SVI ({n_a}×{n_s}) mortality rate difference "
                 f"vs a1s1, adjusted for smoking", fontsize=11, color="#555555")

    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.03, shrink=0.85)
    cbar.set_label("MRD vs a1s1\n(deaths per 100k)", fontsize=10)
    cbar.outline.set_edgecolor("#bbbbbb")
    cbar.ax.tick_params(labelsize=9)

    # size legend (representative |MRD|)
    ticks = [t for t in (5, 10, 20, 40) if t <= maxabs * 1.05][:3]
    if ticks:
        handles = [Line2D([0], [0], marker="o", linestyle="",
                          markerfacecolor="#d2d2d2", markeredgecolor="#888888",
                          markersize=np.sqrt(AREA_MIN + (t / maxabs) * AREA_SPAN) * 0.5,
                          label=f"{t}") for t in ticks]
        leg = ax.legend(handles=handles, title="bubble size = |MRD| (deaths/100k)",
                        loc="upper center", bbox_to_anchor=(0.5, -0.14),
                        frameon=False, ncol=len(ticks), columnspacing=3.0,
                        handletextpad=1.2, fontsize=9, title_fontsize=10)
        ax.add_artist(leg)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f"{grid}_{disease}.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  saved {out}")


def main():
    for grid, d in RESULT_DIRS.items():
        if not d.exists():
            print(f"[skip] {d} not found")
            continue
        print(f"=== {grid} joint ({d}) ===")
        for csv_path in sorted(d.glob("*.csv")):
            make_bubble(csv_path, grid, csv_path.name.split("_")[0])


if __name__ == "__main__":
    main()
