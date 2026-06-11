"""
Composite figure of the SVI trajectory analyses (reads raw results only; no
model is re-fit here).

Inputs
    Data/Processed/SVI/SVI_{LCGA,GBTM}.csv          (individual county values)
    Result/CDC_SVI_Trajectory/{LCGA,GBTM}_trajectories.csv  (class means)
    Result/CDC_SVI_Trajectory/{LCGA,GBTM}_bic.csv           (BIC-vs-K)
    Result/CDC_SVI_Trajectory/LTA_{states,composition,transition}.csv

Layout (2 x 3)
    (a) LCGA mean + IQR band      (b) GBTM mean + IQR band      (c) LTA composition
    (d) Model selection (BIC)     (e) Net change by class       (f) LTA transition

Design follows Howard (2021): mean trajectories over semi-transparent bands,
no occluding spaghetti. The two mixture models (LCGA == GBTM) share a single
BIC panel; panel (e) replaces the redundant second elbow with the within-class
distribution of 2000->2022 change.

Output: Result/CDC_SVI_Trajectory/Trajectory.png
"""

from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RES = Path("Result/CDC_SVI_Trajectory")
LABELS = Path("Data/Processed/SVI")
OUT = RES / "Trajectory.png"

# Diverging tier palette, ordered A->D (lowest->highest SVI):
# A = green (least vulnerable / best), D = orange (most vulnerable / worst);
# muted tones for an academic look.
COLORS = ["#3F7E54", "#9CC196", "#E6C58A", "#D58A4E",
          "#699DCB", "#9E79B8"]

plt.rcParams["font.family"] = "Georgia"
plt.rcParams["font.size"] = 14


def _years(cols):
    svi = [c for c in cols if c.endswith("_SVI_RPL")]
    return svi, np.array([int(c.split("_")[0]) for c in svi])


def mixture_panel(ax, method, title):
    """Model mean trajectory per class over a semi-transparent IQR band."""
    lab = pd.read_csv(LABELS / f"SVI_{method}.csv", dtype={"COUNTY_FIPS": str})
    traj = pd.read_csv(RES / f"{method}_trajectories.csv")
    svi_cols, years = _years(lab.columns)
    classes = sorted(traj["Class"].unique())
    cmap = {c: COLORS[i % len(COLORS)] for i, c in enumerate(classes)}
    cls_col = f"{method}_Class"

    # IQR bands first (background)
    for c in classes:
        Yc = lab.loc[lab[cls_col] == c, svi_cols].to_numpy(float)
        p25 = np.nanpercentile(Yc, 25, axis=0)
        p75 = np.nanpercentile(Yc, 75, axis=0)
        ax.fill_between(years, p25, p75, color=cmap[c], alpha=0.22, lw=0)
    # Model mean lines on top, with the class share labelled at the line end
    for c in classes:
        t = traj[traj["Class"] == c].sort_values("Year")
        pct = t["Pct"].iloc[0]
        ax.plot(t["Year"], t["Mean"], color=cmap[c], lw=3.0, marker="o", ms=7,
                markeredgecolor="white", markeredgewidth=1.1, zorder=5)
        ax.annotate(f"{c.split()[-1]}: {pct:.0f}%",
                    xy=(years[-1], t["Mean"].iloc[-1]),
                    xytext=(-3, -11), textcoords="offset points", va="top",
                    ha="right", fontsize=10, fontweight="bold", color=cmap[c])
    ax.set_xlabel("Year")
    ax.set_ylabel("SVI percentile")
    ax.set_title(title)
    ax.set_xticks(years)
    ax.tick_params(axis="x", rotation=45)
    ax.set_xlim(years.min(), years.max())
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.25)


def bic_panel(ax, title):
    """Single BIC-vs-K panel overlaying LCGA and GBTM (they coincide)."""
    l = pd.read_csv(RES / "LCGA_bic.csv")
    g = pd.read_csv(RES / "GBTM_bic.csv")
    ax.plot(l["K"], l["BIC"], "-o", color="#2F7F4F", ms=7, label="LCGA")
    ax.plot(g["K"], g["BIC"], "--s", color="#E68785", ms=6, mfc="none",
            label="GBTM (identical model)")
    ksel = int(g.loc[g["Selected"] == 1, "K"].iloc[0])
    ax.axvline(ksel, color="#888888", ls=":", lw=1.8)
    ax.annotate(f"elbow K = {ksel}", xy=(ksel, l["BIC"].min()),
                xytext=(ksel + 0.4, l["BIC"].min()), fontsize=11, va="bottom")
    ax.set_title(title)
    ax.set_xlabel("Number of classes (K)")
    ax.set_ylabel("BIC")
    ax.set_xticks(l["K"])
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right")


def net_change_panel(ax, method, title):
    """Within-class distribution of last-minus-first SVI."""
    lab = pd.read_csv(LABELS / f"SVI_{method}.csv", dtype={"COUNTY_FIPS": str})
    svi_cols, years = _years(lab.columns)
    first, last = svi_cols[0], svi_cols[-1]
    delta = lab[last] - lab[first]
    cls_col = f"{method}_Class"
    classes = sorted(lab[cls_col].dropna().unique())
    cmap = {c: COLORS[i % len(COLORS)] for i, c in enumerate(classes)}
    data = [delta[(lab[cls_col] == c)].dropna().to_numpy() for c in classes]

    parts = ax.violinplot(data, showmedians=True, showextrema=False, widths=0.85)
    for body, c in zip(parts["bodies"], classes):
        body.set_facecolor(cmap[c]); body.set_alpha(0.8)
        body.set_edgecolor("#444444"); body.set_linewidth(0.8)
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.5)
    ax.axhline(0, color="#555555", ls="--", lw=1.2)
    ax.set_xticks(range(1, len(classes) + 1))
    ax.set_xticklabels(classes)
    ax.set_xlabel("Class")
    ax.set_ylabel(f"SVI change, {years[-1]} minus {years[0]}")
    ax.set_title(title)
    ax.grid(True, axis="y", alpha=0.25)


def lta_composition_panel(ax, title):
    comp = pd.read_csv(RES / "LTA_composition.csv")
    states = pd.read_csv(RES / "LTA_states.csv")
    lta = pd.read_csv(LABELS / "SVI_LTA.csv", dtype={"COUNTY_FIPS": str})
    counts = lta["LTA_State_2022"].value_counts()
    letters = list(states["State"])
    years = comp["Year"].to_numpy()
    shares = comp[letters].to_numpy().T
    colors = [COLORS[i % len(COLORS)] for i in range(len(letters))]
    # lighter fill, matching the IQR-band shading in panels (a)/(b)
    ax.stackplot(years, shares, colors=colors, alpha=0.22)
    # label each band with: letter (n counties, mean percentile) at 2022
    means = dict(zip(states["State"], states["Mean"]))
    last = shares[:, -1]
    mids = np.cumsum(last) - last / 2
    for s, ymid in zip(letters, mids):
        ax.annotate(f"Category {s} ({int(counts.get(s, 0))} counties)\n"
                    f"Mean percentile: {means[s] * 100:.2f}",
                    xy=(years[-1], ymid), xytext=(-8, 0),
                    textcoords="offset points", ha="right", va="center",
                    fontsize=12, color="black")
    ax.set_xlim(years.min(), years.max())
    ax.set_ylim(0, 1)
    ax.set_xticks(years)
    ax.tick_params(axis="x", rotation=45)
    ax.set_xlabel("Year")
    ax.set_ylabel("Share of counties")
    ax.set_title(title)


def lta_transition_panel(ax, fig, title):
    t = pd.read_csv(RES / "LTA_transition.csv")
    letters = list(t["From"])
    M = t[letters].to_numpy()
    S = len(letters)
    im = ax.imshow(M, cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(S)); ax.set_yticks(range(S))
    ax.set_xticklabels(letters); ax.set_yticklabels(letters)
    ax.set_xlabel("State in 2022")
    ax.set_ylabel("State in 2000")
    ax.set_title(title)
    for i in range(S):
        for j in range(S):
            ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center",
                    color="white" if M[i, j] > 0.5 else "black", fontsize=12)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


def main():
    fig, axes = plt.subplots(2, 3, figsize=(21, 13))
    mixture_panel(axes[0, 0], "LCGA", "(a) LCGA — mean trajectory + IQR band")
    mixture_panel(axes[0, 1], "GBTM", "(b) GBTM — mean trajectory + IQR band")
    lta_composition_panel(axes[0, 2], "(c) LTA — state composition over time")
    bic_panel(axes[1, 0], "(d) Model selection — BIC vs K")
    net_change_panel(axes[1, 1], "GBTM", "(e) Within-class SVI change (2000 to 2022)")
    lta_transition_panel(axes[1, 2], fig, "(f) LTA — transition 2000 to 2022")

    fig.suptitle("County SVI Trajectory Analysis — LCGA / GBTM / LTA (2000–2022)",
                 fontsize=22, y=0.995)

    # One shared legend: the A-D colour means the same vulnerability tier in
    # every panel (Class for LCGA/GBTM/net-change, State for LTA).
    tier_labels = ["A", "B", "C", "D — highest SVI (most vulnerable)"]
    handles = [mpatches.Patch(facecolor=COLORS[i], edgecolor="none",
                              label=tier_labels[i]) for i in range(4)]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=True,
               fontsize=13, bbox_to_anchor=(0.5, 0.0),
               title="Vulnerability Category")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=[0, 0.045, 1, 0.98])
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print(f"Saved panel -> {OUT}")


if __name__ == "__main__":
    main()
