"""
Shared utilities for county SVI trajectory analyses (GBTM / GMM / LTA).

Provides data loading, polynomial design matrices, model-selection helpers,
Nagin-style classification diagnostics, and a common trajectory plot so the
three methods are directly comparable.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

IN_PATH = Path("Data/Processed/SVI/SVI.csv")
RESULT_DIR = Path("Result/CDC_SVI_Trajectory")   # shared raw-result output dir
YEARS = [2000, 2010, 2014, 2016, 2018, 2020, 2022]
SVI_COLS = [f"{y}_SVI_RPL" for y in YEARS]

# Palette ordered low -> high SVI, shared across all methods/figures
CLASS_COLORS = ["#2F7F4F", "#97C889", "#E6EAB8", "#699DCB",
                "#E68785", "#EFC085", "#9E79B8", "#7FB8B0"]

plt.rcParams["font.family"] = "Georgia"
plt.rcParams["font.size"] = 14


def load_svi():
    """Return (df, Y, mask, keep) where Y is (N,T) with NaN->0 under mask."""
    df = pd.read_csv(IN_PATH, dtype={"COUNTY_FIPS": str})
    df["COUNTY_FIPS"] = df["COUNTY_FIPS"].str.zfill(5)
    Y_all = df[SVI_COLS].to_numpy(float)
    mask_all = ~np.isnan(Y_all)
    keep = mask_all.sum(axis=1) > 0
    Y = np.where(mask_all, Y_all, 0.0)
    return df, Y, mask_all, keep


def design_matrix(years, degree):
    """Polynomial design on centred, scaled calendar time (good conditioning)."""
    t = (np.asarray(years, float) - np.mean(years)) / 10.0
    return np.vander(t, degree + 1, increasing=True)  # [1, t, t^2, ...]


def select_elbow(Ks, bics):
    """Kneedle-style elbow on a (monotone) BIC-vs-K curve.

    For mixture models on large N, BIC is often monotone-decreasing with no
    interior minimum, so plain min-BIC degenerates to the largest K. The elbow
    marks diminishing returns -- the standard fallback for a parsimonious,
    interpretable number of trajectory classes.
    """
    Ks = np.asarray(Ks, float)
    bics = np.asarray(bics, float)
    x = (Ks - Ks.min()) / (Ks.max() - Ks.min())
    y = (bics - bics[0]) / (bics[-1] - bics[0])
    return int(Ks[np.argmax(np.abs(y - x))])


def diagnostics(resp):
    """Nagin classification-quality diagnostics from responsibilities (N,K).

    Returns dict with per-class average posterior probability (AvePP, target
    > 0.7), odds of correct classification (OCC, target > 5), mixing weights,
    hard-assigned counts, and overall relative entropy (target > 0.8).
    """
    N, K = resp.shape
    hard = resp.argmax(axis=1)
    pi = resp.mean(axis=0)
    counts = np.bincount(hard, minlength=K)
    avepp = np.array([resp[hard == k, k].mean() if counts[k] else np.nan
                      for k in range(K)])
    occ = (avepp / (1 - avepp)) / (pi / (1 - pi))
    r = np.clip(resp, 1e-12, 1)
    entropy = 1 - (-(r * np.log(r)).sum()) / (N * np.log(K))
    return dict(hard=hard, pi=pi, counts=counts, avepp=avepp,
                occ=occ, entropy=entropy)


def trajectory_plot(ax, years, Y_all, mask_all, full_idx, hard, mu, counts,
                    title, ylabel="SVI overall percentile (RPL_THEMES)"):
    """Faint individual county lines coloured by class + bold class means."""
    K = mu.shape[1]
    years = np.asarray(years)
    for j in range(K):
        col = CLASS_COLORS[j % len(CLASS_COLORS)]
        for i in full_idx[hard == j]:
            m = mask_all[i]
            ax.plot(years[m], Y_all[i][m], color=col, alpha=0.03, lw=0.6)
    for j in range(K):
        col = CLASS_COLORS[j % len(CLASS_COLORS)]
        pct = 100 * counts[j] / counts.sum()
        ax.plot(years, mu[:, j], color=col, lw=3.2, marker="o", ms=8,
                markeredgecolor="white", markeredgewidth=1.2, zorder=5,
                label=f"Class {chr(65 + j)} (n={counts[j]}, {pct:.1f}%)")
    ax.set_xlabel("Year")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(years)
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left", framealpha=0.9, fontsize=11)


def order_classes(mu):
    """Return index order sorting latent classes by mean level (low->high)."""
    return np.argsort(mu.mean(axis=0))


def write_labels(df, keep, resp, mu, method, out_path):
    """Relabel classes low->high, attach labels + posteriors, write CSV.

    Returns (hard_labels_letters, ordered_resp, ordered_mu).
    """
    order = order_classes(mu)
    resp = resp[:, order]
    mu = mu[:, order]
    hard = resp.argmax(axis=1)
    K = resp.shape[1]
    labels = np.array([f"Class {chr(65 + c)}" for c in hard])

    out = df.copy()
    col = f"{method}_Class"
    out[col] = pd.NA
    out.loc[keep, col] = labels
    pos = np.where(keep)[0]
    for j in range(K):
        vals = np.full(len(df), np.nan)
        vals[pos] = resp[:, j]
        out[f"Prob_{chr(65 + j)}"] = vals
    maxp = np.full(len(df), np.nan)
    maxp[pos] = resp.max(axis=1)
    out["Max_Prob"] = maxp

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    return labels, resp, mu, hard


# --------------------------------------------------------------------------- #
# Raw-result export (consumed by Code/Visualization/SVI_Trajectory.py)
# --------------------------------------------------------------------------- #
def export_mixture(method, years, mu, counts, Ks, bics, K_best):
    """Write class-mean trajectories and the BIC-vs-K curve for a mixture model.

    mu (T,K): predicted class mean at each year, ordered low->high.
    """
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    K = mu.shape[1]
    total = counts.sum()
    rows = [(f"Class {chr(65 + j)}", y, mu[ti, j], int(counts[j]),
             100 * counts[j] / total)
            for j in range(K) for ti, y in enumerate(years)]
    pd.DataFrame(rows, columns=["Class", "Year", "Mean", "N", "Pct"]).to_csv(
        RESULT_DIR / f"{method}_trajectories.csv", index=False)
    pd.DataFrame({"K": Ks, "BIC": bics,
                  "Selected": [int(k == K_best) for k in Ks]}).to_csv(
        RESULT_DIR / f"{method}_bic.csv", index=False)


def export_lta(years, mu, sd, comp, Tmat):
    """Write LTA state emissions, per-year composition, and transition matrix."""
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    S = len(mu)
    letters = [chr(65 + s) for s in range(S)]
    pd.DataFrame({"State": letters, "Mean": mu, "SD": sd}).to_csv(
        RESULT_DIR / "LTA_states.csv", index=False)
    comp_df = pd.DataFrame(comp, columns=letters)
    comp_df.insert(0, "Year", years)
    comp_df.to_csv(RESULT_DIR / "LTA_composition.csv", index=False)
    tmat_df = pd.DataFrame(Tmat, columns=letters)
    tmat_df.insert(0, "From", letters)      # rows = state in 2000, cols = 2022
    tmat_df.to_csv(RESULT_DIR / "LTA_transition.csv", index=False)
