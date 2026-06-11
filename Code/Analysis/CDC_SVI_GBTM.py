"""
GBTM — Group-Based Trajectory Modeling (Nagin) of county SVI, 2000-2022.
MAIN ANALYSIS.

Model
    Finite mixture of polynomial growth curves with NO within-group random
    effects: every county in group k follows the same mean trajectory
    X * beta_k and deviates only by i.i.d. normal measurement error (variance
    sigma2_k). This is Nagin's GBTM; under a normal emission it is the same
    estimator as LCGA. Boundary SVI values (0 / 1) occur in only ~0.03% of
    observations, so a censored-normal likelihood is unnecessary.

Estimation
    EM, vectorised over counties, per-county missing-year handling via a mask.
    Number of groups chosen by the BIC elbow over K = 2..K_MAX (BIC is
    monotone here, so min-BIC would just hit the ceiling).

Outputs
    Data/Processed/SVI/SVI_GBTM.csv                  (labels + posteriors)
    Result/CDC_SVI_Trajectory/GBTM_trajectories.csv  (class-mean trajectories)
    Result/CDC_SVI_Trajectory/GBTM_bic.csv           (BIC-vs-K curve)
    Result/CDC_SVI_Trajectory/GBTM_trajectory.png    (standalone figure)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import logsumexp

import CDC_SVI_Trajectory_Utils as U

CSV_OUT = Path("Data/Processed/SVI/SVI_GBTM.csv")
PNG_OUT = U.RESULT_DIR / "GBTM_trajectory.png"

DEGREE = 2
K_MIN, K_MAX = 2, 8
N_INIT = 40
MAX_ITER = 500
TOL = 1e-7
VAR_FLOOR = 1e-4
SEED = 20240603


def e_step(Y, mask, n_obs, X, beta, var, pi):
    mu = X @ beta.T
    resid = Y[:, :, None] - mu[None, :, :]
    sq = np.where(mask[:, :, None], resid ** 2, 0.0).sum(axis=1)
    log_dens = -0.5 * (n_obs[:, None] * (np.log(2 * np.pi) + np.log(var)[None, :])
                       + sq / var[None, :])
    log_unnorm = np.log(pi)[None, :] + log_dens
    ll = logsumexp(log_unnorm, axis=1)
    resp = np.exp(log_unnorm - ll[:, None])
    return resp, ll.sum()


def m_step(Y, mask, n_obs, X, resp, degree):
    K = resp.shape[1]
    p = degree + 1
    pi = resp.mean(axis=0)
    W = mask.T @ resp
    YM = np.where(mask, Y, 0.0)
    SY = (resp.T @ YM).T
    beta = np.zeros((K, p))
    for k in range(K):
        A = (X.T * W[:, k]) @ X
        b = X.T @ SY[:, k]
        beta[k] = np.linalg.solve(A + 1e-8 * np.eye(p), b)
    mu = X @ beta.T
    resid = Y[:, :, None] - mu[None, :, :]
    sq = np.where(mask[:, :, None], resid ** 2, 0.0).sum(axis=1)
    denom = (resp * n_obs[:, None]).sum(axis=0)
    var = np.maximum((resp * sq).sum(axis=0) / denom, VAR_FLOOR)
    return pi, beta, var


def fit(Y, mask, X, K, degree, rng):
    N = Y.shape[0]
    n_obs = mask.sum(axis=1).astype(float)
    labels = rng.integers(0, K, size=N)
    resp = np.zeros((N, K))
    resp[np.arange(N), labels] = 1.0
    pi, beta, var = m_step(Y, mask, n_obs, X, resp, degree)
    prev = -np.inf
    for _ in range(MAX_ITER):
        resp, ll = e_step(Y, mask, n_obs, X, beta, var, pi)
        pi, beta, var = m_step(Y, mask, n_obs, X, resp, degree)
        if ll - prev < TOL:
            break
        prev = ll
    resp, ll = e_step(Y, mask, n_obs, X, beta, var, pi)
    return dict(pi=pi, beta=beta, var=var, resp=resp, loglik=ll, K=K)


def n_params(K, degree):
    return (K - 1) + K * (degree + 1) + K


def main():
    df, Y_all, mask_all, keep = U.load_svi()
    Y = Y_all[keep]
    mask = mask_all[keep]
    X = U.design_matrix(U.YEARS, DEGREE)
    N = Y.shape[0]
    print(f"[GBTM] counties: {N}/{len(df)}")

    rng = np.random.default_rng(SEED)
    rows, best_per_k = [], {}
    for K in range(K_MIN, K_MAX + 1):
        best = None
        for _ in range(N_INIT):
            r = fit(Y, mask, X, K, DEGREE, rng)
            if best is None or r["loglik"] > best["loglik"]:
                best = r
        npar = n_params(K, DEGREE)
        bic = -2 * best["loglik"] + npar * np.log(N)
        best_per_k[K] = best
        rows.append((K, best["loglik"], bic))
        print(f"  K={K}: logL={best['loglik']:.1f}  BIC={bic:.1f}")

    sel = pd.DataFrame(rows, columns=["K", "logL", "BIC"])
    K_best = U.select_elbow(sel["K"].values, sel["BIC"].values)
    print(f"[GBTM] selected K = {K_best} (BIC elbow)")
    model = best_per_k[K_best]

    mu = X @ model["beta"].T
    labels, resp, mu, hard = U.write_labels(df, keep, model["resp"], mu,
                                             "GBTM", CSV_OUT)
    diag = U.diagnostics(resp)
    print(f"[GBTM] entropy={diag['entropy']:.3f}")
    for j in range(K_best):
        print(f"  Class {chr(65+j)}: n={diag['counts'][j]:4d}  "
              f"AvePP={diag['avepp'][j]:.3f}  OCC={diag['occ'][j]:.1f}")
    U.export_mixture("GBTM", U.YEARS, mu, diag["counts"],
                     sel["K"].values, sel["BIC"].values, K_best)
    print(f"[GBTM] labels -> {CSV_OUT}")
    print(f"[GBTM] raw    -> {U.RESULT_DIR}/GBTM_trajectories.csv, GBTM_bic.csv")

    # ---- Figure: trajectories + BIC inset + diagnostics table ---------- #
    fig, ax = plt.subplots(figsize=(11, 7))
    full_idx = np.where(keep)[0]
    U.trajectory_plot(ax, U.YEARS, Y_all, mask_all, full_idx, hard, mu,
                      diag["counts"],
                      f"GBTM — County SVI Trajectories ({K_best} groups, 2000–2022)")

    inset = ax.inset_axes([0.66, 0.06, 0.31, 0.30])
    inset.plot(sel["K"], sel["BIC"], "-o", color="#444444", ms=5)
    inset.axvline(K_best, color="#E68785", ls="--", lw=1.5,
                  label=f"elbow K={K_best}")
    inset.set_title("BIC vs K", fontsize=11)
    inset.set_xlabel("K", fontsize=10)
    inset.set_xticks(sel["K"])
    inset.tick_params(labelsize=9)
    inset.legend(fontsize=8, loc="upper right")

    txt = "AvePP / OCC (Nagin)\n" + "\n".join(
        f"{chr(65+j)}: {diag['avepp'][j]:.2f} / {diag['occ'][j]:.0f}"
        for j in range(K_best))
    ax.text(0.30, 0.62, txt, transform=ax.transAxes, ha="left", va="top",
            fontsize=10, family="monospace",
            bbox=dict(boxstyle="round", fc="white", ec="#999", alpha=0.9))

    PNG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(PNG_OUT, dpi=200, bbox_inches="tight")
    print(f"[GBTM] plot   -> {PNG_OUT}")


if __name__ == "__main__":
    main()
