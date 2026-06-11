"""
Latent Class Growth Analysis (LCGA) of county SVI trajectories, 2000-2022.

LCGA = finite mixture of polynomial growth curves; the mean trajectory depends
only on time and latent class (no within-class random effects). Estimated by EM
(vectorised over counties, per-county missing-year mask). K chosen by BIC elbow.

Outputs
    Data/Processed/SVI/SVI_LCGA.csv               (labels + posteriors)
    Result/CDC_SVI_Trajectory/LCGA_trajectories.csv  (class-mean trajectories)
    Result/CDC_SVI_Trajectory/LCGA_bic.csv           (BIC-vs-K curve)
    Result/CDC_SVI_Trajectory/LCGA_trajectory.png    (standalone figure)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import logsumexp

import CDC_SVI_Trajectory_Utils as U

CSV_OUT = Path("Data/Processed/SVI/SVI_LCGA.csv")
PNG_OUT = U.RESULT_DIR / "LCGA_trajectory.png"

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


def fit_lcga(Y, mask, X, K, degree, rng):
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
    print(f"[LCGA] counties: {N}/{len(df)}")

    rng = np.random.default_rng(SEED)
    Ks, bics, best_per_k = [], [], {}
    for K in range(K_MIN, K_MAX + 1):
        best = None
        for _ in range(N_INIT):
            r = fit_lcga(Y, mask, X, K, DEGREE, rng)
            if best is None or r["loglik"] > best["loglik"]:
                best = r
        bic = -2 * best["loglik"] + n_params(K, DEGREE) * np.log(N)
        Ks.append(K); bics.append(bic); best_per_k[K] = best
        print(f"  K={K}: logL={best['loglik']:.1f}  BIC={bic:.1f}")

    K_best = U.select_elbow(Ks, bics)
    print(f"[LCGA] selected K = {K_best} (BIC elbow)")
    model = best_per_k[K_best]

    mu = X @ model["beta"].T
    labels, resp, mu, hard = U.write_labels(df, keep, model["resp"], mu,
                                            "LCGA", CSV_OUT)
    diag = U.diagnostics(resp)
    counts = diag["counts"]
    print(f"[LCGA] entropy={diag['entropy']:.3f}")

    # ---- export raw results for the composite figure ---- #
    U.export_mixture("LCGA", U.YEARS, mu, counts, Ks, bics, K_best)
    print(f"[LCGA] labels -> {CSV_OUT}")
    print(f"[LCGA] raw    -> {U.RESULT_DIR}/LCGA_trajectories.csv, LCGA_bic.csv")

    # ---- standalone figure ---- #
    fig, ax = plt.subplots(figsize=(11, 7))
    full_idx = np.where(keep)[0]
    U.trajectory_plot(ax, U.YEARS, Y_all, mask_all, full_idx, hard, mu, counts,
                      f"LCGA — County SVI Trajectories ({K_best} classes, 2000–2022)")
    inset = ax.inset_axes([0.66, 0.06, 0.31, 0.30])
    inset.plot(Ks, bics, "-o", color="#444444", ms=5)
    inset.axvline(K_best, color="#E68785", ls="--", lw=1.5, label=f"elbow K={K_best}")
    inset.set_title("BIC vs K", fontsize=11)
    inset.set_xlabel("K", fontsize=10)
    inset.set_xticks(Ks)
    inset.tick_params(labelsize=9)
    inset.legend(fontsize=8, loc="upper right")

    PNG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(PNG_OUT, dpi=200, bbox_inches="tight")
    print(f"[LCGA] plot   -> {PNG_OUT}")


if __name__ == "__main__":
    main()
