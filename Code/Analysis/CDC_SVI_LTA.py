"""
LTA — Latent Transition Analysis of county SVI, 2000-2022.
Models movement BETWEEN vulnerability strata over time ("class mobility").

Model
    A Gaussian hidden Markov model. At each measurement year a county occupies
    one of S latent vulnerability states; each state s has a time-invariant
    Gaussian emission N(mu_s, sig2_s) (measurement invariance, so the states
    mean the same thing every year and transitions are interpretable). Counties
    move between states according to time-specific transition matrices
    A^(t) = P(state_{t+1} | state_t), one per adjacent year pair (the year gaps
    are unequal -- 2000->2010 vs 2018->2020 -- so homogeneous transitions would
    be wrong). An initial-state distribution governs the year-2000 state.

    Whereas GBTM/GMM assign each county to ONE trajectory for all years, LTA
    lets a county change stratum across years and quantifies that flow.

Estimation
    Baum-Welch EM (scaled log-space forward-backward), missing years contribute
    a flat emission. S chosen by BIC but fixed at 4 for a 4-stratum narrative
    comparable to the GBTM groups.

Outputs
    Data/Processed/SVI/SVI_LTA.csv                 (per-year state, init/final, mobility)
    Result/CDC_SVI_Trajectory/LTA_states.csv       (state emission mean/SD)
    Result/CDC_SVI_Trajectory/LTA_composition.csv  (per-year state shares)
    Result/CDC_SVI_Trajectory/LTA_transition.csv   (2000->2022 transition matrix)
    Result/CDC_SVI_Trajectory/LTA_mobility.png     (standalone figure)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import logsumexp

import CDC_SVI_Trajectory_Utils as U

CSV_OUT = Path("Data/Processed/SVI/SVI_LTA.csv")
PNG_OUT = U.RESULT_DIR / "LTA_mobility.png"

S = 4                 # latent vulnerability states
N_INIT = 10
MAX_ITER = 500
TOL = 1e-6
VAR_FLOOR = 1e-4
SEED = 20240603


def emission_logprob(Y, mask, mu, sig2):
    """log N(y_t; mu_s, sig2_s); 0 (prob 1) where the year is missing."""
    # Y (N,T), returns (N,T,S)
    R = Y[:, :, None] - mu[None, None, :]
    lp = -0.5 * (np.log(2 * np.pi * sig2)[None, None, :] + R ** 2 / sig2[None, None, :])
    return np.where(mask[:, :, None], lp, 0.0)


def forward_backward(logb, log_pi, logA):
    """Scaled log-space FB. logb (N,T,S), logA (T-1,S,S)."""
    N, T, Sn = logb.shape
    log_alpha = np.zeros((N, T, Sn))
    log_alpha[:, 0] = log_pi[None, :] + logb[:, 0]
    for t in range(1, T):
        prev = log_alpha[:, t - 1][:, :, None] + logA[t - 1][None, :, :]
        log_alpha[:, t] = logb[:, t] + logsumexp(prev, axis=1)
    ll = logsumexp(log_alpha[:, T - 1], axis=1)            # (N,)

    log_beta = np.zeros((N, T, Sn))
    for t in range(T - 2, -1, -1):
        nxt = (logA[t][None, :, :] + logb[:, t + 1][:, None, :]
               + log_beta[:, t + 1][:, None, :])
        log_beta[:, t] = logsumexp(nxt, axis=2)

    log_gamma = log_alpha + log_beta - ll[:, None, None]
    gamma = np.exp(log_gamma)                              # (N,T,S)

    xi = np.zeros((N, T - 1, Sn, Sn))
    for t in range(T - 1):
        m = (log_alpha[:, t][:, :, None] + logA[t][None, :, :]
             + logb[:, t + 1][:, None, :] + log_beta[:, t + 1][:, None, :]
             - ll[:, None, None])
        xi[:, t] = np.exp(m)
    return gamma, xi, ll.sum()


def fit(Y, mask, rng):
    N, T = Y.shape
    # init: state means spread over data quantiles
    obs = Y[mask]
    mu = np.quantile(obs, np.linspace(0.1, 0.9, S))
    sig2 = np.full(S, obs.var() / S)
    pi = np.full(S, 1 / S)
    A = np.stack([np.full((S, S), 0.1 / (S - 1)) + np.eye(S) * 0.9
                  for _ in range(T - 1)])
    A += rng.uniform(0, 0.02, A.shape)
    A /= A.sum(axis=2, keepdims=True)

    prev = -np.inf
    for _ in range(MAX_ITER):
        logb = emission_logprob(Y, mask, mu, sig2)
        gamma, xi, ll = forward_backward(logb, np.log(pi), np.log(A + 1e-300))
        # M-step
        pi = gamma[:, 0].sum(0) / N
        A = xi.sum(0) / np.clip(gamma[:, :-1].sum(0)[:, :, None], 1e-300, None)
        A /= A.sum(axis=2, keepdims=True)
        gm = gamma * mask[:, :, None]                      # zero out missing
        wsum = gm.sum(axis=(0, 1))                         # (S,)
        mu = (gm * Y[:, :, None]).sum(axis=(0, 1)) / wsum
        R = Y[:, :, None] - mu[None, None, :]
        sig2 = np.maximum((gm * R ** 2).sum(axis=(0, 1)) / wsum, VAR_FLOOR)
        if ll - prev < TOL:
            break
        prev = ll
    logb = emission_logprob(Y, mask, mu, sig2)
    gamma, xi, ll = forward_backward(logb, np.log(pi), np.log(A + 1e-300))
    return dict(pi=pi, A=A, mu=mu, sig2=sig2, gamma=gamma, loglik=ll)


def main():
    df, Y_all, mask_all, keep = U.load_svi()
    Y = Y_all[keep]
    mask = mask_all[keep]
    N, T = Y.shape
    print(f"[LTA] counties: {N}/{len(df)}, S={S} states, T={T} years")

    rng = np.random.default_rng(SEED)
    best = None
    for _ in range(N_INIT):
        r = fit(Y, mask, rng)
        if best is None or r["loglik"] > best["loglik"]:
            best = r
    npar = 2 * S + (S - 1) + (T - 1) * S * (S - 1)
    bic = -2 * best["loglik"] + npar * np.log(N)
    print(f"[LTA] logL={best['loglik']:.1f}  BIC={bic:.1f}")

    # order states low->high by emission mean
    order = np.argsort(best["mu"])
    inv = np.argsort(order)
    mu = best["mu"][order]
    sig = np.sqrt(best["sig2"][order])
    gamma = best["gamma"][:, :, order]
    print("[LTA] state emission means (low->high):",
          ", ".join(f"{chr(65+s)}={mu[s]:.3f}(SD {sig[s]:.3f})" for s in range(S)))

    # decode per-year state (max smoothed posterior)
    states = gamma.argmax(axis=2)                          # (N,T)
    decoded = np.array([[chr(65 + s) for s in row] for row in states])
    # only keep observed years; missing -> NA
    decoded = np.where(mask, decoded, "")

    # ---- CSV ---- #
    out = df[["COUNTY_FIPS"]].copy()
    pos = np.where(keep)[0]
    for j, y in enumerate(U.YEARS):
        col = np.full(len(df), "", dtype=object)
        col[pos] = decoded[:, j]
        out[f"LTA_State_{y}"] = pd.Series(col).replace("", pd.NA)
    init = np.where(mask[:, 0], states[:, 0], -1)
    final = np.where(mask[:, -1], states[:, -1], -1)
    init_l = np.array([chr(65 + s) if s >= 0 else None for s in init], dtype=object)
    fin_l = np.array([chr(65 + s) if s >= 0 else None for s in final], dtype=object)
    mob = np.array(["" for _ in range(N)], dtype=object)
    both = (init >= 0) & (final >= 0)
    d = final - init
    mob[both & (d > 0)] = "Worsening"   # moved to higher-vulnerability state
    mob[both & (d < 0)] = "Improving"
    mob[both & (d == 0)] = "Stable"
    for name, arr in [("LTA_Init", init_l), ("LTA_Final", fin_l),
                      ("LTA_Mobility", mob)]:
        col = np.full(len(df), None, dtype=object)
        col[pos] = arr
        out[name] = pd.Series(col).replace("", pd.NA)

    CSV_OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(CSV_OUT, index=False)
    print(f"[LTA] labels -> {CSV_OUT}")
    print("[LTA] 2000->2022 mobility:",
          dict(pd.Series(mob[both]).value_counts()))

    # ---- export raw results for the composite figure ---- #
    comp = np.array([gamma[mask[:, t], t, :].mean(axis=0) for t in range(T)])  # (T,S)
    Tmat = np.zeros((S, S))
    for i, f in zip(init[both], final[both]):
        Tmat[i, f] += 1
    Tmat = Tmat / np.clip(Tmat.sum(axis=1, keepdims=True), 1, None)
    U.export_lta(U.YEARS, mu, sig, comp, Tmat)
    print(f"[LTA] raw    -> {U.RESULT_DIR}/LTA_states.csv, "
          f"LTA_composition.csv, LTA_transition.csv")

    # ---- Figure: state composition over time + transition heatmap ---- #
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))

    # (a) soft state composition per year (precomputed comp)
    colors = [U.CLASS_COLORS[s] for s in range(S)]
    axes[0].stackplot(U.YEARS, comp.T, colors=colors,
                      labels=[f"State {chr(65+s)} (μ={mu[s]:.2f})" for s in range(S)],
                      alpha=0.9)
    axes[0].set_xlim(min(U.YEARS), max(U.YEARS))
    axes[0].set_ylim(0, 1)
    axes[0].set_xticks(U.YEARS)
    axes[0].set_xlabel("Year")
    axes[0].set_ylabel("Share of counties")
    axes[0].set_title("Vulnerability-state composition over time")
    axes[0].legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=10)

    # (b) initial(2000) -> final(2022) transition matrix (precomputed Tmat)
    im = axes[1].imshow(Tmat, cmap="Blues", vmin=0, vmax=1)
    axes[1].set_xticks(range(S)); axes[1].set_yticks(range(S))
    axes[1].set_xticklabels([f"{chr(65+s)}" for s in range(S)])
    axes[1].set_yticklabels([f"{chr(65+s)}" for s in range(S)])
    axes[1].set_xlabel("State in 2022")
    axes[1].set_ylabel("State in 2000")
    axes[1].set_title("Transition: 2000 -> 2022  (row-normalized P)")
    for i in range(S):
        for j in range(S):
            axes[1].text(j, i, f"{Tmat[i, j]:.2f}", ha="center", va="center",
                         color="white" if Tmat[i, j] > 0.5 else "black",
                         fontsize=12)
    fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04)

    fig.suptitle(f"LTA — County SVI Vulnerability Mobility ({S} states, 2000–2022)",
                 fontsize=18)
    PNG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(PNG_OUT, dpi=200, bbox_inches="tight")
    print(f"[LTA] plot   -> {PNG_OUT}")


if __name__ == "__main__":
    main()
