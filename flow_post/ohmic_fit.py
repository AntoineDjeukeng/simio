#!/usr/bin/env python3
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from config_flow import (
    RESULTS,
    OHMIC_E_MAX_V_NM,
    OHMIC_E_LIST,
    OHMIC_FIT_MODE,
    OHMIC_WEIGHT,
)

MASTER = RESULTS / "per_surface_complete.tsv"
df = pd.read_csv(MASTER, sep="\t")

def _fit_r2(y, yhat):
    ss_res = float(((y - yhat) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum()) if len(y) > 1 else np.nan
    return (1.0 - ss_res / ss_tot) if (ss_tot and ss_tot > 0) else np.nan

def _weighted_r2(y, yhat, w):
    wsum = float(np.sum(w))
    if not np.isfinite(wsum) or wsum <= 0:
        return np.nan
    y_bar_w = float(np.sum(w * y) / wsum)
    ss_res_w = float(np.sum(w * (y - yhat) ** 2))
    ss_tot_w = float(np.sum(w * (y - y_bar_w) ** 2))
    return (1.0 - ss_res_w / ss_tot_w) if ss_tot_w > 0 else np.nan

def linear_fit(x, y, y_sem=None, mode="WLS", weight="I_sem"):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    X = np.column_stack([x, np.ones_like(x)])
    n = len(x)
    dof = n - 2

    mode_up = str(mode).upper()
    weight_low = str(weight).lower()
    use_wls = mode_up == "WLS" and weight_low == "i_sem" and y_sem is not None

    fit_mode_used = "OLS"
    fit_weight_used = "none"
    chi2 = np.nan
    chi2_red = np.nan
    r2_w = np.nan

    if use_wls:
        y_sem = np.asarray(y_sem, float)
        ok = np.isfinite(y_sem) & (y_sem > 0)
        if ok.all():
            w = 1.0 / (y_sem ** 2)
            XTWX = X.T @ (w[:, None] * X)
            XTWy = X.T @ (w * y)
            beta = np.linalg.pinv(XTWX) @ XTWy
            cov = np.linalg.pinv(XTWX)
            yhat = X @ beta
            resid = y - yhat
            chi2 = float(np.sum((resid / y_sem) ** 2))
            chi2_red = float(chi2 / dof) if dof > 0 else np.nan
            r2_w = _weighted_r2(y, yhat, w)
            fit_mode_used = "WLS"
            fit_weight_used = "I_sem"
        else:
            use_wls = False

    if not use_wls:
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        yhat = X @ beta
        resid = y - yhat
        ss_res = float(np.sum(resid ** 2))
        if dof > 0:
            s2 = ss_res / dof
            cov = s2 * np.linalg.pinv(X.T @ X)
        else:
            cov = np.full((2, 2), np.nan)

    a = float(beta[0])
    b = float(beta[1])
    a_se = float(np.sqrt(cov[0, 0])) if np.isfinite(cov[0, 0]) and cov[0, 0] >= 0 else np.nan
    b_se = float(np.sqrt(cov[1, 1])) if np.isfinite(cov[1, 1]) and cov[1, 1] >= 0 else np.nan
    r2 = _fit_r2(y, yhat)

    return {
        "a": a,
        "b": b,
        "a_se": a_se,
        "b_se": b_se,
        "r2": r2,
        "r2_w": float(r2_w),
        "chi2": float(chi2),
        "chi2_red": float(chi2_red),
        "n_points": int(n),
        "dof": int(dof),
        "fit_mode_used": fit_mode_used,
        "fit_weight_used": fit_weight_used,
    }

def in_ohmic_regime(df_in):
    g = df_in[df_in["E_V_nm"] > 0].copy()
    if OHMIC_E_LIST:
        target = np.asarray([float(v) for v in OHMIC_E_LIST], float)
        mask = np.any(np.isclose(g["E_V_nm"].to_numpy(float)[:, None], target[None, :], atol=1e-12), axis=1)
        g = g[mask]
    elif OHMIC_E_MAX_V_NM is not None:
        g = g[g["E_V_nm"] <= OHMIC_E_MAX_V_NM]
    return g

df_nz = in_ohmic_regime(df)

rows = []
for (H,L,ch), g in df_nz.groupby(["H_nm","L_nm","charge"], sort=True):
    g = g.sort_values("E_V_nm")
    x = g["E_V_nm"].to_numpy(float)  # V/nm
    y = g["I_mean_A"].to_numpy(float)
    if len(x) < 2:
        continue
    y_sem = g["I_sem_A"].to_numpy(float) if "I_sem_A" in g.columns else None
    fit = linear_fit(x, y, y_sem=y_sem, mode=OHMIC_FIT_MODE, weight=OHMIC_WEIGHT)

    # conductivity from slope a:
    # a = dI/dE with E in V/nm
    # convert to SI:
    # dI/dE_SI = a * 1e-9  (since 1 V/nm = 1e9 V/m)
    S_m2 = float(g["S_nm2"].iloc[0]) * 1e-18
    sigma_fit = (fit["a"] * 1e-9) / S_m2  # S/m
    sigma_fit_se = (fit["a_se"] * 1e-9) / S_m2 if np.isfinite(fit["a_se"]) else np.nan

    rows.append({
        "H_nm": H, "L_nm": L, "charge": ch,
        "fit_mode": fit["fit_mode_used"],
        "fit_weight": fit["fit_weight_used"],
        "n_points": fit["n_points"],
        "dof": fit["dof"],
        "e_fields_used": json.dumps(g["E_V_nm"].tolist()),
        "a_dIdE_A_per_Vnm": float(fit["a"]),
        "a_se_A_per_Vnm": float(fit["a_se"]),
        "b_A": float(fit["b"]),
        "b_se_A": float(fit["b_se"]),
        "R2": float(fit["r2"]),
        "R2_weighted": float(fit["r2_w"]),
        "chi2": float(fit["chi2"]),
        "chi2_red": float(fit["chi2_red"]),
        "S_nm2": float(g["S_nm2"].iloc[0]),
        "sigma_fit_S_m": float(sigma_fit),
        "sigma_fit_se_S_m": float(sigma_fit_se),
    })

if not rows:
    raise SystemExit("No valid groups for Ohmic fit after applying field filters.")

out = pd.DataFrame(rows).sort_values(["charge","H_nm","L_nm"])
out_tsv = RESULTS / "ohmic_fit.tsv"
out.to_csv(out_tsv, sep="\t", index=False)
print("Wrote:", out_tsv)

# composite plot: sigma_fit vs L, all surfaces together per height
charge_order = ["neg", "neutral", "pos"]
for H, gH in out.groupby("H_nm", sort=True):
    plt.figure(figsize=(7.8, 5.0))
    present = list(gH["charge"].unique())
    ordered = [c for c in charge_order if c in present]
    for c in sorted(present):
        if c not in ordered:
            ordered.append(c)

    for ch in ordered:
        g = gH[gH["charge"] == ch].sort_values("L_nm")
        if g.empty:
            continue
        yerr = g["sigma_fit_se_S_m"] if "sigma_fit_se_S_m" in g.columns else None
        plt.errorbar(g["L_nm"], g["sigma_fit_S_m"], yerr=yerr, marker="o", linewidth=2, capsize=3, label=ch)

    plt.xlabel("L (nm)")
    plt.ylabel("σ_fit (S/m)")
    plt.title(f"σ_fit vs L | H={H:.1f} nm")
    plt.legend(title="surface")
    plt.tight_layout()
    pdf = RESULTS / f"sigma_fit_vs_L_H{H:.1f}.pdf"
    plt.savefig(pdf, bbox_inches="tight")
    plt.close()
    print("Saved:", pdf)
