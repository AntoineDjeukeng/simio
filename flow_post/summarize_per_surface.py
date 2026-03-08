#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from config_flow import (
    RESULTS,
    FIELD_TO_E,
    E_ORDER,
    CHARGES,
    OHMIC_E_LIST,
    OHMIC_E_MAX_V_NM,
    SIGMA_E_MAX_V_NM,
)

MASTER = RESULTS / "master_transport_partial.tsv"
FIELD_SET = set(FIELD_TO_E.keys())

df = pd.read_csv(MASTER, sep="\t")
df = df[df["charge"].isin(CHARGES)].copy()
df["charge"] = df["charge"].replace({"positive":"pos", "negative":"neg"})

def has_all_fields(g):
    return set(g["field"].unique()) >= FIELD_SET

ok = []
for (H,L,ch), g in df.groupby(["H_nm","L_nm","charge"], sort=False):
    if has_all_fields(g):
        ok.append(g)
if not ok:
    raise SystemExit("No (H,L,charge) has all 4 fields yet.")

complete = pd.concat(ok, ignore_index=True)
complete = complete.sort_values(["charge","H_nm","L_nm","E_V_nm"])

out1 = RESULTS / "per_surface_complete.tsv"
complete.to_csv(out1, sep="\t", index=False)
print("Wrote:", out1)

# wide table
pivot_values = [
    c for c in [
        "I_mean_A", "I_sem_A",
        "I_NA_mean_A", "I_NA_sem_A",
        "I_CL_mean_A", "I_CL_sem_A",
        "sigma_mean_S_m", "sigma_sem_S_m",
        "sigma_NA_mean_S_m", "sigma_NA_sem_S_m",
        "sigma_CL_mean_S_m", "sigma_CL_sem_S_m",
    ]
    if c in complete.columns
]
wide = complete.pivot_table(
    index=["H_nm","L_nm","charge"],
    columns="E_V_nm",
    values=pivot_values,
    aggfunc="first"
).reset_index()
wide.columns = [f"{a}_{b}" if b != "" else str(a) for (a,b) in wide.columns.to_flat_index()]

out2 = RESULTS / "per_surface_complete_wide.tsv"
wide.to_csv(out2, sep="\t", index=False)
print("Wrote:", out2)
print("Completed surfaces:", wide.shape[0])

def in_linear_regime(df_in):
    g = df_in[df_in["E_V_nm"] > 0].copy()
    if OHMIC_E_LIST:
        vals = np.asarray([float(v) for v in OHMIC_E_LIST], float)
        mask = np.any(np.isclose(g["E_V_nm"].to_numpy(float)[:, None], vals[None, :], atol=1e-12), axis=1)
        g = g[mask]
    elif OHMIC_E_MAX_V_NM is not None:
        g = g[g["E_V_nm"] <= OHMIC_E_MAX_V_NM]
    return g

complete_linear = in_linear_regime(complete)
out_linear = RESULTS / "per_surface_complete_linear.tsv"
complete_linear.to_csv(out_linear, sep="\t", index=False)
print("Wrote:", out_linear)

if not complete_linear.empty:
    wide_linear = complete_linear.pivot_table(
        index=["H_nm","L_nm","charge"],
        columns="E_V_nm",
        values=pivot_values,
        aggfunc="first"
    ).reset_index()
    wide_linear.columns = [f"{a}_{b}" if b != "" else str(a) for (a,b) in wide_linear.columns.to_flat_index()]
    out_linear_wide = RESULTS / "per_surface_complete_linear_wide.tsv"
    wide_linear.to_csv(out_linear_wide, sep="\t", index=False)
    print("Wrote:", out_linear_wide)

def _as_ordered(g):
    gg = g[g["E_V_nm"].isin(E_ORDER)].copy()
    gg["E_V_nm"] = pd.Categorical(gg["E_V_nm"], categories=E_ORDER, ordered=True)
    return gg.sort_values("E_V_nm")

def _ohmic_fit_subset(g):
    out = g[g["E_V_nm"] > 0].copy()
    if OHMIC_E_LIST:
        vals = np.asarray([float(v) for v in OHMIC_E_LIST], float)
        mask = np.any(np.isclose(out["E_V_nm"].to_numpy(float)[:, None], vals[None, :], atol=1e-12), axis=1)
        out = out[mask]
    elif OHMIC_E_MAX_V_NM is not None:
        out = out[out["E_V_nm"] <= OHMIC_E_MAX_V_NM]
    return out

def _charge_order(values):
    order = ["neg", "neutral", "pos"]
    present = [c for c in order if c in values]
    for c in sorted(values):
        if c not in present:
            present.append(c)
    return present

def _plot_I_vs_E_all_L(gH, H):
    ls = sorted(gH["L_nm"].unique())
    chs = _charge_order(gH["charge"].unique())
    colors = {"neg": "tab:blue", "neutral": "tab:gray", "pos": "tab:red"}
    fig, axes = plt.subplots(1, len(ls), figsize=(5.2 * len(ls), 4.6), sharey=True)
    if len(ls) == 1:
        axes = [axes]
    for ax, L in zip(axes, ls):
        gL = gH[gH["L_nm"] == L]
        for ch in chs:
            g = _as_ordered(gL[gL["charge"] == ch])
            if g.empty:
                continue
            c = colors.get(ch, None)
            x = g["E_V_nm"].to_numpy(float)
            y = g["I_mean_A"].to_numpy(float)
            yerr = g["I_sem_A"].to_numpy(float) if "I_sem_A" in g.columns else None
            ax.errorbar(x, y, yerr=yerr, marker="o", capsize=3, linewidth=1.8, color=c, label=ch)

            gfit = _ohmic_fit_subset(g)
            if len(gfit) >= 2:
                xf = gfit["E_V_nm"].to_numpy(float)
                yf = gfit["I_mean_A"].to_numpy(float)
                a, b = np.polyfit(xf, yf, 1)
                xline = np.array([float(np.nanmin(x)), float(np.nanmax(x))], float)
                yline = a * xline + b
                ax.plot(xline, yline, linestyle="--", linewidth=1.2, color=c, alpha=0.9)
        ax.axhline(0.0, linestyle="--", linewidth=0.8, color="black")
        ax.set_title(f"L={int(L)} nm")
        ax.set_xlabel("E (V/nm)")
    axes[0].set_ylabel("I (A)")
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles, labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.94),
            ncol=min(len(labels), 4),
            frameon=False,
        )
    fig.suptitle(f"I vs E (points) + linear-regime fit (dashed) | H={H:.1f} nm", y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out = RESULTS / f"I_vs_E_all_L_H{H:.1f}.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", out)

def _plot_I_components_all(gH, H):
    if not {"I_NA_mean_A", "I_CL_mean_A", "I_NA_sem_A", "I_CL_sem_A"} <= set(gH.columns):
        return
    ls = sorted(gH["L_nm"].unique())
    chs = _charge_order(gH["charge"].unique())
    nrows, ncols = len(ls), len(chs)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.8 * ncols, 3.8 * nrows), sharex=True, sharey=True)
    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = np.array([axes])
    elif ncols == 1:
        axes = np.array([[ax] for ax in axes])

    for i, L in enumerate(ls):
        for j, ch in enumerate(chs):
            ax = axes[i, j]
            g = _as_ordered(gH[(gH["L_nm"] == L) & (gH["charge"] == ch)])
            if g.empty:
                ax.set_visible(False)
                continue
            x = g["E_V_nm"].to_numpy(float)
            ax.errorbar(x, g["I_mean_A"].to_numpy(float), yerr=g["I_sem_A"].to_numpy(float),
                        marker="o", capsize=3, linewidth=1.6, label="I_total")
            ax.errorbar(x, g["I_NA_mean_A"].to_numpy(float), yerr=g["I_NA_sem_A"].to_numpy(float),
                        marker="s", capsize=3, linewidth=1.6, label="I_NA")
            ax.errorbar(x, g["I_CL_mean_A"].to_numpy(float), yerr=g["I_CL_sem_A"].to_numpy(float),
                        marker="^", capsize=3, linewidth=1.6, label="I_CL")
            ax.axhline(0.0, linestyle="--", linewidth=0.8, color="black")
            ax.set_title(f"L={int(L)} nm | {ch}")
            if i == nrows - 1:
                ax.set_xlabel("E (V/nm)")
            if j == 0:
                ax.set_ylabel("I (A)")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles, labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.94),
            ncol=3,
            frameon=False,
        )
    fig.suptitle(f"I components vs E | H={H:.1f} nm", y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out = RESULTS / f"I_components_all_H{H:.1f}.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", out)

def _plot_sigma_vs_E_all_L(gH, H):
    if "sigma_mean_S_m" not in gH.columns:
        return
    ls = sorted(gH["L_nm"].unique())
    chs = _charge_order(gH["charge"].unique())
    colors = {"neg": "tab:blue", "neutral": "tab:gray", "pos": "tab:red"}
    fig, axes = plt.subplots(1, len(ls), figsize=(5.2 * len(ls), 4.6), sharey=True)
    if len(ls) == 1:
        axes = [axes]
    for ax, L in zip(axes, ls):
        gL = gH[gH["L_nm"] == L]
        for ch in chs:
            g = _as_ordered(gL[gL["charge"] == ch])
            g = g[g["E_V_nm"] > 0]
            if SIGMA_E_MAX_V_NM is not None:
                g = g[g["E_V_nm"] <= SIGMA_E_MAX_V_NM]
            g = g[np.isfinite(g["sigma_mean_S_m"].to_numpy(float))]
            if g.empty:
                continue
            c = colors.get(ch, None)
            x = g["E_V_nm"].to_numpy(float)
            y = g["sigma_mean_S_m"].to_numpy(float)
            yerr = g["sigma_sem_S_m"].to_numpy(float) if "sigma_sem_S_m" in g.columns else None
            ax.errorbar(x, y, yerr=yerr, marker="o", capsize=3, linewidth=1.8, color=c, label=ch)
        ax.set_title(f"L={int(L)} nm")
        ax.set_xlabel("E (V/nm)")
    axes[0].set_ylabel("σ (S/m)")
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles, labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.94),
            ncol=min(len(labels), 4),
            frameon=False,
        )
    fig.suptitle(f"σ vs E (linear regime) | H={H:.1f} nm", y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out = RESULTS / f"sigma_vs_E_all_L_H{H:.1f}.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", out)

for H, gH in complete.groupby("H_nm", sort=True):
    _plot_I_vs_E_all_L(gH, H)
    _plot_I_components_all(gH, H)
    _plot_sigma_vs_E_all_L(gH, H)
