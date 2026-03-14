#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from config_flow import PATTERN, I_START, I_END

def step_hold_interp(x_new, x, y):
    idx = np.searchsorted(x, x_new, side="right") - 1
    idx = np.clip(idx, 0, len(y) - 1)
    return y[idx]

def robust_dt_ps(t):
    d = np.diff(t)
    d = d[d > 0]
    return float(np.median(d)) if d.size else float("nan")

def make_grid_ps(t_min, t_max, dt_ps):
    n = int(np.floor((t_max - t_min) / dt_ps)) + 1
    return t_min + np.arange(n) * dt_ps

def mean_sem(Y):
    m = Y.mean(axis=1)
    if Y.shape[1] > 1:
        sd = Y.std(axis=1, ddof=1)
        se = sd / np.sqrt(Y.shape[1])
    else:
        se = np.zeros_like(m)
    return m, se

def read_replica(path: Path, baseline_zero=True):
    df = pd.read_csv(path, sep=r"\s+", comment="#", header=None, engine="python")
    df = df.iloc[:, :4].copy()
    df.columns = ["time_ps", "SOL", "NA", "CL"]  # project truth
    df = df.sort_values("time_ps").reset_index(drop=True)
    t = df["time_ps"].to_numpy(float)
    SOL = df["SOL"].to_numpy(float)
    NA  = df["NA"].to_numpy(float)
    CL  = df["CL"].to_numpy(float)
    if baseline_zero:
        SOL -= SOL[0]; NA -= NA[0]; CL -= CL[0]
    return t, SOL, NA, CL

def reduce_compile_dir(compile_dir: Path, baseline_zero=True, plot_replicas=True):
    indir = compile_dir
    outdir = indir / "reduced"
    outdir.mkdir(parents=True, exist_ok=True)

    reps, used = [], []
    for i in range(I_START, I_END + 1):
        f = indir / PATTERN.format(i=i)
        if not f.exists():
            continue
        reps.append(read_replica(f, baseline_zero=baseline_zero))
        used.append(str(f))
    if not reps:
        raise RuntimeError(f"No replicas found in {indir}")

    # common time grid
    dts = [robust_dt_ps(t) for (t, *_r) in reps]
    dts = [dt for dt in dts if np.isfinite(dt) and dt > 0]
    dt_ps = float(np.median(dts))

    t_min = max(t[0] for (t, *_r) in reps)
    t_max = min(t[-1] for (t, *_r) in reps)
    grid_ps = make_grid_ps(t_min, t_max, dt_ps)

    SOL = np.vstack([step_hold_interp(grid_ps, t, sol) for (t, sol, na, cl) in reps]).T
    NA  = np.vstack([step_hold_interp(grid_ps, t, na)  for (t, sol, na, cl) in reps]).T
    CL  = np.vstack([step_hold_interp(grid_ps, t, cl)  for (t, sol, na, cl) in reps]).T
    qN  = NA - CL

    sol_m, sol_se = mean_sem(SOL)
    na_m,  na_se  = mean_sem(NA)
    cl_m,  cl_se  = mean_sem(CL)
    qN_m,  qN_se  = mean_sem(qN)

    out = pd.DataFrame({
        "time_ps": grid_ps,
        "mean_SOL": sol_m, "sem_SOL": sol_se,
        "mean_NA":  na_m,  "sem_NA":  na_se,
        "mean_CL":  cl_m,  "sem_CL":  cl_se,
        "mean_qN_e": qN_m, "sem_qN_e": qN_se,
        "R": SOL.shape[1],
        "dt_ps": dt_ps,
    })
    out_tsv = outdir / "reduced_counts.tsv"
    out.to_csv(out_tsv, sep="\t", index=False, float_format="%.8e")

    meta = dict(
        compile_dir=str(indir),
        outdir=str(outdir),
        baseline_zero=baseline_zero,
        R=int(SOL.shape[1]),
        dt_ps=dt_ps,
        t_min_ps=float(grid_ps[0]),
        t_max_ps=float(grid_ps[-1]),
        used_files=used,
        interp="step-hold previous sample",
        columns_in="time_ps SOL NA CL",
        columns_out=list(out.columns),
    )
    (outdir / "reduced_counts_meta.json").write_text(json.dumps(meta, indent=2))

    # plots
    x = grid_ps * 1e-3
    def plot_one(name, Y, mean, sem):
        plt.figure(figsize=(7.2, 4.8))
        if plot_replicas:
            cmap = plt.cm.tab20
            for j in range(Y.shape[1]):
                plt.plot(x, Y[:, j], color=cmap(j % 20), alpha=0.5, linewidth=1)
        plt.plot(x, mean, color="black", linewidth=2.5, label=f"mean {name}")
        plt.fill_between(x, mean-sem, mean+sem, alpha=0.25, label="±SEM")
        plt.axhline(0, linestyle="--", linewidth=0.8)
        plt.xlabel("time (ns)")
        plt.ylabel(f"{name} cumulative net crossings")
        plt.legend()
        plt.tight_layout()
        pdf = outdir / f"{name}_mean_sem.pdf"
        plt.savefig(pdf, bbox_inches="tight")
        plt.close()
        return pdf

    plot_one("SOL", SOL, sol_m, sol_se)
    plot_one("NA",  NA,  na_m,  na_se)
    plot_one("CL",  CL,  cl_m,  cl_se)
    plot_one("qN_e", qN, qN_m, qN_se)

    return out_tsv

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("compile_dir", type=str, help=".../FIELD_XX/compile")
    ap.add_argument("--no-replica-lines", action="store_true")
    args = ap.parse_args()
    reduce_compile_dir(Path(args.compile_dir), plot_replicas=not args.no_replica_lines)
    print("OK")
