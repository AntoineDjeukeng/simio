#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from config_flow import (
    I_END,
    I_START,
    PATTERN,
    RUNS_TRANSPORT_INPUT_GLOB,
    TRANSPORT_REDUCTION_DIRNAME,
)
from transport_io import infer_field_dir, read_replica, resolve_replica_files


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


def reduce_input_root(
    input_root: Path,
    outdir: Path | None = None,
    baseline_zero=True,
    plot_replicas=True,
    replica_files=None,
    replica_glob=RUNS_TRANSPORT_INPUT_GLOB,
    pattern=PATTERN,
    i_start=I_START,
    i_end=I_END,
):
    input_root = Path(input_root)
    if outdir is None:
        outdir = infer_field_dir(input_root) / TRANSPORT_REDUCTION_DIRNAME
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    files = resolve_replica_files(
        input_root,
        i_start=i_start,
        i_end=i_end,
        pattern=pattern,
        replica_files=replica_files,
        replica_glob=replica_glob,
    )
    reps, used = [], []
    for f in files:
        reps.append(read_replica(f, baseline_zero=baseline_zero))
        used.append(str(f))
    if not reps:
        raise RuntimeError(f"No replicas found in {input_root}")

    # common time grid
    dts = [robust_dt_ps(t) for (t, *_r) in reps]
    dts = [dt for dt in dts if np.isfinite(dt) and dt > 0]
    dt_ps = float(np.median(dts))

    t_min = max(t[0] for (t, *_r) in reps)
    t_max = min(t[-1] for (t, *_r) in reps)
    grid_ps = make_grid_ps(t_min, t_max, dt_ps)

    SOL = np.vstack([step_hold_interp(grid_ps, t, sol) for (t, sol, na, cl) in reps]).T
    NA = np.vstack([step_hold_interp(grid_ps, t, na) for (t, sol, na, cl) in reps]).T
    CL = np.vstack([step_hold_interp(grid_ps, t, cl) for (t, sol, na, cl) in reps]).T
    qN = NA - CL

    sol_m, sol_se = mean_sem(SOL)
    na_m, na_se = mean_sem(NA)
    cl_m, cl_se = mean_sem(CL)
    qN_m, qN_se = mean_sem(qN)

    out = pd.DataFrame(
        {
            "time_ps": grid_ps,
            "mean_SOL": sol_m,
            "sem_SOL": sol_se,
            "mean_NA": na_m,
            "sem_NA": na_se,
            "mean_CL": cl_m,
            "sem_CL": cl_se,
            "mean_qN_e": qN_m,
            "sem_qN_e": qN_se,
            "R": SOL.shape[1],
            "dt_ps": dt_ps,
        }
    )
    out_tsv = outdir / "reduced_counts.tsv"
    out.to_csv(out_tsv, sep="\t", index=False, float_format="%.8e")

    meta = dict(
        input_root=str(input_root),
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
        plt.fill_between(x, mean - sem, mean + sem, alpha=0.25, label="±SEM")
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
    plot_one("NA", NA, na_m, na_se)
    plot_one("CL", CL, cl_m, cl_se)
    plot_one("qN_e", qN, qN_m, qN_se)

    return out_tsv


# Backward compatibility helper for old callers.
def reduce_compile_dir(compile_dir: Path, baseline_zero=True, plot_replicas=True):
    return reduce_input_root(
        compile_dir,
        baseline_zero=baseline_zero,
        plot_replicas=plot_replicas,
        replica_glob=None,
        pattern=PATTERN,
    )


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument(
        "input_root",
        type=str,
        help="Input root (normal: .../FIELD_XX/runs)",
    )
    ap.add_argument(
        "--outdir",
        type=str,
        default=None,
        help="Output directory (default: FIELD_XX/transport_reduction)",
    )
    ap.add_argument("--no-replica-lines", action="store_true")
    ap.add_argument(
        "--replica-file",
        action="append",
        default=[],
        help="Explicit replica input file. Repeat for multiple.",
    )
    ap.add_argument(
        "--replica-glob",
        type=str,
        default=RUNS_TRANSPORT_INPUT_GLOB,
        help=f"Glob relative to input_root (default: {RUNS_TRANSPORT_INPUT_GLOB})",
    )
    ap.add_argument(
        "--legacy-indexed-counts",
        action="store_true",
        help="Use legacy indexed files (count_XX.dat) via --pattern and index range.",
    )
    ap.add_argument(
        "--pattern",
        type=str,
        default=PATTERN,
        help="Legacy pattern for indexed replica files.",
    )
    ap.add_argument("--i-start", type=int, default=I_START)
    ap.add_argument("--i-end", type=int, default=I_END)
    args = ap.parse_args()

    replica_glob = None if args.legacy_indexed_counts else args.replica_glob
    reduce_input_root(
        Path(args.input_root),
        outdir=Path(args.outdir) if args.outdir else None,
        plot_replicas=not args.no_replica_lines,
        replica_files=[Path(p) for p in args.replica_file],
        replica_glob=replica_glob,
        pattern=args.pattern,
        i_start=args.i_start,
        i_end=args.i_end,
    )
    print("OK")

