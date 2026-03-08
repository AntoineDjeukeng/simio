#!/usr/bin/env python3
from pathlib import Path
import json
import re
import numpy as np
import pandas as pd

from config_flow import (
    FIELD_TO_E,
    FIT_LAST_NS_DEFAULT,
    I_END,
    I_START,
    LY_NM,
    PATTERN,
    RUNS_TRANSPORT_INPUT_GLOB,
    SIGMA_E_MAX_V_NM,
    TRANSPORT_FITS_DIRNAME,
)
from transport_io import infer_field_dir, read_replica, resolve_replica_files

E_CHARGE_C = 1.602176634e-19  # C
NM_TO_M = 1e-9


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


def fit_slope(t_ns, y):
    a, _b = np.polyfit(t_ns, y, 1)
    return float(a)


def current_from_slope(slope_per_ns, slope_sem_per_ns, charge_sign):
    i_mean = charge_sign * E_CHARGE_C * slope_per_ns * 1e9
    i_sem = charge_sign * E_CHARGE_C * slope_sem_per_ns * 1e9
    return float(i_mean), float(abs(i_sem))


def conductivity_from_current(i_mean_a, i_sem_a, e_v_nm, area_m2):
    if e_v_nm <= 0:
        return float("nan"), float("nan")
    if SIGMA_E_MAX_V_NM is not None and e_v_nm > SIGMA_E_MAX_V_NM:
        return float("nan"), float("nan")
    e_v_m = e_v_nm * 1e9
    sigma_mean = i_mean_a / (e_v_m * area_m2)
    sigma_sem = i_sem_a / (e_v_m * area_m2)
    return float(sigma_mean), float(sigma_sem)


def parse_condition(input_root: Path):
    parts = input_root.resolve().parts
    field = next((p for p in parts if p.startswith("FIELD_")), None)
    charge = next((p for p in parts if p in ("neutral", "pos", "neg", "positive", "negative")), None)
    L_nm = None
    H_nm = None
    for p in parts:
        if re.fullmatch(r"L_\d+", p):
            L_nm = float(p.split("_")[1])
        if re.fullmatch(r"single_H_\d+", p):
            H_nm = float(p.split("_")[2]) / 10.0
    if field is None or charge is None or L_nm is None or H_nm is None:
        raise ValueError(f"Cannot parse condition from path: {input_root}")
    if field not in FIELD_TO_E:
        raise ValueError(f"Unknown field tag {field}")
    return dict(field=field, E_V_nm=float(FIELD_TO_E[field]), charge=charge, L_nm=L_nm, H_nm=H_nm)


def load_reps(
    input_root: Path,
    replica_files=None,
    replica_glob=RUNS_TRANSPORT_INPUT_GLOB,
    pattern=PATTERN,
    i_start=I_START,
    i_end=I_END,
):
    reps, used = [], []
    files = resolve_replica_files(
        input_root,
        i_start=i_start,
        i_end=i_end,
        pattern=pattern,
        replica_files=replica_files,
        replica_glob=replica_glob,
    )
    for f in files:
        reps.append(read_replica(f, baseline_zero=True))
        used.append(str(f))
    if not reps:
        raise RuntimeError(f"No replicas found in {input_root}")
    return reps, used


def align(reps):
    dts = []
    for (t, *_r) in reps:
        dt = robust_dt_ps(t)
        if np.isfinite(dt) and dt > 0:
            dts.append(dt)
    if not dts:
        raise RuntimeError("Cannot estimate dt.")
    dt_ps = float(np.median(dts))
    t_min = max(t[0] for (t, *_r) in reps)
    t_max = min(t[-1] for (t, *_r) in reps)
    grid_ps = make_grid_ps(t_min, t_max, dt_ps)
    SOL = np.vstack([step_hold_interp(grid_ps, t, sol) for (t, sol, na, cl) in reps]).T
    NA = np.vstack([step_hold_interp(grid_ps, t, na) for (t, sol, na, cl) in reps]).T
    CL = np.vstack([step_hold_interp(grid_ps, t, cl) for (t, sol, na, cl) in reps]).T
    return grid_ps, dt_ps, SOL, NA, CL


def slope_stats(t_ns, Y, fit_last_ns):
    t1 = t_ns[-1]
    t0 = max(t_ns[0], t1 - fit_last_ns)
    mask = t_ns >= t0
    n_points = int(mask.sum())
    if n_points < 5:
        raise RuntimeError("Not enough points in fit window.")
    slopes = np.array([fit_slope(t_ns[mask], Y[mask, j]) for j in range(Y.shape[1])], float)
    mean = float(slopes.mean())
    sd = float(slopes.std(ddof=1)) if slopes.size > 1 else 0.0
    sem = float(sd / np.sqrt(slopes.size)) if slopes.size > 1 else 0.0
    return dict(
        mean=mean,
        sem=sem,
        n=int(slopes.size),
        n_points=n_points,
        t0=float(t0),
        t1=float(t1),
        slopes=slopes,
    )


def process_input_root(
    input_root: Path,
    fit_last_ns: float,
    replica_files=None,
    replica_glob=RUNS_TRANSPORT_INPUT_GLOB,
    pattern=PATTERN,
    i_start=I_START,
    i_end=I_END,
):
    if LY_NM <= 0:
        raise SystemExit("Set LY_NM > 0 in config_flow.py to compute conductivity.")

    cond = parse_condition(input_root)
    reps, used = load_reps(
        input_root,
        replica_files=replica_files,
        replica_glob=replica_glob,
        pattern=pattern,
        i_start=i_start,
        i_end=i_end,
    )
    grid_ps, dt_ps, SOL, NA, CL = align(reps)
    t_ns = grid_ps * 1e-3

    qN = NA - CL  # e
    st_qn = slope_stats(t_ns, qN, fit_last_ns)
    st_na = slope_stats(t_ns, NA, fit_last_ns)
    st_cl = slope_stats(t_ns, CL, fit_last_ns)

    slope_qn_mean = st_qn["mean"]
    slope_qn_sem = st_qn["sem"]
    slope_na_mean = st_na["mean"]
    slope_na_sem = st_na["sem"]
    slope_cl_mean = st_cl["mean"]
    slope_cl_sem = st_cl["sem"]

    # Currents from per-species slopes
    I_mean_A, I_sem_A = current_from_slope(slope_qn_mean, slope_qn_sem, charge_sign=+1.0)
    I_NA_mean_A, I_NA_sem_A = current_from_slope(slope_na_mean, slope_na_sem, charge_sign=+1.0)
    I_CL_mean_A, I_CL_sem_A = current_from_slope(slope_cl_mean, slope_cl_sem, charge_sign=-1.0)

    # Geometry
    S_nm2 = cond["H_nm"] * LY_NM
    S_m2 = S_nm2 * (NM_TO_M**2)

    # Current density
    j_mean_A_m2 = I_mean_A / S_m2
    j_sem_A_m2 = I_sem_A / S_m2
    j_NA_mean_A_m2 = I_NA_mean_A / S_m2
    j_NA_sem_A_m2 = I_NA_sem_A / S_m2
    j_CL_mean_A_m2 = I_CL_mean_A / S_m2
    j_CL_sem_A_m2 = I_CL_sem_A / S_m2

    # Conductivity (undefined at E=0 and optionally restricted to linear regime)
    E_V_nm = cond["E_V_nm"]
    sigma_mean_S_m, sigma_sem_S_m = conductivity_from_current(I_mean_A, I_sem_A, E_V_nm, S_m2)
    sigma_NA_mean_S_m, sigma_NA_sem_S_m = conductivity_from_current(
        I_NA_mean_A, I_NA_sem_A, E_V_nm, S_m2
    )
    sigma_CL_mean_S_m, sigma_CL_sem_S_m = conductivity_from_current(
        I_CL_mean_A, I_CL_sem_A, E_V_nm, S_m2
    )

    sigma_in_linear_regime = int(
        (E_V_nm > 0) and (SIGMA_E_MAX_V_NM is None or E_V_nm <= SIGMA_E_MAX_V_NM)
    )

    row = {
        **cond,
        "Ly_nm": float(LY_NM),
        "S_nm2": float(S_nm2),
        "R": int(len(reps)),
        "dt_ps": float(dt_ps),
        "fit_last_ns": float(fit_last_ns),
        "fit_t0_ns": float(st_qn["t0"]),
        "fit_t1_ns": float(st_qn["t1"]),
        "fit_n_points": int(st_qn["n_points"]),
        "slope_qN_mean_e_per_ns": float(slope_qn_mean),
        "slope_qN_sem_e_per_ns": float(slope_qn_sem),
        "slope_NA_mean_per_ns": float(slope_na_mean),
        "slope_NA_sem_per_ns": float(slope_na_sem),
        "slope_CL_mean_per_ns": float(slope_cl_mean),
        "slope_CL_sem_per_ns": float(slope_cl_sem),
        "slopes_qN_rep_e_per_ns": json.dumps(st_qn["slopes"].tolist()),
        "slopes_NA_rep_per_ns": json.dumps(st_na["slopes"].tolist()),
        "slopes_CL_rep_per_ns": json.dumps(st_cl["slopes"].tolist()),
        "I_mean_A": float(I_mean_A),
        "I_sem_A": float(I_sem_A),
        "I_NA_mean_A": float(I_NA_mean_A),
        "I_NA_sem_A": float(I_NA_sem_A),
        "I_CL_mean_A": float(I_CL_mean_A),
        "I_CL_sem_A": float(I_CL_sem_A),
        "j_mean_A_m2": float(j_mean_A_m2),
        "j_sem_A_m2": float(j_sem_A_m2),
        "j_NA_mean_A_m2": float(j_NA_mean_A_m2),
        "j_NA_sem_A_m2": float(j_NA_sem_A_m2),
        "j_CL_mean_A_m2": float(j_CL_mean_A_m2),
        "j_CL_sem_A_m2": float(j_CL_sem_A_m2),
        "sigma_mean_S_m": float(sigma_mean_S_m),
        "sigma_sem_S_m": float(sigma_sem_S_m),
        "sigma_NA_mean_S_m": float(sigma_NA_mean_S_m),
        "sigma_NA_sem_S_m": float(sigma_NA_sem_S_m),
        "sigma_CL_mean_S_m": float(sigma_CL_mean_S_m),
        "sigma_CL_sem_S_m": float(sigma_CL_sem_S_m),
        "sigma_in_linear_regime": int(sigma_in_linear_regime),
        "used_files_n": int(len(used)),
    }
    return row


# Backward compatibility helper for old callers.
def process_compile_dir(compile_dir: Path, fit_last_ns: float):
    return process_input_root(
        compile_dir,
        fit_last_ns,
        replica_glob=None,
        pattern=PATTERN,
    )


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("input_root", type=str, help="Input root (normal: .../FIELD_XX/runs)")
    ap.add_argument("--fit-last-ns", type=float, default=FIT_LAST_NS_DEFAULT)
    ap.add_argument(
        "--out",
        type=str,
        default=None,
        help="Output TSV (default: FIELD_XX/transport_fits/transport.tsv)",
    )
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
    ap.add_argument("--pattern", type=str, default=PATTERN, help="Legacy pattern for indexed replica files.")
    ap.add_argument("--i-start", type=int, default=I_START)
    ap.add_argument("--i-end", type=int, default=I_END)
    args = ap.parse_args()

    input_root = Path(args.input_root)
    replica_glob = None if args.legacy_indexed_counts else args.replica_glob
    row = process_input_root(
        input_root,
        args.fit_last_ns,
        replica_files=[Path(p) for p in args.replica_file],
        replica_glob=replica_glob,
        pattern=args.pattern,
        i_start=args.i_start,
        i_end=args.i_end,
    )

    if args.out:
        out = Path(args.out)
    else:
        out = infer_field_dir(input_root) / TRANSPORT_FITS_DIRNAME / "transport.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row]).to_csv(out, sep="\t", index=False)
    print("Wrote:", out)

