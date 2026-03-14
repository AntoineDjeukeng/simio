#!/usr/bin/env python3
from pathlib import Path
import json
import re
import numpy as np
import pandas as pd
from config_flow import (
    FIELD_TO_E,
    LY_NM,
    FIT_LAST_NS_DEFAULT,
    PATTERN,
    I_START,
    I_END,
    SIGMA_E_MAX_V_NM,
)

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
    a, b = np.polyfit(t_ns, y, 1)
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

def parse_condition(compile_dir: Path):
    parts = compile_dir.resolve().parts
    field = next((p for p in parts if p.startswith("FIELD_")), None)
    charge = next((p for p in parts if p in ("neutral","pos","neg","positive","negative")), None)
    L_nm = None
    H_nm = None
    for p in parts:
        if re.fullmatch(r"L_\d+", p):
            L_nm = float(p.split("_")[1])
        if re.fullmatch(r"single_H_\d+", p):
            H_nm = float(p.split("_")[2]) / 10.0
    if field is None or charge is None or L_nm is None or H_nm is None:
        raise ValueError(f"Cannot parse condition from path: {compile_dir}")
    if field not in FIELD_TO_E:
        raise ValueError(f"Unknown field tag {field}")
    return dict(field=field, E_V_nm=float(FIELD_TO_E[field]), charge=charge, L_nm=L_nm, H_nm=H_nm)

def read_replica(path: Path, baseline_zero=True):
    df = pd.read_csv(path, sep=r"\s+", comment="#", header=None, engine="python")
    df = df.iloc[:, :4].copy()
    df.columns = ["time_ps","SOL","NA","CL"]
    df = df.sort_values("time_ps").reset_index(drop=True)
    t = df["time_ps"].to_numpy(float)
    SOL = df["SOL"].to_numpy(float)
    NA  = df["NA"].to_numpy(float)
    CL  = df["CL"].to_numpy(float)
    if baseline_zero:
        SOL -= SOL[0]; NA -= NA[0]; CL -= CL[0]
    return t, SOL, NA, CL

def load_reps(compile_dir: Path):
    reps, used = [], []
    for i in range(I_START, I_END + 1):
        f = compile_dir / PATTERN.format(i=i)
        if not f.exists():
            continue
        reps.append(read_replica(f, baseline_zero=True))
        used.append(str(f))
    if not reps:
        raise RuntimeError(f"No replicas found in {compile_dir}")
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
    NA  = np.vstack([step_hold_interp(grid_ps, t, na)  for (t, sol, na, cl) in reps]).T
    CL  = np.vstack([step_hold_interp(grid_ps, t, cl)  for (t, sol, na, cl) in reps]).T
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

def process(compile_dir: Path, fit_last_ns: float):
    if LY_NM <= 0:
        raise SystemExit("Set LY_NM > 0 in config_flow.py to compute conductivity.")

    cond = parse_condition(compile_dir)
    reps, used = load_reps(compile_dir)
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
    j_sem_A_m2  = I_sem_A  / S_m2
    j_NA_mean_A_m2 = I_NA_mean_A / S_m2
    j_NA_sem_A_m2  = I_NA_sem_A  / S_m2
    j_CL_mean_A_m2 = I_CL_mean_A / S_m2
    j_CL_sem_A_m2  = I_CL_sem_A  / S_m2

    # Conductivity (undefined at E=0 and optionally restricted to linear regime)
    E_V_nm = cond["E_V_nm"]
    sigma_mean_S_m, sigma_sem_S_m = conductivity_from_current(I_mean_A, I_sem_A, E_V_nm, S_m2)
    sigma_NA_mean_S_m, sigma_NA_sem_S_m = conductivity_from_current(I_NA_mean_A, I_NA_sem_A, E_V_nm, S_m2)
    sigma_CL_mean_S_m, sigma_CL_sem_S_m = conductivity_from_current(I_CL_mean_A, I_CL_sem_A, E_V_nm, S_m2)

    sigma_in_linear_regime = int(
        (E_V_nm > 0)
        and (SIGMA_E_MAX_V_NM is None or E_V_nm <= SIGMA_E_MAX_V_NM)
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

def process_compile_dir(compile_dir: Path, fit_last_ns: float):
    return process(compile_dir, fit_last_ns)

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("compile_dir", type=str, help=".../FIELD_XX/compile")
    ap.add_argument("--fit-last-ns", type=float, default=FIT_LAST_NS_DEFAULT)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    compile_dir = Path(args.compile_dir)
    row = process(compile_dir, args.fit_last_ns)

    out = Path(args.out) if args.out else (compile_dir / "reduced" / "transport.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row]).to_csv(out, sep="\t", index=False)
    print("Wrote:", out)
