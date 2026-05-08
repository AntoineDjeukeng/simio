#!/usr/bin/env python3
import argparse
import csv
import json
import math
import os
import re
import sys
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


SPECIES = [
    ("SOL", "water"),
    ("NA", "na"),
    ("CL", "cl"),
]

PMF_IONS = ["NA", "CL"]
PMF_COLUMNS = ["deltaF_center", "barrier_mean_center", "barrier_asym_center"]
DEFAULT_PMF_CSV = "/home/antoine/CDI/plots/New_pmf_60/pmf_summary_compact.csv"

FIELD_TO_E = {
    "FIELD_00": 0.0,
    "FIELD_01": 0.004,
    "FIELD_02": 0.01,
    "FIELD_03": 0.03,
}

REGIONS = [
    "res_left",
    "mouth_left",
    "center",
    "mouth_right",
    "res_right",
]

COORD_COLUMNS = [
    ("IBC", "IBC_mean"),
    ("IBA", "IBA_mean"),
    ("BNC", "BNC_mean"),
    ("BNA", "BNA_mean"),
    ("BWW", "BWW_mean"),
    ("HBww_don", "HBww_don_mean"),
    ("HBww_acc", "HBww_acc_mean"),
    ("HBww_tot", "HBww_tot_mean"),
    ("HBwcl_don", "HBwcl_don_mean"),
]

NACL_CLUSTER_COLUMNS = [
    ("cluster_count", "nacl_cluster_count_mean"),
    ("ions_in_clusters", "ions_in_nacl_clusters_mean"),
    ("na_in_clusters", "na_in_nacl_clusters_mean"),
    ("cl_in_clusters", "cl_in_nacl_clusters_mean"),
    ("cluster_size", "nacl_cluster_size_mean"),
]

UPSERT_KEYS = ["H", "L", "surface_charge_q", "charge", "field", "replica"]


def read_csv_dicts(path: str) -> List[Dict[str, str]]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def to_float(v: str) -> float:
    if v is None or v == "":
        return math.nan
    return float(v)


def mean(values: Iterable[float]) -> float:
    xs = [v for v in values if not math.isnan(v)]
    if not xs:
        return math.nan
    return sum(xs) / float(len(xs))


def safe_ratio(num: object, den: object) -> float:
    try:
        n = float(num)
        d = float(den)
    except (TypeError, ValueError):
        return math.nan
    if math.isnan(n) or math.isnan(d) or d == 0.0:
        return math.nan
    return n / d


def same_key(a: object, b: object) -> bool:
    sa = str(a).strip()
    sb = str(b).strip()
    try:
        return int(float(sa)) == int(float(sb))
    except (TypeError, ValueError):
        return sa == sb


def charge_strength_from_path(path: str) -> int:
    parts = os.path.abspath(path).split(os.sep)
    for p in parts:
        m = re.fullmatch(r"L[_-][0-9.]+(?:[_-]([24])q)?", p, re.IGNORECASE)
        if m:
            return int(m.group(1)) if m.group(1) else 1
    return 1


def signed_surface_charge(charge: object, strength: int) -> int:
    c = str(charge).strip().lower()
    if c in ("negative", "neg", "-"):
        return -strength
    if c in ("positive", "pos", "+"):
        return strength
    if c in ("neutral", "neu", "0"):
        return 0
    return 0


def pmf_variant_q_from_surface_charge(surface_charge_q: object) -> Tuple[str, int]:
    try:
        qabs = abs(int(float(str(surface_charge_q).strip())))
    except (TypeError, ValueError):
        return "base", 0
    if qabs == 2:
        return "2q", 2
    if qabs == 4:
        return "4q", 4
    return "base", 0


def linear_slope(xs: Sequence[float], ys: Sequence[float]) -> float:
    good = [(x, y) for x, y in zip(xs, ys) if not (math.isnan(x) or math.isnan(y))]
    if len(good) < 2:
        return math.nan
    xbar = sum(x for x, _ in good) / float(len(good))
    ybar = sum(y for _, y in good) / float(len(good))
    den = sum((x - xbar) * (x - xbar) for x, _ in good)
    if den == 0.0:
        return math.nan
    num = sum((x - xbar) * (y - ybar) for x, y in good)
    return num / den


def load_json(path: str) -> Dict:
    with open(path) as f:
        return json.load(f)


def resolve_path(path: str, base_dir: str) -> str:
    if os.path.isabs(path):
        return path
    direct = os.path.abspath(path)
    if os.path.exists(direct):
        return direct
    return os.path.abspath(os.path.join(base_dir, path))


def geometry_from_config(config_path: str) -> Dict[str, float]:
    cfg = load_json(config_path)
    config_dir = os.path.dirname(os.path.abspath(config_path))

    topo_path = cfg.get("topology_json") or cfg.get("mini_topology_json") or cfg.get("setup_json")
    if topo_path:
        topo_abs = resolve_path(str(topo_path), config_dir)
        if os.path.exists(topo_abs):
            topo = load_json(topo_abs)
            ch = topo.get("channel")
            if isinstance(ch, dict) and "min" in ch and "max" in ch:
                out = {
                    "x_left_carbon_nm": float(ch["min"][0]),
                    "x_right_carbon_nm": float(ch["max"][0]),
                    "channel_y_width_nm": float(ch["max"][1]) - float(ch["min"][1]),
                    "channel_height_nm": float(ch["max"][2]) - float(ch["min"][2]),
                    "box_lx_nm": math.nan,
                    "box_ly_nm": math.nan,
                    "box_lz_nm": math.nan,
                }
                gro_box = topo.get("gro_box")
                if isinstance(gro_box, list):
                    if len(gro_box) > 0:
                        out["box_lx_nm"] = float(gro_box[0])
                    if len(gro_box) > 1:
                        out["box_ly_nm"] = float(gro_box[1])
                    if len(gro_box) > 2:
                        out["box_lz_nm"] = float(gro_box[2])
                return out

    if "xmin" in cfg and "xmax" in cfg:
        out = {
            "x_left_carbon_nm": float(cfg["xmin"]),
            "x_right_carbon_nm": float(cfg["xmax"]),
            "channel_y_width_nm": math.nan,
            "channel_height_nm": math.nan,
            "box_lx_nm": math.nan,
            "box_ly_nm": math.nan,
            "box_lz_nm": math.nan,
        }
        for key in ("box_lx", "box_lx_nm", "lx", "Lx"):
            if key in cfg:
                out["box_lx_nm"] = float(cfg[key])
                break
        return out

    raise ValueError(
        "Could not determine carbon/channel x edges. Provide topology_json with channel.min/max "
        "or xmin/xmax in the config."
    )


def infer_box_lx_from_density(results_dir: str) -> float:
    path = os.path.join(results_dir, "density_x.csv")
    if not os.path.exists(path):
        return math.nan

    rows = read_csv_dicts(path)
    xs = sorted(to_float(row.get("x_center_nm", "")) for row in rows)
    xs = [x for x in xs if not math.isnan(x)]
    if len(xs) < 2:
        return math.nan

    dxs = [b - a for a, b in zip(xs, xs[1:]) if b > a]
    if not dxs:
        return math.nan
    dx = min(dxs)
    return xs[-1] + 0.5 * dx


def field_to_e_v_per_nm(field: object) -> object:
    text = str(field).strip()
    if text == "":
        return ""
    key = text.upper()
    if re.fullmatch(r"\d+", key):
        key = f"FIELD_{int(key):02d}"
    elif not key.startswith("FIELD_"):
        key = f"FIELD_{key}"
    return FIELD_TO_E.get(key, "NaN")


def region_for_x(x: float, x_left: float, x_right: float, mouth_width_nm: float,
                 box_lx_nm: float, reservoir_fraction: float) -> Optional[str]:
    if x_left <= x < x_right:
        return "center"
    if x_left - mouth_width_nm <= x < x_left:
        return "mouth_left"
    if x_right <= x < x_right + mouth_width_nm:
        return "mouth_right"
    if not math.isnan(box_lx_nm) and x < reservoir_fraction * box_lx_nm:
        return "res_left"
    if not math.isnan(box_lx_nm) and x >= (1.0 - reservoir_fraction) * box_lx_nm:
        return "res_right"
    return None


def add_density_features(out: Dict[str, object], results_dir: str, x_left: float, x_right: float,
                         mouth_width_nm: float, box_lx_nm: float, reservoir_fraction: float) -> None:
    path = os.path.join(results_dir, "density_x.csv")
    if not os.path.exists(path):
        return
    rows = read_csv_dicts(path)
    buckets: Dict[Tuple[str, str], List[float]] = {}
    for row in rows:
        x = to_float(row.get("x_center_nm", ""))
        region = region_for_x(x, x_left, x_right, mouth_width_nm, box_lx_nm, reservoir_fraction)
        if region is None:
            continue
        for label, internal in SPECIES:
            col = f"rho_{internal}_mean"
            if col in row:
                buckets.setdefault((label, region), []).append(to_float(row[col]))

    for label, _ in SPECIES:
        for region in REGIONS:
            out[f"{label}_{region}_mean"] = mean(buckets.get((label, region), []))

        res_left = out.get(f"{label}_res_left_mean", math.nan)
        res_right = out.get(f"{label}_res_right_mean", math.nan)
        mouth_left = out.get(f"{label}_mouth_left_mean", math.nan)
        mouth_right = out.get(f"{label}_mouth_right_mean", math.nan)
        center = out.get(f"{label}_center_mean", math.nan)
        res_avg = (float(res_left) + float(res_right)) / 2.0

        out[f"{label}_res_asym"] = float(res_right) - float(res_left)
        out[f"{label}_mouth_asym"] = float(mouth_right) - float(mouth_left)
        out[f"{label}_center_enrichment"] = float(center) - res_avg


def add_profile_region_features(out: Dict[str, object], path: str, prefix: str,
                                columns: Sequence[Tuple[str, str]], x_left: float, x_right: float,
                                mouth_width_nm: float, box_lx_nm: float,
                                reservoir_fraction: float) -> None:
    if not os.path.exists(path):
        return

    rows = read_csv_dicts(path)
    buckets: Dict[Tuple[str, str], List[float]] = {}
    for row in rows:
        x = to_float(row.get("x_center_nm", ""))
        region = region_for_x(x, x_left, x_right, mouth_width_nm, box_lx_nm, reservoir_fraction)
        if region is None:
            continue
        for feature, col in columns:
            if col in row:
                buckets.setdefault((feature, region), []).append(to_float(row[col]))

    for feature, _ in columns:
        for region in REGIONS:
            out[f"{prefix}_{feature}_{region}_mean"] = mean(buckets.get((feature, region), []))


def add_coord_features(out: Dict[str, object], results_dir: str, x_left: float, x_right: float,
                       mouth_width_nm: float, box_lx_nm: float, reservoir_fraction: float) -> None:
    add_profile_region_features(out, os.path.join(results_dir, "coord_x.csv"), "COORD",
                                COORD_COLUMNS, x_left, x_right, mouth_width_nm,
                                box_lx_nm, reservoir_fraction)


def add_nacl_cluster_features(out: Dict[str, object], results_dir: str, x_left: float,
                              x_right: float, mouth_width_nm: float, box_lx_nm: float,
                              reservoir_fraction: float) -> None:
    add_profile_region_features(out, os.path.join(results_dir, "coord_x.csv"), "NACL",
                                NACL_CLUSTER_COLUMNS, x_left, x_right, mouth_width_nm,
                                box_lx_nm, reservoir_fraction)


def pmf_l_from_replica_l(value: object) -> str:
    try:
        return str(int(float(str(value).strip())) * 10)
    except (TypeError, ValueError):
        return str(value).strip()


def pmf_state_from_charge(value: object) -> str:
    return {
        "negative": "neg",
        "neutral": "neu",
        "positive": "pos",
        "neg": "neg",
        "neu": "neu",
        "pos": "pos",
    }.get(str(value).strip().lower(), str(value).strip())


def init_pmf_features(out: Dict[str, object]) -> None:
    for ion in PMF_IONS:
        for col in PMF_COLUMNS:
            out[f"PMF_{ion}_{col}"] = math.nan


def add_pmf_features(out: Dict[str, object], pmf_csv: str) -> None:
    init_pmf_features(out)
    if not pmf_csv or not os.path.exists(pmf_csv):
        return

    rows = read_csv_dicts(pmf_csv)
    if not rows:
        return

    target_h = out.get("H", "")
    target_l = pmf_l_from_replica_l(out.get("L", ""))
    target_state = pmf_state_from_charge(out.get("charge", ""))
    target_variant, target_q = pmf_variant_q_from_surface_charge(out.get("surface_charge_q", 0))

    for ion in PMF_IONS:
        match = None
        for row in rows:
            if str(row.get("ion", "")).strip().upper() != ion:
                continue
            if "variant" in row and str(row.get("variant", "")).strip() not in ("", target_variant):
                continue
            if "q" in row and not same_key(row.get("q", ""), target_q):
                continue
            if not same_key(row.get("H", ""), target_h):
                continue
            if not same_key(row.get("L", ""), target_l):
                continue
            if str(row.get("state", "")).strip().lower() != target_state:
                continue
            match = row
            break

        if match is None:
            continue
        for col in PMF_COLUMNS:
            out[f"PMF_{ion}_{col}"] = to_float(match.get(col, ""))


def cumulative(values: Sequence[float]) -> List[float]:
    out: List[float] = []
    s = 0.0
    for v in values:
        s += v
        out.append(s)
    return out


def warn_gating_mismatch(path: str, col: str, max_delta: float) -> None:
    print(f"warning: {path}: {col} differs from cumulative per-frame counts; max delta={max_delta:g}",
          file=sys.stderr)


def validate_cumulative_column(path: str, rows: List[Dict[str, str]], event_col: str,
                               cum_col: str, tolerance: float = 1e-9) -> None:
    if cum_col not in rows[0]:
        return
    built = cumulative([to_float(r[event_col]) for r in rows])
    reported = [to_float(r[cum_col]) for r in rows]
    max_delta = max(abs(a - b) for a, b in zip(built, reported))
    if max_delta > tolerance:
        warn_gating_mismatch(path, cum_col, max_delta)


def add_gating_features(out: Dict[str, object], results_dir: str, tail_ns: float,
                        validate_gating: bool = False) -> None:
    path = os.path.join(results_dir, "gating_flux.csv")
    if not os.path.exists(path):
        return
    rows = read_csv_dicts(path)
    if not rows:
        return

    time_ns = [to_float(r.get("time_ps", "")) / 1000.0 for r in rows]
    t_end = max(t for t in time_ns if not math.isnan(t))
    tail_start = t_end - tail_ns
    tail_idx = [i for i, t in enumerate(time_ns) if not math.isnan(t) and t >= tail_start]

    for label, internal in SPECIES:
        left_col = f"center_{internal}_left"
        right_col = f"center_{internal}_right"
        cum_left_col = f"center_{internal}_cum_left"
        cum_right_col = f"center_{internal}_cum_right"

        if left_col not in rows[0] or right_col not in rows[0]:
            continue

        if validate_gating:
            validate_cumulative_column(path, rows, left_col, cum_left_col)
            validate_cumulative_column(path, rows, right_col, cum_right_col)

        if cum_left_col in rows[0]:
            left_cum = [to_float(r[cum_left_col]) for r in rows]
        else:
            left_cum = cumulative([to_float(r[left_col]) for r in rows])

        if cum_right_col in rows[0]:
            right_cum = [to_float(r[cum_right_col]) for r in rows]
        else:
            right_cum = cumulative([to_float(r[right_col]) for r in rows])
        net_cum = [r - l for r, l in zip(right_cum, left_cum)]

        out[f"{label}_cum_left_final"] = left_cum[-1]
        out[f"{label}_cum_right_final"] = right_cum[-1]
        out[f"{label}_net_final"] = net_cum[-1]
        out[f"{label}_net_slope"] = linear_slope([time_ns[i] for i in tail_idx],
                                                 [net_cum[i] for i in tail_idx])


def infer_metadata_from_path(path: str) -> Dict[str, str]:
    parts = os.path.abspath(path).split(os.sep)
    meta: Dict[str, str] = {}
    for p in parts:
        m = re.search(r"(?:^|[_-])H[_-]?([0-9.]+)(?:$|[_-])", p, re.IGNORECASE)
        if m:
            meta["H"] = m.group(1)
        m = re.search(r"(?:^|[_-])L[_-]?([0-9.]+)(?:$|[_-])", p, re.IGNORECASE)
        if m:
            meta["L"] = m.group(1)
        m = re.fullmatch(r"FIELD[_-]?([0-9]+)", p, re.IGNORECASE)
        if m:
            meta["field"] = m.group(1)
        m = re.fullmatch(r"rep[_-]?([0-9]+)", p, re.IGNORECASE)
        if m:
            meta["replica"] = m.group(1)
        if p.lower() in ("neutral", "negative", "positive"):
            meta["charge"] = p.lower()

    strength = charge_strength_from_path(path)
    meta["surface_charge_q"] = str(signed_surface_charge(meta.get("charge", ""), strength))

    basename = os.path.basename(path)
    m = re.search(r"_(\d+)\.(?:xtc|trr|gro|pdb)$", basename, re.IGNORECASE)
    if m:
        meta["replica"] = m.group(1)
    return meta


def infer_metadata_from_config(config_path: str) -> Dict[str, str]:
    cfg = load_json(config_path)
    config_dir = os.path.dirname(os.path.abspath(config_path))
    xtc_path = cfg.get("xtc_path") or cfg.get("trajectory") or cfg.get("traj_path")
    if not xtc_path:
        return {}
    return infer_metadata_from_path(resolve_path(str(xtc_path), config_dir))


def ordered_fields(row: Dict[str, object]) -> List[str]:
    fields: List[str] = ["H", "L", "surface_charge_q", "charge", "field", "replica", "E_V_per_nm"]

    fields += [
        "x_left_carbon_nm",
        "x_right_carbon_nm",
        "channel_length_nm",
        "channel_center_nm",
        "channel_y_width_nm",
        "channel_height_nm",
        "box_lx_nm",
        "box_ly_nm",
        "box_lz_nm",
        "mouth_width_nm",
        "reservoir_fraction",
        "reservoir_width_nm",
        "slope_tail_ns",
    ]

    for ion in PMF_IONS:
        for col in PMF_COLUMNS:
            fields.append(f"PMF_{ion}_{col}")

    for label, _ in SPECIES:
        for region in REGIONS:
            fields.append(f"{label}_{region}_mean")

    for label, _ in SPECIES:
        fields += [
            f"{label}_cum_left_final",
            f"{label}_cum_right_final",
            f"{label}_net_final",
            f"{label}_net_slope",
        ]

    for label, _ in SPECIES:
        fields += [
            f"{label}_res_asym",
            f"{label}_mouth_asym",
            f"{label}_center_enrichment",
        ]

    for feature, _ in COORD_COLUMNS:
        for region in REGIONS:
            fields.append(f"COORD_{feature}_{region}_mean")

    for feature, _ in NACL_CLUSTER_COLUMNS:
        for region in REGIONS:
            fields.append(f"NACL_{feature}_{region}_mean")

    return fields


def clean_output_row(fields: Sequence[str], row: Dict[str, object]) -> Dict[str, object]:
    clean_row: Dict[str, object] = {}
    for k in fields:
        v = row.get(k, "NaN")
        if isinstance(v, float) and math.isnan(v):
            v = "NaN"
        clean_row[k] = v
    return clean_row


def upsert_key(row: Dict[str, object]) -> Optional[Tuple[str, ...]]:
    values = tuple(str(row.get(k, "")).strip() for k in UPSERT_KEYS)
    if any(v == "" or v == "NaN" for v in values):
        return None
    return values


def write_or_upsert_row(path: str, row: Dict[str, object]) -> str:
    fields = ordered_fields(row)
    clean_row = clean_output_row(fields, row)
    key = upsert_key(clean_row)
    action = "written"

    rows: List[Dict[str, str]] = []
    if os.path.exists(path) and os.path.getsize(path) > 0:
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            rows = list(reader)

    if rows and key is not None:
        replaced = False
        new_rows: List[Dict[str, object]] = []
        for old in rows:
            old_key = upsert_key(old)
            if old_key == key:
                if not replaced:
                    new_rows.append(clean_row)
                    replaced = True
                    action = "updated"
                continue
            new_rows.append(clean_output_row(fields, old))
        if not replaced:
            new_rows.append(clean_row)
            action = "appended"
        rows = new_rows
    elif rows:
        rows = [clean_output_row(fields, old) for old in rows]
        rows.append(clean_row)
        action = "appended"
    else:
        rows = [clean_row]

    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore", restval="NaN")
        w.writeheader()
        w.writerows(rows)
    return action


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Export one compact structure/transport/association feature row for one replica."
    )
    ap.add_argument("--results-dir", required=True, help="Directory containing density_x.csv and gating_flux.csv")
    ap.add_argument("--config", required=True, help="Run config JSON; carbon/channel edges are read from it")
    ap.add_argument("--out", required=True, help="Output CSV path")
    ap.add_argument("--mouth-width-nm", type=float, default=0.62)
    ap.add_argument("--reservoir-fraction", type=float, default=0.12,
                    help="Reservoir slab width on each box edge, as a fraction of Lx")
    ap.add_argument("--slope-tail-ns", type=float, default=45.0)
    ap.add_argument("--validate-gating", action="store_true",
                    help="Check reported cumulative gating columns against sums of per-frame counts.")
    ap.add_argument("--pmf-csv", default=DEFAULT_PMF_CSV,
                    help="Optional PMF summary CSV. Adds static PMF_NA_* and PMF_CL_* descriptors when matched.")
    ap.add_argument("--H", dest="H", default=None)
    ap.add_argument("--L", dest="L", default=None)
    ap.add_argument("--charge", default=None)
    ap.add_argument("--field", default=None)
    ap.add_argument("--replica", default=None)
    ap.add_argument("--E-V-per-nm", dest="E_V_per_nm", default=None)
    args = ap.parse_args()

    results_dir = os.path.abspath(args.results_dir)
    geom = geometry_from_config(args.config)
    x_left = geom["x_left_carbon_nm"]
    x_right = geom["x_right_carbon_nm"]
    if math.isnan(geom["box_lx_nm"]):
        geom["box_lx_nm"] = infer_box_lx_from_density(results_dir)

    row: Dict[str, object] = {}
    row.update(infer_metadata_from_path(results_dir))
    row.update(infer_metadata_from_config(args.config))
    for key in ("H", "L", "charge", "field", "replica"):
        value = getattr(args, key)
        row[key] = value if value is not None else row.get(key, "")
    strength = charge_strength_from_path(results_dir)
    row["surface_charge_q"] = signed_surface_charge(row.get("charge", ""), strength)
    row["E_V_per_nm"] = args.E_V_per_nm if args.E_V_per_nm is not None else field_to_e_v_per_nm(row.get("field", ""))

    row["x_left_carbon_nm"] = geom["x_left_carbon_nm"]
    row["x_right_carbon_nm"] = geom["x_right_carbon_nm"]
    row["channel_length_nm"] = x_right - x_left
    row["channel_center_nm"] = 0.5 * (x_left + x_right)
    row["channel_y_width_nm"] = geom["channel_y_width_nm"]
    row["channel_height_nm"] = geom["channel_height_nm"]
    row["box_lx_nm"] = geom["box_lx_nm"]
    row["box_ly_nm"] = geom["box_ly_nm"]
    row["box_lz_nm"] = geom["box_lz_nm"]
    row["mouth_width_nm"] = args.mouth_width_nm
    row["reservoir_fraction"] = args.reservoir_fraction
    row["reservoir_width_nm"] = args.reservoir_fraction * geom["box_lx_nm"]
    row["slope_tail_ns"] = args.slope_tail_ns
    add_pmf_features(row, args.pmf_csv)

    add_density_features(row, results_dir, x_left, x_right, args.mouth_width_nm,
                         geom["box_lx_nm"], args.reservoir_fraction)
    add_gating_features(row, results_dir, args.slope_tail_ns, args.validate_gating)
    add_coord_features(row, results_dir, x_left, x_right, args.mouth_width_nm,
                       geom["box_lx_nm"], args.reservoir_fraction)
    add_nacl_cluster_features(row, results_dir, x_left, x_right, args.mouth_width_nm,
                              geom["box_lx_nm"], args.reservoir_fraction)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    action = write_or_upsert_row(args.out, row)
    print(f"{action}: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
