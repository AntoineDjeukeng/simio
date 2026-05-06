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

REGIONS = [
    "center",
    "res_left",
    "res_right",
    "mouth_left",
    "mouth_right",
]

CLUSTER_COLUMNS = [
    ("cluster_count", "nacl_cluster_count_mean"),
    ("ions_in_clusters", "ions_in_nacl_clusters_mean"),
    ("na_in_clusters", "na_in_nacl_clusters_mean"),
    ("cl_in_clusters", "cl_in_nacl_clusters_mean"),
    ("cluster_size", "nacl_cluster_size_mean"),
]


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


def channel_edges_from_config(config_path: str) -> Tuple[float, float]:
    cfg = load_json(config_path)
    config_dir = os.path.dirname(os.path.abspath(config_path))

    topo_path = cfg.get("topology_json") or cfg.get("mini_topology_json") or cfg.get("setup_json")
    if topo_path:
        topo_abs = resolve_path(str(topo_path), config_dir)
        if os.path.exists(topo_abs):
            topo = load_json(topo_abs)
            ch = topo.get("channel")
            if isinstance(ch, dict) and "min" in ch and "max" in ch:
                return float(ch["min"][0]), float(ch["max"][0])

    if "xmin" in cfg and "xmax" in cfg:
        return float(cfg["xmin"]), float(cfg["xmax"])

    raise ValueError(
        "Could not determine carbon/channel x edges. Provide topology_json with channel.min/max "
        "or xmin/xmax in the config."
    )


def region_for_x(x: float, x_left: float, x_right: float, mouth_width_nm: float) -> Optional[str]:
    if x_left <= x < x_right:
        return "center"
    if x_left - mouth_width_nm <= x < x_left:
        return "mouth_left"
    if x_right <= x < x_right + mouth_width_nm:
        return "mouth_right"
    if x < x_left - mouth_width_nm:
        return "res_left"
    if x >= x_right + mouth_width_nm:
        return "res_right"
    return None


def add_density_features(out: Dict[str, object], results_dir: str, x_left: float, x_right: float,
                         mouth_width_nm: float) -> None:
    path = os.path.join(results_dir, "density_x.csv")
    rows = read_csv_dicts(path)
    buckets: Dict[Tuple[str, str], List[float]] = {}
    for row in rows:
        x = to_float(row.get("x_center_nm", ""))
        region = region_for_x(x, x_left, x_right, mouth_width_nm)
        if region is None:
            continue
        for label, internal in SPECIES:
            col = f"rho_{internal}_mean"
            if col in row:
                buckets.setdefault((label, region), []).append(to_float(row[col]))

    for label, _ in SPECIES:
        for region in REGIONS:
            out[f"{label}_{region}_mean"] = mean(buckets.get((label, region), []))


def add_cluster_features(out: Dict[str, object], results_dir: str, x_left: float, x_right: float,
                         mouth_width_nm: float) -> None:
    path = os.path.join(results_dir, "nacl_cluster_x.csv")
    if not os.path.exists(path):
        return

    rows = read_csv_dicts(path)
    buckets: Dict[Tuple[str, str], List[float]] = {}
    for row in rows:
        x = to_float(row.get("x_center_nm", ""))
        region = region_for_x(x, x_left, x_right, mouth_width_nm)
        if region is None:
            continue
        for feature, col in CLUSTER_COLUMNS:
            if col in row:
                buckets.setdefault((feature, region), []).append(to_float(row[col]))

    for feature, _ in CLUSTER_COLUMNS:
        for region in REGIONS:
            out[f"NACL_{feature}_{region}_mean"] = mean(buckets.get((feature, region), []))


def cumulative(values: Sequence[float]) -> List[float]:
    out: List[float] = []
    s = 0.0
    for v in values:
        s += v
        out.append(s)
    return out


def add_gating_features(out: Dict[str, object], results_dir: str, tail_ns: float) -> None:
    path = os.path.join(results_dir, "gating_flux.csv")
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
        cum_right_col = f"center_{internal}_cum_right"

        if left_col not in rows[0] or right_col not in rows[0]:
            continue

        left_events = [to_float(r[left_col]) for r in rows]
        left_cum = cumulative(left_events)
        if cum_right_col in rows[0]:
            right_cum = [to_float(r[cum_right_col]) for r in rows]
        else:
            right_cum = cumulative([to_float(r[right_col]) for r in rows])
        net_cum = [r - l for r, l in zip(right_cum, left_cum)]

        out[f"center_{label}_cum_left_final"] = left_cum[-1]
        out[f"center_{label}_cum_right_final"] = right_cum[-1]
        out[f"{label}_net_slope"] = linear_slope([time_ns[i] for i in tail_idx],
                                                 [net_cum[i] for i in tail_idx])
        out[f"{label}_left_slope"] = linear_slope([time_ns[i] for i in tail_idx],
                                                  [left_cum[i] for i in tail_idx])
        out[f"{label}_right_slope"] = linear_slope([time_ns[i] for i in tail_idx],
                                                   [right_cum[i] for i in tail_idx])


def infer_metadata_from_path(path: str) -> Dict[str, str]:
    parts = os.path.abspath(path).split(os.sep)
    meta: Dict[str, str] = {}
    for p in parts:
        m = re.fullmatch(r"H[_-]?([0-9.]+)", p, re.IGNORECASE)
        if m:
            meta["H"] = m.group(1)
        m = re.fullmatch(r"L[_-]?([0-9.]+)", p, re.IGNORECASE)
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
    return meta


def ordered_fields(row: Dict[str, object]) -> List[str]:
    first = ["H", "L", "charge", "field", "replica"]
    rest = [k for k in row.keys() if k not in first]
    return first + rest


def write_one_row(path: str, row: Dict[str, object]) -> None:
    fields = ordered_fields(row)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerow(row)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Export one compact structure/transport/cluster feature row for one replica."
    )
    ap.add_argument("--results-dir", required=True, help="Directory containing density_x.csv and gating_flux.csv")
    ap.add_argument("--config", required=True, help="Run config JSON; carbon/channel edges are read from it")
    ap.add_argument("--out", required=True, help="Output CSV path")
    ap.add_argument("--mouth-width-nm", type=float, default=0.62)
    ap.add_argument("--slope-tail-ns", type=float, default=45.0)
    ap.add_argument("--H", dest="H", default=None)
    ap.add_argument("--L", dest="L", default=None)
    ap.add_argument("--charge", default=None)
    ap.add_argument("--field", default=None)
    ap.add_argument("--replica", default=None)
    args = ap.parse_args()

    results_dir = os.path.abspath(args.results_dir)
    x_left, x_right = channel_edges_from_config(args.config)

    row: Dict[str, object] = {}
    row.update(infer_metadata_from_path(results_dir))
    for key in ("H", "L", "charge", "field", "replica"):
        value = getattr(args, key)
        row[key] = value if value is not None else row.get(key, "")

    row["carbon_left_x_nm"] = x_left
    row["carbon_right_x_nm"] = x_right
    row["mouth_width_nm"] = args.mouth_width_nm
    row["slope_tail_ns"] = args.slope_tail_ns

    add_density_features(row, results_dir, x_left, x_right, args.mouth_width_nm)
    add_gating_features(row, results_dir, args.slope_tail_ns)
    add_cluster_features(row, results_dir, x_left, x_right, args.mouth_width_nm)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    write_one_row(args.out, row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
