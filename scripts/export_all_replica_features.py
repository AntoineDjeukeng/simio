#!/usr/bin/env python3
import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Iterable, List


def parse_csv_list(value: str) -> List[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def l_dirs_for(root: Path, h: str, length: str) -> List[Path]:
    h_dir = root / f"single_H_{h}"
    if not h_dir.is_dir():
        return []
    out: List[Path] = []
    for path in h_dir.glob(f"L_{length}*"):
        if not path.is_dir():
            continue
        name = path.name
        if name == f"L_{length}" or name in (f"L_{length}_2q", f"L_{length}_4q"):
            out.append(path)
    return sorted(out)


def discover_replica_dirs(root: Path, hs: Iterable[str], ls: Iterable[str],
                          charges: Iterable[str], fields: Iterable[str],
                          reps: Iterable[str], run_dir: str) -> List[Path]:
    out: List[Path] = []
    for h in hs:
        for length in ls:
            for l_dir in l_dirs_for(root, h, length):
                for charge in charges:
                    for field in fields:
                        for rep in reps:
                            rep_dir = l_dir / charge / field / run_dir / f"rep_{rep}"
                            if rep_dir.is_dir():
                                out.append(rep_dir)
    return out


def topo_for_replica(rep_dir: Path) -> Path:
    field_dir = rep_dir.parent.parent
    topo_setup = field_dir / "topo_setup.json"
    if topo_setup.is_file():
        return topo_setup
    charge_dir = field_dir.parent
    return charge_dir / "mini_topology.json"


def has_required_inputs(rep_dir: Path) -> bool:
    return (rep_dir / "density_x.csv").is_file() and (rep_dir / "gating_flux.csv").is_file()


def run_one(exporter: Path, rep_dir: Path, out_csv: Path, mouth_width_nm: float,
            reservoir_fraction: float, slope_tail_ns: float, strict: bool,
            validate_gating: bool, pmf_csv: str) -> bool:
    topo = topo_for_replica(rep_dir)
    if not topo.is_file():
        print(f"skip missing topology: {rep_dir}", file=sys.stderr)
        return False
    if not has_required_inputs(rep_dir):
        print(f"skip missing density_x.csv or gating_flux.csv: {rep_dir}", file=sys.stderr)
        return False

    cfg = {
        "topology_json": str(topo),
        "use_channel_bounds_from_topology": True,
    }

    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        json.dump(cfg, f)
        config_path = f.name

    try:
        cmd = [
            sys.executable,
            str(exporter),
            "--results-dir",
            str(rep_dir),
            "--config",
            config_path,
            "--out",
            str(out_csv),
            "--mouth-width-nm",
            str(mouth_width_nm),
            "--reservoir-fraction",
            str(reservoir_fraction),
            "--slope-tail-ns",
            str(slope_tail_ns),
            "--pmf-csv",
            pmf_csv,
        ]
        if validate_gating:
            cmd.append("--validate-gating")
        proc = subprocess.run(cmd)
        if proc.returncode != 0:
            print(f"skip exporter failed ({proc.returncode}): {rep_dir}", file=sys.stderr)
            if strict:
                raise subprocess.CalledProcessError(proc.returncode, cmd)
            return False
        return True
    finally:
        try:
            os.unlink(config_path)
        except OSError:
            pass


def main() -> int:
    here = Path(__file__).parent
    ap = argparse.ArgumentParser(
        description="Export/upsert compact feature rows for many FLOW_CDI replica result directories."
    )
    ap.add_argument("--root", default="/home/antoine/FLOW_CDI",
                    help="Root containing single_H_*/L_*/charge/FIELD_*/RUN_DIR/rep_*")
    ap.add_argument("--out", required=True, help="One master compact CSV to create/update")
    ap.add_argument("--H", dest="hs", default="7,9")
    ap.add_argument("--L", dest="ls", default="1,2,6")
    ap.add_argument("--charges", default="negative,neutral,positive")
    ap.add_argument("--fields", default="FIELD_00,FIELD_01,FIELD_02,FIELD_03")
    ap.add_argument("--reps", default="01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20")
    ap.add_argument("--run-dir", default="new_runs",
                    help="Replica output directory name under each FIELD_* directory. Use 'runs' for old outputs.")
    ap.add_argument("--mouth-width-nm", type=float, default=0.62)
    ap.add_argument("--reservoir-fraction", type=float, default=0.12)
    ap.add_argument("--slope-tail-ns", type=float, default=45.0)
    ap.add_argument("--exporter", default=None,
                    help="Path to export_replica_features.py. Defaults to the same directory as this script.")
    ap.add_argument("--strict", action="store_true",
                    help="Stop at the first failed replica instead of skipping it.")
    ap.add_argument("--validate-gating", action="store_true",
                    help="Check reported cumulative gating columns against sums of per-frame counts.")
    ap.add_argument("--pmf-csv", default="/home/antoine/CDI/plots/New_pmf_60/pmf_summary_full.csv",
                    help="Optional PMF summary CSV passed to export_replica_features.py.")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    root = Path(args.root)
    out_csv = Path(args.out)
    exporter = Path(args.exporter) if args.exporter else here / "export_replica_features.py"
    if not exporter.is_file():
        print(f"error: exporter not found: {exporter}", file=sys.stderr)
        print("Use --exporter /path/to/export_replica_features.py", file=sys.stderr)
        return 2

    rep_dirs = discover_replica_dirs(
        root,
        parse_csv_list(args.hs),
        parse_csv_list(args.ls),
        parse_csv_list(args.charges),
        parse_csv_list(args.fields),
        parse_csv_list(args.reps),
        args.run_dir,
    )

    print(f"found {len(rep_dirs)} replica directories")
    if args.dry_run:
        for d in rep_dirs:
            status = "ok" if has_required_inputs(d) and topo_for_replica(d).is_file() else "missing inputs"
            print(f"{status}: {d}")
        return 0

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    n_ok = 0
    for rep_dir in rep_dirs:
        if run_one(exporter, rep_dir, out_csv, args.mouth_width_nm,
                   args.reservoir_fraction, args.slope_tail_ns, args.strict,
                   args.validate_gating, args.pmf_csv):
            n_ok += 1

    print(f"exported {n_ok}/{len(rep_dirs)} rows into {out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
