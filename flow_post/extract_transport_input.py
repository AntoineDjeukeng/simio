#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from config_flow import (
    COUNT_PATTERN,
    GATING_FLUX_FILENAME,
    I_END,
    I_START,
    RUN_REPLICA_DIR_PATTERN,
    RUNS_DIRNAME,
    RUNS_TRANSPORT_INPUT_RELATIVE,
)

CENTER_COLUMNS = (
    "time_ps",
    "center_water_left",
    "center_water_right",
    "center_na_left",
    "center_na_right",
    "center_cl_left",
    "center_cl_right",
)


def _pd():
    import pandas as pd

    return pd


def _load_center_plane_transport(gating_flux_path: Path):
    pd = _pd()
    df = pd.read_csv(gating_flux_path, sep=None, engine="python")
    missing = [c for c in CENTER_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"{gating_flux_path} missing required center-plane columns: {missing}"
        )

    time_ps = pd.to_numeric(df["time_ps"], errors="coerce")
    water_net = pd.to_numeric(df["center_water_right"], errors="coerce") - pd.to_numeric(
        df["center_water_left"], errors="coerce"
    )
    na_net = pd.to_numeric(df["center_na_right"], errors="coerce") - pd.to_numeric(
        df["center_na_left"], errors="coerce"
    )
    cl_net = pd.to_numeric(df["center_cl_right"], errors="coerce") - pd.to_numeric(
        df["center_cl_left"], errors="coerce"
    )

    out = pd.DataFrame(
        {
            "time_ps": time_ps,
            "SOL": water_net.cumsum(),
            "NA": na_net.cumsum(),
            "CL": cl_net.cumsum(),
        }
    )
    out = out.dropna(subset=["time_ps", "SOL", "NA", "CL"])
    out = out.sort_values("time_ps").reset_index(drop=True)
    if out.empty:
        raise ValueError(f"{gating_flux_path} produced empty transport_input.tsv")
    return out


def _inventory_gating_flux(path: Path):
    try:
        df = _pd().read_csv(path, sep=None, engine="python")
        cols = list(df.columns)
        missing = [c for c in CENTER_COLUMNS if c not in cols]
        return {"readable": True, "columns": cols, "missing_center_columns": missing}
    except Exception as exc:
        return {"readable": False, "error": str(exc)}


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Extract transport_input.tsv from runs/rep_XX/gating_flux.csv using center-plane crossings:\n"
            "SOL = cumsum(center_water_right - center_water_left)\n"
            "NA  = cumsum(center_na_right - center_na_left)\n"
            "CL  = cumsum(center_cl_right - center_cl_left)"
        )
    )
    ap.add_argument("field_dir", type=str, help=".../FIELD_XX directory")
    ap.add_argument("--i-start", type=int, default=I_START)
    ap.add_argument("--i-end", type=int, default=I_END)
    ap.add_argument(
        "--gating-file",
        type=str,
        default=GATING_FLUX_FILENAME,
        help="Per-replica gating flux filename inside runs/rep_XX/",
    )
    ap.add_argument(
        "--transport-name",
        type=str,
        default=RUNS_TRANSPORT_INPUT_RELATIVE,
        help="Output filename inside runs/rep_XX/",
    )
    ap.add_argument(
        "--inventory-only",
        action="store_true",
        help="Report gating_flux presence/columns per replica and exit.",
    )
    ap.add_argument(
        "--legacy-count-outdir",
        type=str,
        default=None,
        help="Optional legacy compatibility output directory for count_XX.dat (secondary mode).",
    )
    args = ap.parse_args()

    field_dir = Path(args.field_dir).resolve()
    runs_dir = field_dir / RUNS_DIRNAME
    legacy_outdir = Path(args.legacy_count_outdir).resolve() if args.legacy_count_outdir else None

    manifest = {
        "field_dir": str(field_dir),
        "runs_dir": str(runs_dir),
        "gating_file": args.gating_file,
        "transport_name": args.transport_name,
        "legacy_count_outdir": str(legacy_outdir) if legacy_outdir else None,
        "replicas": [],
    }
    if legacy_outdir:
        legacy_outdir.mkdir(parents=True, exist_ok=True)

    if args.inventory_only:
        for i in range(args.i_start, args.i_end + 1):
            rep_dir = runs_dir / RUN_REPLICA_DIR_PATTERN.format(i=i)
            rec = {"replica": i, "rep_dir": str(rep_dir), "exists": rep_dir.exists()}
            if rep_dir.exists():
                gf = rep_dir / args.gating_file
                rec["gating_flux"] = str(gf)
                rec["gating_flux_exists"] = gf.exists()
                if gf.exists():
                    rec.update(_inventory_gating_flux(gf))
            manifest["replicas"].append(rec)
        print(json.dumps(manifest, indent=2))
        return

    ok, failed = 0, 0
    for i in range(args.i_start, args.i_end + 1):
        rep_dir = runs_dir / RUN_REPLICA_DIR_PATTERN.format(i=i)
        rec = {"replica": i, "rep_dir": str(rep_dir)}
        if not rep_dir.exists():
            rec["status"] = "missing_replica_dir"
            manifest["replicas"].append(rec)
            failed += 1
            continue

        gating_flux_path = rep_dir / args.gating_file
        if not gating_flux_path.exists():
            rec["status"] = "missing_gating_flux"
            rec["gating_flux"] = str(gating_flux_path)
            manifest["replicas"].append(rec)
            failed += 1
            continue

        try:
            transport_df = _load_center_plane_transport(gating_flux_path)
            out_path = rep_dir / args.transport_name
            out_path.parent.mkdir(parents=True, exist_ok=True)
            transport_df.to_csv(out_path, sep="\t", index=False, float_format="%.8e")

            rec.update(
                {
                    "status": "ok",
                    "gating_flux": str(gating_flux_path),
                    "transport_out": str(out_path),
                    "rows": int(len(transport_df)),
                }
            )

            if legacy_outdir is not None:
                legacy_path = legacy_outdir / COUNT_PATTERN.format(i=i)
                transport_df.to_csv(
                    legacy_path,
                    sep=" ",
                    index=False,
                    header=False,
                    float_format="%.8e",
                )
                rec["legacy_count_out"] = str(legacy_path)

            ok += 1
        except Exception as exc:
            rec.update({"status": "error", "error": str(exc)})
            failed += 1
        manifest["replicas"].append(rec)

    manifest_path = field_dir / "transport_input_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"Wrote manifest: {manifest_path}")
    print(f"Replicas converted: {ok}")
    print(f"Replicas missing/failed: {failed}")


if __name__ == "__main__":
    main()

