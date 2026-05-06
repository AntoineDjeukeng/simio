#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Export the subset of axial_profile_x.csv needed for the standard "
            "3-panel x-profile plot, in a form that is easy to use from gnuplot."
        )
    )
    parser.add_argument("input_csv", type=Path, help="Path to axial_profile_x.csv")
    parser.add_argument(
        "--out-prefix",
        type=Path,
        default=Path("axial_profile_plot"),
        help="Output prefix for <prefix>.csv, <prefix>_vars.gp, <prefix>_example.gp",
    )
    parser.add_argument(
        "--xmin",
        type=float,
        default=-4.2,
        help="Lower x limit for exported rows and gnuplot example",
    )
    parser.add_argument(
        "--xmax",
        type=float,
        default=4.2,
        help="Upper x limit for exported rows and gnuplot example",
    )
    return parser.parse_args()


FIELDNAMES = [
    "x_nm",
    "in_channel",
    "rho_ow",
    "rho_hw_half",
    "rho_na_x10",
    "rho_cl_x10",
    "muz_fold",
    "mux",
    "ac",
    "na_bond",
    "dn",
    "cl_bond",
    "ac_plus_na",
    "dn_plus_cl",
]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        raise SystemExit(f"No data rows found in {path}")
    return rows


def as_float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def fmt(value: float) -> str:
    return f"{value:.10g}"


def build_export_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    export_rows: list[dict[str, str]] = []
    for row in rows:
        export_rows.append(
            {
                "x_nm": fmt(as_float(row, "x_relative_nm")),
                "in_channel": str(int(float(row["in_channel"]))),
                "rho_ow": fmt(as_float(row, "rho_ow_mean")),
                "rho_hw_half": fmt(0.5 * as_float(row, "rho_hw_mean")),
                "rho_na_x10": fmt(10.0 * as_float(row, "rho_na_mean")),
                "rho_cl_x10": fmt(10.0 * as_float(row, "rho_cl_mean")),
                "muz_fold": fmt(as_float(row, "muz_fold_mean")),
                "mux": fmt(as_float(row, "mux_mean")),
                "ac": fmt(as_float(row, "ac_mean")),
                "na_bond": fmt(as_float(row, "na_bond_mean")),
                "dn": fmt(as_float(row, "dn_mean")),
                "cl_bond": fmt(as_float(row, "cl_bond_mean")),
                "ac_plus_na": fmt(as_float(row, "ac_plus_na_mean")),
                "dn_plus_cl": fmt(as_float(row, "dn_plus_cl_mean")),
            }
        )
    return export_rows


def write_export_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def write_vars_file(path: Path, data_file: Path, xl: float, xr: float, xmin: float, xmax: float) -> None:
    text = "\n".join(
        [
            f"datafile = '{data_file.name}'",
            f"xl = {xl:.12g}",
            f"xr = {xr:.12g}",
            f"xmin = {xmin:.12g}",
            f"xmax = {xmax:.12g}",
            "",
            "# Column map for the exported CSV:",
            "# 1  x_nm",
            "# 2  in_channel",
            "# 3  rho_ow",
            "# 4  rho_hw_half",
            "# 5  rho_na_x10",
            "# 6  rho_cl_x10",
            "# 7  muz_fold",
            "# 8  mux",
            "# 9  ac",
            "# 10 na_bond",
            "# 11 dn",
            "# 12 cl_bond",
            "# 13 ac_plus_na",
            "# 14 dn_plus_cl",
            "",
            "col_x = 1",
            "col_in_channel = 2",
            "col_rho_ow = 3",
            "col_rho_hw_half = 4",
            "col_rho_na_x10 = 5",
            "col_rho_cl_x10 = 6",
            "col_muz_fold = 7",
            "col_mux = 8",
            "col_ac = 9",
            "col_na_bond = 10",
            "col_dn = 11",
            "col_cl_bond = 12",
            "col_ac_plus_na = 13",
            "col_dn_plus_cl = 14",
            "",
        ]
    )
    path.write_text(text)


def write_example_file(path: Path, vars_file: Path) -> None:
    text = "\n".join(
        [
            "set datafile separator ','",
            f"load '{vars_file.name}'",
            "",
            "set terminal qt size 1000,800",
            "set multiplot layout 3,1 rowsfirst",
            "set xrange [xmin:xmax]",
            "set grid",
            "",
            "unset key",
            "set ylabel 'Density'",
            "set arrow 1 from xl, graph 0 to xl, graph 1 nohead dt 2 lc rgb 'black'",
            "set arrow 2 from xr, graph 0 to xr, graph 1 nohead dt 2 lc rgb 'black'",
            "plot datafile using col_x:col_rho_ow with lines lw 2 title 'O_W', \\",
            "     datafile using col_x:col_rho_hw_half with lines lw 2 title 'H_W / 2', \\",
            "     datafile using col_x:col_rho_na_x10 with lines lw 2 title '10*Na+', \\",
            "     datafile using col_x:col_rho_cl_x10 with lines lw 2 title '10*Cl-'",
            "",
            "unset arrow 1",
            "unset arrow 2",
            "set ylabel 'Dipole'",
            "set arrow 1 from xl, graph 0 to xl, graph 1 nohead dt 2 lc rgb 'black'",
            "set arrow 2 from xr, graph 0 to xr, graph 1 nohead dt 2 lc rgb 'black'",
            "plot datafile using col_x:col_muz_fold with lines lw 2 title 'muz_fold', \\",
            "     datafile using col_x:col_mux with lines lw 2 title 'mux'",
            "",
            "unset arrow 1",
            "unset arrow 2",
            "set ylabel 'Bonds / water'",
            "set xlabel 'x from pore center (nm)'",
            "set arrow 1 from xl, graph 0 to xl, graph 1 nohead dt 2 lc rgb 'black'",
            "set arrow 2 from xr, graph 0 to xr, graph 1 nohead dt 2 lc rgb 'black'",
            "plot datafile using col_x:col_ac with lines lw 2 title 'Ac', \\",
            "     datafile using col_x:col_na_bond with lines lw 2 title 'Na+', \\",
            "     datafile using col_x:col_dn with lines lw 2 title 'Dn', \\",
            "     datafile using col_x:col_cl_bond with lines lw 2 title 'Cl-', \\",
            "     datafile using col_x:col_ac_plus_na with lines dt 2 lw 2 lc rgb 'black' title 'Ac+Na+', \\",
            "     datafile using col_x:col_dn_plus_cl with lines dt 2 lw 2 lc rgb 'gray40' title 'Dn+Cl-'",
            "",
            "unset multiplot",
            "",
        ]
    )
    path.write_text(text)


def main() -> int:
    args = parse_args()

    rows = read_rows(args.input_csv)
    rows.sort(key=lambda row: as_float(row, "x_relative_nm"))
    xl = as_float(rows[0], "channel_edge_left_relative_nm")
    xr = as_float(rows[0], "channel_edge_right_relative_nm")

    rows = [
        row
        for row in rows
        if args.xmin <= as_float(row, "x_relative_nm") <= args.xmax
    ]

    if not rows:
        raise SystemExit("No rows remain after applying the x-range filter.")

    out_prefix = args.out_prefix
    out_dir = out_prefix.parent if out_prefix.parent != Path("") else Path(".")
    out_dir.mkdir(parents=True, exist_ok=True)

    data_path = out_prefix.with_suffix(".csv")
    vars_path = out_dir / f"{out_prefix.name}_vars.gp"
    example_path = out_dir / f"{out_prefix.name}_example.gp"

    write_export_csv(data_path, build_export_rows(rows))
    write_vars_file(vars_path, data_path, xl, xr, args.xmin, args.xmax)
    write_example_file(example_path, vars_path)

    print(f"Wrote {data_path}")
    print(f"Wrote {vars_path}")
    print(f"Wrote {example_path}")
    print(f"Channel edges: xl={xl:.10g} xr={xr:.10g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
