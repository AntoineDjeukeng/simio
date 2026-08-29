#!/usr/bin/env python3
import argparse
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd


def parse_csv_list(value: str) -> List[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def average_cluster_count(csv_paths: Iterable[Path]) -> Tuple[Optional[pd.DataFrame], int]:
    frames = []
    for path in csv_paths:
        if not path.is_file():
            continue
        df = pd.read_csv(path)
        if "x_center_nm" not in df.columns or "nacl_cluster_count_mean" not in df.columns:
            continue
        frames.append(df[["x_center_nm", "nacl_cluster_count_mean"]])

    if not frames:
        return None, 0

    all_data = pd.concat(frames, ignore_index=True)
    mean_data = (
        all_data.groupby("x_center_nm", as_index=False)["nacl_cluster_count_mean"]
        .mean()
        .sort_values("x_center_nm")
    )
    return mean_data, len(frames)


def replica_paths(root: Path, h: str, length: str, charge: str, field: str,
                  run_dir: str, reps: Iterable[str]) -> List[Path]:
    field_dir = root / f"single_H_{h}" / f"L_{length}" / charge / field
    return [field_dir / run_dir / f"rep_{rep}" / "coord_x.csv" for rep in reps]


def plot_one_h(root: Path, h: str, lengths: List[str], charges: List[str],
               fields: List[str], reps: List[str], run_dir: str, out_dir: Path) -> Path:
    fig, axes = plt.subplots(
        len(charges),
        len(lengths),
        figsize=(4.4 * len(lengths), 3.2 * len(charges)),
        sharey=True,
        squeeze=False,
    )

    colors = {
        "FIELD_00": "tab:blue",
        "FIELD_01": "tab:orange",
        "FIELD_02": "tab:green",
        "FIELD_03": "tab:red",
    }

    for row, charge in enumerate(charges):
        for col, length in enumerate(lengths):
            ax = axes[row][col]
            n_lines = 0
            for field in fields:
                paths = replica_paths(root, h, length, charge, field, run_dir, reps)
                data, n_rep = average_cluster_count(paths)
                if data is None:
                    continue
                ax.plot(
                    data["x_center_nm"],
                    data["nacl_cluster_count_mean"],
                    label=f"{field} (n={n_rep})",
                    color=colors.get(field),
                )
                n_lines += 1

            ax.set_title(f"{charge}, L={length}")
            ax.set_xlabel("x (nm)")
            if col == 0:
                ax.set_ylabel("NaCl cluster count")
            ax.grid(True, alpha=0.25)
            if n_lines == 0:
                ax.text(0.5, 0.5, "no coord_x data", ha="center", va="center",
                        transform=ax.transAxes)
            else:
                ax.legend(fontsize=8)

    fig.suptitle(f"NaCl cluster count vs x, H={h}", y=0.995)
    fig.tight_layout()

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"nacl_cluster_count_grid_H{h}.png"
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Plot NaCl cluster count averaged over replicas as a charge x L grid."
    )
    ap.add_argument("--root", default="/data/antoine/Flow_CDI",
                    help="Root containing single_H_*/L_*/charge/FIELD_*/new_runs/rep_*/coord_x.csv")
    ap.add_argument("--H", dest="hs", default="7,9")
    ap.add_argument("--L", dest="lengths", default="1,2,6")
    ap.add_argument("--charges", default="positive,negative,neutral")
    ap.add_argument("--fields", default="FIELD_00,FIELD_01,FIELD_02,FIELD_03")
    ap.add_argument("--reps", default="01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20")
    ap.add_argument("--run-dir", default="new_runs",
                    help="Replica result directory name under each FIELD_* directory.")
    ap.add_argument("--out-dir", default=None,
                    help="Output directory. Defaults to ROOT/plots.")
    args = ap.parse_args()

    root = Path(args.root)
    out_dir = Path(args.out_dir) if args.out_dir else root / "plots"

    for h in parse_csv_list(args.hs):
        out_path = plot_one_h(
            root=root,
            h=h,
            lengths=parse_csv_list(args.lengths),
            charges=parse_csv_list(args.charges),
            fields=parse_csv_list(args.fields),
            reps=parse_csv_list(args.reps),
            run_dir=args.run_dir,
            out_dir=out_dir,
        )
        print(f"wrote: {out_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
