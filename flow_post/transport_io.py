from __future__ import annotations

from pathlib import Path
import re
from typing import Iterable, Sequence

import pandas as pd


REQUIRED_COLUMNS = ("time_ps", "SOL", "NA", "CL")


def _replica_sort_key(path: Path):
    m = re.search(r"(?:count_|rep_)(\d+)", path.name)
    if m:
        return (0, int(m.group(1)), path.name)
    return (1, path.name)


def infer_field_dir(path: Path) -> Path:
    p = Path(path).resolve()
    for cand in (p, *p.parents):
        if cand.name.startswith("FIELD_"):
            return cand
    return p


def resolve_replica_files(
    input_dir: Path,
    i_start: int,
    i_end: int,
    pattern: str,
    replica_files: Sequence[Path] | None = None,
    replica_glob: str | None = None,
) -> list[Path]:
    input_dir = Path(input_dir)

    if replica_files:
        files = [Path(p) for p in replica_files]
    elif replica_glob:
        files = sorted(input_dir.glob(replica_glob), key=_replica_sort_key)
    elif pattern:
        files = []
        for i in range(i_start, i_end + 1):
            path = input_dir / pattern.format(i=i)
            if path.exists():
                files.append(path)
    else:
        files = []

    files = [p for p in files if p.exists()]
    if not files:
        raise RuntimeError(
            "No replica files found. "
            "Provide --replica-file or --replica-glob (normal path: rep_*/transport_input.tsv)."
        )
    return sorted(files, key=_replica_sort_key)


def _read_legacy_dat(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", comment="#", header=None, engine="python")
    if df.shape[1] < 4:
        raise ValueError(f"{path} has {df.shape[1]} columns, expected at least 4")
    df = df.iloc[:, :4].copy()
    df.columns = list(REQUIRED_COLUMNS)
    return df


def _lower_map(columns: Iterable[str]) -> dict[str, str]:
    return {str(c).strip().lower(): str(c) for c in columns}


def _read_normalized_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    cmap = _lower_map(df.columns)
    missing = [c for c in REQUIRED_COLUMNS if c.lower() not in cmap]
    if missing:
        raise ValueError(
            f"{path} does not provide normalized columns {REQUIRED_COLUMNS}. "
            f"Missing: {missing}. Use extract_transport_input.py first."
        )
    out = pd.DataFrame(
        {
            "time_ps": df[cmap["time_ps"]],
            "SOL": df[cmap["sol"]],
            "NA": df[cmap["na"]],
            "CL": df[cmap["cl"]],
        }
    )
    return out


def read_replica(path: Path, baseline_zero: bool = True):
    path = Path(path)
    if path.suffix.lower() in {".tsv", ".csv"}:
        df = _read_normalized_table(path)
    elif path.suffix.lower() == ".dat":
        df = _read_legacy_dat(path)
    else:
        # Try normalized parser first for unknown extensions, then legacy fallback.
        try:
            df = _read_normalized_table(path)
        except Exception:
            df = _read_legacy_dat(path)

    df = df.copy()
    for col in REQUIRED_COLUMNS:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=list(REQUIRED_COLUMNS))
    df = df.sort_values("time_ps").reset_index(drop=True)
    if df.empty:
        raise ValueError(f"{path} produced an empty trajectory after numeric filtering")

    t = df["time_ps"].to_numpy(float)
    sol = df["SOL"].to_numpy(float)
    na = df["NA"].to_numpy(float)
    cl = df["CL"].to_numpy(float)
    if baseline_zero:
        sol = sol - sol[0]
        na = na - na[0]
        cl = cl - cl[0]
    return t, sol, na, cl
