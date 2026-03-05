#!/usr/bin/env python3
import csv
import math
import os
import sys
from typing import Dict, List, Tuple

EPS_ABS_DEFAULT = 1e-12
EPS_REL_DEFAULT = 1e-9

def is_float(s: str) -> bool:
    try:
        float(s)
        return True
    except Exception:
        return False

def list_files(root: str) -> List[str]:
    out = []
    for dirpath, _, filenames in os.walk(root):
        for fn in filenames:
            rel = os.path.relpath(os.path.join(dirpath, fn), root)
            out.append(rel)
    return sorted(out)

def read_csv(path: str) -> Tuple[List[str], List[List[str]]]:
    with open(path, newline="") as f:
        r = csv.reader(f)
        rows = list(r)
    if not rows:
        return [], []
    header = rows[0]
    data = rows[1:]
    return header, data

def compare_csv(a_path: str, b_path: str, eps_abs: float, eps_rel: float) -> List[str]:
    errs = []
    ah, ad = read_csv(a_path)
    bh, bd = read_csv(b_path)

    if ah != bh:
        errs.append(f"HEADER_MISMATCH: {a_path} vs {b_path}\n  A={ah}\n  B={bh}")
        return errs

    if len(ad) != len(bd):
        errs.append(f"ROW_COUNT_MISMATCH: {a_path} ({len(ad)}) vs {b_path} ({len(bd)})")
        # keep going on common prefix

    nrows = min(len(ad), len(bd))
    ncols = min(len(ah), len(bh))

    for i in range(nrows):
        ar = ad[i]
        br = bd[i]
        if len(ar) != len(br):
            errs.append(f"COL_COUNT_MISMATCH row={i+1}: {a_path} ({len(ar)}) vs {b_path} ({len(br)})")
            ncols = min(len(ar), len(br), ncols)
        for j in range(ncols):
            av = ar[j]
            bv = br[j]
            if av == bv:
                continue
            # numeric compare if both parse
            if is_float(av) and is_float(bv):
                af = float(av); bf = float(bv)
                diff = abs(af - bf)
                tol = max(eps_abs, eps_rel * max(abs(af), abs(bf), 1.0))
                if math.isnan(af) and math.isnan(bf):
                    continue
                if diff > tol:
                    errs.append(
                        f"VAL_MISMATCH {os.path.basename(a_path)} row={i+1} col={j+1} ({ah[j]}): "
                        f"A={af} B={bf} diff={diff} tol={tol}"
                    )
            else:
                errs.append(
                    f"STR_MISMATCH {os.path.basename(a_path)} row={i+1} col={j+1} ({ah[j]}): A='{av}' B='{bv}'"
                )
    return errs

def main():
    if len(sys.argv) < 3:
        print("Usage: compare_outputs.py <dirA> <dirB> [eps_abs=1e-12] [eps_rel=1e-9]", file=sys.stderr)
        return 2
    a = sys.argv[1]
    b = sys.argv[2]
    eps_abs = float(sys.argv[3]) if len(sys.argv) >= 4 else EPS_ABS_DEFAULT
    eps_rel = float(sys.argv[4]) if len(sys.argv) >= 5 else EPS_REL_DEFAULT

    afiles = list_files(a)
    bfiles = list_files(b)

    if afiles != bfiles:
        print("FILE_SET_MISMATCH")
        only_a = sorted(set(afiles) - set(bfiles))
        only_b = sorted(set(bfiles) - set(afiles))
        if only_a:
            print("Only in A:")
            for f in only_a[:200]:
                print("  ", f)
        if only_b:
            print("Only in B:")
            for f in only_b[:200]:
                print("  ", f)
        return 1

    # ignore logs by default
    csv_files = [f for f in afiles if f.lower().endswith(".csv")]
    non_csv = [f for f in afiles if not f.lower().endswith(".csv") and f != "run.log"]

    # non-csv must match byte-for-byte (except run.log)
    for f in non_csv:
        ap = os.path.join(a, f)
        bp = os.path.join(b, f)
        with open(ap, "rb") as fa, open(bp, "rb") as fb:
            if fa.read() != fb.read():
                print(f"NONCSV_BINARY_MISMATCH: {f}")
                return 1

    all_errs: List[str] = []
    for f in csv_files:
        ap = os.path.join(a, f)
        bp = os.path.join(b, f)
        all_errs.extend(compare_csv(ap, bp, eps_abs, eps_rel))

    if all_errs:
        print(f"CSV_MISMATCHES: {len(all_errs)}")
        for e in all_errs[:200]:
            print(e)
        if len(all_errs) > 200:
            print(f"... truncated, {len(all_errs)-200} more")
        return 1

    print("OK: outputs match (CSV numeric compare + non-CSV binary compare; run.log ignored)")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
