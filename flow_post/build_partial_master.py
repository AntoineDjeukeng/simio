#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
from config_flow import ROOT, RESULTS

RESULTS.mkdir(parents=True, exist_ok=True)

files = sorted(ROOT.rglob("FIELD_*/transport_fits/transport.tsv"))
if not files:
    files = sorted(ROOT.rglob("FIELD_*/compile/reduced/transport.tsv"))
if not files:
    raise SystemExit("No transport.tsv found yet. Run fit_transport.py first.")

dfs = []
for f in files:
    try:
        df = pd.read_csv(f, sep="\t")
        if not df.empty:
            dfs.append(df)
    except Exception as e:
        print("SKIP", f, "->", e)

master = pd.concat(dfs, ignore_index=True)
sort_cols = [c for c in ["charge","H_nm","L_nm","E_V_nm","field"] if c in master.columns]
if sort_cols:
    master = master.sort_values(sort_cols, kind="mergesort")

out = RESULTS / "master_transport_partial.tsv"
master.to_csv(out, sep="\t", index=False)
print("Wrote:", out)
print("Rows:", len(master))
print("Unique conditions (H,L,charge):", master.groupby(["H_nm","L_nm","charge"]).ngroups)
