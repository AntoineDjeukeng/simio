#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

from fit_transport import process_compile_dir, process_input_root  # same folder import

ROOT = Path("/data/antoine/Flow_CDI")

# heights/lengths you specified (acts as a filter)
# H_DIRS = ["single_H_7", "single_H_9"]
H_DIRS = ["single_H_9"]
L_DIRS = [ "L_2", "L_6"]
# L_DIRS = ["L_1", "L_2", "L_6"]
# CHARGE_DIRS = ["neutral", "positive", "negative"]
CHARGE_DIRS = ["neutral"]
FIELD_DIRS = ["FIELD_00", "FIELD_01", "FIELD_02", "FIELD_03"]

FIT_LAST_NS = 45.0

rows = []
for h in H_DIRS:
    for l in L_DIRS:
        for ch in CHARGE_DIRS:
            base = ROOT / h / l / ch
            if not base.exists():
                continue
            for fld in FIELD_DIRS:
                input_root = base / fld / "runs"
                try:
                    if input_root.exists():
                        row = process_input_root(input_root, fit_last_ns=FIT_LAST_NS)
                    else:
                        compile_dir = base / fld / "compile"
                        if not compile_dir.exists():
                            continue
                        row = process_compile_dir(compile_dir, fit_last_ns=FIT_LAST_NS)
                    rows.append(row)
                except Exception as e:
                    print("SKIP", input_root, "->", e)

df = pd.DataFrame(rows)
out = ROOT / "RESULTS" / "master_transport.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, sep="\t", index=False)
print("Wrote:", out)
print("Rows:", len(df))
