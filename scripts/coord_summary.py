import sys
import pandas as pd
import numpy as np

csv_path = sys.argv[1]
df = pd.read_csv(csv_path)

n = len(df)
i0 = int(0.2 * n)
i1 = int(0.8 * n)

res = pd.concat([df.iloc[:i0], df.iloc[i1:]])
mid = df.iloc[i0:i1]

def mean_sem(s: pd.Series):
    x = s.to_numpy(dtype=float)
    if len(x) <= 1:
        return float(x.mean()) if len(x) else 0.0, 0.0
    return float(x.mean()), float(x.std(ddof=1) / np.sqrt(len(x)))

for name, sub in [("res", res), ("mid", mid)]:
    hbww, _ = mean_sem(sub["HBww_tot_mean"])
    hbwcl, _ = mean_sem(sub["HBwcl_don_mean"])
    bww,  _ = mean_sem(sub["BWW_mean"])
    print(name, "HBww_tot", hbww, "HBwcl_don", hbwcl, "BWW", bww)
