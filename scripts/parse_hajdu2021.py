"""Parse Table 'tab:binprop' (parameters of the 87 bulge RR Lyrae LTTE binary candidates) from the arXiv LaTeX source of
Hajdu et al. 2021, ApJ 915, 50 (arXiv:2105.03750; source downloaded from https://arxiv.org/e-print/2105.03750 and kept at
data/external/hajdu2021_arxiv2105.03750.tex). Output: data/external/hajdu2021_binprop.csv.

Usage: python scripts/parse_hajdu2021.py
"""
import re
from pathlib import Path

import pandas as pd

src = Path("data/external/hajdu2021_arxiv2105.03750.tex").read_text()
i = src.index("\\label{tab:binprop}")
body = src[src.index("\\startdata", i) + len("\\startdata"):src.index("\\enddata", i)]
cols = ["ID", "P", "Porb", "Porb_err", "a1sini_au", "a1sini_err", "e", "e_err", "omega", "omega_err", "T0", "T0_err",
        "beta", "beta_err", "K1", "K1_err", "fm", "fm_err", "Msmin", "quality"]
rows = []
for line in re.split(r"\\\\", body):
    line = line.strip()
    if not line or line.startswith("%"):
        continue
    parts = [p.strip() for p in line.split("&")]
    if len(parts) != 12:
        continue
    vals = [parts[0], parts[1]]
    for p in parts[2:10]:
        a = [x.strip() for x in p.split("$\\pm$")]
        vals += a if len(a) == 2 else [a[0], ""]
    vals += [parts[10], parts[11]]
    rows.append(vals)
d = pd.DataFrame(rows, columns=cols)
for c in cols[1:-1]:
    d[c] = pd.to_numeric(d[c], errors="coerce")
d.to_csv("data/external/hajdu2021_binprop.csv", index=False)
print(len(d), "rows; quality:", d.quality.value_counts().to_dict())
print(d.describe().loc[["min", "50%", "max"], ["Porb", "a1sini_au", "e", "K1", "fm", "Msmin"]].round(3))
