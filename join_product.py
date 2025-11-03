#!/usr/bin/env python3
"""
add_tpm_and_product.py  counts.tsv  augustus.gff  > counts_plus.tsv
  - counts.tsv : featureCounts 出力 (先頭に '# Program:' 行)
  - augustus.gff : CDS 行に product= が付く GFF3
出力 : counts.tsv に 各サンプルTPM列 + Product 列 を追加
"""

import sys, re
import pandas as pd
from pathlib import Path
from math import isclose

# ----------- 引数チェック -----------
if len(sys.argv) != 3:
    sys.exit("usage: python add_tpm_and_product.py counts.tsv augustus.gff > out.tsv")

counts_file, gff_file = map(Path, sys.argv[1:])
if not counts_file.exists() or not gff_file.exists():
    sys.exit("input file not found")

# ----------- 1) DataFrame 読み込み (skip '# Program:' 行) -----------
df = pd.read_csv(counts_file, sep="\t", comment="#", skiprows=1)

# サンプル列を取得 (Geneid〜Length の 6 列を除く残り全部)
sample_cols = df.columns[6:]
if sample_cols.empty:
    sys.exit("No sample column detected in counts file")

# ----------- 2) TPM 計算 -----------
length_kb = df["Length"] / 1000.0

for col in sample_cols:
    rpk  = df[col] / length_kb
    scaling = rpk.sum() / 1e6            # ΣRPK / 10^6
    if isclose(scaling, 0.0):
        df[f"TPM_{col}"] = 0
    else:
        df[f"TPM_{col}"] = (rpk / scaling).round(3)

# ----------- 3) transcriptID → product 辞書 -----------
prod_dict = {}
pat = re.compile(r"Parent=([^;]+).*?product=([^;]+)", re.I)

with gff_file.open(newline="") as fh:
    for ln in fh:
        if ln.startswith("#"):
            continue
        cols = ln.rstrip("\r\n").split("\t")
        if len(cols) != 9 or cols[2] != "CDS":
            continue
        m = pat.search(cols[8])
        if m:
            tid, prod = m.groups()
            prod_dict.setdefault(tid, prod)   # 最初のみ採用

df["Product"] = df["Geneid"].map(prod_dict).fillna("NA")

# ----------- 4) TSV 書き出し （元の '# Program:' 行を残す） -----------
with counts_file.open() as fh_in:
    first_line = fh_in.readline().rstrip("\n")

print(first_line)                   # '# Program:featureCounts …'
df.to_csv(sys.stdout, sep="\t", index=False)