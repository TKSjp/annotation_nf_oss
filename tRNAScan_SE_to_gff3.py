#!/usr/bin/env python3
"""
tRNAScan_SE_to_gff3.py
  tRNAscan-SE の tabular 出力 (default) → GFF3 へ

使い方:
  python tRNAScan_SE_to_gff3.py tRNAs.tab > tRNAs.gff3
"""
import sys
import re

def parse_line(line):
    """tab 行を GFF 要素へパース"""
    cols = line.rstrip().split('\t')
    if len(cols) < 12:
        return None  # 期待列数でなければスキップ

    seqid   = cols[0]
    start   = int(cols[1])
    end     = int(cols[2])
    name    = cols[3]                     # 例 ctg001.tRNA80-AlaTGC
    score   = cols[4]
    strand  = cols[5] if cols[5] in ('+','-') else '+'

    # tRNA80-AlaTGC → 抗コドン TGC、アミノ酸 Ala を抽出
    m = re.search(r'-(\w{3})([ACGTU]{3})$', name)
    product = f"tRNA-{m.group(1)}" if m else "tRNA"

    attrs = [
        f"ID={name}",
        f"Name={name}",
        f"Product={product}",
        f"anticodon={m.group(2) if m else 'NNN'}",
        f"Score={score}"
    ]
    gff_fields = [
        seqid, "tRNAscan-SE", "tRNA",
        str(start+1), str(end), # GFF3 は 1-based 座標, tRNAscan-SE は 0-based 座標
        score, strand, ".",  ";".join(attrs)
    ]
    return "\t".join(gff_fields)

def main(tab):
    print("##gff-version 3")
    with open(tab) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            gff = parse_line(line)
            if gff:
                print(gff)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: tRNAScan_SE_to_gff3.py <tRNAs.tab> > out.gff3")
    main(sys.argv[1])
