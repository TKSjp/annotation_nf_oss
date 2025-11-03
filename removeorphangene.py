#!/usr/bin/env python3
# remove_orphan_gene.py  in.gff  > out.gff
import sys, re

if len(sys.argv) != 2:
    sys.exit("usage: python remove_orphan_gene.py input.gff > filtered.gff")

gff_in = sys.argv[1]

gene_id_re = re.compile(r'ID=([^;]+)')
parent_re  = re.compile(r'Parent=([^;]+)')

# -------------------------------------------------
# 1) まず gene と 子フィーチャの対応を調べる
gene_lines   = []             # (line, gene_id)
parent_found = set()          # gene_id を Parent= に持つ行

with open(gff_in) as fh:
    for line in fh:
        if line.startswith('#') or not line.strip():
            continue
        fields = line.rstrip('\n').split('\t')
        ftype  = fields[2]
        quals  = fields[8]

        if ftype == 'gene':
            m = gene_id_re.search(quals)
            if m:
                gene_lines.append((line, m.group(1)))
        else:
            m = parent_re.search(quals)
            if m:
                parent_found.add(m.group(1))

# orphan gene の集合
orphan = {gid for _, gid in gene_lines if gid not in parent_found}

# -------------------------------------------------
# 2) 再読み込みして orphan gene と子孫行をスキップ
with open(gff_in) as fh:
    for line in fh:
        if not orphan:
            sys.stdout.write(line); continue

        if line.startswith('#') or not line.strip():
            sys.stdout.write(line); continue

        fields = line.rstrip('\n').split('\t')
        quals  = fields[8]

        # gene 行か？ → その ID が orphan ならスキップ
        if fields[2] == 'gene':
            m = gene_id_re.search(quals)
            if m and m.group(1) in orphan:
                continue  # drop gene
        # それ以外の行で Parent が orphan gene ならスキップ
        m = parent_re.search(quals)
        if m and m.group(1) in orphan:
            continue

        sys.stdout.write(line)