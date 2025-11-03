#!/usr/bin/env python3
"""
tbl+GFF -> GFF (CDS に product を付加)
  1) DIAMOND 12 列 TSV から (mRNA_ID -> product) を取得
  2) Helixer GFF3 を読み、CDS 行だけ
       Parent=<mRNA_ID> が一致すれば
       product=<description> を追加
使い方:
  python add_product_cds.py helixer.gff diamond.tsv > helixer_product.gff
"""

import re, sys

# ------------------------------------------------------------------
def parse_diamond(tsv_path):
    """DIAMOND TSV -> {mRNA_ID: product-description}"""
    id2prod = {}
    with open(tsv_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 13:
                continue
            qid   = cols[0]       # mRNA ID (例: ctg001_000001.1)
            descr = cols[12]

            # 'UniRef90_XXXX <product> n= ... Tax= ...' から product 部分だけ抜く
            m = re.match(r'\S+\s+(.+?)(?:\s+n=|\s+Tax=|\s+RepID=|$)', descr)
            product = m.group(1) if m else descr
            id2prod[qid] = product
    return id2prod

# ------------------------------------------------------------------
def gff_add_product_cds(gff_path, id2prod, out=sys.stdout):
    parent_re = re.compile(r'Parent=([^;]+)')
    for line in open(gff_path):
        if line.startswith('#') or not line.strip():
            print(line, end='', file=out)
            continue

        fields = line.rstrip().split('\t')
        if len(fields) != 9:
            print(line, end='', file=out); continue

        feature = fields[2]
        attrs   = fields[8]

        # --- 追記対象は CDS のみ ---
        if feature == 'CDS':
            par_match = parent_re.search(attrs)
            if par_match:
                mrna_id = par_match.group(1)
                if mrna_id in id2prod and 'product=' not in attrs:
                    attrs += f'product={id2prod[mrna_id]}' #;product= helixer, product= braker
                    fields[8] = attrs
        print('\t'.join(fields), file=out)

# ------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(f'usage: {sys.argv[0]} helixer.gff diamond.tsv > out.gff')
    gff_in, dia_tsv = sys.argv[1:]
    id2p = parse_diamond(dia_tsv)
    gff_add_product_cds(gff_in, id2p)
    