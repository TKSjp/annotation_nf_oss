#!/usr/bin/env python3
"""
add_hypothetical_product.py
  GFF3 の CDS で product が欠落している行に
  ;product=hypothetical protein を追記

使い方:
  python add_hypothetical_product.py in.gff > out.gff
"""

import sys

def patch_gff(gff_path, out=sys.stdout):
    with open(gff_path) as fh:
        for line in fh:
            # コメントや空行はそのまま
            if line.startswith('#') or not line.strip():
                print(line, end='', file=out)
                continue

            cols = line.rstrip('\n').split('\t')
            if len(cols) != 9:
                print(line, end='', file=out)
                continue

            feature = cols[2]
            attrs   = cols[8]

            if feature == 'CDS' and 'product=' not in attrs:
                # 末尾にセミコロンが無ければ追加
                if not attrs.endswith(';'):
                    attrs += ';'
                attrs += 'product=hypothetical protein'
                cols[8] = attrs

            print('\t'.join(cols), file=out)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(f'usage: {sys.argv[0]} <in.gff> > out.gff')
    patch_gff(sys.argv[1])
    