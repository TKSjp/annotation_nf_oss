#!/usr/bin/env python3
# fix_cmscan_gff.py  input.gff > fixed.gff
import sys, re

MAP = {
    '5S_rRNA'          : ('rRNA', '5S ribosomal RNA',   None),
    'LSU_rRNA_eukarya' : ('rRNA', '28S ribosomal RNA',  None),
    'SSU_rRNA_eukarya' : ('rRNA', '18S ribosomal RNA',  None),
}

fh = open(sys.argv[1]) if len(sys.argv) == 2 else sys.stdin

for line in fh:
    if line.startswith('#') or not line.strip():
        sys.stdout.write(line)
        continue

    cols = line.rstrip('\n').split('\t')
    if len(cols) != 9:
        sys.stdout.write(line)
        continue

    source, ftype = cols[1], cols[2]

    if source == 'cmscan':
        # ---------- MAP にあるタイプ ----------
        if ftype in MAP:
            new_type, product, ncrna_class = MAP[ftype]
            attrs = cols[8]
            if 'product=' not in attrs:
                attrs += (';' if attrs and not attrs.endswith(';') else '') + f'product={product}'
            if new_type == 'ncRNA' and ncrna_class and 'ncrna_class=' not in attrs:
                attrs += f';ncrna_class={ncrna_class}'
            cols[2], cols[8] = new_type, attrs
            sys.stdout.write('\t'.join(cols) + '\n')

        # ---------- MAP にないタイプ ----------
        else:
            # feature type を ncRNA に変更
            original = ftype                      # 元のタイプを保存
            cols[2] = 'ncRNA'                     # 第三カラムを書き換え
            attrs = cols[8]
            # product に元タイプを記録
            if 'product=' not in attrs:
                attrs += (';' if attrs and not attrs.endswith(';') else '') + f'product={original}'
            cols[8] = attrs
            sys.stdout.write('\t'.join(cols) + '\n')

    else:
        # cmscan 以外はそのまま
        sys.stdout.write(line)

if fh is not sys.stdin:
    fh.close()