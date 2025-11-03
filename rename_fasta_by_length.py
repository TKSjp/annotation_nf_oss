#!/usr/bin/env python3
"""
rename_fasta_by_length.py

FASTA 中の配列を長さの降順に並べ替え、
Chr01, Chr02, … のように連番でリネームして出力します。

使用例:
  python rename_fasta_by_length.py --in input.fa --out renamed.fa --width 3 --prefix Scaffold
"""
import argparse
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def rename_by_length(in_fa: Path, out_fa: Path, prefix: str = "Chr", width: int = None) -> None:
    # 全レコードを読み込む
    records = list(SeqIO.parse(str(in_fa), "fasta"))
    if not records:
        raise ValueError(f"{in_fa} に FASTA レコードが見つかりません")

    # 長さで降順ソート
    records.sort(key=lambda r: len(r.seq), reverse=True)

    # 連番の桁数を決める
    count = len(records)
    # width が指定されていない場合、自動計算（最低 2 桁）
    width = width or max(2, len(str(count)))

    # 新しい ID として prefix + 0-padded index を付与
    renamed = []
    for idx, rec in enumerate(records, start=1):
        new_id = f"{prefix}{idx:0{width}d}"
        # description を空にする
        renamed.append(SeqRecord(rec.seq, id=new_id, description=""))

    # 出力
    SeqIO.write(renamed, str(out_fa), "fasta")
    print(f"✓ 出力: {out_fa} （{count} 配列をリネーム, 桁幅: {width}）")

def main():
    p = argparse.ArgumentParser(
        description="FASTA を長い順にソートして連番でリネームします"
    )
    p.add_argument("--in",  "-i", dest="in_fa",  required=True, help="入力 FASTA ファイル")
    p.add_argument("--out", "-o", dest="out_fa", required=True, help="出力 FASTA ファイル")
    p.add_argument("--prefix", "-p", default="Chr",
                   help="連番前の接頭辞（デフォルト: Chr）")
    p.add_argument("--width", "-w", type=int, default=None,
                   help="連番の桁数（省略時は自動決定, 最低2桁）")
    args = p.parse_args()

    rename_by_length(Path(args.in_fa), Path(args.out_fa), args.prefix, args.width)

if __name__ == "__main__":
    main()
