#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import argparse
import regex as re
from Bio import SeqIO
from Bio.Seq import Seq

def parse_args():
    p = argparse.ArgumentParser(
        description="ITS 프라이머 트리밍 (Indel 허용, 방향 자동 보정)"
    )
    p.add_argument("input_fasta", help="ITS FASTA")
    p.add_argument("primer_fasta", help="프라이머 FASTA (정확히 2개 서열)")
    p.add_argument("--max-mm", type=int, default=3, help="허용 에러(Indel 포함) 최대값 (기본=3)")
    p.add_argument("-o", "--out", default=None, help="출력 FASTA (기본: 입력명*_trimmed.fasta)")
    return p.parse_args()

def read_primers(primer_file):
    primers = [str(rec.seq).upper() for rec in SeqIO.parse(primer_file, "fasta")]
    if len(primers) != 2:
        print("❌ primer.fasta 파일에는 2개의 프라이머 서열이 있어야 합니다.")
        sys.exit(1)
    fwd = primers[0]
    rev = primers[1]
    rev_rc = str(Seq(rev).reverse_complement())
    return fwd, rev_rc 

def find_primer(seq, primer, max_mismatch=3, start_pos=0):
    """
    fuzzy 매칭으로 primer 검색.
    {e<=N}: Substitution + Insertion + Deletion 모두 포함하여 N개 이하 에러 허용
    BESTMATCH: 가장 에러가 적은 최적의 위치를 찾음
    """
    # 수정됨: s(치환) -> e(모든 에러)
    pat = f"({re.escape(primer)}){{e<={max_mismatch}}}"
    return re.search(pat, seq, flags=re.IGNORECASE | re.BESTMATCH, pos=start_pos)

def try_trim_on_sequence(seq, fwd_primer, rev_primer_rc, max_mismatch):
    # 1. Forward Primer 찾기
    fwd = find_primer(seq, fwd_primer, max_mismatch=max_mismatch, start_pos=0)
    if not fwd:
        return None, False
    
    # Trim Start: Fwd Primer가 끝난 지점부터
    start = fwd.end()

    # 2. Reverse Primer (RC) 찾기 (Fwd 이후 구간에서)
    rev = find_primer(seq, rev_primer_rc, max_mismatch=max_mismatch, start_pos=start)
    if not rev:
        return None, False
    
    # Trim End: Rev Primer가 시작된 지점까지
    end = rev.start()

    if start < end:
        return seq[start:end], True
    return None, False

def main():
    args = parse_args()
    
    # 출력 파일명 설정
    if args.out:
        output_fasta = args.out
    else:
        base = args.input_fasta.rsplit(".", 1)[0]
        output_fasta = f"{base}_trimmed.fasta"

    fwd_primer, rev_primer_rc = read_primers(args.primer_fasta)

    trimmed_count = 0
    total = 0
    records_to_write = []

    print(f" 분석 시작: {args.input_fasta}")

    for record in SeqIO.parse(args.input_fasta, "fasta"):
        total += 1
        seq = str(record.seq).upper()

        # 1) 정방향 시도
        trimmed_seq, ok = try_trim_on_sequence(seq, fwd_primer, rev_primer_rc, args.max_mm)
        if ok:
            record.seq = Seq(trimmed_seq)
            record.description += " direction=Fwd" # 로깅용
            records_to_write.append(record)
            trimmed_count += 1
            continue

        # 2) 역방향 시도 (Sequence 자체를 RC 변환)
        seq_rc = str(Seq(seq).reverse_complement())
        trimmed_rc, ok_rc = try_trim_on_sequence(seq_rc, fwd_primer, rev_primer_rc, args.max_mm)
        if ok_rc:
            # 수정됨: 다시 뒤집지 않고, 정방향으로 보정된 서열 그대로 저장
            record.seq = Seq(trimmed_rc)
            record.description += " direction=Rev_converted" # 로깅용
            records_to_write.append(record)
            trimmed_count += 1
            continue

    if records_to_write:
        SeqIO.write(records_to_write, output_fasta, "fasta")
    
    print("-" * 40)
    print(f"✅ 완료! 저장 파일: {output_fasta}")
    print(f" {trimmed_count} / {total} reads trimmed ({trimmed_count/total*100:.1f}%)")

if __name__ == "__main__":
    main()
