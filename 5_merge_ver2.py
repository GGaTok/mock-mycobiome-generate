#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
import glob
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def merge_pair(r1_file, r2_file, out_merged, out_stats):
    """
    _1.fastq와 _2.fastq 파일을 입력받아 병합을 수행하는 함수
    """
    if not os.path.exists(r1_file) or not os.path.exists(r2_file):
        print(f"❌ 파일 누락: {r1_file} 또는 {r2_file}")
        return

    fwd_reads = list(SeqIO.parse(r1_file, "fastq"))
    rev_reads = list(SeqIO.parse(r2_file, "fastq"))

    if len(fwd_reads) != len(rev_reads):
        print(f"⚠️ 경고: Read 개수가 다릅니다 ({r1_file}). Index 순서대로 매칭합니다.")

    merged_records = []
    merged_count = 0
    linked_count = 0

    for fwd, rev in zip(fwd_reads, rev_reads):
        # 서열 준비
        fwd_seq = str(fwd.seq)
        rev_seq_rc = str(rev.seq.reverse_complement())

        # Quality Score 준비 (Phred)
        fwd_qual = list(fwd.letter_annotations.get("phred_quality", [30] * len(fwd_seq)))
        rev_qual = list(rev.letter_annotations.get("phred_quality", [30] * len(rev.seq)))
        rev_qual_rc = list(reversed(rev_qual))

        # Overlap 계산 (단순 Suffix-Prefix 비교)
        max_overlap = min(len(fwd_seq), len(rev_seq_rc))
        overlap_len = 0
        
        # 최소 10bp 이상 겹치는지 확인
        for i in range(10, max_overlap + 1): 
            if fwd_seq[-i:] == rev_seq_rc[:i]:
                overlap_len = i

        # 병합 수행
        if overlap_len > 0:
            # 겹치는 부분 병합
            merged_seq = fwd_seq + rev_seq_rc[overlap_len:]
            merged_qual = fwd_qual + rev_qual_rc[overlap_len:]
            merged_count += 1
        else:
            # 겹치지 않으면 단순 연결
            merged_seq = fwd_seq + rev_seq_rc
            merged_qual = fwd_qual + rev_qual_rc
            linked_count += 1

        # ID 설정 (파일명에서 _1 제거 후 _merged 추가)
        # 예: Sample_A_1 -> Sample_A_merged
        merged_id = fwd.id.replace("_1", "_merged")
        
        merged_rec = SeqRecord(
            Seq(merged_seq),
            id=merged_id,
            description=""
        )
        merged_rec.letter_annotations["phred_quality"] = merged_qual
        merged_records.append(merged_rec)

    # 저장
    SeqIO.write(merged_records, out_merged, "fastq")

    # 통계 작성
    total = len(fwd_reads)
    if total == 0: total = 1
    
    with open(out_stats, "w") as f:
        f.write(f"Source: {os.path.basename(r1_file)} + {os.path.basename(r2_file)}\n")
        f.write(f"Total pairs: {total}\n")
        f.write(f"Merged (Overlapped): {merged_count} ({merged_count/total*100:.2f}%)\n")
        f.write(f"Linked (No Overlap): {linked_count} ({linked_count/total*100:.2f}%)\n")

    print(f"   ✅ [완료] {os.path.basename(out_merged)} (Merged: {merged_count/total*100:.1f}%)")


def main():
    parser = argparse.ArgumentParser(description="Paired-End 리드(_1.fastq, _2.fastq)를 일괄 병합합니다.")
    parser.add_argument("target", help="병합할 파일이 있는 '폴더 경로' 또는 '특정 파일명(_1.fastq)'")
    args = parser.parse_args()

    target_files = []

    # 1. 폴더인 경우 -> 해당 폴더의 모든 *_1.fastq 검색
    if os.path.isdir(args.target):
        print(f" 폴더 '{args.target}' 내의 모든 _1.fastq 파일을 검색합니다...")
        target_files = sorted(glob.glob(os.path.join(args.target, "*_1.fastq")))
    
    # 2. 파일인 경우
    elif os.path.isfile(args.target):
        target_files = [args.target]
    
    # 3. 패턴인 경우
    else:
        # 입력값이 "*_1.fastq" 형태로 들어오지 않았을 경우를 대비해 처리
        pattern = args.target if args.target.endswith("_1.fastq") else args.target + "*_1.fastq"
        target_files = sorted(glob.glob(pattern))

    if not target_files:
        print("❌ 처리할 파일(*_1.fastq)을 찾지 못했습니다.")
        sys.exit(1)

    print(f" 총 {len(target_files)}개의 Paired-End 세트를 병합합니다.\n")

    for r1_file in target_files:
        # 파일명 치환: _1.fastq -> _2.fastq
        if "_1.fastq" in r1_file:
            r2_file = r1_file.replace("_1.fastq", "_2.fastq")
            out_merged = r1_file.replace("_1.fastq", "_assembled.fastq")
            out_stats = r1_file.replace("_1.fastq", "_merging_stats.txt")
        else:
            print(f"⚠️ 건너뜀: 파일명에 '_1.fastq'가 포함되어 있지 않습니다 ({r1_file})")
            continue

        if not os.path.exists(r2_file):
            print(f"⚠️ 짝 파일(_2.fastq) 없음: {r2_file} (R1만 존재함)")
            continue

        merge_pair(r1_file, r2_file, out_merged, out_stats)

    print("\n 모든 병합 작업이 종료되었습니다.")

if __name__ == "__main__":
    main()
