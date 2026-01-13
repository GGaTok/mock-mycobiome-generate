#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import random
import sys
from Bio import SeqIO

def sample_fasta(input_fasta, n=10):
    records = list(SeqIO.parse(input_fasta, "fasta"))

    if len(records) <= n:
        print(f"입력된 FASTA에 {len(records)}개 시퀀스만 존재 — 전부 복사합니다.")
        sampled_records = records
    else:
        sampled_records = random.sample(records, n)

    output_fasta = input_fasta.replace(".fasta", "_sampled.fasta")

    SeqIO.write(sampled_records, output_fasta, "fasta")
    print(f" {len(sampled_records)}개 시퀀스를 {output_fasta}에 저장했습니다.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("사용법: python sample_fasta.py <ITS_result.fasta> [샘플 개수]")
        sys.exit(1)

    input_fasta = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 10

    sample_fasta(input_fasta, n)

