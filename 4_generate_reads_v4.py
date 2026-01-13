#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import random
import numpy as np
import multiprocessing
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import write as SeqIO_write
from pathlib import Path

# ==============================================================================
# 1. MiSeq 품질/오류 모델링 함수 (기존 동일)
# ==============================================================================
MISEQ_MEAN_Q = np.array(
    [40] * 10 + list(np.linspace(40, 32, 200)) + list(np.linspace(32, 28, 40)) + [28] * 100
)
MISEQ_STD_DEV = np.array(
    [1] * 10 + list(np.linspace(1, 3, 200)) + list(np.linspace(3, 4, 40)) + [4] * 100
)

def get_miseq_qual_scores_kde_style(length):
    scores = []
    for i in range(length):
        pos_idx = min(i, len(MISEQ_MEAN_Q) - 1)
        mean_q = MISEQ_MEAN_Q[pos_idx]
        std_dev = MISEQ_STD_DEV[pos_idx]
        q_score = np.random.normal(loc=mean_q, scale=std_dev)
        q_score = int(np.round(np.clip(q_score, 2, 41)))
        scores.append(q_score)
    return scores

def apply_sequencing_error(sequence, qual_scores):
    bases = ['A', 'T', 'C', 'G']
    seq_list = list(sequence)
    for i in range(len(seq_list)):
        q = qual_scores[i]
        p_error = 10**(-q / 10)
        if random.random() < p_error:
            true_base = seq_list[i]
            error_base = random.choice([b for b in bases if b != true_base])
            seq_list[i] = error_base
    return "".join(seq_list)

# ==============================================================================
# 2. 단일 샘플 생성 함수 (수정됨: 다중 Region 처리)
# ==============================================================================
def generate_single_sample(args):
    """
    하나의 샘플에 대해 '고정된 군집 비율'을 생성한 뒤,
    입력된 모든 Region(V3V4, V1V3 등)에 대해 리드를 생성합니다.
    """
    # region_data_map: { "filename": { "ID": SeqObject, ... }, ... }
    region_data_map, common_ids, sample_idx, reads_per_sample, read_len, alpha = args

    # 시드 설정 (샘플별 고유 시드)
    random.seed(os.getpid() + sample_idx)
    np.random.seed(os.getpid() + sample_idx)

    num_organisms = len(common_ids)
    
    # ---------------------------------------------------------
    # [핵심 변경] 군집 비율(Composition)을 먼저 한 번만 정의
    # ---------------------------------------------------------
    # 이 비율은 모든 Region(V3-V4, V1-V3 등)에 공통으로 적용됨
    proportions = np.random.dirichlet([alpha] * num_organisms)
    proportions = np.maximum(proportions, 0.05) # 0이 되는 것 방지
    proportions /= proportions.sum()
    
    # ID 별 비율 매핑: { "Bacteria_A": 0.2, "Bacteria_B": 0.1, ... }
    id_to_prop = {id_: prop for id_, prop in zip(common_ids, proportions)}

    # Summary 파일 작성 (이 샘플의 정답 비율)
    output_summary = f"Sample_{sample_idx}_Composition.txt"
    with open(output_summary, "w") as f:
        f.write("Taxa_ID\tGlobal_Percentage\n")
        for id_, prop in id_to_prop.items():
            f.write(f"{id_}\t{prop*100:.4f}%\n")

    # ---------------------------------------------------------
    # 각 Region 파일(V3V4, V1V3...) 별로 리드 생성
    # ---------------------------------------------------------
    for region_name, seq_dict in region_data_map.items():
        
        # 해당 Region에서의 리드 수 계산 (비율 * 총 리드 수)
        # 반올림 오차 보정을 위해 정수 변환 후 차이 보정
        reads_counts = []
        target_ids = []
        
        current_total = 0
        for id_ in common_ids:
            if id_ not in seq_dict: continue # 혹시 모를 예외 처리
            
            count = int(id_to_prop[id_] * reads_per_sample)
            reads_counts.append(count)
            target_ids.append(id_)
            current_total += count
            
        # 남는 리드(반올림 오차)는 첫 번째 균주에 더해줌
        diff = reads_per_sample - current_total
        if reads_counts:
            reads_counts[0] += diff

        reads_all_1 = []
        reads_all_2 = []
        read_counter = 1

        # 실제 리드 생성 루프
        for id_, n_reads in zip(target_ids, reads_counts):
            if n_reads <= 0: continue
            
            # 해당 Region의 서열 가져오기
            rec = seq_dict[id_]
            seq_str = str(rec.seq)
            seq_len = len(seq_str)

            # 서열 길이에 따른 처리 (기존 로직 유지)
            if seq_len < read_len:
                true_seq_1 = seq_str
                true_seq_2 = str(Seq(seq_str).reverse_complement())
            else:
                true_seq_1 = seq_str[:read_len]
                true_seq_2 = str(Seq(seq_str[-read_len:]).reverse_complement())

            for _ in range(n_reads):
                # R1
                len_1 = len(true_seq_1)
                q1 = get_miseq_qual_scores_kde_style(len_1)
                obs_seq1 = apply_sequencing_error(true_seq_1, q1)
                
                h1 = f"{id_}_{read_counter:06d}_1"
                reads_all_1.append(SeqRecord(Seq(obs_seq1), id=h1, description="", letter_annotations={"phred_quality": q1}))

                # R2
                len_2 = len(true_seq_2)
                q2 = get_miseq_qual_scores_kde_style(len_2)
                obs_seq2 = apply_sequencing_error(true_seq_2, q2)
                
                h2 = f"{id_}_{read_counter:06d}_2"
                reads_all_2.append(SeqRecord(Seq(obs_seq2), id=h2, description="", letter_annotations={"phred_quality": q2}))

                read_counter += 1

        # 파일명 생성: Sample_1_V3V4_R1.fastq
        # region_name은 파일명(확장자제외)을 사용
        out_r1 = f"Sample_{sample_idx}_{region_name}_1.fastq"
        out_r2 = f"Sample_{sample_idx}_{region_name}_2.fastq"
        
        SeqIO_write(reads_all_1, out_r1, "fastq")
        SeqIO_write(reads_all_2, out_r2, "fastq")
        
        print(f"   [Region: {region_name}] Saved -> {out_r1}, {out_r2}")

    print(f"✅ Sample {sample_idx} Done (Composition Saved).")

# ==============================================================================
# 3. 메인 프로세스
# ==============================================================================
def main_process(fasta_files, num_samples, reads_per_sample, read_len, num_cores, alpha):
    
    # 1. 모든 FASTA 파일 로드 및 ID 매핑
    print(" FASTA 파일 로딩 중...")
    
    region_data_map = {} # {'file_prefix': {ID: SeqRecord}}
    common_ids = set()
    first_file = True

    for fpath in fasta_files:
        path_obj = Path(fpath)
        region_name = path_obj.stem # 파일명에서 확장자 뺀 부분 (예: v3v4)
        
        record_dict = SeqIO.to_dict(SeqIO.parse(fpath, "fasta"))
        current_ids = set(record_dict.keys())
        
        if not current_ids:
            print(f"❌ Error: {fpath} 에 시퀀스가 없습니다.")
            sys.exit(1)

        if first_file:
            common_ids = current_ids
            first_file = False
        else:
            # 교집합(Intersection)만 남김: 모든 파일에 다 있는 균주만 시뮬레이션 가능
            common_ids = common_ids.intersection(current_ids)
            
        region_data_map[region_name] = record_dict
        print(f"   - {region_name}: {len(current_ids)} seqs loaded.")

    common_ids = sorted(list(common_ids)) # 리스트로 변환 및 정렬 (순서 고정)
    
    if not common_ids:
        print("❌ Error: 입력된 모든 파일에 공통으로 존재하는 ID가 하나도 없습니다.")
        print("   각 파일의 Fasta Header(ID)가 동일한 균주를 가리키는지 확인하세요.")
        sys.exit(1)
        
    print(f"✨ 공통으로 시뮬레이션할 균주(ID) 수: {len(common_ids)}개")

    # 2. 멀티프로세싱 인자 준비
    # region_data_map이 꽤 클 수 있지만, 리눅스 fork 방식에서는 카피 온 라이트로 효율적 처리됨
    # 윈도우에서는 오버헤드가 있을 수 있으나 일반적으로 수천 개 수준은 문제 없음
    args_list = []
    for i in range(1, num_samples + 1):
        args_list.append((region_data_map, common_ids, i, reads_per_sample, read_len, alpha))

    print(f" {num_samples}개 샘플에 대해 시뮬레이션 시작 (Cores: {num_cores})...")
    
    with multiprocessing.Pool(processes=num_cores) as pool:
        pool.map(generate_single_sample, args_list)

    print(" 모든 작업이 완료되었습니다.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="동일한 군집 비율을 유지하며 여러 Primer Region(Input Files)에 대한 MiSeq 데이터를 생성합니다."
    )

    # 수정됨: 여러 개의 파일을 입력받음
    parser.add_argument(
        "--input_fastas", 
        nargs='+', 
        required=True,
        help="사용할 Trimmed FASTA 파일들 (예: v3v4.fasta v1v3.fasta)"
    )

    parser.add_argument("-n", "--num_samples", type=int, required=True, help="샘플 수")
    parser.add_argument("-r", "--reads_per_sample", type=int, required=True, help="샘플/Region 당 리드 수")
    parser.add_argument("-l", "--read_length", type=int, required=True, help="리드 길이 (bp)")
    parser.add_argument("-c", "--num_cores", type=int, required=True, help="CPU 코어 수")
    parser.add_argument("-a", "--alpha", type=float, default=0.5, help="Dirichlet Alpha (기본 0.5)")

    args = parser.parse_args()
    
    main_process(
        args.input_fastas,
        args.num_samples,
        args.reads_per_sample,
        args.read_length,
        args.num_cores,
        args.alpha
    )
