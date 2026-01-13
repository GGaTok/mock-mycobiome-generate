import os
import glob
import subprocess
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def get_args():
    parser = argparse.ArgumentParser(description="BLAST based sequence extractor with orientation fix.")
    
    parser.add_argument('--input_dir', required=True, help='Path to genome files')
    parser.add_argument('--output_dir', required=True, help='Path to save extracted sequences')
    parser.add_argument('--primer_db', required=True, help='Path to blast database')
    parser.add_argument('--missing_log', required=True, help='Path to save missing list (txt)')
    parser.add_argument('--rules', nargs='+', required=True, 
                        help='Extraction rules (e.g., NS1:+3000 LN7:-3000)')
    
    return parser.parse_args()

def run_blast(query_file, db_path):
    """
    지놈(Query) vs 프라이머(DB) BLAST 수행
    """
    cmd = [
        "blastn",
        "-query", query_file,
        "-db", db_path,
        "-outfmt", "6 qseqid sseqid qstart qend sstrand", # sstrand: plus/minus 확인용
        "-task", "blastn-short", 
        "-max_target_seqs", "1",
        "-word_size", "16"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout

def extract_sequence(fasta_path, hit_info, marker_name, offset, output_path):
    """
    offset: 양수면 Start 기준 Downstream, 음수면 End 기준 Upstream
    sstrand가 minus일 경우, 좌표 계산을 뒤집고 결과 서열을 Reverse Complement 함.
    """
    # 지놈 로딩
    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    target_id = hit_info['qseqid']
    
    if target_id not in record_dict:
        return False
    
    seq_record = record_dict[target_id]
    full_seq = seq_record.seq
    genome_len = len(full_seq)

    # BLAST 좌표 파싱 (1-based -> 0-based 변환)
    # qstart, qend는 strand와 상관없이 물리적 위치의 min, max로 정규화
    raw_qstart = int(hit_info['qstart'])
    raw_qend = int(hit_info['qend'])
    
    p_start = min(raw_qstart, raw_qend) - 1 # 물리적 시작 (0-based)
    p_end = max(raw_qstart, raw_qend)       # 물리적 끝
    
    strand = hit_info['sstrand'] # 'plus' or 'minus'
    
    extract_start = 0
    extract_end = 0
    final_seq = ""

    # =========================================================
    # Case 1: PLUS STRAND (정방향 매칭)
    # =========================================================
    if strand == 'plus':
        if offset > 0: 
            # NS1 (+3000): 매칭 시작점(p_start)부터 뒤로
            extract_start = p_start
            extract_end = min(genome_len, p_start + offset)
        else: 
            # LN7 (-3000): 매칭 끝점(p_end)부터 앞으로 (offset이 음수임)
            extract_end = p_end
            extract_start = max(0, extract_end + offset)
            
        final_seq = full_seq[extract_start:extract_end]

    # =========================================================
    # Case 2: MINUS STRAND (역방향 매칭) -> 좌표 뒤집기 필요
    # =========================================================
    else:
        # 역방향일 때는 물리적 p_end가 마커의 "Start(5')"이고,
        # 물리적 p_start가 마커의 "End(3')"임.
        
        if offset > 0:
            # NS1 (+3000): 마커 Start(p_end)에서 Downstream(물리적 앞쪽)으로 이동
            extract_end = p_end
            extract_start = max(0, p_end - offset)
        else:
            # LN7 (-3000): 마커 End(p_start)에서 Upstream(물리적 뒤쪽)으로 이동 (offset 음수)
            extract_start = p_start
            extract_end = min(genome_len, p_start - offset) # -(-3000) = +3000
            
        # 추출 후 Reverse Complement 수행
        temp_seq = full_seq[extract_start:extract_end]
        final_seq = temp_seq.reverse_complement()

    # 결과 저장
    # 헤더 정보에 방향성(strand) 정보 추가
    output_header = f"{target_id}_extracted_{marker_name}_{offset}_{strand}"
    
    with open(output_path, "w") as out_f:
        out_f.write(f">{output_header}\n{str(final_seq)}\n")
    
    return True

def main():
    args = get_args()
    
    # 출력 디렉토리 생성
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 규칙 파싱
    parsed_rules = []
    for rule in args.rules:
        marker, val = rule.split(':')
        parsed_rules.append((marker, int(val)))

    # 지놈 파일 목록 가져오기
    genome_files = glob.glob(os.path.join(args.input_dir, "*.fasta")) + \
                   glob.glob(os.path.join(args.input_dir, "*.fna"))
    
    missing_list = []
    processed_count = 0

    print(f"Total genomes found: {len(genome_files)}")
    print(f"Rules: {parsed_rules}")

    for genome_file in genome_files:
        filename = os.path.basename(genome_file)
        
        # 1. BLAST 실행
        blast_out = run_blast(genome_file, args.primer_db)
        
        # 2. 결과 파싱
        hits = {} 
        if blast_out:
            lines = blast_out.strip().split('\n')
            for line in lines:
                cols = line.split('\t')
                sseqid = cols[1] # Marker Name
                # 가장 높은 점수(첫번째)만 저장
                if sseqid not in hits:
                    hits[sseqid] = {
                        'qseqid': cols[0],
                        'qstart': cols[2],
                        'qend': cols[3],
                        'sstrand': cols[4] # strand 정보 추가
                    }

        # 3. 규칙 적용 및 추출
        extracted = False
        output_filename = os.path.join(args.output_dir, f"{os.path.splitext(filename)[0]}_region.fasta")
        
        for marker_name, offset in parsed_rules:
            if marker_name in hits:
                success = extract_sequence(genome_file, hits[marker_name], marker_name, offset, output_filename)
                if success:
                    extracted = True
                    processed_count += 1
                    print(f"[{filename}] Found {marker_name} ({hits[marker_name]['sstrand']}) -> Extracted.")
                    break 
        
        # 매칭 실패 시 리스트에 추가
        if not extracted:
            print(f"[{filename}] Target Missing.")
            missing_list.append(filename)

    # 4. 누락된 파일 목록 저장 (TXT 파일)
    with open(args.missing_log, "w") as f:
        if missing_list:
            for fname in missing_list:
                f.write(f"{fname}\n")
        else:
            f.write("No missing targets.\n")
    
    print("-" * 30)
    print(f"Processing complete.")
    print(f"Extracted: {processed_count}/{len(genome_files)}")
    print(f"Missing log saved to: {args.missing_log}")

if __name__ == "__main__":
    main()
