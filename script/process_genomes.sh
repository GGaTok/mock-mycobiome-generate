#!/bin/bash

# 1. 결과 저장 폴더 생성
OUTPUT_DIR="genome_fasta"
mkdir -p "$OUTPUT_DIR"

echo "=== Genome Extraction & Renaming (Space -> Underscore) Start ==="

# 2. 현재 디렉터리의 모든 *_ref.zip 파일 탐색
for zipfile in *_ref.zip; do
    
    # 파일 존재 확인
    [ -e "$zipfile" ] || continue

    # 3. 파일명 처리
    # (1) _ref.zip 제거하여 원본 이름 추출 (예: "Escherichia coli")
    raw_name=$(basename "$zipfile" "_ref.zip")
    
    # (2) 공백(띄어쓰기)을 언더바(_)로 치환 (예: "Escherichia_coli")
    # ${변수//찾을문자열/바꿀문자열} 문법 사용
    clean_name=${raw_name// /_}
    
    echo "Processing: '$raw_name' -> Saving as: '$clean_name'"

    # 4. 임시 폴더 생성
    temp_dir="temp_extract_${clean_name}"
    mkdir -p "$temp_dir"

    # 5. 압축 해제 (-q: 조용히)
    unzip -q "$zipfile" -d "$temp_dir"

    # 6. .fna 파일 찾기 (하위 폴더 깊숙이 있는 파일 탐색)
    found_fasta=$(find "$temp_dir" -type f -name "*.fna" | head -n 1)

    if [ -n "$found_fasta" ]; then
        # 7. 파일 복사 및 이름 변경 (공백이 제거된 clean_name 사용)
        cp "$found_fasta" "${OUTPUT_DIR}/${clean_name}.fna"
        
        echo " -> [OK] Created: ${OUTPUT_DIR}/${clean_name}.fna"
    else
        echo " -> [FAIL] Fasta file not found inside $zipfile"
    fi

    # 8. 임시 폴더 삭제
    rm -rf "$temp_dir"

done

echo "=== All Tasks Finished ==="
echo "결과물은 './$OUTPUT_DIR' 폴더에 저장되었습니다."
