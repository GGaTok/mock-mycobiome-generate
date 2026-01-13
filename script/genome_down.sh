#!/bin/bash

input_file="pretreatment_list.txt"
fail_log="download_failures.txt"

# 1. 실패 목록 파일 초기화 (기존 내용 지우기)
> "$fail_log"

echo "=== Genome Download Start ==="

# 2. 리스트를 한 줄씩 읽으며 다운로드 시도
while IFS= read -r line || [ -n "$line" ]; do
    # 빈 줄이 있으면 건너뜀
    [[ -z "$line" ]] && continue

    echo "Downloading: $line ..."

    # 3. datasets 명령어 실행
    # if ! 명령어; then ... fi 구문은 명령어가 에러를 뱉으면(실패하면) 내부를 실행합니다.
    if ! datasets download genome taxon "$line" --reference --include genome,gff3,protein --filename "${line}_ref.zip"; then
        
        # 실패 시 실행되는 부분
        echo " -> [FAIL] $line 다운로드 실패!"
        echo "$line" >> "$fail_log"
        
    else
        # 성공 시 실행되는 부분
        echo " -> [OK] $line 다운로드 성공."
    fi

done < "$input_file"

echo "=== Download Finished ==="

# 4. 결과 요약
if [ -s "$fail_log" ]; then
    count=$(wc -l < "$fail_log")
    echo "총 ${count}개의 지놈 다운로드에 실패했습니다."
    echo "실패 목록 파일: $fail_log"
else
    echo "모든 지놈이 성공적으로 다운로드되었습니다."
    rm "$fail_log" # 실패한 게 없으면 로그 파일 삭제
fi
