#!/bin/bash

# 텍스트 파일 이름 정의
LIST_FILE="pretreatment_list.txt"

# 파일이 존재하는지 확인
if [ ! -f "$LIST_FILE" ]; then
    echo "Error: $LIST_FILE 파일을 찾을 수 없습니다."
    exit 1
fi

# 한 줄씩 읽으며 작업 수행
while IFS= read -r taxon; do
    # 빈 줄은 건너뜀
    [ -z "$taxon" ] && continue

    # 파일명 조합 (예: Saccharomyces cerevisiae + _ref.zip)
    zip_filename="${taxon}_ref.zip"

    # 압축 파일이 실제로 존재하는지 확인
    if [ -f "$zip_filename" ]; then
        echo "Processing: $zip_filename..."
        
        # 압축 해제 (-d 옵션으로 별도 폴더에 풀고 싶다면 아래 주석 참고)
        # unzip -o "$zip_filename" -d "${taxon}_ref" 
        unzip -o "$zip_filename"
        
        echo "Done: $zip_filename"
    else
        echo "Warning: $zip_filename 파일이 현재 폴더에 없습니다."
    fi

done < "$LIST_FILE"
