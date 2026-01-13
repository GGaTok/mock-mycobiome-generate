#!/bin/bash

# ==============================================================================
# 사용법 (Usage)
# bash mock_ITS_analysis.sh <Forward_Marker.fasta> [Reverse_Marker.fasta]
# 예: bash mock_ITS_analysis.sh NS1.fasta LN7.fasta
# 작업하고자 하는 dir 안에 pretreatment_list.txt 파일이 필요하며, 내가 찾고자 하는 균의 scientific name이 들어있어야함
# ex)
# Candida albicans
# Candida glabrata


# 인자 확인
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <Forward_Marker.fasta> [Reverse_Marker.fasta]"
    echo "Error: At least one marker file is required."
    exit 1
fi

FWD_FILE=$1
REV_FILE=$2

# ==============================================================================
# 설정 (Configuration)
# ==============================================================================
BASE_DIR="."
SCRIPT_DIR="$BASE_DIR/script"
GENOME_DIR="./genome_fasta"
OUTPUT_DIR="./"
LOG_FILE="./log.txt"

# 임시 통합 마커 파일 및 DB 경로
COMBINED_MARKERS="./combined_markers.fasta"
MARKER_DB="./markers_db"
PYTHON_EXTRACTOR="$SCRIPT_DIR/extract_generic.py"

echo "========================================================" | tee -a $LOG_FILE
echo "Starting Dynamic ITS Analysis: $(date)" | tee -a $LOG_FILE
echo "========================================================" | tee -a $LOG_FILE

# ==============================================================================
# Step 1: 마커 병합 및 규칙(Rules) 자동 생성
# ==============================================================================
echo "[Step 1] Processing Marker Arguments..." | tee -a $LOG_FILE

# 1. Forward Marker 처리 (첫 번째 인자) -> Offset +3000
# FASTA 헤더(>)에서 ID만 추출 (첫 번째 단어)
FWD_ID=$(grep "^>" "$FWD_FILE" | head -n1 | sed 's/>//' | awk '{print $1}')
RULES="${FWD_ID}:+3000"

echo "  -> Forward Marker: $FWD_ID (from $FWD_FILE)" | tee -a $LOG_FILE
cat "$FWD_FILE" > "$COMBINED_MARKERS"

# 2. Reverse Marker 처리 (두 번째 인자, 있을 경우) -> Offset -3000
if [ ! -z "$REV_FILE" ]; then
    REV_ID=$(grep "^>" "$REV_FILE" | head -n1 | sed 's/>//' | awk '{print $1}')
    RULES="$RULES ${REV_ID}:-3000"
    
    echo "  -> Reverse Marker: $REV_ID (from $REV_FILE)" | tee -a $LOG_FILE
    cat "$REV_FILE" >> "$COMBINED_MARKERS"
fi

echo "  -> Extraction Rules: $RULES" | tee -a $LOG_FILE


# ==============================================================================
# Step 2: 전처리 스크립트 실행 (지놈 다운로드 등)
# ==============================================================================
echo "[Step 2] Running Genome Pre-processing..." | tee -a $LOG_FILE

bash $SCRIPT_DIR/genome_down.sh
bash $SCRIPT_DIR/unzip_genomes.sh
bash $SCRIPT_DIR/process_genomes.sh


# ==============================================================================
# Step 3: BLAST DB 생성
# ==============================================================================
echo "[Step 3] Building BLAST Database..." | tee -a $LOG_FILE

makeblastdb \
    -in "$COMBINED_MARKERS" \
    -dbtype nucl \
    -out "$MARKER_DB" 
   

echo "  -> DB created at: $MARKER_DB" | tee -a $LOG_FILE


# ==============================================================================
# Step 4: Python 스크립트 실행 (동적 규칙 적용)
# ==============================================================================
echo "[Step 4] Extracting Sequences..." | tee -a $LOG_FILE

python $PYTHON_EXTRACTOR \
    --input_dir "$GENOME_DIR" \
    --output_dir "$OUTPUT_DIR" \
    --primer_db "$MARKER_DB" \
    --missing_log "$SCRIPT_DIR/missing_targets.txt" \
    --rules $RULES

if [ $? -eq 0 ]; then
    echo "  -> Success! Check output in: $OUTPUT_DIR" | tee -a $LOG_FILE
else
    echo "  -> Extraction Failed." | tee -a $LOG_FILE
    exit 1
fi


echo "[Step 5] Clean up..."


rm *.zip*
mkdir out_fasta
mv *_region.fasta out_fasta
rm -r md5sum.txt ncbi_dataset/ README.md *db.* combined_markers.fasta
