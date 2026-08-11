#!/bin/bash
# Description: Align paired-end FASTQ reads to human reference (hg38) using BWA-MEM
# Usage: bash 01_bwa_alignment.sh [sample_list.txt]

set -euo pipefail

# Configurable environment variables 
DATA_DIR="${DATA_DIR:-./data}"
OUT_DIR="${OUT_DIR:-${DATA_DIR}/results/bam}"
REF="${REF_FASTA:-./ref/Homo_sapiens_assembly38.fa}"
THREADS="${THREADS:-4}"
JOBS="${JOBS:-12}"

mkdir -p "$OUT_DIR" "${DATA_DIR}/logs"

SAMPLE_LIST="${1:-}"

process_sample() {
    local sample_dir="$1"
    local sample_id=$(basename "$sample_dir")

    local tumour_1=$(ls "${sample_dir}"/*_TUMOUR_1.fq.gz 2>/dev/null | head -n1)
    local tumour_2=$(ls "${sample_dir}"/*_TUMOUR_2.fq.gz 2>/dev/null | head -n1)
    local normal_1=$(ls "${sample_dir}"/*_NORMAL_1.fq.gz 2>/dev/null | head -n1)
    local normal_2=$(ls "${sample_dir}"/*_NORMAL_2.fq.gz 2>/dev/null | head -n1)

    if [[ -z "$tumour_1" || -z "$normal_1" ]]; then
        echo "⚠️ Skipping $sample_id: FASTQ files not found."
        return 0
    fi

    echo "⚡ Processing alignment for: $sample_id"

    # Align Tumour
    local RG_T="@RG\tID:${sample_id}_TUMOUR\tSM:${sample_id}_TUMOUR\tPL:ILLUMINA"
    bwa mem -t "$THREADS" -K 100000000 -R "$RG_T" "$REF" "$tumour_1" "$tumour_2" \
        | samtools view -Shb - > "${OUT_DIR}/${sample_id}_TUMOUR.bam"

    # Align Normal
    local RG_N="@RG\tID:${sample_id}_NORMAL\tSM:${sample_id}_NORMAL\tPL:ILLUMINA"
    bwa mem -t "$THREADS" -K 100000000 -R "$RG_N" "$REF" "$normal_1" "$normal_2" \
        | samtools view -Shb - > "${OUT_DIR}/${sample_id}_NORMAL.bam"

    echo "✅ Finished: $sample_id"
}

export -f process_sample
export OUT_DIR REF THREADS

if [[ -n "$SAMPLE_LIST" && -f "$SAMPLE_LIST" ]]; then
    echo "Processing specified sample list: $SAMPLE_LIST"
    cat "$SAMPLE_LIST" | parallel -j "$JOBS" 'process_sample "'"${DATA_DIR}"'/{}"'
else
    echo "Processing all samples under: $DATA_DIR"
    find "$DATA_DIR" -maxdepth 1 -type d -name "DDN*" | parallel -j "$JOBS" process_sample {}
fi
