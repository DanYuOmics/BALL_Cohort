#!/bin/bash
# Description: Mark and remove PCR duplicates using GATK MarkDuplicatesSpark

set -euo pipefail

BASE_DIR="${BASE_DIR:-./data}"
BAM_DIR="${BAM_DIR:-$BASE_DIR/results/bam}"
DEDUP_DIR="${DEDUP_DIR:-$BASE_DIR/results/dedup_bam}"
LOG_DIR="${LOG_DIR:-$BASE_DIR/logs/markduplicates}"
TMP_DIR="${TMP_DIR:-/tmp}"
MAX_JOBS="${MAX_JOBS:-10}"

mkdir -p "$DEDUP_DIR" "$LOG_DIR" "$TMP_DIR"

export BAM_DIR DEDUP_DIR LOG_DIR TMP_DIR

run_markdup() {
    local file="$1"
    local sample_name=$(basename "$file" .bam)
    local out_bam="$DEDUP_DIR/${sample_name}.sort.dedup.bam"
    local log_file="$LOG_DIR/${sample_name}.log"

    echo "🔁 Running MarkDuplicatesSpark for: $sample_name"

    gatk MarkDuplicatesSpark \
        -I "$file" \
        -O "$out_bam" \
        -M "$DEDUP_DIR/${sample_name}.sort.dedup.metrics" \
        --conf "spark.executor.memory=16g" \
        --conf "spark.driver.memory=16g" \
        --remove-sequencing-duplicates true \
        --tmp-dir "$TMP_DIR" \
        > "$log_file" 2>&1

    echo "✅ Finished MarkDuplicates: $sample_name"
}

export -f run_markdup

find "$BAM_DIR" -name "*.bam" | parallel -j "$MAX_JOBS" --eta run_markdup {}
