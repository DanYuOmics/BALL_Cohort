#!/bin/bash
# Description: High-throughput MiXCR alignment, assembly, and clone export for V(D)J repertoire profiling

set -euo pipefail

DATA_DIR="${1:-./data}"
THREADS="${THREADS:-16}"
MAX_JOBS="${MAX_JOBS:-8}"
LOGFILE="${DATA_DIR}/mixcr_parallel.log"

echo "🚀 Starting MiXCR parallel processing at $(date)" > "$LOGFILE"

job_list=$(find "$DATA_DIR" -maxdepth 2 -name "*_clean_1.fq.gz" | sed 's/_clean_1.fq.gz//')

run_mixcr() {
    local prefix="$1"
    local r1="${prefix}_clean_1.fq.gz"
    local r2="${prefix}_clean_2.fq.gz"
    local out="${prefix}_mixcr"
    mkdir -p "$out"

    echo "Running MiXCR pipeline on: $prefix" >> "$LOGFILE"

    # Step 1: Align RNA-seq VDJ reads
    mixcr align -p rna-seq -s hsa --report "$out/alignment_report.log" \
                "$r1" "$r2" "$out/alignments.vdjca" --threads "$THREADS" >> "$LOGFILE" 2>&1

    # Step 2: Assemble clonotypes
    mixcr assemble --report "$out/assembly_report.log" \
                   "$out/alignments.vdjca" "$out/clones.clns" >> "$LOGFILE" 2>&1

    # Step 3: Export clones for all TCR/BCR chains
    mixcr exportClones --chains IGH,IGK,IGL,TRA,TRB,TRD,TRG \
                       "$out/clones.clns" "$out/clones.tsv" >> "$LOGFILE" 2>&1
}

export -f run_mixcr
export LOGFILE THREADS

echo "$job_list" | parallel -j "$MAX_JOBS" run_mixcr {}

echo "🎉 MiXCR processing completed at $(date)" >> "$LOGFILE"
