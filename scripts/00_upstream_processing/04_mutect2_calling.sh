#!/bin/bash
# Description: Somatic variant calling using GATK Mutect2 in paired tumour-normal mode

set -euo pipefail

BASE_DIR="${BASE_DIR:-./data}"
DEDUP_DIR="${DEDUP_DIR:-$BASE_DIR/results/dedup_bam}"
OUT_DIR="${OUT_DIR:-$BASE_DIR/mutect2_output}"
LOG_DIR="${LOG_DIR:-$BASE_DIR/logs/mutect2}"

REF="${REF_FASTA:-./ref/Homo_sapiens_assembly38.fa}"
PON="${PON_VCF:-./ref/1000g_pon.hg38.vcf.gz}"
GERMLINE="${GERMLINE_VCF:-./ref/af-only-gnomad.hg38.vcf.gz}"
TMP_DIR="${TMP_DIR:-/tmp}"
JOBS="${JOBS:-8}"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP_DIR"

run_mutect2_single() {
    local tumour_bam="$1"
    local sample_id=$(basename "$tumour_bam" | sed -E 's/_TUMOUR\.sort\.dedup\.bam$//')
    local normal_bam="${DEDUP_DIR}/${sample_id}_NORMAL.sort.dedup.bam"
    local vcf_out="${OUT_DIR}/${sample_id}.vcf.gz"
    local log_file="${LOG_DIR}/${sample_id}.log"

    if [[ ! -f "$normal_bam" ]]; then
        echo "⚠️ Error: Matched normal BAM not found for $sample_id: $normal_bam"
        return 1
    fi

    exec > "$log_file" 2>&1
    echo "[$(date)] Starting Mutect2 for sample: $sample_id"

    # Extract Normal sample name from BAM header
    local normal_sample=$(samtools view -H "$normal_bam" | grep '^@RG' | head -n1 | sed -E 's/.*SM:([^\t]+).*/\1/')

    gatk Mutect2 \
        -R "$REF" \
        -I "$tumour_bam" \
        -I "$normal_bam" \
        -normal "$normal_sample" \
        --panel-of-normals "$PON" \
        --germline-resource "$GERMLINE" \
        --native-pair-hmm-threads 2 \
        --tmp-dir "$TMP_DIR" \
        -O "$vcf_out"

    echo "[$(date)] Finished Mutect2 for: $sample_id"
}

export -f run_mutect2_single
export DEDUP_DIR OUT_DIR LOG_DIR REF PON GERMLINE TMP_DIR

# Find Tumour BAMs that haven't been processed yet
unprocessed_bams=$(comm -23 \
  <(ls "$DEDUP_DIR"/*_TUMOUR.sort.dedup.bam | sort) \
  <(ls "$OUT_DIR"/*.vcf.gz 2>/dev/null | sed "s/.*\///; s/\.vcf\.gz$//" | sed "s|^|${DEDUP_DIR}/|; s|$|_TUMOUR.sort.dedup.bam|" | sort))

if [[ -z "$unprocessed_bams" ]]; then
    echo "🎉 All samples have already been processed by Mutect2."
else
    echo "$unprocessed_bams" | parallel -j "$JOBS" run_mutect2_single {}
fi
