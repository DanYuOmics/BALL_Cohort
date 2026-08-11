#!/bin/bash
# Description: Filter raw somatic mutations using GATK FilterMutectCalls

set -euo pipefail

OUTDIR="${OUTDIR:-./data/mutect2_output}"
REF="${REF_FASTA:-./ref/Homo_sapiens_assembly38.fa}"
JOBS="${JOBS:-16}"

filter_single() {
    local vcf_in="$1"
    local sample=$(basename "$vcf_in" .vcf.gz)
    local vcf_out="${OUTDIR}/${sample}.filtered.vcf.gz"

    echo "🔍 Filtering: $sample..."
    gatk FilterMutectCalls \
        -V "$vcf_in" \
        -R "$REF" \
        -O "$vcf_out"
    echo "✅ Finished: $sample.filtered.vcf.gz"
}

export -f filter_single
export OUTDIR REF

find "$OUTDIR" -name "*.vcf.gz" ! -name "*.filtered.vcf.gz" | parallel -j "$JOBS" filter_single {}
