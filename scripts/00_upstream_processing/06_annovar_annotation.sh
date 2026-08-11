#!/bin/bash
# Description: Annotate filtered somatic variants using ANNOVAR

set -euo pipefail

VCF_DIR="${VCF_DIR:-./data/mutect2_output}"
ANNOVAR_PATH="${ANNOVAR_DIR:-./tools/annovar}"
HUMANDB="${ANNOVAR_PATH}/humandb/hg38"
BUILDVER="hg38"
PROTOCOL="refGene,cytoBand,exac03,avsnp147,dbnsfp30a,clinvar_20210501"
OPERATION="g,r,f,f,f,f"
OUTDIR="${VCF_DIR}/annovar_output"
JOBS="${JOBS:-16}"

mkdir -p "$OUTDIR"

annotate_one() {
    local file="$1"
    local sample=$(basename "$file" .filtered.vcf.gz)
    local outprefix="${OUTDIR}/${sample}"

    echo "🧬 Annotating: $sample"

    perl "${ANNOVAR_PATH}/table_annovar.pl" "$file" "$HUMANDB" \
        -buildver "$BUILDVER" \
        -out "$outprefix" \
        -remove \
        -protocol "$PROTOCOL" \
        -operation "$OPERATION" \
        -nastring . \
        -vcfinput \
        -polish

    echo "✅ Done: ${outprefix}.${BUILDVER}_multianno.vcf"
}

export -f annotate_one
export ANNOVAR_PATH HUMANDB BUILDVER PROTOCOL OPERATION OUTDIR

find "$VCF_DIR" -name "*.filtered.vcf.gz" | parallel -j "$JOBS" annotate_one {}

echo "🎉 All annotations completed. Results in: $OUTDIR"
