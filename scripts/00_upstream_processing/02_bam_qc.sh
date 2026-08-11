#!/bin/bash
# Description: Rapid integrity check for generated BAM files using samtools

set -euo pipefail

BAM_DIR="${1:-./data/results/bam}"
FAILED=0

echo "🔍 Checking BAM integrity in: $BAM_DIR"

for bam in "$BAM_DIR"/*.bam; do
    [[ -e "$bam" ]] || continue
    if samtools quickcheck "$bam"; then
        echo "✅ OK: $bam"
    else
        echo "❌ ERROR: Integrity check failed for $bam"
        ((FAILED++))
    fi
done

if [ "$FAILED" -eq 0 ]; then
    echo "🎉 All BAM files passed integrity check!"
else
    echo "⚠️ $FAILED BAM file(s) failed integrity check."
    exit 1
fi
