#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="../"
TRIMMED_DIR="trimmed"
THREADS=32

mkdir -p "${TRIMMED_DIR}"
for R1 in ${INPUT_DIR}/*_1.fq.gz; do
    # Extract the sample name by removing the _1.fq.gz suffix
    SAMPLE=$(basename "${R1}" _1.fq.gz)
    R2="${INPUT_DIR}/${SAMPLE}_2.fq.gz"

    if [[ -f "${R2}" ]]; then
        echo ">>> Running fastp on ${SAMPLE}"
        fastp \
            -i "${R1}" \
            -I "${R2}" \
            -o "${TRIMMED_DIR}/${SAMPLE}_R1.trimmed.fastq.gz" \
            -O "${TRIMMED_DIR}/${SAMPLE}_R2.trimmed.fastq.gz" \
            --detect_adapter_for_pe \
            --thread ${THREADS} \
            --html "${TRIMMED_DIR}/${SAMPLE}_fastp.html" \
            --json "${TRIMMED_DIR}/${SAMPLE}_fastp.json"
    else
        echo ">>> Skipping ${SAMPLE}: R2 not found"
    fi
done
