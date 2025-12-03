#!/bin/bash
set -euo pipefail
# Number of threads to use
THREADS=8

# Input directory with raw BAMs (from alignment)
INPUT_DIR="/projects/health_sciences/bms/biochemistry/kenny_group/Roseanna/Paua_sprng2024_sumr2025/Josh_test/Paua_salmon_mapping/all_bams"

# Output directory for intermediate and final files
OUTDIR="/projects/health_sciences/bms/biochemistry/kenny_group/Roseanna/Paua_sprng2024_sumr2025/Josh_test/dedup_bams"

mkdir "$OUTDIR"

find "${INPUT_DIR}" -name "*.sorted.bam" | while read bamfile; do
    filename=$(basename "$bamfile")
    sample="${filename%.sorted.bam}"

    echo "==> Processing $sample"

    ##### Step 1: Name-sort for fixmate #####
    echo " → Name-sorting"
    samtools sort -@ ${THREADS} -n "$bamfile" -o "${OUTDIR}/${sample}.namesort.bam"

    ##### Step 2: Fixmate #####
    echo " → Running fixmate"
    samtools fixmate -m "${OUTDIR}/${sample}.namesort.bam" "${OUTDIR}/${sample}.fixmate.bam"

    ##### Step 3: Coordinate-sort again #####
    echo " → Coordinate-sorting after fixmate"
    samtools sort -@ ${THREADS} "${OUTDIR}/${sample}.fixmate.bam" -o "${OUTDIR}/${sample}.fixmate.sorted.bam"

    ##### Step 4: Add read groups #####
    echo " → Adding read groups"
    picard AddOrReplaceReadGroups \
        I="${OUTDIR}/${sample}.fixmate.sorted.bam" \
        O="${OUTDIR}/${sample}.rg.bam" \
        RGID="$sample" \
        RGLB="lib1" \
        RGPL="ILLUMINA" \
        RGPU="unit1" \
        RGSM="$sample" \
        VALIDATION_STRINGENCY=SILENT

    ##### Step 5: Mark duplicates #####
    echo " → Removing duplicates"
    picard MarkDuplicates \
        I="${OUTDIR}/${sample}.rg.bam" \
        O="${OUTDIR}/${sample}.rmdup.bam" \
        M="${OUTDIR}/${sample}.rmdup.metrics.txt" \
        REMOVE_DUPLICATES=true \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=SILENT

    ##### Step 6: Index BAM #####
    echo " → Indexing final BAM"
    samtools index -@ ${THREADS} "${OUTDIR}/${sample}.rmdup.bam"

    ##### Step 7: Mapping stats #####
    echo " → Running flagstat"
    samtools flagstat -@ ${THREADS} "${OUTDIR}/${sample}.rmdup.bam" > "${OUTDIR}/${sample}.flagstat.txt"

    echo " ✅ Done: $sample"
    echo
done

echo "🏁 All BAMs processed with fixmate + dedup."
