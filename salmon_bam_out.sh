# Define paths
INDEX_PREFIX="/projects/health_sciences/bms/biochemistry/kenny_group/Roseanna/Paua_sprng2024_sumr2025/Josh_test/Paua_salmon_index"
DATA_DIR="/projects/health_sciences/bms/biochemistry/kenny_group/Roseanna/Paua_sprng2024_sumr2025/Josh_test/trimmed"
OUT_DIR="/projects/health_sciences/bms/biochemistry/kenny_group/Roseanna/Paua_sprng2024_sumr2025/Josh_test/Paua_salmon_mapping"
mkdir -p $OUT_DIR

# Process each sample
for R1 in ${DATA_DIR}/*_R1.trimmed.fastq.gz; do
    R2=${R1/_R1/_R2}
    # e.g. V350293938_L01_AP303B_76_R1.trimmed.fastq.gz -> L01_AP303B_76
    STEM=$(basename "$R1" | sed -E 's/_R[12]\.trimmed\.fastq\.gz$//')
    SAMPLE=$(awk -F'_' '{print $2"_"$3"_"$4}' <<< "$STEM")
    SAMPLE_OUT_DIR="${OUT_DIR}/${SAMPLE}"
    mkdir -p $SAMPLE_OUT_DIR

    echo "Quantifying $SAMPLE with Salmon..."

    # Run Salmon and write mappings to file
    salmon quant \
        -i $INDEX_PREFIX \
        -l A \
        -1 $R1 -2 $R2 \
        -p 8 \
        --validateMappings \
        --writeMappings "${SAMPLE_OUT_DIR}/${SAMPLE}.sam" \
        -o $SAMPLE_OUT_DIR

    # Convert SAM to sorted BAM
    samtools view -@ 8 -bS "${SAMPLE_OUT_DIR}/${SAMPLE}.sam" | \
    samtools sort -@ 8 -o "${SAMPLE_OUT_DIR}/${SAMPLE}.sorted.bam"

    # Index the BAM file
    samtools index -@ 8 "${SAMPLE_OUT_DIR}/${SAMPLE}.sorted.bam"

    # Generate mapping statistics
    samtools flagstat -@ 8 "${SAMPLE_OUT_DIR}/${SAMPLE}.sorted.bam" > "${SAMPLE_OUT_DIR}/${SAMPLE}.flagstat.txt"

    # (Optional) Clean up SAM to save space
    rm "${SAMPLE_OUT_DIR}/${SAMPLE}.sam"
done

echo "Salmon quantification (with BAM + stats) complete!"
