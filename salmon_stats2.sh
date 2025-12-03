#!/usr/bin/env bash

R1_DIR="/projects/health_sciences/bms/biochemistry/kenny_group/1.Lincoln_Genomics_16th_Spetember_2025/Paua_salmon_mapping_round2"   # <-- change this

echo -e "sample\tnum_processed\tnum_mapped\tpercent_mapped" > round2_mapping.tsv

for meta in "$R1_DIR"/*/aux_info/meta_info.json "$R1_DIR"/*/auxinfo/meta_info.json; do
    [ -e "$meta" ] || continue

    # Go two levels up to get the true sample name
    sample=$(basename "$(dirname "$(dirname "$meta")")")

    num_processed=$(jq '.num_processed' "$meta")
    num_mapped=$(jq '.num_mapped' "$meta")
    percent_mapped=$(jq '.percent_mapped' "$meta")

    echo -e "${sample}\t${num_processed}\t${num_mapped}\t${percent_mapped}" >> round2_mapping.tsv
done

echo "Wrote round2_mapping.tsv"
