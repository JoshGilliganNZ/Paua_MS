awk -F'\t' '$0 ~ /GO:0009408|GO:0034605|GO:0009266/' emapper_geneBest.annotations.tsv  > by_heat_only.tsv
