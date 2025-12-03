awk -F'\t' 'BEGIN{IGNORECASE=1}
{
  go=$0
  if (go ~ /GO:0009408|GO:0034605|GO:0009266|GO:0006986|GO:0006457|GO:0061077|GO:0031072|GO:0051087|GO:0033554|GO:0034620/)
    print
}' emapper_geneBest.annotations.tsv  > by_GO_heat.tsv
