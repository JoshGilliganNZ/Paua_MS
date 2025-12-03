#!/usr/bin/env bash
set -euo pipefail

# Usage message
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <emapper_tsv> [output_prefix]" >&2
  exit 1
fi

IN="$1"
PREFIX="${2:-emapper_geneBest}"

# Output file names
OUT_ANN="${PREFIX}.annotations.tsv"       # filtered annotations
OUT_MAP="${PREFIX}.retained_mapping.tsv"  # gene_id -> full query ID
OUT_STATS="${PREFIX}.stats.txt"           # summary stats

# Choose how to read (gz-safe)
read_cmd="cat"
[[ "$IN" =~ \.gz$ ]] && read_cmd="zcat"

# Count non-comment lines BEFORE filtering
before_lines=$($read_cmd "$IN" | awk '!/^#/ {c++} END{print c}')

# Main AWK: pick best per gene based on SCORE (desc), tie-break by EVALUE (asc)
$read_cmd "$IN" | awk -v OFS="\t" -v MAP="$OUT_MAP" -v STATS="$OUT_STATS" '
BEGIN{
  FS = OFS = "\t"
}
# Skip comments (lines starting with #)
/^#/ { next }

{
  full_q = $1          # full query ID
  evalue = $3 + 0
  score  = $4 + 0

  # Extract gene ID (strip everything from _iNN onward)
  gene = full_q
  sub(/_i[0-9]+.*/, "", gene)

  # Isoform and translation numbers
  iso = (match(full_q, /_i([0-9]+)/, m) ? m[1] : "NA")
  trn = (match(full_q, /\.p([0-9]+)/, n) ? n[1] : "NA")

  # First encounter for this gene
  if (!(gene in best_score)) {
    best_score[gene] = score
    best_evalue[gene] = evalue
    best_query[gene] = full_q
    best_line[gene] = $0
    best_iso[gene] = iso
    best_trn[gene] = trn
  } else {
    # Compare to stored best
    take = 0
    if (score > best_score[gene]) take = 1
    else if (score == best_score[gene] && evalue < best_evalue[gene]) take = 1

    if (take) {
      best_score[gene] = score
      best_evalue[gene] = evalue
      best_query[gene] = full_q
      best_line[gene] = $0
      best_iso[gene] = iso
      best_trn[gene] = trn
    }
  }

  total_rows++
}

END{
  # Write mapping (gene_id -> full query ID)
  print "#gene_id", "retained_query_id" > MAP
  kept = 0
  for (g in best_query) {
    print g, best_query[g] >> MAP
    kept++
  }
  close(MAP)

  # Collect counts of isoforms/translations picked
  for (g in best_query) {
    isoCount[best_iso[g]]++
    trnCount[best_trn[g]]++
  }

  # Write stats
  print "Total input (non-comment) lines:\t" total_rows > STATS
  print "Genes retained (unique gene IDs):\t" kept >> STATS
  print "" >> STATS
  print "[Isoform chosen among retained genes]" >> STATS
  for (iso in isoCount) print "i" iso "\t" isoCount[iso] >> STATS
  print "" >> STATS
  print "[Translation chosen among retained genes]" >> STATS
  for (t in trnCount) print "p" t "\t" trnCount[t] >> STATS
  close(STATS)

  # Emit final annotation file, replacing col1 with gene ID
  for (g in best_line) {
    nf = split(best_line[g], F, "\t")
    F[1] = g
    line = F[1]
    for (i = 2; i <= nf; i++) line = line OFS F[i]
    print line
  }
}
' | sort -k1,1 > "$OUT_ANN"

# Count non-comment lines AFTER filtering
after_lines=$(grep -vc '^#' "$OUT_ANN" || true)

# Append before/after counts to stats
{
  echo ""
  echo "Line counts (non-comment):"
  echo "Before: $before_lines"
  echo "After : $after_lines"
} >> "$OUT_STATS"

# Final messages
echo "Done."
echo "Filtered annotations : $OUT_ANN"
echo "Retained mapping     : $OUT_MAP"
echo "Stats                : $OUT_STATS"
