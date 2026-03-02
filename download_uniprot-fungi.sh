mkdir -p refdb

BASE="https://rest.uniprot.org/uniprotkb/search"
QUERY="(taxonomy_id:4751)"                 # or "(taxonomy_id:4751 AND reviewed:true)"
OUTGZ="refdb/uniprotkb_saccharomycetales.fasta.gz"

: > "$OUTGZ"

# First request uses -G + --data-urlencode to build the query URL safely.
hdr=$(mktemp)
curl -sS -D "$hdr" -G "$BASE" \
  --data-urlencode "query=$QUERY" \
  --data-urlencode "format=fasta" \
  --data-urlencode "size=500" \
  --data-urlencode "compressed=true" \
  >> "$OUTGZ"

# Next pages: follow rel="next" from the Link header (contains cursor and size)
next=$(grep -i '^link:' "$hdr" | sed -n 's/.*<\([^>]*\)>; rel="next".*/\1/p')
rm "$hdr"

while [ -n "$next" ]; do
  hdr=$(mktemp)
  curl -sS -D "$hdr" "$next" >> "$OUTGZ"
  next=$(grep -i '^link:' "$hdr" | sed -n 's/.*<\([^>]*\)>; rel="next".*/\1/p')
  rm "$hdr"
done

# Decompress for DIAMOND
gunzip -c "$OUTGZ" > refdb/uniprotkb_fungi.fasta

# Quick sanity check: count sequences
grep -c '^>' refdb/uniprotkb_fungi.fasta
