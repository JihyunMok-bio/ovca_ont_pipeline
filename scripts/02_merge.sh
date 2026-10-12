#!/usr/bin/env bash
# Custom code 2: merge the per-chromosome collapse outputs into one catalog.
#  0. Collapse is given the whole reference GTF, so each chromosome's output also contains the annotated
#     isoforms of all other chromosomes -> from each output keep only the isoforms on that chromosome
#     (column 1 = chromosome); this leaves exactly one copy of every isoform.
#  1. GTF: concatenate in chromosome order and rename IDs to <chr>_<ID> so that IDs are unique across chromosomes
#  2. Remove malformed IDs that are neither read names (UUID) nor transcript IDs (e.g. 4-1-0)
#     -> Merge_OvCa_cleaned.gtf (input of SQANTI3)
#  3. BED / FASTA: the isoforms kept in the cleaned GTF, with the same <chr>_<ID> names (FASTA = input of quantify)
#     ID = the part of the isoform name before the first "_" (read name, ENST..., or ENST...-1)
set -euo pipefail
source "$(dirname "$0")/config.sh"

C=$B/01_collapse/Collapse
M=$B/02_merge
OUT=$M/Merge_OvCa
mkdir -p "$M"

log "1. GTF merge (from each chromosome's output, only the isoforms on that chromosome)"
> "$OUT.gtf"
for CHR in "${CHRS[@]}"; do
  FILE=$C/${CHR}.isoforms.gtf
  [ -s "$FILE" ] || { log "  [SKIP] $CHR missing"; continue; }
  awk -v chr="$CHR" 'BEGIN{OFS="\t"}
  $1 != chr { next }
  $0 ~ /transcript_id/ {
      match($0, /transcript_id "([^"]+)"/, arr)
      old_id = arr[1]; new_id = chr "_" old_id
      gsub("transcript_id \"" old_id "\"", "transcript_id \"" new_id "\"")
      gsub("gene_id \"" old_id "\"",       "gene_id \""       new_id "\"")
  }
  { print }' "$FILE" >> "$OUT.gtf"
done
log "  transcripts: $(grep -c $'\ttranscript\t' "$OUT.gtf")  (sum over the per-chromosome outputs: $(cat "$C"/*.isoforms.gtf | grep -c $'\ttranscript\t'))"

log "2. remove malformed IDs"
grep -oP 'transcript_id "[^"]+"' "$OUT.gtf" | sed 's/transcript_id "//;s/"//' | sort -u \
  | grep -P '^(chr[^_]+_)?\d+-\d+(-\d+)?$' > "$M/bad_ids.txt" || true
log "  bad IDs: $(wc -l < "$M/bad_ids.txt")"
if [ -s "$M/bad_ids.txt" ]; then
  grep -vF -f "$M/bad_ids.txt" "$OUT.gtf" > "${OUT}_cleaned.gtf"
else
  cp "$OUT.gtf" "${OUT}_cleaned.gtf"
fi
log "  cleaned transcripts: $(grep -c $'\ttranscript\t' "${OUT}_cleaned.gtf")"

log "3. BED / FASTA restricted to cleaned GTF"
KEEP=$M/keep_uuid.txt
grep -P '\ttranscript\t' "${OUT}_cleaned.gtf" | grep -oP 'transcript_id "[^"]+"' \
  | sed 's/transcript_id "//;s/"//' | sort -u > "$KEEP"       # <chr>_ID
> "$OUT.isoforms.bed"; > "$OUT.isoforms.fa"
for CHR in "${CHRS[@]}"; do
  BED_FILE=$C/${CHR}.isoforms.bed; FA_FILE=$C/${CHR}.isoforms.fa
  [ -s "$BED_FILE" ] || continue
  awk -v chr="$CHR" -v keep="$KEEP" '
  BEGIN { while ((getline line < keep) > 0) uuid_set[line] = 1 }
  $1 == chr { split($4, a, "_"); id = chr "_" a[1]; if (id in uuid_set) { $4 = id; print } }
  ' OFS='\t' "$BED_FILE" >> "$OUT.isoforms.bed"
  # the FASTA has no chromosome field, so the names of this chromosome's isoforms are read from the BED first
  awk -v chr="$CHR" -v keep="$KEEP" '
  BEGIN { while ((getline line < keep) > 0) uuid_set[line] = 1; p = 0 }
  FNR == NR { if ($1 == chr) own[$4] = 1; next }
  /^>/ { name = substr($0, 2); split(name, a, "_"); id = chr "_" a[1]
         p = ((name in own) && (id in uuid_set)); if (p) print ">" id; next }
  p { print }
  ' "$BED_FILE" "$FA_FILE" >> "$OUT.isoforms.fa"
done
rm -f "$KEEP"
log "  BED lines: $(wc -l < "$OUT.isoforms.bed")  FASTA entries: $(grep -c '^>' "$OUT.isoforms.fa")"
log "done"

# sanity checks: the GTF and the FASTA must hold the same number of isoforms, and no name may occur twice
ng=$(grep -c $'\ttranscript\t' "${OUT}_cleaned.gtf"); nf=$(grep -c '^>' "$OUT.isoforms.fa")
[ "$ng" -eq "$nf" ] || { log "ERROR: cleaned GTF $ng != FASTA $nf"; exit 1; }
nu=$(grep '^>' "$OUT.isoforms.fa" | sort -u | wc -l)
[ "$nf" -eq "$nu" ] || { log "ERROR: duplicated FASTA names ($nf entries, $nu unique)"; exit 1; }
log "check OK: cleaned GTF = FASTA = $ng isoforms, no duplicated names"
