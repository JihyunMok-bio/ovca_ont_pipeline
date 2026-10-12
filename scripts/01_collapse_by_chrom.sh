#!/usr/bin/env bash
# Custom code 1: FLAIR collapse run per chromosome, with the alignment to the annotated transcripts computed once
# (FLAIR 3.0.0b1; options --stringent --check_splice --generate_map --annotation_reliant generate)
#
#   Step 1  Pool the corrected BED files of all samples and split them by chromosome (chr1-22, X, Y, M).
#           Only the BED is split: every collapse run is given the whole reference GTF and all reads.
#           (If the GTF were cut to one chromosome, reads from other chromosomes would be assigned to
#           similar transcripts on that chromosome.)
#   Step 2  Align all reads to the annotated transcripts and assign the reads to them. This is the first
#           step of `flair collapse --annotation_reliant generate`; it depends on the reads, the GTF and the
#           alignment options but not on the query BED, so it is the same for every chromosome and is
#           computed once, with the commands FLAIR runs internally:
#             minimap2 -a -t T -N 4 --MD annotated_transcripts.fa <reads>
#               | filter_transcriptome_align.py --sam - -o <counts> -t T --quality 0 --generate_map <map> \
#                   --stringent --check_splice -i annotated_transcripts.bed
#           Reads are processed in pieces that run in parallel and are merged afterwards (the assignment
#           of a read does not depend on other reads).
#   Step 3  Run the rest of `flair collapse` for each chromosome through a copy of FLAIR's flair_collapse.py
#           with one added option, --annot_shared, which reads the result of step 2 instead of repeating
#           the alignment. flair_collapse_shared_annot.diff is the whole change; it is applied here to the
#           installed FLAIR with `patch`. Every chromosome's output also lists the annotated isoforms of
#           the other chromosomes; 02_merge.sh keeps the isoforms located on the chromosome itself.
#
# On a test subset (4 million reads, three chromosomes) the isoforms.bed/gtf/fa files were byte-identical to
# those of plain `flair collapse` with the same options.
#
# Input : set in config.sh, one file per sample
#         $BED_DIR/<sample>_all_corrected.bed   corrected reads (from `flair correct`)
#         $FASTQ_DIR/<sample>.fastq.gz          reads
# Output: $B/01_collapse/Collapse/<chr>.isoforms.{bed,gtf,fa}, <chr>.isoform.read.map.txt
# Finished parts are skipped, so an interrupted run can be restarted with the same command.
# Run it from a directory where temporary folders may be created (FLAIR writes them to the current directory).
# With chromosome names as arguments, step 3 is run for those only (e.g. bash 01_collapse_by_chrom.sh chr21).
set -euo pipefail
S=$(dirname "$(readlink -f "$0")")
source "$S/config.sh"

export PATH=$FLAIR_ENV:$PATH
FLAIR_PKG=$(python -c 'import flair, os; print(os.path.dirname(flair.__file__))')
export PATH=$FLAIR_PKG:$PATH            # filter_transcriptome_align.py
export TMPDIR=$B/tmp/collapse
O=$B/01_collapse/beds                   # step 1
A=$B/01_collapse/annot; P=$A/parts; ANNOT=$A/shared      # step 2
DRIVER=$B/01_collapse/flair_collapse_shared_annot.py    # step 3
mkdir -p "$O" "$P" "$TMPDIR" "$B/01_collapse/Collapse" "$B/logs/collapse"
for s in "${SAMPLES[@]}"; do
  for f in "$BED_DIR/${s}_all_corrected.bed" "$FASTQ_DIR/$s.fastq.gz"; do
    [ -s "$f" ] || { log "input not found: $f"; exit 1; }
  done
done

# ---------------------------------------------------------------- Step 1: pool and split by chromosome
if [ "$(ls "$O"/chr*_total_corrected.bed 2>/dev/null | wc -l)" -eq "${#CHRS[@]}" ] && [ -s "$O/total_corrected.bed" ]; then
  log "step 1 already done -> skip"
else
  log "step 1: merge samples -> $O/total_corrected.bed"
  for s in "${SAMPLES[@]}"; do cat "$BED_DIR/${s}_all_corrected.bed"; done > "$O/total_corrected.bed"
  log "step 1: split by primary chromosome -> $O/<chr>_total_corrected.bed"
  rm -f "$O"/chr*_total_corrected.bed
  awk -F'\t' -v d="$O" '$1 ~ /^chr([0-9]+|X|Y|M)$/ { print > (d "/" $1 "_total_corrected.bed") }' "$O/total_corrected.bed"
  for c in "${CHRS[@]}"; do
    f=$O/${c}_total_corrected.bed
    printf '%-6s %12s\n' "$c" "$( [ -e "$f" ] && wc -l < "$f" || echo 0 )"
  done
fi

# ---------------------------------------------------------------- Step 2: alignment to the annotated transcripts, once
align_piece() {   # stdin = FASTQ, $1 = output prefix
  set -o pipefail
  minimap2 -a -t "$ANNOT_THREADS" -N 4 --MD "$A/annotated_transcripts.fa" - 2> "$1.mm2.log" \
    | filter_transcriptome_align.py --sam - -o "$1.counts" -t "$ANNOT_THREADS" --quality 0 \
        --generate_map "$1.map" --stringent --check_splice -i "$A/annotated_transcripts.bed" > "$1.filter.log" 2>&1 \
    && touch "$1.done"
}
align_sample() {  # $1 = sample name, $2 = fastq.gz
  s=$1; fq=$2
  n=$(( ($(stat -L -c %s "$fq") + ANNOT_PIECE_GB * 1000000000 - 1) / (ANNOT_PIECE_GB * 1000000000) ))
  if [ "$(ls "$P/$s".part*.done 2>/dev/null | wc -l)" -eq "$n" ]; then echo "[skip] $s"; return 0; fi
  rm -f "$P/$s".part*
  echo "[$(date '+%F %T')] start $s ($n parts)"
  # each FASTQ record (4 lines) is joined into one line and dealt round-robin to n pieces; every piece runs align_piece
  pigz -dc "$fq" | paste -d $'\001' - - - - \
    | SHELL=/bin/bash split -n r/"$n" --numeric-suffixes=1 -a 2 --filter='tr "\001" "\n" | align_piece $FILE' - "$P/$s.part" \
    && [ "$(ls "$P/$s".part*.done 2>/dev/null | wc -l)" -eq "$n" ] \
    && echo "[$(date '+%F %T')] done  $s" \
    || { echo "[$(date '+%F %T')] FAIL  $s (see $P/$s.part*.log)"; return 1; }
}
export -f align_piece align_sample; export A P ANNOT_THREADS ANNOT_PIECE_GB

if [ -e "$ANNOT.done" ]; then
  log "step 2 already done ($ANNOT.*) -> skip"
else
  if [ ! -e "$A/annotated_transcripts.done" ]; then
    log "step 2: annotated transcripts BED/FASTA (same functions as flair collapse)"
    python - "$A" "$GTF_FLAIR" "$REF" <<'PY'
import sys
from flair.gtf_to_bed import gtf_to_bed
from flair.bed_to_sequence import bed_to_sequence
a, gtf, ref = sys.argv[1:]
gtf_to_bed(a + '/annotated_transcripts.bed', gtf, include_gene=True)
bed_to_sequence(query=a + '/annotated_transcripts.bed', genome=ref, outfilename=a + '/annotated_transcripts.fa')
PY
    touch "$A/annotated_transcripts.done"
  fi
  log "step 2: align all samples at once, pieces of at most ${ANNOT_PIECE_GB} GB (fastq.gz), ${ANNOT_THREADS} threads per piece"
  # filter_transcriptome_align.py creates its temporary directory in the current directory
  ( cd "$P" && for s in "${SAMPLES[@]}"; do printf '%s\t%s\n' "$s" "$FASTQ_DIR/$s.fastq.gz"; done \
      | xargs -P "${#SAMPLES[@]}" -L 1 bash -c 'align_sample "$0" "$1"' )

  log "step 2: merge parts -> $ANNOT.*"
  python - "$P" "$ANNOT" <<'PY'
import sys, glob, re
P, out = sys.argv[1:]
counts, reads, nseq = {}, {}, {}
for f in sorted(glob.glob(P + '/*.part*.counts')):
    for l in open(f):
        t, n = l.rstrip('\n').split('\t'); counts[t] = counts.get(t, 0) + int(n)
for f in sorted(glob.glob(P + '/*.part*.map')):
    for l in open(f):
        t, r = l.rstrip('\n').split('\t'); reads.setdefault(t, []).append(r)
assert set(counts) == set(reads), 'transcripts differ between counts and map'
with open(out + '.annotated_transcripts.alignment.counts', 'w') as fc, \
     open(out + '.annotated_transcripts.isoform.read.map.txt', 'w') as fm:
    for t in sorted(counts):
        r = ','.join(reads[t])
        assert r.count(',') + 1 == counts[t], f'{t}: count {counts[t]} != reads {r.count(",") + 1}'
        fc.write(f'{t}\t{counts[t]}\n'); fm.write(f'{t}\t{r}\n')
for f in sorted(glob.glob(P + '/*.part*.mm2.log')):     # reads seen by minimap2, per sample
    s = re.sub(r'\.part\d+\.mm2\.log$', '', f.split('/')[-1])
    nseq[s] = nseq.get(s, 0) + sum(int(m.group(1)) for m in re.finditer(r'mapped (\d+) sequences', open(f).read()))
for s in nseq: print(f'  {s}: reads {nseq[s]:,}')
print(f'  all reads {sum(nseq.values()):,} | reads assigned to annotated transcripts {sum(counts.values()):,} | '
      f'transcripts with reads {len(counts):,} (3 or more: {sum(v >= 3 for v in counts.values()):,})')
PY
  touch "$ANNOT.done"
fi

# ---------------------------------------------------------------- Step 3: collapse per chromosome
if [ ! -s "$DRIVER" ]; then
  log "step 3: create $DRIVER (installed flair_collapse.py + flair_collapse_shared_annot.diff)"
  patch -s -o "$DRIVER" "$FLAIR_PKG/flair_collapse.py" "$S/flair_collapse_shared_annot.diff"
fi

READ_LIST=$(for s in "${SAMPLES[@]}"; do echo "$FASTQ_DIR/$s.fastq.gz"; done | paste -sd, -)
export B REF GTF_FLAIR READ_LIST COLLAPSE_THREADS TMPDIR ANNOT DRIVER

# Number of chromosomes run at once: each job holds the names of all reads assigned to annotated transcripts
# in memory (about 150 bytes per read); the estimated total is kept below 40% of the available memory.
nreads=$(awk -F'\t' '{s+=$2} END{printf "%d", s}' "$ANNOT.annotated_transcripts.alignment.counts")
per_gb=$(( nreads * 150 / 1000000000 + 6 ))
avail_gb=$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo)
JOBS=$COLLAPSE_SHARED_JOBS
max_jobs=$(( avail_gb * 40 / 100 / per_gb ))
[ "$max_jobs" -lt "$JOBS" ] && JOBS=$max_jobs
[ "$JOBS" -lt 2 ] && JOBS=2
log "step 3: assigned reads $nreads, about ${per_gb} GB per job, ${avail_gb} GB available -> $JOBS jobs at once"

run_one() {
  c=$1
  out=$B/01_collapse/Collapse/$c
  bed=$B/01_collapse/beds/${c}_total_corrected.bed
  if [ -s "$out.isoforms.gtf" ]; then echo "[skip] $c"; return 0; fi
  [ -s "$bed" ] || { echo "[empty] $c"; return 0; }
  echo "[$(date '+%F %T')] start $c"
  rm -f "$out".*          # intermediate files left by an interrupted run
  mkdir -p "$TMPDIR/$c"
  python "$DRIVER" -g "$REF" --gtf "$GTF_FLAIR" -q "$bed" -r "$READ_LIST" \
      --stringent --check_splice --generate_map --annotation_reliant generate \
      --annot_shared "$ANNOT" \
      --temp_dir "$TMPDIR/$c" \
      --output "$out" --thread "$COLLAPSE_THREADS" \
      > "$B/logs/collapse/$c.log" 2>&1 \
    && echo "[$(date '+%F %T')] done  $c" \
    || { echo "[$(date '+%F %T')] FAIL  $c (see logs/collapse/$c.log)"; return 1; }
}
export -f run_one

if [ $# -gt 0 ]; then
  printf '%s\n' "$@" | xargs -P "$JOBS" -I{} bash -c 'run_one {}'
  log "done: $*"; exit 0
fi
# Largest chromosomes (by BED size) first. Failed chromosomes are tried once more (finished ones are skipped).
for pass in 1 2; do
  ls -S "$B"/01_collapse/beds/chr*_total_corrected.bed \
    | xargs -n1 basename | sed 's/_total_corrected.bed//' \
    | xargs -P "$JOBS" -I{} bash -c 'run_one {}' || true
  n=$(ls "$B"/01_collapse/Collapse/*.isoforms.gtf 2>/dev/null | wc -l)
  [ "$n" -eq "${#CHRS[@]}" ] && break
  log "pass $pass: $n / ${#CHRS[@]} finished"
done

n=$(ls "$B"/01_collapse/Collapse/*.isoforms.gtf 2>/dev/null | wc -l)
log "collapse finished: $n / ${#CHRS[@]} chromosomes"
[ "$n" -eq "${#CHRS[@]}" ] || { log "some chromosomes failed: check logs/collapse/ and run again"; exit 1; }
