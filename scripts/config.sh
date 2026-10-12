# Settings for 01_collapse_by_chrom.sh and 02_merge.sh. Edit the paths in this file; nothing else needs changing.
B=/path/to/project                     # working directory: outputs are written to $B/01_collapse and $B/02_merge

# Inputs, one file of each kind per sample
SAMPLES=(P1-LO P1-LF P1-RO P1-RF P2-LOC P2-OMC P2-ROC)
FASTQ_DIR=/path/to/fastq               # reads:           $FASTQ_DIR/<sample>.fastq.gz
BED_DIR=/path/to/corrected_bed         # corrected reads: $BED_DIR/<sample>_all_corrected.bed  (from `flair correct`)

REF=/path/to/GRCh38.p14.genome.fa      # GENCODE GRCh38.p14 genome sequence
# Reference annotation for FLAIR collapse (and SQANTI3): GENCODE v47, primary chromosomes only (chr1-22, X, Y, M).
# It is the GENCODE file used for `flair correct`, with the rows on patches, haplotypes and scaffolds removed:
#   awk -F'\t' '/^#/ || $1 ~ /^chr([0-9]+|X|Y|M)$/' gencode.v47.chr_patch_hapl_scaff.annotation.gtf \
#       > gencode.v47.primary_chr.annotation.gtf        (412,034 -> 385,659 transcripts)
GTF_FLAIR=/path/to/gencode.v47.primary_chr.annotation.gtf

FLAIR_ENV=/path/to/envs/flair3.0.0b1/bin      # FLAIR 3.0.0b1 (minimap2 2.30, samtools 1.22.1)

CHRS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12
      chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)

ANNOT_PIECE_GB=25        # step 2: largest read piece (compressed FASTQ size); sets the number of pieces per sample
ANNOT_THREADS=6          # step 2: threads per piece (read assignment in FLAIR is bound to one thread, more does not help)
COLLAPSE_SHARED_JOBS=13  # step 3: chromosomes run at the same time
COLLAPSE_THREADS=8       # step 3: threads per chromosome

log() { echo "[$(date '+%F %T')] $*"; }
