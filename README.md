# Ovarian cancer ONT long-read RNA-seq: processing scripts

This repository contains the custom code used to process the Oxford Nanopore
full-length cDNA data of seven tissue specimens (ovary, fallopian tube and
omentum) from two ovarian cancer patients, as described in the accompanying
Data Descriptor.

Dorado, SQANTI3 and the FLAIR align, correct and quantify modules were run as
distributed, with the versions and parameters listed under Workflow. FLAIR
collapse was run with the options listed under Workflow through a patched copy
of `flair_collapse.py`. Three steps used custom code, which is in `scripts/`:

1. **Isoform collapse per chromosome** (`01_collapse_by_chrom.sh`,
   `flair_collapse_shared_annot.diff`): FLAIR collapse run chromosome by
   chromosome.
2. **Catalog merge** (`02_merge.sh`): the per-chromosome catalogs combined
   into one.
3. **Coding subset** (`03_coding_subset.py`): the isoforms of known genes that
   SQANTI3 predicts to be protein-coding, selected from the SQANTI3-filtered
   catalog.

## Data availability

- Raw reads (FASTQ): NCBI SRA SRP744628 (BioProject PRJNA1542798)
- Processed data (isoform catalog, SQANTI3 classification, count and TPM matrices): NCBI GEO GSE350874

To run the scripts, download the FASTQ files from SRA, run `flair align` and
`flair correct` as listed under Workflow, and set the locations of the FASTQ
files and the corrected BED files in `scripts/config.sh`. Isoform IDs are
derived from read names, so a rerun reproduces the IDs in GEO only if the
original read names are kept when the FASTQ files are downloaded.

## Contents

```
├── environment.yml
├── LICENSE
├── LICENSE-FLAIR
└── scripts/
    ├── config.sh                        paths, reference files, thread settings
    ├── 01_collapse_by_chrom.sh          pool and split the corrected reads, align all reads to the
    │                                    annotated transcripts once, FLAIR collapse per chromosome
    ├── flair_collapse_shared_annot.diff the change made to FLAIR's flair_collapse.py (applied by 01)
    ├── 02_merge.sh                      merge the per-chromosome catalogs
    └── 03_coding_subset.py              known-gene, protein-coding subset
```

## Workflow

| Step | Code | Tool (version) and main parameters | Output |
|---|---|---|---|
| Basecalling | Dorado | Dorado 0.9.6, HAC, `--trim --estimate-poly-a` | FASTQ per sample |
| Alignment, splice-site correction | FLAIR | FLAIR 2.0.0 `align`, `correct` (minimap2 2.24 `-ax splice --secondary=no`, MAPQ ≥ 1) | corrected BED per sample |
| Pool and split | custom, `01` step 1 | awk | pooled corrected BED, split into chr1–22, X, Y, M |
| Alignment to annotated transcripts | custom, `01` step 2 | minimap2 2.30 `-a -N 4 --MD`; FLAIR 3.0.0b1 `filter_transcriptome_align.py --stringent --check_splice` | read-to-transcript assignments |
| Isoform collapse | custom, `01` step 3 | FLAIR 3.0.0b1 `collapse --stringent --check_splice --generate_map --annotation_reliant generate`, per chromosome | 25 per-chromosome catalogs |
| Catalog merge | custom, `02` | awk / grep | merged catalog (GTF, BED, FASTA) |
| Classification, filtering | SQANTI3 | SQANTI3 5.2.2 `sqanti3_qc.py` (`--force_id_ignore`, refTSS v3.1 CAGE peaks, polyA motif list, `--isoAnnotLite`); `sqanti3_filter.py rules`, default rule set | classification table, filtered catalog |
| Quantification | FLAIR | FLAIR 3.0.0b1 `quantify --tpm` | count and TPM matrices |
| Coding subset | custom, `03` | Python (pandas): filter-passing FSM, ISM, NIC, NNC and fusion isoforms of known genes with SQANTI3 `coding` = coding | protein-coding subset |

References: GRCh38.p14 genome and GENCODE v47 annotation. `flair correct` used
the comprehensive annotation file (`chr_patch_hapl_scaff`); collapse and
SQANTI3 used the same file restricted to chr1–22, X, Y and M (the command is
given in `config.sh`). SQANTI3 used the refTSS v3.1 CAGE peaks and the polyA
motif list distributed with SQANTI3.

## Running the custom code

1. Create the environment: `conda env create -f environment.yml`.
2. Edit `scripts/config.sh`: the working directory, the reference files, the
   `bin` folder of the environment (`FLAIR_ENV`), and the directories that hold
   the FASTQ files (`<sample>.fastq.gz`) and the corrected BED files from
   `flair correct` (`<sample>_all_corrected.bed`). Adjust the thread and job
   settings to the machine if needed.
3. `bash scripts/01_collapse_by_chrom.sh`, then `bash scripts/02_merge.sh`.
   The merged catalog is the input of SQANTI3 (`Merge_OvCa_cleaned.gtf`) and
   of `flair quantify` (`Merge_OvCa.isoforms.fa`).
4. After the SQANTI3 rules filter:
   `python scripts/03_coding_subset.py <filter_result_RulesFilter_result_classification.txt>`

An interrupted run of `01_collapse_by_chrom.sh` can be restarted with the same
command: finished steps, samples and chromosomes are skipped.

## Notes

- **Collapse per chromosome.** The pooled corrected BED of all samples is too
  large for a single `flair collapse` run, so it is split by chromosome. Only
  the BED is split: every run is given the whole reference GTF and all reads,
  so that reads are not assigned to similar transcripts on another chromosome.
- **Alignment to the annotated transcripts is computed once.** This first step
  of `flair collapse` does not depend on the query BED, so step 2 of
  `01_collapse_by_chrom.sh` computes it once and step 3 runs the rest of
  `flair collapse` for each chromosome through a copy of FLAIR's
  `flair_collapse.py` that reads the result
  (`flair_collapse_shared_annot.diff` is the whole change; the script applies
  it to the installed FLAIR with `patch`). On a test subset (4 million reads,
  three chromosomes), the isoforms.bed/gtf/fa files were byte-identical to
  those of unmodified `flair collapse` run with the same options on the same
  per-chromosome BED files.
- **Merge.** `02_merge.sh` keeps, from each per-chromosome output, the isoforms
  located on that chromosome, prefixes isoform IDs with the chromosome, removes
  the 82 isoforms whose IDs are neither read names nor transcript IDs (for
  example, 4-1-0), and checks that the GTF, the BED and the FASTA hold the
  same isoforms without duplicated names.

## Software

Versions of the main tools are given in the Workflow table and in
`environment.yml`. Others: samtools 1.22.1, bedtools 2.31.1 and pigz 2.4;
Python 3.12 (FLAIR 3.0.0b1 environment); pandas 3.0.1 (coding subset). The
shell scripts use GNU awk, grep (`-P`), sed, coreutils and patch, and read
`/proc/meminfo` (Linux).

## License

The code in this repository is released under the MIT License (`LICENSE`).
`scripts/flair_collapse_shared_annot.diff` contains lines of FLAIR source code,
which remain under the FLAIR license (BSD 3-Clause, `LICENSE-FLAIR`;
https://github.com/BrooksLabUCSC/flair).
