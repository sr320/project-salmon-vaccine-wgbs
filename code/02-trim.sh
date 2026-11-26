#!/usr/bin/env bash
# Trim WGBS reads (adaptase-style library) based on FastQC review.
#   pass 1: remove Illumina adapters, 2-colour poly-G/quality trim (q20),
#           clip 15 bp from R2 5' end (G-rich adaptase tail)
#   pass 2: clip 15 bp from R1 3' end (complement of the tail, left behind
#           after adapter removal on short inserts); drop pairs < 20 bp
# Usage: code/02-trim.sh <in_dir> <out_dir> <library> [threads]
set -euo pipefail
in=$1; out=$2; lib=$3; t=${4:-8}
CUTADAPT=/home/shared/cutadapt-venv/bin/cutadapt
mkdir -p "$out"
$CUTADAPT -j "$t" --interleaved \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --nextseq-trim 20 -U 15 \
  "$in/${lib}_1.fastq.gz" "$in/${lib}_2.fastq.gz" \
  2> "$out/${lib}.pass1.log" |
$CUTADAPT -j "$t" --interleaved -u -15 -m 20 \
  -o "$out/${lib}_1.trimmed.fastq.gz" -p "$out/${lib}_2.trimmed.fastq.gz" - \
  > "$out/${lib}.pass2.log" 2>&1
