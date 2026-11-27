#!/usr/bin/env bash
# Bismark alignment of trimmed WGBS reads to Atlantic salmon Ssal_v3.1
# (NCBI RefSeq GCF_905237065.1). Directional libraries, paired-end.
#
# Steps:
#   0. wait for trimming (code/02-trim.sh) to finish
#   1. build the Bismark (bowtie2) bisulfite index, if not already built
#   2. align each library in turn
#   3. bismark2report / bismark2summary
#
# Usage: code/03-bismark.sh            (run from anywhere)
# Log:   output/03-bismark/03-bismark.log
set -euo pipefail

REPO=/home/shared/8TB_HDD_03/sr320/github/project-salmon-vaccine-wgbs
BISMARK=/home/shared/Bismark-0.24.0
BOWTIE2=/home/shared/bowtie2-2.4.4-linux-x86_64
SAMTOOLS=/usr/bin
export PATH="$BISMARK:$BOWTIE2:$SAMTOOLS:/usr/local/bin:/usr/bin:/bin"

GENOME_FA=$REPO/data/genome/GCF_905237065.1_Ssal_v3.1_genomic.fna.gz
GENOME_DIR=$REPO/data/genome/bismark
TRIM=$REPO/output/02-trim
OUT=$REPO/output/03-bismark
LIBS="SsParrT0_1 SsParrT0_2 SsParrVac_1 SsParrVac_2"

# Per-library parallelism: --parallel N runs N Bismark instances, each
# running 2 bowtie2 processes (directional) with -p threads, so roughly
# N * 2 * P cores plus overhead. 3 * 2 * 2 = 12 bowtie2 threads.
# raven is a shared machine: keep this modest and run at the lowest CPU
# priority. (The Nov 2026 alignment used PARALLEL=8, about 5 h per library;
# expect roughly twice that now.)
renice -n 19 -p $$ >/dev/null
PARALLEL=${PARALLEL:-3}
BT2_THREADS=${BT2_THREADS:-2}

# Maximum insert (fragment) length for a valid pair. Bismark's default of
# 500 bp discarded 10-17% of pairs in SsParrT0_2 and both Vac libraries,
# whose fragments average ~400 bp with a tail past 500 (checked on 100k-pair
# subsets: 1000 bp raised SsParrT0_2 unique mapping from 61.7% to 72.6%).
MAXINS=1000

mkdir -p "$OUT"
exec >> "$OUT/03-bismark.log" 2>&1
echo "== start $(date)"

# One-shot: remove this job's cron entry so it never fires again.
( crontab -l 2>/dev/null | grep -v '# salmon-wgbs-bismark-oneshot' ) | crontab - || true

# 0. Wait for trimming: every library's pass-2 cutadapt log must be complete.
for i in $(seq 1 144); do   # up to 12 h, checking every 5 min
  done_n=0
  for l in $LIBS; do
    grep -q "Total written" "$TRIM/$l.pass2.log" 2>/dev/null && done_n=$((done_n+1))
  done
  [ "$done_n" -eq 4 ] && break
  echo "waiting for trimming ($done_n/4 done) $(date)"
  sleep 300
done
[ "$done_n" -eq 4 ] || { echo "trimming not finished after 12 h; giving up"; exit 1; }
echo "trimming complete $(date)"

# 1. Genome index (Bismark needs a .fa/.fa.gz extension).
if [ ! -d "$GENOME_DIR/Bisulfite_Genome" ]; then
  mkdir -p "$GENOME_DIR"
  ln -sf "$GENOME_FA" "$GENOME_DIR/Ssal_v3.1.fa.gz"
  echo "building bisulfite index $(date)"
  bismark_genome_preparation --verbose --parallel 4 \
    --path_to_aligner "$BOWTIE2" "$GENOME_DIR"
  echo "index built $(date)"
fi

# 2. Align.
for l in $LIBS; do
  if [ -s "$OUT/${l}_1.trimmed_bismark_bt2_pe.bam" ]; then
    echo "$l already aligned, skipping"; continue
  fi
  echo "aligning $l $(date)"
  bismark --genome "$GENOME_DIR" \
    -1 "$TRIM/${l}_1.trimmed.fastq.gz" -2 "$TRIM/${l}_2.trimmed.fastq.gz" \
    --parallel "$PARALLEL" -p "$BT2_THREADS" --maxins "$MAXINS" \
    --path_to_bowtie2 "$BOWTIE2" --samtools_path "$SAMTOOLS" \
    --temp_dir "$OUT/tmp" -o "$OUT"
  echo "done $l $(date)"
done

# 3. Reports.
cd "$OUT"
bismark2report || true
bismark2summary || true
echo "== finished $(date)"
