#!/usr/bin/env bash
# Deduplicate Bismark alignments and extract methylation calls.
#
# Run in two passes:
#
#   1. code/04-dedup-methylation.sh mbias
#        Deduplicates each BAM, then writes only the M-bias table
#        (per-position methylation along R1 and R2). Use it to decide how many
#        bases to ignore at each read end.
#
#   2. IGNORE_R1=2 IGNORE_R2=2 code/04-dedup-methylation.sh extract
#        Full extraction with the chosen ignore settings: per-context calls,
#        bedGraph/coverage, genome-wide CpG cytosine report, a strand-merged CpG
#        coverage file for DMR tools, and an HTML report per library.
#
# Ignore settings (bp; default 0). Trimming already removed 15 bp from the
# R2 5' end and the R1 3' end, so any residual bias should be small:
#   IGNORE_R1     5' end of R1        IGNORE_3P_R1  3' end of R1
#   IGNORE_R2     5' end of R2        IGNORE_3P_R2  3' end of R2
#
# Deduplication is skipped for a library whose deduplicated BAM already
# exists, and extraction is skipped when its merged CpG coverage file
# exists, so the script can be re-run after an interruption.
#
# Output: output/04-methylation/   Log: output/04-methylation/04-<mode>.log
set -euo pipefail

MODE=${1:-}
case "$MODE" in
  mbias|extract) ;;
  *) echo "usage: $0 mbias|extract" >&2; exit 2 ;;
esac

REPO=/home/shared/8TB_HDD_03/sr320/github/project-salmon-vaccine-wgbs
BISMARK=/home/shared/Bismark-0.24.0
SAMTOOLS=/usr/bin
export PATH="$BISMARK:$SAMTOOLS:/usr/local/bin:/usr/bin:/bin"

GENOME_DIR=$REPO/data/genome/bismark
ALN=$REPO/output/03-bismark
OUT=$REPO/output/04-methylation
DEDUP=$OUT/dedup
LIBS="SsParrT0_1 SsParrT0_2 SsParrVac_1 SsParrVac_2"

# raven is a shared machine: run at the lowest CPU priority and keep the
# core count modest. The extractor with --parallel N uses about 3N cores;
# libraries run one at a time.
renice -n 19 -p $$ >/dev/null
PARALLEL=${PARALLEL:-3}
BUFFER=${BUFFER:-20G}

IGNORE_R1=${IGNORE_R1:-0}
IGNORE_R2=${IGNORE_R2:-0}
IGNORE_3P_R1=${IGNORE_3P_R1:-0}
IGNORE_3P_R2=${IGNORE_3P_R2:-0}

mkdir -p "$DEDUP"
exec > >(tee -a "$OUT/04-$MODE.log") 2>&1
echo "== start $MODE $(date)"
echo "ignore: R1 5'=$IGNORE_R1 3'=$IGNORE_3P_R1 | R2 5'=$IGNORE_R2 3'=$IGNORE_3P_R2"

# Alignment must be complete.
grep -q "== finished" "$ALN/03-bismark.log" 2>/dev/null \
  || { echo "Bismark alignment (03-bismark.sh) has not finished; stopping."; exit 1; }
for l in $LIBS; do
  [ -s "$ALN/${l}_1.trimmed_bismark_bt2_pe.bam" ] \
    || { echo "missing alignment for $l; stopping."; exit 1; }
done

ignore_opts=(
  --ignore "$IGNORE_R1" --ignore_r2 "$IGNORE_R2"
  --ignore_3prime "$IGNORE_3P_R1" --ignore_3prime_r2 "$IGNORE_3P_R2"
)

for l in $LIBS; do
  bam=$ALN/${l}_1.trimmed_bismark_bt2_pe.bam
  ddbam=$DEDUP/${l}_1.trimmed_bismark_bt2_pe.deduplicated.bam

  # 1. Remove PCR duplicates (pairs with identical start/end positions).
  if [ ! -s "$ddbam" ]; then
    echo "deduplicating $l $(date)"
    deduplicate_bismark --paired --bam --output_dir "$DEDUP" "$bam"
  fi

  if [ "$MODE" = mbias ]; then
    # 2a. M-bias table only.
    mkdir -p "$OUT/mbias"
    echo "M-bias $l $(date)"
    bismark_methylation_extractor --mbias_only --paired-end --no_overlap \
      --parallel "$PARALLEL" --samtools_path "$SAMTOOLS" \
      -o "$OUT/mbias" "$ddbam"
    continue
  fi

  # 2b. Full extraction.
  report=$OUT/${l}.CpG.CpG_report.merged_CpG_evidence.cov.gz
  if [ -s "$report" ]; then
    echo "$l already extracted, skipping"; continue
  fi
  echo "extracting $l $(date)"
  bismark_methylation_extractor --paired-end --no_overlap --comprehensive \
    "${ignore_opts[@]}" \
    --gzip --bedGraph --buffer_size "$BUFFER" \
    --parallel "$PARALLEL" --samtools_path "$SAMTOOLS" \
    -o "$OUT" "$ddbam"

  # 3. Genome-wide CpG report (<lib>.CpG.CpG_report.txt.gz) and the
  #    strand-merged CpG coverage (one row per CpG dinucleotide) used by DSS.
  #    Run here rather than through the extractor's --cytosine_report, which
  #    in Bismark 0.24.0 cannot merge strands and would write the same
  #    genome-wide report a second time (~1 h and 435 MB per library).
  cov=$OUT/${l}_1.trimmed_bismark_bt2_pe.deduplicated.bismark.cov.gz
  coverage2cytosine --merge_CpG --gzip --genome_folder "$GENOME_DIR" \
    --dir "$OUT" -o "${l}.CpG" "$cov"

  # 4. Per-library HTML report (alignment, dedup, splitting, M-bias).
  bismark2report --dir "$OUT" \
    --alignment_report "$ALN/${l}_1.trimmed_bismark_bt2_PE_report.txt" \
    --dedup_report "$DEDUP/${l}_1.trimmed_bismark_bt2_pe.deduplication_report.txt" \
    --splitting_report "$OUT/${l}_1.trimmed_bismark_bt2_pe.deduplicated_splitting_report.txt" \
    --mbias_report "$OUT/${l}_1.trimmed_bismark_bt2_pe.deduplicated.M-bias.txt" \
    || echo "bismark2report failed for $l (continuing)"
  echo "done $l $(date)"
done

echo "== finished $MODE $(date)"
