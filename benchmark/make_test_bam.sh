#!/bin/bash
# Test BAM for the benchmark: NA12878 rel7 reads (40x, minimap2, hg19) aligned to one chromosome.
# Run on the server node that holds the full BAM (n14):
#   bash make_test_bam.sh /mnt/sda/home/andrey/ont/rel7_hg19/rel7_on_hg19_minimap2_mapont_sorted.40.bam chr22 8
set -euo pipefail

FULL_BAM=$1
CHROM=${2:-chr22}
THREADS=${3:-8}
OUT=NA12878_rel7_hg19_40x_$CHROM.bam

samtools view -@ "$THREADS" -b -o "$OUT" "$FULL_BAM" "$CHROM"
samtools index "$OUT"
samtools idxstats "$OUT" | awk -v c="$CHROM" '$1 == c {print c ": " $3 " alignments"}'
ls -lh "$OUT" "$OUT.bai"
md5sum "$OUT" > "$OUT.md5"
cat "$OUT.md5"
