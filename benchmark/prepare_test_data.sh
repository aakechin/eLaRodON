#!/bin/bash
# Stage the fixed public chr22 dataset for an image build; never copy a user's reference.
set -euo pipefail
BENCH=$(dirname "$(readlink -f "$0")")
usage() {
    echo 'Usage: bash benchmark/prepare_test_data.sh /path/to/NA12878_rel7_hg19_40x_chr22.bam'
    echo '   or: bash benchmark/prepare_test_data.sh --url BASE_URL'
    echo 'BASE_URL must host the BAM and BAM.bai with their original names (GitHub Release or Zenodo).'
}
[ "$#" -gt 0 ] || { usage >&2; exit 1; }
[ "$1" != --help ] || { usage; exit 0; }
name=NA12878_rel7_hg19_40x_chr22.bam
mkdir -p "$BENCH/test-data"
if [ "$1" = --url ]; then
    [ "$#" -eq 2 ] || { usage >&2; exit 1; }
    for file in "$name" "$name.bai"; do
        curl --fail --location --retry 3 "$2/$file" -o "$BENCH/test-data/$file.download"
        mv "$BENCH/test-data/$file.download" "$BENCH/test-data/$file"
    done
else
    [ "$#" -eq 1 ] || { usage >&2; exit 1; }
    bam=$(readlink -f "$1")
    [ -f "$bam" ] && [ -f "$bam.bai" ] || { echo 'BAM and its .bam.bai index are required' >&2; exit 1; }
    for suffix in '' .bai; do
        target="$BENCH/test-data/$name$suffix"
        [ "$bam$suffix" = "$target" ] || cp "$bam$suffix" "$target"
    done
fi
(cd "$BENCH/test-data" && sha256sum -c SHA256SUMS)
echo 'Test BAM and index verified; Docker will bundle them in /opt/benchmark/data/.'
