#!/bin/bash
# Host-side launcher: users supply only a reference and an output directory.
set -euo pipefail
usage() {
    echo 'Usage: bash benchmark/run.sh --ref /path/to/hg19.fa --out /path/to/results [--image IMAGE] [runner options]'
    echo 'Runner options: --tools "LIST", --steps "LIST", --threads N, --chrom chr22, --min-length N'
}
REF=
OUT=
IMAGE=elarodon-benchmark:check
args=()
while [ "$#" -gt 0 ]; do
    case "$1" in
        --ref|--out|--image|--tools|--steps|--threads|--chrom|--min-length)
            [ "$#" -ge 2 ] || { echo "Missing value for $1" >&2; exit 1; }
            case "$1" in
                --ref) REF=$2 ;;
                --out) OUT=$2 ;;
                --image) IMAGE=$2 ;;
                *) args+=("$1" "$2") ;;
            esac
            shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done
[ -f "$REF" ] || { echo "Reference FASTA not found: $REF" >&2; exit 1; }
REF=$(readlink -f "$REF")
[ -f "$REF.fai" ] || { echo "Reference index not found: $REF.fai; create it with samtools faidx." >&2; exit 1; }
[ -n "$OUT" ] || { echo '--out is required' >&2; exit 1; }
mkdir -p "$OUT"
OUT=$(readlink -f "$OUT")
# --mount rejects missing directories rather than silently creating empty input mounts.
exec docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp \
    --mount "type=bind,source=$(dirname "$REF"),target=/reference,readonly" \
    --mount "type=bind,source=$OUT,target=/data" \
    "$IMAGE" run_benchmark.sh --ref "/reference/$(basename "$REF")" --out /data "${args[@]}"
