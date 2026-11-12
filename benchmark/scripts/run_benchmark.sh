#!/bin/bash
# Runs eLaRodON and five other SV callers on a BAM file and compares their calls with the NA12878
# reference sets, as in the article (Figure 5B) and in the answer to the reviewers (Truvari).
set -uo pipefail

usage() {
    cat <<EOF
Usage: run_benchmark.sh --bam <sorted indexed BAM> --ref <hg19 FASTA with .fai> [options]

Options:
  --out DIR        output directory (default: results)
  --chrom NAME     chromosome compared with the reference sets (default: chr22)
  --threads N      threads for each tool (default: 4)
  --tools "LIST"   tools to run (default: "eLaRodON Sniffles2 SVIM NanoVar NanoSV Severus")
  --min-length L   minimal SV length for the evaluate step: 50, 100, 200 or 500 (default: 200);
                   the golden set built with the same threshold is used (truth/na12878_golden_set_min<L>.vcf)
  --nanosv-bed BED positions for the NanoSV coverage distribution (default: the hg19 BED bundled with NanoSV)
  --steps "LIST"   steps to run (default: "call evaluate truvari")
                     call     - run the tools, VCF files go to <out>/vcf/<tool>.vcf
                     evaluate - precision-recall and ROC curves against the golden set (Figure 5B rules)
                     truvari  - Truvari against the full SvABA set of at least 200 bp

Outputs:
  <out>/vcf/<tool>.vcf                  calls of each tool
  <out>/logs/<tool>.log, <tool>.time    tool log, run time and peak memory (/usr/bin/time -v)
  <out>/evaluation/metrics.tsv          calls, found records, recall, AP, best F1 (and its QUAL, precision,
                                        recall), ROC AUC for each tool
  <out>/evaluation/roc_curves.png       ROC curves (AUC in the legend)
  <out>/evaluation/pr_curves.png        precision-recall curves (AP and best F1 in the legend, best F1 marked)
  <out>/evaluation/pr_points.tsv        precision and recall at each QUAL threshold
  <out>/truvari/truvari_summary.tsv     Truvari recall, precision and F1 for each tool
EOF
}

BENCH=$(dirname "$(readlink -f "$0")")
TRUTH_GOLDEN=$BENCH/truth/na12878_golden_set.vcf
TRUTH_FULL=$BENCH/truth/svaba_full.typed.vcf.gz
ELARODON=/opt/eLaRodON/main.py
OUT=results
CHROM=chr22
THREADS=4
TOOLS="eLaRodON Sniffles2 SVIM NanoVar NanoSV Severus"
STEPS="call evaluate truvari"
BAM=
REF=
NANOSV_BED=
MIN_LENGTH=200

while [ $# -gt 0 ]; do
    case $1 in
        --bam) BAM=$2; shift 2 ;;
        --ref) REF=$2; shift 2 ;;
        --out) OUT=$2; shift 2 ;;
        --chrom) CHROM=$2; shift 2 ;;
        --threads) THREADS=$2; shift 2 ;;
        --tools) TOOLS=$2; shift 2 ;;
        --steps) STEPS=$2; shift 2 ;;
        --nanosv-bed) NANOSV_BED=$(readlink -f "$2"); shift 2 ;;
        --min-length) MIN_LENGTH=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "unknown option: $1"; usage; exit 1 ;;
    esac
done
[ -n "$REF" ] && [ -f "$REF.fai" ] || { echo "reference FASTA with .fai is required (--ref)"; exit 1; }
OUT=$(readlink -f "$OUT")
mkdir -p "$OUT/vcf" "$OUT/logs"

in_env() { micromamba run -n "$@"; }

run_tool() {
    local tool=$1 d=$OUT/runs/$1
    rm -rf "$d"
    mkdir -p "$d"
    local timed="/usr/bin/time -v -o $OUT/logs/$tool.time"
    echo "[$(date '+%F %T')] $tool"
    case $tool in
    eLaRodON)
        $timed micromamba run -n elarodon python3 "$ELARODON" -dir "$d" -bam "$BAM" -ref "$REF" \
            -vcfanno "$(in_env elarodon which vcfanno)" -th "$THREADS"
        vcf=$(ls "$d"/*_all_LGRS.vcf 2>/dev/null | head -1) ;;
    Sniffles2)
        # as in the article run on NA12878: tandem repeat BED of Sniffles2 for hg19 (chr names), no --reference
        $timed micromamba run -n sniffles sniffles -i "$BAM" -v "$d/sniffles.vcf" --tandem-repeats /opt/benchmark/human_hg19_chr.trf.bed --threads "$THREADS"
        vcf=$d/sniffles.vcf ;;
    SVIM)
        $timed micromamba run -n svim svim alignment "$d" "$BAM" "$REF" --sample NA12878
        vcf=$d/variants.vcf ;;
    NanoVar)
        $timed micromamba run -n nanovar nanovar -t "$THREADS" "$BAM" "$REF" "$d"
        vcf=$(ls "$d"/*.nanovar.pass.vcf 2>/dev/null | head -1) ;;
    NanoSV)
        local bed=$NANOSV_BED
        if [ -z "$bed" ]; then
            # the BED file bundled with NanoSV names chromosomes "1", "2", ...; the BAM uses "chr1", ...
            local pkg
            # "import nanosv" starts its argument parser, so the package is located without importing it
            pkg=$(in_env nanosv python -c "import importlib.util; print(importlib.util.find_spec('nanosv').submodule_search_locations[0])")
            bed=$d/human_hg19_chr.bed
            sed 's/^/chr/' "$pkg/bedfiles/human_hg19.bed" > "$bed"
        fi
        $timed micromamba run -n nanosv NanoSV "$BAM" -b "$bed" \
            -s "$(in_env nanosv which samtools)" -t "$THREADS" -o "$d/nanosv.vcf"
        vcf=$d/nanosv.vcf ;;
    Severus)
        $timed micromamba run -n severus severus --target-bam "$BAM" --out-dir "$d" -t "$THREADS"
        vcf=$d/all_SVs/severus_all.vcf ;;
    *) echo "unknown tool: $tool"; return 1 ;;
    esac
    local status
    status=$(awk -F': ' '/Exit status/ {print $2}' "$OUT/logs/$tool.time" 2>/dev/null)
    if [ "${status:-1}" != 0 ]; then
        echo "[$(date '+%F %T')] $tool: failed (exit status ${status:-unknown}), see $OUT/logs/$tool.log"
        rm -f "$OUT/vcf/$tool.vcf"
    elif [ -n "$vcf" ] && [ -s "$vcf" ]; then
        cp "$vcf" "$OUT/vcf/$tool.vcf"
        echo "[$(date '+%F %T')] $tool: $(grep -vc '^#' "$vcf") records -> $OUT/vcf/$tool.vcf"
    else
        echo "[$(date '+%F %T')] $tool: no VCF, see $OUT/logs/$tool.log"
    fi
}

if [[ " $STEPS " == *" call "* ]]; then
    [ -n "$BAM" ] && [ -f "$BAM" ] || { echo "BAM file is required for the call step (--bam)"; exit 1; }
    BAM=$(readlink -f "$BAM")
    REF=$(readlink -f "$REF")
    for tool in $TOOLS; do
        run_tool "$tool" > "$OUT/logs/$tool.log" 2>&1
        tail -1 "$OUT/logs/$tool.log"
    done
fi

if [[ " $STEPS " == *" evaluate "* ]]; then
    args=()
    for tool in $TOOLS; do
        [ -s "$OUT/vcf/$tool.vcf" ] && args+=(--vcf "$tool=$OUT/vcf/$tool.vcf")
    done
    truth=$TRUTH_GOLDEN
    [ "$MIN_LENGTH" = 200 ] || truth=$BENCH/truth/na12878_golden_set_min$MIN_LENGTH.vcf
    [ -f "$truth" ] || { echo "no golden set for --min-length $MIN_LENGTH: $truth"; exit 1; }
    in_env bench python3 "$BENCH/evaluate.py" --truth "$truth" --chrom "$CHROM" --min-length "$MIN_LENGTH" \
        --out "$OUT/evaluation" "${args[@]}"
fi

if [[ " $STEPS " == *" truvari "* ]]; then
    d=$OUT/truvari
    rm -rf "$d"
    mkdir -p "$d"
    awk -v c="$CHROM" '$1 == c {print c "\t0\t" $2}' "$REF.fai" > "$d/region.bed"
    [ -s "$d/region.bed" ] || { echo "$CHROM is not in $REF.fai"; exit 1; }
    for tool in $TOOLS; do
        vcf=$OUT/vcf/$tool.vcf
        [ -s "$vcf" ] || continue
        if [ "$tool" = eLaRodON ]; then
            # the same normalization as for the comparison in the answer to the reviewers
            in_env bench python3 "$BENCH/sv_to_std_vcf.py" elarodon -i "$vcf" -o "$d/$tool.vcf.gz" \
                -r "$REF" --mode typed --trl ins_intra --pass-only > "$d/$tool.prepare.log" 2>&1
        else
            # NanoSV writes ##cmdline before ##fileformat, which bcftools does not accept
            { grep -m1 '^##fileformat' "$vcf"; grep -v '^##fileformat' "$vcf"; } > "$d/$tool.input.vcf"
            in_env bench bcftools sort -Oz -o "$d/$tool.vcf.gz" "$d/$tool.input.vcf" > "$d/$tool.prepare.log" 2>&1 \
                && in_env bench tabix -f -p vcf "$d/$tool.vcf.gz"
            rm -f "$d/$tool.input.vcf"
        fi
        if [ ! -s "$d/$tool.vcf.gz" ]; then
            echo "$tool: VCF could not be prepared for Truvari, see $d/$tool.prepare.log"
            continue
        fi
        in_env bench truvari bench -b "$TRUTH_FULL" -c "$d/$tool.vcf.gz" -o "$d/$tool" \
            --includebed "$d/region.bed" --passonly --pctseq 0 --sizemin 200 --sizefilt 200 \
            --sizemax 300000000 > "$d/$tool.truvari.log" 2>&1 \
            || echo "$tool: Truvari failed, see $d/$tool.truvari.log"
    done
    in_env bench python3 "$BENCH/truvari_summary.py" "$d" > "$d/truvari_summary.tsv"
    cat "$d/truvari_summary.tsv"
fi
