#!/bin/bash
# Runs eLaRodON and five other SV callers on a BAM file and compares their calls with the NA12878
# reference sets, as in the article (Figure 5B) and in the answer to the reviewers (Truvari).
set -euo pipefail

usage() {
    cat <<EOF
Usage: run_benchmark.sh --ref <hg19 FASTA with .fai> [options]

The bundled NA12878 chr22 BAM is used unless --bam is supplied.

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
                     truvari  - Truvari against BOTH the golden set and full SvABA (at least 200 bp); ignores FILTER

Outputs:
  <out>/vcf/<tool>.vcf                  calls of each tool
  <out>/logs/<tool>.log, <tool>.time    tool log, run time and peak memory (/usr/bin/time -v)
  <out>/evaluation/metrics.tsv          calls, found records, recall, AP, best F1 (and its QUAL, precision,
                                        recall), ROC AUC for each tool
  <out>/evaluation/roc_curves.png       ROC curves (AUC in the legend)
  <out>/evaluation/pr_curves.png        precision-recall curves (AP and best F1 in the legend, best F1 marked)
  <out>/evaluation/pr_points.tsv        precision and recall at each QUAL threshold
  <out>/truvari/golden_set/truvari_summary.tsv   Truvari against the golden set
  <out>/truvari/full_svaba/truvari_summary.tsv   Truvari against full SvABA
  <out>/truvari/run_manifest.json              inputs, hashes, FILTER policy and parameters
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
BAM=/opt/benchmark/data/NA12878_rel7_hg19_40x_chr22.bam
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
[ -n "$REF" ] && [ -f "$REF" ] && [ -f "$REF.fai" ] || { echo "reference FASTA with .fai is required (--ref)"; exit 1; }
case "$MIN_LENGTH" in 50|100|200|500) ;; *) echo "--min-length must be 50, 100, 200 or 500" >&2; exit 1 ;; esac
if [[ " $STEPS " == *" truvari "* ]] && [ "$MIN_LENGTH" != 200 ]; then
    echo "Truvari's bundled full SvABA baseline is >=200 bp; use --min-length 200 or --steps 'call evaluate'." >&2
    exit 1
fi
for tool in $TOOLS; do
    case "$tool" in eLaRodON|Sniffles2|SVIM|NanoVar|NanoSV|Severus) ;; *) echo "Unknown tool: $tool" >&2; exit 1 ;; esac
done
for step in $STEPS; do
    case "$step" in call|evaluate|truvari) ;; *) echo "Unknown step: $step" >&2; exit 1 ;; esac
done
[ -n "$TOOLS" ] && [ -n "$STEPS" ] || { echo "Select at least one tool and step" >&2; exit 1; }
REF=$(readlink -f "$REF")
OUT=$(realpath -m "$OUT")
if [[ " $STEPS " == *" truvari "* ]] && [ -e "$OUT/truvari" ]; then
    echo "Output already exists: $OUT/truvari; use a new --out directory to preserve previous results." >&2
    exit 1
fi
mkdir -p "$OUT/vcf" "$OUT/logs"

in_env() { micromamba run -n "$@"; }

run_tool() {
    local tool=$1 d=$OUT/runs/$1 vcf=
    rm -rf "$d" || return 1
    rm -f "$OUT/vcf/$tool.vcf" "$OUT/vcf/$tool.source.txt" || return 1
    mkdir -p "$d" || return 1
    local -a timed=(/usr/bin/time -v -o "$OUT/logs/$tool.time")
    echo "[$(date '+%F %T')] $tool"
    case $tool in
    eLaRodON)
        "${timed[@]}" micromamba run -n elarodon python3 "$ELARODON" -dir "$d" -bam "$BAM" -ref "$REF" \
            -vcfanno "$(in_env elarodon which vcfanno)" -th "$THREADS"
        vcf=$(ls "$d"/*_all_LGRS.vcf 2>/dev/null | head -1) ;;
    Sniffles2)
        # as in the article run on NA12878: tandem repeat BED of Sniffles2 for hg19 (chr names), no --reference
        "${timed[@]}" micromamba run -n sniffles sniffles -i "$BAM" -v "$d/sniffles.vcf" --tandem-repeats /opt/benchmark/human_hg19_chr.trf.bed --threads "$THREADS"
        vcf=$d/sniffles.vcf ;;
    SVIM)
        "${timed[@]}" micromamba run -n svim svim alignment "$d" "$BAM" "$REF" --sample NA12878
        vcf=$d/variants.vcf ;;
    NanoVar)
        "${timed[@]}" micromamba run -n nanovar nanovar --debug -t "$THREADS" "$BAM" "$REF" "$d"
        vcf=$(ls "$d"/*.nanovar.total.vcf 2>/dev/null | head -1) ;;
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
        "${timed[@]}" micromamba run -n nanosv NanoSV "$BAM" -b "$bed" \
            -s "$(in_env nanosv which samtools)" -t "$THREADS" -o "$d/nanosv.vcf"
        vcf=$d/nanosv.vcf ;;
    Severus)
        "${timed[@]}" micromamba run -n severus severus --target-bam "$BAM" --out-dir "$d" -t "$THREADS"
        vcf=$d/all_SVs/severus_all.vcf ;;
    *) echo "unknown tool: $tool"; return 1 ;;
    esac
    local status
    status=$(awk -F': ' '/Exit status/ {print $2}' "$OUT/logs/$tool.time" 2>/dev/null)
    if [ "${status:-1}" != 0 ]; then
        echo "[$(date '+%F %T')] $tool: failed (exit status ${status:-unknown}), see $OUT/logs/$tool.log"
        rm -f "$OUT/vcf/$tool.vcf" "$OUT/vcf/$tool.source.txt"
        return 1
    elif [ "$tool" = SVIM ] && grep -q '^Traceback (most recent call last):' "$OUT/logs/$tool.log"; then
        echo "[$(date '+%F %T')] SVIM: failed with a Python traceback, see $OUT/logs/SVIM.log"
        rm -f "$OUT/vcf/$tool.vcf" "$OUT/vcf/$tool.source.txt"
        return 1
    elif [ -n "$vcf" ] && [ -s "$vcf" ]; then
        cp "$vcf" "$OUT/vcf/$tool.vcf.tmp" || return 1
        mv "$OUT/vcf/$tool.vcf.tmp" "$OUT/vcf/$tool.vcf" || return 1
        printf '%s\n' "$vcf" > "$OUT/vcf/$tool.source.txt" || return 1
        echo "[$(date '+%F %T')] $tool: $(grep -vc '^#' "$vcf") records -> $OUT/vcf/$tool.vcf"
    else
        echo "[$(date '+%F %T')] $tool: no VCF, see $OUT/logs/$tool.log"
        return 1
    fi
}

if [[ " $STEPS " == *" call "* ]]; then
    [ -n "$BAM" ] && [ -f "$BAM" ] || { echo "BAM file is required for the call step (--bam)"; exit 1; }
    [ -f "$BAM.bai" ] || [ -f "${BAM%.bam}.bai" ] || { echo "BAM index missing: $BAM.bai" >&2; exit 1; }
    BAM=$(readlink -f "$BAM")
    REF=$(readlink -f "$REF")
    for tool in $TOOLS; do
        if run_tool "$tool" > "$OUT/logs/$tool.log" 2>&1; then :; else
            tail -20 "$OUT/logs/$tool.log" >&2
            exit 1
        fi
        tail -1 "$OUT/logs/$tool.log"
    done
fi

if [[ " $STEPS " == *" evaluate "* ]] || [[ " $STEPS " == *" truvari "* ]]; then
    for tool in $TOOLS; do
        [ -s "$OUT/vcf/$tool.vcf" ] || { echo "Missing VCF for $tool: $OUT/vcf/$tool.vcf" >&2; exit 1; }
        if [ "$tool" = NanoVar ]; then
            grep -q '\.nanovar\.total\.vcf$' "$OUT/vcf/NanoVar.source.txt" 2>/dev/null || {
                echo "NanoVar requires its complete total.vcf; rerun --tools NanoVar --steps call with this runner before evaluation." >&2
                exit 1
            }
        fi
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
    [ ! -e "$d" ] || { echo "Output already exists: $d; use a new --out directory to preserve previous results." >&2; exit 1; }
    mkdir -p "$d/prepared" "$d/truth"
    awk -v c="$CHROM" '$1 == c {print c "\t0\t" $2}' "$REF.fai" > "$d/region.bed"
    [ -s "$d/region.bed" ] || { echo "$CHROM is not in $REF.fai" >&2; exit 1; }

    # The golden-set source is raw SvABA (including headerless records and chromosome names without chr).
    # Keep all selected events; do not cluster away distinct golden-set entries.
    in_env bench python3 "$BENCH/sv_to_std_vcf.py" svaba -i "$TRUTH_GOLDEN" \
        -o "$d/truth/golden_set.typed.vcf.gz" -r "$REF" --mode typed --no-collapse \
        > "$d/truth/golden_set.prepare.log" 2>&1

    # Prepare each caller ONCE; exactly these same VCFs are used for both baselines.
    for tool in $TOOLS; do
        vcf=$OUT/vcf/$tool.vcf
        if [ "$tool" = eLaRodON ]; then
            in_env bench python3 "$BENCH/sv_to_std_vcf.py" elarodon -i "$vcf" \
                -o "$d/prepared/$tool.vcf.gz" -r "$REF" --mode typed --trl ins_intra \
                > "$d/prepared/$tool.prepare.log" 2>&1
        else
            # NanoSV puts ##cmdline before ##fileformat. Preserve all FILTER values.
            { grep -m1 '^##fileformat' "$vcf"; grep -v '^##fileformat' "$vcf"; } > "$d/prepared/$tool.input.vcf"
            in_env bench bcftools sort -Oz -o "$d/prepared/$tool.vcf.gz" "$d/prepared/$tool.input.vcf" \
                > "$d/prepared/$tool.prepare.log" 2>&1
            in_env bench tabix -f -p vcf "$d/prepared/$tool.vcf.gz"
            rm "$d/prepared/$tool.input.vcf"
        fi
    done
    in_env bench python3 "$BENCH/benchmark_manifest.py" "$d" "$REF" "$CHROM" \
        "$TRUTH_GOLDEN" "$TRUTH_FULL" $TOOLS
    for baseline in golden_set full_svaba; do
        truth=$TRUTH_FULL
        [ "$baseline" != golden_set ] || truth=$d/truth/golden_set.typed.vcf.gz
        mkdir -p "$d/$baseline"
        for tool in $TOOLS; do
            echo "Truvari: $baseline / $tool (all FILTER values)"
            if ! in_env bench truvari bench -b "$truth" -c "$d/prepared/$tool.vcf.gz" \
                -o "$d/$baseline/$tool" --includebed "$d/region.bed" \
                --pctseq 0 --pctsize 0.7 --refdist 500 --bnddist 100 \
                --sizemin 200 --sizefilt 200 --sizemax 300000000 \
                > "$d/$baseline/$tool.truvari.log" 2>&1; then
                echo "Truvari failed: $d/$baseline/$tool.truvari.log" >&2
                exit 1
            fi
        done
        in_env bench python3 "$BENCH/truvari_summary.py" "$d/$baseline" > "$d/$baseline/truvari_summary.tsv"
        echo "Truvari results: $baseline"
        cat "$d/$baseline/truvari_summary.tsv"
    done
fi
