# NA12878 chr22 example output

These are the small summary files from the complete 12 November 2026 benchmark run
with all six callers. `metrics.tsv` and the PNGs come from the Figure 5B matching
script; the two `truvari_*.tsv` files use Truvari 5.4.0 and the corresponding
golden-set or full-SvABA baseline. Do not compare their precision or F1 values as
though the two methods used identical matching rules.

`run_manifest.json` records input checksums, the source commit
`b1918e3763eff74223ec31eae43e77896d48a03e`, parameters, and the common
all-FILTER policy. The raw BAM, reference FASTA, caller VCFs and full logs are not
included here. The BAM and BAI will be published as GitHub Release assets; the
reference comes from the user's own hg19 installation.

To reproduce the run, follow [the benchmark guide](../../README.md).
