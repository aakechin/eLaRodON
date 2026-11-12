# eLaRodON benchmark on NA12878 chromosome 22

This directory reproduces the comparison of eLaRodON with Sniffles2, SVIM, NanoVar, NanoSV and Severus
on the NA12878 cell line (Figure 5B of the article) on a small test dataset: NA12878 reads (ONT, rel7, 40x,
minimap2, hg19) aligned to chromosome 22.

The Dockerfile installs every tool in the version used in the article, each in its own conda environment:

| Tool | Version | Command (as in the article) |
|---|---|---|
| eLaRodON | ec95ede (caller fixes; release tag pending) | `python3 main.py -dir OUT -bam BAM -ref REF -vcfanno vcfanno -th N` |
| Sniffles2 | 2.7.1 | `sniffles -i BAM -v OUT.vcf --tandem-repeats human_hg19_chr.trf.bed --threads N` |
| SVIM | 2.0.0 | `svim alignment OUT BAM REF --sample NA12878` |
| NanoVar | 1.8.3 | `nanovar -t N BAM REF OUT` |
| NanoSV | 1.2.4 | `NanoSV BAM -b human_hg19_chr.bed -s samtools -t N -o OUT.vcf` |
| Severus | 1.6 | `severus --target-bam BAM --out-dir OUT -t N` |
| Truvari | 5.4.0 | `truvari bench --passonly --pctseq 0 --sizemin 200 --sizefilt 200 --sizemax 300000000` |

Exact package versions of each environment are listed in `envs/*.txt` (conda explicit lists; an environment
can be recreated with `conda create -n NAME --file envs/NAME.txt`).

On one M. tuberculosis genome from the article (ERR12245887), the image gives the same number of VCF records as
the runs of the article for SVIM (1219), NanoVar (164), NanoSV (272) and Severus (85). The six tools take about
2 minutes on 4 threads.

## 1. Get the image

```bash
# Build from the published pre-release commit; run from this benchmark directory.
docker build --build-arg ELARODON_VERSION=fad19a6 -t elarodon-benchmark:revision .
```

## 2. Get the data

```bash
mkdir -p data && cd data
# Place NA12878_rel7_hg19_40x_chr22.bam and its .bam.bai index here.
# Download URLs will be added when the dataset is deposited.
# Alternatively, extract the test BAM from the full NA12878 BAM using make_test_bam.sh.
# hg19 reference (UCSC, chr names), about 3 GB unpacked
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
docker run --rm -v "$PWD":/data elarodon-benchmark:revision micromamba run -n elarodon samtools faidx /data/hg19.fa
cd ..
```

The whole reference is needed even for one chromosome: split reads of chromosome 22 have supplementary
alignments on other chromosomes, and the tools read the reference around them.

## 3. Run

The run needs at least 16 GB of memory: eLaRodON takes up to 9.3 GB and NanoVar up to 11 GB on the test BAM.

```bash
docker run --rm -v "$PWD/data":/data elarodon-benchmark:revision \
    run_benchmark.sh --bam /data/NA12878_rel7_hg19_40x_chr22.bam --ref /data/hg19.fa \
    --chrom chr22 --threads 4 --out /data/results
```

Steps (`--steps`, all by default):

* `call` runs the six tools. VCF files: `results/vcf/<tool>.vcf`; logs, run time and peak memory:
  `results/logs/<tool>.log`, `results/logs/<tool>.time`.
* `evaluate` compares the calls with the golden set by the rules of Figure 5B and draws the curves:
  `results/evaluation/metrics.tsv` (calls, found records, recall, AP, best F1 with its QUAL, precision and
  recall, ROC AUC), `results/evaluation/roc_curves.png` (AUC in the legend), `results/evaluation/pr_curves.png`
  (AP and best F1 in the legend, the best F1 point marked), `results/evaluation/pr_points.tsv`.
  `--min-length 50|100|500` repeats the evaluation with another minimal SV length and the golden set built with
  the same threshold.
* `truvari` compares the calls with the full SvABA set of at least 200 bp in Truvari:
  `results/truvari/truvari_summary.tsv`.

`--tools "eLaRodON Sniffles2"` runs only some of the tools. `--steps "evaluate truvari"` evaluates
VCF files that are already in `results/vcf/` (for example, calls of another version of a tool).

## Reference sets (`scripts/truth/`)

* `na12878_golden_set.vcf` — the reference set: 2363 SvABA SVs (Wala et al., 2018) of at least 200 bp and
  junctions between chromosomes (no length) that were detected by at least one of the six tools, one record per
  SV; 27 of them are on chromosome 22. `na12878_golden_set_min50.vcf`, `_min100.vcf`, `_min500.vcf`: the same
  with minimal lengths of 50, 100 and 500 bp (5518, 3101 and 1127 records). `na12878_golden_set_2336.vcf` and
  `na12878_golden_set_article_2181.vcf`: earlier versions, before the junctions between chromosomes were added
  and before the fixes described in the answer to Reviewer 2, comment Q1.
* `svaba_full.typed.vcf.gz` — the full SvABA call set (SV and indel files) of at least 200 bp converted to the
  standard form with `scripts/sv_to_std_vcf.py`: 4049 SVs, 4125 VCF records; it does not depend on the compared
  tools.

## Matching rules of Figure 5B (`scripts/evaluate.py`)

* calls shorter than 200 bp and eLaRodON INVTD and BND_INVTD are not compared; FILTER is not used;
  junctions between chromosomes have no length and are always compared;
* a call matches a golden-set record of a compatible type if both breakpoints (for insertions, the position
  and the length) differ by at most 10% of the call length (100 bp for calls without a length and for junctions
  between chromosomes);
* DUP is treated as TD; BND records are typed by ALT, a BND whose mate is on another chromosome is BND_TRL
  (for Sniffles2 and NanoVar the mate is taken from ALT); an eLaRodON TRL between chromosomes (a fragment of
  another chromosome inserted at POS) has two junctions, and the call matches if either of them does;
* if several calls match one record, only the first one is counted;
* precision and recall are calculated at every QUAL threshold; AP is the sum of precision values weighted by
  the increase in recall; in the ROC curve, golden-set records missed by a tool get QUAL = -1.

With the genome-wide VCF files (the latest versions of the tools, full 40x BAM), `evaluate.py` gives the
values of Table S5 (eLaRodON: 2228 of 2363 records found, AP 0.315, best F1 0.517, ROC AUC 0.871).

For Truvari, eLaRodON calls are converted with `scripts/sv_to_std_vcf.py --mode typed --trl ins_intra
--pass-only`: a TRL within a chromosome is written as an insertion, a TRL between chromosomes as two junctions
(BND pairs), as SvABA and the other tools write them.

## Expected result on chromosome 22

Calls of the genome-wide runs (full 40x BAM) with a breakpoint on chromosome 22, evaluated with
`evaluate.py --chrom chr22` (`results_chr22_30oct/expected_from_genome_vcfs/`). Reference set: 27 records on
chromosome 22.

| Tool | Calls compared | Found | AP | Best F1 | ROC AUC |
|---|---|---|---|---|---|
| eLaRodON | 11263 | 25 | 0.165 | 0.359 | 0.829 |
| NanoVar | 176 | 18 | 0.110 | 0.200 | 0.352 |
| SVIM | 2558 | 21 | 0.105 | 0.206 | 0.750 |
| Severus | 161 | 19 | 0.084 | 0.203 | 0.354 |
| Sniffles2 | 166 | 18 | 0.073 | 0.200 | 0.338 |
| NanoSV | 323 | 21 | 0.065 | 0.181 | 0.485 |

## Result of the run on the chromosome 22 BAM

Run of 29.10.2026 on a laptop with 6 CPUs and 14 GB of memory, `--threads 4` (`call` step:
`results_chr22_29oct/`); `evaluate` and `truvari` repeated on the same VCF files on 30.10.2026 with the current
scripts (`results_chr22_30oct/`). Full SvABA set in Truvari: 62 records on chromosome 22.

| Tool | Calls compared | Found | AP | Best F1 | ROC AUC | Truvari recall | Truvari precision | Time, m:s | Peak memory, MB |
|---|---|---|---|---|---|---|---|---|---|
| eLaRodON | 5359 | 23 | 0.143 | 0.326 | 0.759 | 0.452 | 0.009 | 2:11 | 9115 |
| SVIM | 2198 | 21 | 0.133 | 0.217 | 0.747 | 0.323 | 0.010 | 0:41 | 1273 |
| Severus | 156 | 19 | 0.086 | 0.208 | 0.352 | 0.323 | 0.127 | 0:28 | 555 |
| Sniffles2 | 155 | 18 | 0.079 | 0.212 | 0.338 | 0.290 | 0.122 | 0:14 | 156 |
| NanoSV | 289 | 21 | 0.069 | 0.192 | 0.472 | – | – | 2:48 | 392 |
| NanoVar | – | – | – | – | – | – | – | killed | 10 957 |

The five tools and the evaluation take about 10 minutes. NanoVar does not finish on a 14 GB computer: at the step
"Correcting DUP and detecting TE" it aligns reads to the whole reference with minimap2, which uses 11 GB. A
computer with at least 16 GB of memory is needed for the full run.

The numbers of calls on the chromosome 22 BAM differ from the genome-wide runs (Sniffles2 155 vs. 166, Severus
156 vs. 161, eLaRodON 5359 vs. 11263), mostly by junctions between chromosomes: reads with the primary alignment
on another chromosome are present in this BAM only through their supplementary alignments. With 27 reference records,
AP depends on a few calls, so the test checks reproducibility, not the ranking of the tools.

## Notes

* NanoSV marks all records with the MapQual filter (median mapping quality below 80; minimap2 reports at most
  60), so Truvari with `--passonly` keeps none of them. The golden-set evaluation does not use FILTER.
* `run_benchmark.sh` finds the BED file of NanoSV with `importlib.util.find_spec`, because `import nanosv`
  starts its argument parser.
* The test BAM contains only the alignments located on chromosome 22. Reads with the primary alignment on
  another chromosome are present only through their supplementary alignments on chromosome 22.

## Additional analysis scripts

The repository-level `scripts/` directory contains the reference-set construction and table scripts.
These use explicit input paths and do not require the authors' directory layout.

```bash
python3 scripts/make_chr1_table.py chr1.xlsx --vcf eLaRodON=/path/to/elarodon.vcf
python3 scripts/make_na12878_reference_set.py 200 reference_results --inputs inputs.json
python3 scripts/mtub_roc_pr.py mtub_results --vcf-root=/path/to/mtub/vcf
```

Run these commands from the repository root. `make_chr1_table.py` accepts `--truth FILE`,
`--chrom NAME`, and repeated `--vcf TOOL=FILE` arguments. Its default reference is the bundled
2363-event set. The reference-set builder requires directories for the keys `svaba_path`,
`elarodon_path`, `svim_path`, `sniffles_path`, `nanovar_path`, `nanosv_path`, and `severus_path`
in a JSON object. Relative paths are resolved against the JSON file location. Keep each caller's
VCFs in its own directory and keep the original SvABA filenames (`*.svaba.sv.vcf` and
`*.svaba.indel.vcf`). Input VCFs for these additional scripts must be uncompressed.

For M. tuberculosis, the VCF layout is `VCF_ROOT/CALLER/SAMPLE/*.vcf`, with caller directories
`elarodon`, `nanosv`, `nanovar`, `sniffles`, `svim`, and `severus`; use `--csv-root=DIR`
to override the bundled curated CSVs. `--mosaic-root=DIR` supplies optional mosaic runs,
and `--lowsupport --lowsupport-root=DIR` supplies the additional low-support runs.
The VCFs used for the full study are external inputs; the repository does not bundle those calls.

The pre-release build is not a published v1.1.0 release. Final image identifiers, dataset download
links and DOI will be recorded after deposition and container validation.
