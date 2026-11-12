# eLaRodON benchmark on NA12878 chromosome 22

This directory reproduces the comparison of eLaRodON with Sniffles2, SVIM, NanoVar, NanoSV and Severus
on the NA12878 cell line (Figure 5B of the article) on a small test dataset: NA12878 reads (ONT, rel7, 40x,
minimap2, hg19) aligned to chromosome 22.

## Files and availability

The `v1.1.0` source archive contains this Dockerfile, all comparison scripts, the NA12878
reference sets and the curated M. tuberculosis CSVs. It does **not** contain the 540 MiB
chr22 BAM or its `.bam.bai` index. Those files are intended as separate assets of the
[`v1.1.0` GitHub Release](https://github.com/aakechin/eLaRodON/releases/tag/v1.1.0).
Check that both files appear under Assets before using the download command below;
the tag by itself does not include them. The hg19 reference FASTA is external and is
mounted read-only at run time. This guide builds the Docker image locally.

Readers who only need the program can use the [root README](../README.md#quick-start)
without downloading the BAM, reference sets or Docker image.

The Dockerfile installs every tool in the version used in the article, each in its own conda environment:

| Tool | Version | Benchmark command |
|---|---|---|
| eLaRodON | GitHub tag `v1.1.0` (image pins source commit b1918e3, including caller fixes 747bd4a) | `python3 main.py -dir OUT -bam BAM -ref REF -vcfanno vcfanno -th N` |
| Sniffles2 | 2.7.1 | `sniffles -i BAM -v OUT.vcf --tandem-repeats human_hg19_chr.trf.bed --threads N` |
| SVIM | 2.0.0 | `svim alignment OUT BAM REF --sample NA12878` |
| NanoVar | 1.8.3 | `nanovar --debug -t N BAM REF OUT` |
| NanoSV | 1.2.4 | `NanoSV BAM -b human_hg19_chr.bed -s samtools -t N -o OUT.vcf` |
| Severus | 1.6 | `severus --target-bam BAM --out-dir OUT -t N` |
| Truvari | 5.4.0 | `truvari bench --pctseq 0 --pctsize 0.7 --refdist 500 --bnddist 100 --sizemin 200 --sizefilt 200 --sizemax 300000000` |

The original environment snapshots are in `envs/*.txt`. Resolved packages for each actual
image build are saved under `/opt/benchmark/environments/*.explicit.txt` inside that image.
Use those build-specific files when recreating an environment.

The test run below targets NA12878 chr22. For the purpose and limitations of the bundled
M. tuberculosis CSV reference annotations, see the final section of this README.

## 1. Build the image

Download the tagged source and, once the Release assets are published, stage the BAM and
BAI. `prepare_test_data.sh` checks both files against `test-data/SHA256SUMS`. Run from
the repository root:

```bash
git clone --depth 1 --branch v1.1.0 https://github.com/aakechin/eLaRodON.git
cd eLaRodON
bash benchmark/prepare_test_data.sh \
    --url https://github.com/aakechin/eLaRodON/releases/download/v1.1.0
docker build -t elarodon-benchmark:check benchmark
```

The download command will work only after both assets appear under the GitHub Release.
If you already have the BAM and its `.bam.bai` index, use
`bash benchmark/prepare_test_data.sh /path/to/NA12878_rel7_hg19_40x_chr22.bam`
instead. The Docker build verifies their checksums again.
The 540 MiB BAM is excluded from ordinary Git history. The image includes it at
`/opt/benchmark/data/NA12878_rel7_hg19_40x_chr22.bam` together with its BAI index.

The Dockerfile currently pins eLaRodON to source commit
`b1918e3763eff74223ec31eae43e77896d48a03e`; override `--build-arg ELARODON_VERSION=COMMIT_OR_TAG`
when testing another published revision. The local benchmark scripts come from this build context.
The SVIM environment uses Matplotlib 3.8 because SVIM 2.0.0 calls `Legend.legendHandles`,
which Matplotlib removed in 3.9; a Python traceback in the SVIM log stops the benchmark.
Resolved conda packages and the eLaRodON commit are recorded inside the image under
`/opt/benchmark/environments/`.

## 2. Run with your existing reference

Use the complete UCSC hg19 FASTA (chromosome names with `chr`) matching the BAM, with its
`.fai` index alongside it. The reference is not bundled or copied into the image. If you
need to download one, UCSC provides
[`hg19.fa.gz`](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz);
uncompress it and run `samtools faidx hg19.fa`. For an exact comparison with the recorded
run, check the FASTA checksum against `truvari/run_manifest.json` in the example results.
From the repository root:

```bash
bash benchmark/run.sh \
    --ref /absolute/path/to/hg19.fa \
    --out ./benchmark/data/results_check
```

The launcher mounts the reference directory read-only and the result directory writable.
It uses the bundled chr22 BAM and runs all six callers, the original ROC/PR evaluation and
both Truvari comparisons. `--image IMAGE` selects a different image; default:
`elarodon-benchmark:check`. At least 16 GB memory is needed (actual availability to Docker
matters). Results belong to your user account. Use a new result directory for a new run.

Equivalent direct Docker command (replace the reference directory and filename):

```bash
mkdir -p benchmark/data/results_check
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp \
    -v /absolute/path/to/reference_directory:/reference:ro \
    -v "$PWD/benchmark/data/results_check":/data \
    elarodon-benchmark:check run_benchmark.sh \
    --ref /reference/hg19.fa --out /data --chrom chr22 --threads 4
```

## 3. Inspect results

Paths below are relative to the directory given by `--out`:

| File | Purpose |
|---|---|
| `evaluation/pr_curves.png` | Precision–recall curves, AP and best F1 (Figure 5B matching rules) |
| `evaluation/roc_curves.png` | ROC curves and AUC (Figure 5B matching rules) |
| `evaluation/metrics.tsv` | Calls, matches, recall, AP, best F1 and its QUAL threshold, ROC AUC |
| `evaluation/pr_points.tsv` | Precision and recall at each QUAL threshold |
| `truvari/golden_set/truvari_summary.tsv` | Truvari against the normalized golden set |
| `truvari/full_svaba/truvari_summary.tsv` | Truvari against the full SvABA baseline |
| `truvari/run_manifest.json` | Input hashes, reference index, parameters, versions and FILTER policy |
| `vcf/TOOL.vcf` | Caller output used for evaluation (original caller FILTER values retained) |
| `logs/TOOL.log`, `logs/TOOL.time` | Diagnostics, runtime and peak memory |

Open PNG files in an image viewer and TSV files in a spreadsheet application with tab as
the delimiter. Raw VCF record counts printed after calling include short events and may
include paired breakends; `evaluation/metrics.tsv` counts calls after its size/type/region
selection. Those counts need not be equal. Truvari applies its own record/event rules;
its counts can also differ from the original Figure 5B evaluation.

## 4. Two Truvari comparisons with one common policy

Both comparisons use exactly the same prepared caller VCFs in `truvari/prepared/`, the same
chr22 region and the same parameters: `--pctseq 0 --pctsize 0.7 --refdist 500 --bnddist 100
--sizemin 200 --sizefilt 200 --sizemax 300000000`. Caller FILTER values are not used to
exclude records, for any caller. eLaRodON normalization also does not select PASS records.
NanoVar runs with `--debug` to retain its complete `*.nanovar.total.vcf`, which is used
instead of `*.nanovar.pass.vcf`. An old NanoVar PASS-only output cannot be silently reused.
Caller-internal discovery thresholds still differ; this policy removes VCF FILTER selection,
not the algorithms' own discovery criteria.

The golden-set source is normalized with the SvABA converter in typed mode and
`--no-collapse`, retaining distinct selected events. On chr22 it produces 27 VCF records;
the full pre-normalized SvABA baseline has 62 records on chr22. These are different
baselines, so report the two result tables separately. The golden set depends on discovery
by at least one compared caller; the full SvABA baseline does not use that selection.
Truvari does not reproduce the Figure 5B tolerance rules even when the golden set is the same.

Normalization retains the existing handling of eLaRodON linked records and INVTD and maps
intrachromosomal TRL to INS, interchromosomal TRL to paired BND records. Original eLaRodON
FILTER annotations remain available in `vcf/eLaRodON.vcf`; its normalized records have
placeholder PASS annotations. Neither Truvari run uses those annotations for selection.

`--steps "evaluate truvari"` reuses existing complete VCFs in the output directory.
The runner fails if any selected caller VCF is missing or any comparison fails. Existing
`truvari/` output is preserved: use a new output directory, copying the desired `vcf/`
inputs there if reusing calls. `--tools "eLaRodON Sniffles2"` selects fewer callers.
`--min-length 50|100|500 --steps "call evaluate"` runs the original evaluation at another
length threshold. The paired Truvari stage requires 200 bp because the bundled full SvABA
baseline is constructed at that threshold.

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

Historical evaluation of the genome-wide VCF files (full 40x BAM) gave the
values of Table S5 (eLaRodON: 2228 of 2363 records found, AP 0.315, best F1 0.517, ROC AUC 0.871).

For Truvari, eLaRodON calls are converted with `scripts/sv_to_std_vcf.py --mode typed --trl ins_intra`: a TRL within a chromosome is written as an insertion, a TRL between chromosomes as two junctions
(BND pairs), as SvABA and the other tools write them.

## Historical genome-wide calls restricted to chromosome 22

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

## Observed run on 12 November 2026 (all caller FILTER values)

The completed chr22 run produced both Truvari tables. Its summary tables, plots and
input manifest are bundled in [`example-results/chr22-v3/`](example-results/chr22-v3/).
The source was commit `b1918e3` from the `v1.1.0` history; the benchmark scripts
were those committed in `47a657a`. The SVIM log ended with `Done.`.

| Tool | Golden TP / 27 | Golden precision | Full SvABA TP / 62 | Full SvABA precision |
|---|---|---|---|---|
| NanoSV | 21 | 0.046 | 22 | 0.048 |
| NanoVar | 22 | 0.103 | 22 | 0.103 |
| SVIM | 20 | 0.010 | 20 | 0.010 |
| Severus | 20 | 0.127 | 20 | 0.127 |
| Sniffles2 | 18 | 0.116 | 18 | 0.116 |
| eLaRodON | 26 | 0.005 | 30 | 0.006 |

These precision values count calls unmatched to the given reference as FP; an
unmatched call is not independently proven false. The Figure 5B-style evaluation
uses different matching rules: eLaRodON finds 24/27, AP 0.1706 and best F1 0.3684
in this run. These chr22 results demonstrate the pipeline; the full-genome values
reported in the article come from a separate analysis. Record the image ID before
publishing a prebuilt image.

## Historical result of the run on the chromosome 22 BAM (PASS-only Truvari)

These values predate the all-FILTER comparison and the latest eLaRodON run. They are
retained as historical observations, not acceptance targets for the current runner.

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

* In the historical run, all NanoSV records carried MapQual (median mapping quality below 80),
  so PASS-only Truvari excluded them all. The current runner ignores FILTER for every caller
  in both Truvari comparisons, as the original golden-set evaluation already did.
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

The `v1.1.0` source tag identifies the code. Release assets are required to distribute
the BAM/BAI. A prebuilt image and a DOI are optional separate publications; add their
verified links here if they become available.

## M. tuberculosis data: purpose and limits

`data/M_tuberculosis/curated_LGRs/` contains the manually curated comparison reference:
38 sample CSVs and 2420 entries. Filenames identify ENA read accessions. The columns are
`POS` (reference coordinate), `SVTYPE` (rearrangement category) and `SVLEN` (length;
compound events can encode two lengths separated by `_`). Open a CSV in a spreadsheet
to inspect a sample's reference events, their types, positions and lengths. These are
reference annotations, not the output of the six callers and not validation of a new sample.

The container includes these CSVs under `/opt/eLaRodON/data/M_tuberculosis/curated_LGRs/`
and the analysis script `/opt/eLaRodON/scripts/mtub_roc_pr.py`. It does not include the
large bacterial BAM/VCF inputs. The CSVs alone cannot reproduce ROC or PR curves; those
also require the caller VCFs in the layout described above. The script produces
`combined_ROC_curves_mtub_200.png`, `metrics_table_mtub.csv`, `pr_metrics_mtub.tsv`,
`pr_points_mtub.tsv` and per-caller `aggregated_df_*.csv`; it does not currently produce
a PR PNG. `run_benchmark.sh` is the NA12878/hg19 demo and must not be used to evaluate
M. tuberculosis against the human reference sets.
