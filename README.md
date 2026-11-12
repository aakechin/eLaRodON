# eLaRodON

eLaRodON detects genomic rearrangements from Oxford Nanopore read alignments.
It combines CIGAR-derived insertions and deletions with split-read junctions,
classifies rearrangements and reports sequence and genomic context in a VCF.

## Choose what to download

| Purpose | Files | Instructions |
|---|---|---|
| Analyze your own BAM with the current program | Install [eLaRodON 1.1.1 from PyPI](https://pypi.org/project/elarodon/1.1.1/) | [Quick start](#quick-start) |
| Obtain the current program source | [v1.1.1 ZIP](https://github.com/aakechin/eLaRodON/archive/refs/tags/v1.1.1.zip) or clone the repository | [Run from source](#run-from-source) |
| Reproduce the recorded NA12878 chr22 benchmark | [v1.1.0 source](https://github.com/aakechin/eLaRodON/archive/refs/tags/v1.1.0.zip), BAM, BAI and checksums from [Release v1.1.0](https://github.com/aakechin/eLaRodON/releases/tag/v1.1.0), plus a matching hg19 FASTA | [Benchmark guide](benchmark/README.md) |
| Inspect curated reference annotations | Small reference files already included in `benchmark/scripts/truth/` and `data/M_tuberculosis/curated_LGRs/` | [Reference data and limits](benchmark/README.md#reference-sets-scriptstruth) |

The pip package contains the program. Benchmark scripts and curated data are
in the GitHub source archive. The 540 MiB BAM and its index are separate
Release assets; GitHub's automatic source ZIP/TAR files do not contain them.
The reference FASTA is supplied separately. The benchmark Docker image is
built locally from `benchmark/Dockerfile`.

## Quick start

Use Linux with Python 3.10, 3.11 or 3.12. Create a virtual environment on a
filesystem that supports symbolic links (for example, your home directory):

```bash
python3 -m venv ~/envs/elarodon-1.1.1
source ~/envs/elarodon-1.1.1/bin/activate
python -m pip install elarodon==1.1.1
elarodon --help
```

Install the external tools listed below, then provide a coordinate-sorted,
indexed BAM and its matching indexed reference FASTA:

```bash
elarodon \
    -dir /absolute/path/to/results \
    -bam /absolute/path/to/sample.bam \
    -ref /absolute/path/to/reference.fa \
    -vcfanno /absolute/path/to/vcfanno \
    -th 4
```

The final output is `sample_all_LGRS.vcf` in the result directory. Use a new
result directory for each run. Use absolute paths without spaces: annotation
commands pass paths to external tools through the shell.

`python -m elarodon` accepts the same arguments as `elarodon`.

## Dependencies and inputs

Python dependencies are installed by pip and pinned in
[requirements.txt](requirements.txt). Install these programs separately and
make them available on `PATH`:

- [minimap2](https://github.com/lh3/minimap2)
- [samtools](https://www.htslib.org/) and HTSlib's `bgzip` and `tabix`
- [bedtools](https://github.com/arq5x/bedtools2)
- [vcfanno](https://github.com/brentp/vcfanno/releases)

The BAM must be coordinate-sorted and indexed. The FASTA must match the
assembly and chromosome names used for alignment and have a `.fai` index:

```bash
samtools sort -o sample.sorted.bam sample.bam
samtools index sample.sorted.bam
samtools faidx reference.fa
```

An optional annotation BED contains tab-separated chromosome, start and end
coordinates, with an optional fourth column naming the element. Pass it with
`-bed /absolute/path/to/repeats.bed`. Omitting the BED uses a temporary placeholder
annotation file; it does not provide biological repeat annotations.

## Run from source

To obtain the source corresponding to the installed PyPI program:

```bash
git clone --depth 1 --branch v1.1.1 https://github.com/aakechin/eLaRodON.git
cd eLaRodON
python3 -m venv ~/envs/elarodon-source-1.1.1
source ~/envs/elarodon-source-1.1.1/bin/activate
python -m pip install -r requirements.txt
python main.py --help
```

Replace `elarodon` with `python main.py` in the run command above. To install
this source as a package instead, use `python -m pip install .`.

## Arguments

The spelling below follows the actual command-line parser. Run `elarodon --help`
for the full list, including explicit intermediate filenames.

| Argument | Default | Purpose |
|---|---|---|
| `-dir`, `--workdir` | required | Result directory |
| `-bam`, `--bam-file` | required | Sorted, indexed BAM |
| `-ref`, `--ref-genome` | required | Matching indexed reference FASTA |
| `-vcfanno`, `--vcf_anno` | required | Path to the vcfanno executable |
| `-bed`, `--bed_file` | omitted | Genomic annotation BED |
| `-th`, `--threads` | 4 | Number of threads/processes |
| `-len`, `--minimal-length` | 50 | Minimum insertion/deletion length in bp |
| `-clip`, `--minimal-clipped-length` | 100 | Minimum clipped-part length in bp |
| `-dist`, `--dist-to-join-trl` | 1000 | Maximum distance used to join translocation coordinates |
| `-join`, `--maximal-distance-join` | 30 | Maximum breakpoint distance for junction merging |
| `-dvlen`, `--div-length` | 0 | Region length for chromosome subdivision; 0 uses a whole chromosome |
| `-cont`, `--continue` | all | Starting stage: `all`, `bam`, `join` or `def` |
| `-out`, `--out_vcf` | automatic | Output VCF filename |
| `-nrt_ins`, `--not_remove_trash_align` | false | Keep temporary insertion-alignment files |
| `-nrt_anno`, `--not_remove_trash_anno` | false | Keep temporary annotation files |

`--continue join` requires the junction CSVs from the BAM stage and runs
joining, insertion alignment and annotation. `--continue def` requires the
mapped-insertion and junction tables and runs classification and annotation.
All four required arguments remain required when continuing a run. For example,
to repeat annotation after an earlier completed run:

```bash
elarodon -cont def \
    -dir /absolute/path/to/previous_results \
    -bam /absolute/path/to/sample.bam \
    -ref /absolute/path/to/reference.fa \
    -vcfanno /absolute/path/to/vcfanno \
    -bed /absolute/path/to/new_annotation.bed
```

## Analysis and output

The program reads primary alignments and their supplementary-alignment tags,
collects split-read junctions and CIGAR-derived insertions/deletions, and groups
nearby calls. It aligns insertion sequences to the reference with minimap2,
then classifies rearrangements and annotates their boundaries with vcfanno.

The result directory contains:

| File | Content |
|---|---|
| `*_all_LGRS.vcf` | Final rearrangements and annotations |
| `elarodon.log` | Processing log and paths |
| `supplementary/*.fusions.csv`, `*.insertions.csv` | Initial junction and insertion tables |
| `supplementary/*LRs_join100*.csv`, `*INS_join100*.csv` | Joined calls and insertion-mapping intermediates; `join100` is part of the default filename |

The VCF header defines the INFO and FORMAT fields, including rearrangement type,
coordinates, supporting-read counts, mapping quality and sequence/context annotations.
Some compound rearrangements use eLaRodON-specific representations. The
[benchmark guide](benchmark/README.md) describes normalization before comparison
with other callers. A single supporting read can produce a call; interpret
support, QUAL and annotations together when selecting results for a study.

## Versions and benchmark reproducibility

Version **1.1.1** is available on PyPI and includes the fix that retains the
last junction during merging, considers it for merging with a neighboring call,
and handles progress reporting for an empty junction table. Installed-wheel
and source runs produced the same DEL/INS records in a small synthetic test.

The recorded chr22 benchmark uses **v1.1.0**, including its earlier caller fixes.
The additional 1.1.1 junction fix can change calls; the recorded benchmark has
not been recomputed with 1.1.1. Use the v1.1.0 source and the pinned Docker build
to reproduce those tables. See [benchmark/README.md](benchmark/README.md) for
asset download commands, tool versions, comparison rules, example outputs and
the limitations of the chromosome-only dataset.

## Troubleshooting

- Missing executable: install the external tool and check its location with
  `command -v minimap2`, `command -v bedtools` or `command -v vcfanno`.
- Missing BAM/FASTA index: run `samtools index` or `samtools faidx` as above.
- Reference mismatch: use the same assembly and chromosome naming as the BAM.
- Permission or path errors: choose a writable result directory and paths without spaces.
- Insufficient memory: reduce parallelism or use a machine with more available memory.

## Citation and license

Please cite the eLaRodON article, *eLaRodON: identification of large genomic
rearrangements in Oxford Nanopore sequencing data*, and record the software
version used in your analysis.

[MIT License](LICENSE) covers the original software and documentation.
Third-party datasets retain their source terms.
