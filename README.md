# eLaRodON

![Python 3.10–3.12](https://img.shields.io/badge/Python-3.10--3.12-blue)
[![PyPI](https://img.shields.io/pypi/v/elarodon)](https://pypi.org/project/elarodon/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green)](LICENSE)

eLaRodON detects genomic rearrangements from Oxford Nanopore read alignments.
It combines CIGAR-derived insertions and deletions with split-read junctions,
classifies rearrangements and reports sequence and genomic context in a VCF.

## Contents

- [Choose what to download](#choose-what-to-download)
- [Quick start](#quick-start)
- [Dependencies and inputs](#dependencies-and-inputs)
- [Run from source](#run-from-source)
- [Arguments](#arguments)
- [Algorithm and workflow](#algorithm-and-workflow)
- [Analysis and output](#analysis-and-output)
- [Reading the VCF](#reading-the-vcf)
- [Quality and the new sequence pattern](#quality-and-the-new-sequence-pattern)
- [Versions and benchmark reproducibility](#versions-and-benchmark-reproducibility)
- [Troubleshooting](#troubleshooting)
- [Citation and license](#citation-and-license)

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

You can also name intermediate files explicitly:

| Argument | Purpose |
|---|---|
| `-in`, `--input-files` | Pattern selecting the junction CSV files to join |
| `-lrs`, `--output-lrs` | Joined rearrangement CSV filename |
| `-ins`, `--output-ins` | Joined insertion CSV filename |
| `-n`, `--name` | Name of the insertion-alignment run |
| `-sam`, `--sam_file` | Existing insertion-alignment SAM filename |
| `-juncf`, `--junc_file` | Rearrangement table after insertion alignment |
| `-insf`, `--ins_file` | Insertion table after insertion alignment |

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

## Algorithm and workflow

[![](https://mermaid.ink/img/pako:eNpVkFFrgzAUhf9KuM-26NQa8zBodYXC6KAbe5gWCXptHZpIjNu62v--aDfK7lNOzndykpwhlwUCg4Pi7ZG8xKkgZpbJ0_blcZfzuka1J7PZ_bDjn2TU3UBWybusxER0-2tgNTHrqtaosCDGGEiU8Lo6iAaFzjbb53_kRnSodCXFyF2NaDKWQkjNNRYDiZMCy0pgpk8tZrlCs5195GV2q41_awWviex12-uBPCSv0ZqUVY17sMy7qgKYVj1a0KBq-CjhPMZT0EdsMAVmlqaK97VOIRUXE2u5eJOy-Usq2R-OwEped0b1bWGuElfcfNoNQVGgimQvNDAnmI4AdoYvYG7gzD0npAubuiH13MC34ATMC-bUufMdx_M8O3ToxYLvqdOe08C3zQSh6_shpYvLD7JJhTs?type=png)](https://mermaid.live/edit#pako:eNpVkFFrgzAUhf9KuM-26NQa8zBodYXC6KAbe5gWCXptHZpIjNu62v--aDfK7lNOzndykpwhlwUCg4Pi7ZG8xKkgZpbJ0_blcZfzuka1J7PZ_bDjn2TU3UBWybusxER0-2tgNTHrqtaosCDGGEiU8Lo6iAaFzjbb53_kRnSodCXFyF2NaDKWQkjNNRYDiZMCy0pgpk8tZrlCs5195GV2q41_awWviex12-uBPCSv0ZqUVY17sMy7qgKYVj1a0KBq-CjhPMZT0EdsMAVmlqaK97VOIRUXE2u5eJOy-Usq2R-OwEped0b1bWGuElfcfNoNQVGgimQvNDAnmI4AdoYvYG7gzD0npAubuiH13MC34ATMC-bUufMdx_M8O3ToxYLvqdOe08C3zQSh6_shpYvLD7JJhTs)

The workflow diagram and its editable Mermaid source are included in the
repository. The four modules are run by `main.py`:

| Module | Input | Work performed | Output |
|---|---|---|---|
| [ONTLRcaller.py](ONTLRcaller.py) | Sorted BAM | Reads primary alignments, CIGAR operations and supplementary-alignment tags; records candidate junctions and insertions | Per-region junction and insertion CSVs |
| [joinONTLRs.py](joinONTLRs.py) | Detection CSVs | Groups compatible junctions and insertions and aggregates read evidence | Joined rearrangement and insertion tables |
| [alignment_INS.py](alignment_INS.py) | Joined calls and reference FASTA | Maps insertion sequences with minimap2 and recognizes nearby reference-derived sequence, including tandem duplications | Insertion alignments and updated calls |
| [define_type_create_vcf_LRs.py](define_type_create_vcf_LRs.py) | Updated calls, reference and optional annotation BED | Assigns types, links related records, calculates junction features and annotates genomic context with vcfanno | Annotated VCF |

### Detection and read support

Split-read junctions describe how adjacent fragments of the same read align to
the reference: their chromosomes, positions, strands and the sides adjoining
the junction. CIGAR insertions and deletions are collected separately. The
minimum insertion/deletion length defaults to 50 bp; the article's 200-bp
benchmark cutoff is an analysis choice, rather than this detection default.

The split-read analysis uses supplementary alignments recorded in the primary
alignment's `SA` tag. A read must have at least 100 clipped bases by default;
short supplementary fragments are excluded. Compatible junctions from different
reads are grouped using their orientations and coordinates, with a default
breakpoint distance of 30 bp. CIGAR insertions use their own position and
length criteria. Each read is counted once per grouped record.

### One-boundary and two-boundary rearrangements

A read that spans both junctions can resolve a rearranged fragment and its two
boundaries. A read covering only one junction gives a boundary representation,
often named `BND_*`. This distinction concerns the alignment pattern: a BND call
can have several supporting reads, and a two-boundary call can have only one.

Related one-boundary and two-boundary records can remain separate and be linked
by `SBLR`, `SBTRL`, `SBINV` or `SBTD`. Those links indicate possible alternative
representations or a second boundary; they are not additional independent
experimental confirmations.

### Junction and genomic context

For each rearrangement, eLaRodON records overlapping or unaligned sequence
between read fragments, microhomology and homeology at the breakpoints, and
read-support features. With an annotation BED, it also reports repeats or other
named genomic elements containing or flanking each boundary. The meaning of an
annotation depends on the BED supplied by the user.

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

## Reading the VCF

The VCF combines the standard columns with eLaRodON-specific tags. Its header
lists the fields emitted by the installed version. The tables below explain
the main groups; tags can be absent or contain a missing value when the
corresponding evidence is unavailable.

### Coordinates and rearrangement types

| Field | Meaning |
|---|---|
| `SVTYPE` | Native eLaRodON rearrangement type |
| `CHROM2`, `END` | Chromosome and coordinate of the other boundary |
| `SVLEN` | Insertion/deletion length or an estimated rearranged-fragment length |
| `TDRN` | Tandem repeat count for duplication calls |
| `J1`, `J2` | Side of each read fragment adjoining the junction: L, R or ND |
| `S1`, `S2` | Strand orientation of the fragments |
| `D1`, `D2` | Distances used to characterize the junction pattern |
| `RL1`, `RL2` | Read-fragment lengths |

| Type | Interpretation |
|---|---|
| `DEL`, `BND_DEL` | Deletion or a deletion-like boundary |
| `INS` | Insertion relative to the reference |
| `TD`, `BND_TD` | Tandem duplication or a tandem-duplication-like boundary |
| `INV`, `BND_INV` | Inversion or an inversion-like boundary |
| `INVTD`, `BND_INVTD` | Inverted tandem duplication or its boundary representation |
| `TRL`, `BND_TRL` | Translocation or a junction connecting the reported genomic regions |

Native records can describe different evidence patterns for the same event.
For comparisons with other callers, use an explicit normalization and matching
procedure, as in the [benchmark scripts](benchmark/scripts/).

### Read support and sequence features

| Field | Meaning |
|---|---|
| `SR` | Number of supporting reads |
| `MQ` | Mapping-quality measure of the supporting fragments |
| `DP`, `VAF` | Read depth and estimated variant allele frequency |
| `MH`, `MHS` | Microhomology length and sequence |
| `HOM`, `HOMS` | Homeology length and sequence |
| `MUTM`, `MUTV` | Median and variance of read-fragment mapping-error measures |
| `ISN`, `ISNMFS` | Count and stored length measure for overlapping read fragments |
| `NSN`, `NSNMFS` | Count and stored length measure for unaligned sequence between fragments |
| `SLN` | Number of reads whose fragments meet without an overlap or intervening sequence |
| `CN`, `IVN` | CIGAR and inversion evidence counts |
| `NSP` | Homopolymer-tract count component of the new sequence pattern score |
| `MNF`, `PFF`, `PLF` | Rearrangement-fragment number and first/last-fragment summary measures |
| `SBLR`, `SBTRL`, `SBINV`, `SBTD` | Links to a possible second boundary or related rearrangement representation |

The sample `FORMAT` fields include `GT`, `DP`, `AD` and `VF`. Read these together
with the INFO annotations and the actual supporting alignments when interpreting
a candidate. eLaRodON does not establish whether an event is somatic or germline.

### Annotations from the BED file

| Field | Meaning |
|---|---|
| `LERN`, `RERN` | Names of the nearest external elements at the left/right boundary |
| `LIRN`, `RIRN` | Names of the nearest internal elements |
| `LCRN`, `RCRN` | Names of elements containing the boundary |
| `LERD`, `RERD` | Distances to the nearest external elements |
| `LIRD`, `RIRD` | Distances to the nearest internal elements |

For example, a RepeatMasker BED makes these repeat annotations; another BED can
supply another class of genomic features. Without a biological annotation BED,
these fields do not supply a repeat classification.

## Quality and the new sequence pattern

`QUAL` combines read support, mapping and junction features. `FILTER=PASS` means
that eLaRodON assigned a rearrangement type; `FAIL` means no type was assigned.
A PASS label alone does not establish experimental validity or a somatic origin.

The new sequence pattern concerns unaligned sequence between adjacent read
fragments. For a sequence `s`, count runs of the same nucleotide that are at
least four bases long, then subtract one per 100 bases of sequence:

```text
score(s) = number of homopolymer tracts of length >= 4 - length(s) / 100
```

For example, `ACGTTTTTGCAAAAGT` is 16 bases long and contains two such tracts
(`TTTTT` and `AAAA`), giving `2 - 16/100 = 1.84`. For aggregated calls, the stored
score is calculated as `NSP - NSNMFS/100`; `NSP` alone is the count component,
not the final length-adjusted score. If this score is above zero, QUAL is divided
by two. The feature is an annotation of a pattern associated with artifacts,
not a test that proves a particular rearrangement false.

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

For repeated failures, include the command, software version and relevant
`elarodon.log` messages in a [GitHub issue](https://github.com/aakechin/eLaRodON/issues).
For `Too many open files`, check the process file-descriptor limit with
`ulimit -n`; on systems where it is permitted, `ulimit -n 4096` raises the soft
limit for the current shell.

## Citation and license

Please cite the eLaRodON article, *eLaRodON: identification of large genomic
rearrangements in Oxford Nanopore sequencing data*, and record the software
version used in your analysis.

The [Zenodo software record](https://doi.org/10.5281/zenodo.19177352) provides
an archival citation. Record the version of the archive actually used, alongside
the GitHub tag or commit and the installed package version.

[MIT License](LICENSE) covers the original software and documentation.
Third-party datasets retain their source terms.
