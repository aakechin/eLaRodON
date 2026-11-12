# eLaRodON

eLaRodON detects and characterizes genomic rearrangements from Oxford Nanopore
read alignments. It reports rearrangement calls in VCF and intermediate tables.

Version 1.1.1 packages the analysis code from the revised GitHub benchmark
release v1.1.0 and fixes handling of the last junction during merging. A single
junction is now retained, the final junction can be merged with a nearby call,
and writing an empty junction table no longer divides by zero. These changes
can affect calls relative to v1.1.0; the published v1.1.0 benchmark results were
not recomputed for this package. The earlier PyPI release 0.1.3 predates the
revision corrections.

## Installation

Use Linux and Python 3.10, 3.11 or 3.12. In a Python virtual environment, run:

```bash
python -m pip install elarodon==1.1.1
elarodon --help
```

Python dependencies are installed by pip. Install these external command-line
tools separately and make them available on `PATH`:

- minimap2
- samtools
- bgzip and tabix from HTSlib
- bedtools
- vcfanno

The [repository documentation](https://github.com/aakechin/eLaRodON#readme)
lists the tool versions and the input and output formats.

## Run on your own data

Provide a coordinate-sorted, indexed BAM and the matching indexed reference
FASTA. Use paths without spaces for the external annotation commands.

```bash
elarodon \
    -dir /path/to/results \
    -bam /path/to/sample.bam \
    -ref /path/to/reference.fa \
    -vcfanno /path/to/vcfanno \
    -th 4
```

The default final output is `sample_all_LGRS.vcf` in the results directory.
`python -m elarodon` accepts the same arguments as `elarodon`.

## Benchmark data and source code

The pip package contains the program, license and package metadata. Reference
genomes, BAM files, curated benchmark data and benchmark scripts are distributed
separately through the
[GitHub repository](https://github.com/aakechin/eLaRodON) and the
[v1.1.0 benchmark release](https://github.com/aakechin/eLaRodON/releases/tag/v1.1.0).
See the [benchmark guide](https://github.com/aakechin/eLaRodON/tree/v1.1.0/benchmark)
for the test dataset and Docker workflow.

## License

MIT. Copyright Andrey Kechin and Regina Mikheeva.
