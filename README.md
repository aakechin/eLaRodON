# eLaRodON

A comprehensive pipeline for detection and analysis of large genomic rearrangements from Oxford Nanopore Technologies (ONT) sequencing data.

![Python](https://img.shields.io/badge/Python-3.7%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Detailed Usage](#detailed-usage)
- [Output Files](#output-files)
- [Example](#example)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

# Overview

The ONTLRcaller pipeline performs detection and characterization of large genomic rearrangements through four integrated modules:

[![](https://mermaid.ink/img/pako:eNpVkFFrgzAUhf9KuM-26NQa8zBodYXC6KAbe5gWCXptHZpIjNu62v--aDfK7lNOzndykpwhlwUCg4Pi7ZG8xKkgZpbJ0_blcZfzuka1J7PZ_bDjn2TU3UBWybusxER0-2tgNTHrqtaosCDGGEiU8Lo6iAaFzjbb53_kRnSodCXFyF2NaDKWQkjNNRYDiZMCy0pgpk8tZrlCs5195GV2q41_awWviex12-uBPCSv0ZqUVY17sMy7qgKYVj1a0KBq-CjhPMZT0EdsMAVmlqaK97VOIRUXE2u5eJOy-Usq2R-OwEped0b1bWGuElfcfNoNQVGgimQvNDAnmI4AdoYvYG7gzD0npAubuiH13MC34ATMC-bUufMdx_M8O3ToxYLvqdOe08C3zQSh6_shpYvLD7JJhTs?type=png)](https://mermaid.live/edit#pako:eNpVkFFrgzAUhf9KuM-26NQa8zBodYXC6KAbe5gWCXptHZpIjNu62v--aDfK7lNOzndykpwhlwUCg4Pi7ZG8xKkgZpbJ0_blcZfzuka1J7PZ_bDjn2TU3UBWybusxER0-2tgNTHrqtaosCDGGEiU8Lo6iAaFzjbb53_kRnSodCXFyF2NaDKWQkjNNRYDiZMCy0pgpk8tZrlCs5195GV2q41_awWviex12-uBPCSv0ZqUVY17sMy7qgKYVj1a0KBq-CjhPMZT0EdsMAVmlqaK97VOIRUXE2u5eJOy-Usq2R-OwEped0b1bWGuElfcfNoNQVGgimQvNDAnmI4AdoYvYG7gzD0npAubuiH13MC34ATMC-bUufMdx_M8O3ToxYLvqdOe08C3zQSh6_shpYvLD7JJhTs)

## Detailed description of the algorithm

The eLaRodON algorithm is a specialized computational pipeline designed for comprehensive detection of large genomic rearrangements (LGRs) from Oxford Nanopore sequencing data. Unlike conventional tools developed primarily for germline variants, eLaRodON incorporates several innovative features specifically optimized for identifying somatic LGRs, including those supported by single reads.
Input Processing

The algorithm begins by processing aligned sequencing data in BAM format. It performs chromosome-by-chromosome analysis to optimize memory usage, with an option to focus on specific genomic regions of interest. The tool only considers primary alignments containing complete mapping information to ensure analysis quality.

### Core Detection Mechanism

#### 1. Split-read Analysis:

- Identifies reads with fragments mapped to different genomic locations
- Detects strand changes and structural variants >50bp through CIGAR tag analysis
- Extracts all insertions and deletions from primary alignments

#### 2. Junction Characterization:

##### Records all junction sites with their genomic features in two separate files:

- fusions.csv: Contains genome region junctions
- insertions.csv: Stores CIGAR-derived insertions

##### For fusion events, combines split reads corresponding to:

- Translocation boundaries (TRL)
- Inversion breakpoints (INV)
- Tandem duplication junctions (TD)

### Variant Classification

The algorithm employs a sophisticated classification system that:

#### Merges Similar Events:

- Combines fusions and insertions across genomic regions
- Uses precise coordinates, strand orientation, and junction characteristics
- Optional merging via companion script for mechanistic studies

#### Insertion Sequence Analysis:

- Maps CIGAR-derived insertion sequences to reference genome using minimap2
- Reclassifies as tandem duplications when sequences map near original positions

#### Structural Annotation:

- Determines LGR types using strand orientation and junction characteristics
- Evaluates four key genomic features for each rearrangement:

    - Proximity to repeat sequences and mobile elements (via vcfanno)
    - Presence of 2-5 nucleotide microhomology
    - ≥80% sequence similarity over ≤35 nucleotide homeology between breakpoints
    - Presence of novel inserted sequences not matching either breakpoint region

### Quality Assessment

eLaRodON incorporates several quality control metrics:

#### New Sequence Pattern (NSP) Scoring:

    score=∑x ∈[An, Tn,Gn,Cn]1L(x)−L(s)/100score=∑x ∈[An​, Tn​,Gn​,Cn​]​1L(x)​−L(s)/100

    - Identifies homopolymeric tracts (≥4bp) in novel junction sequences
    - Helps discriminate true rearrangements from artifacts

### Output Generation

The final output includes:

#### 1. Comprehensive VCF File:

- Special tags for LGRs sharing coordinates but differing in variant type
- Annotations for:
    - Intersection sequences (ISMFS tag)
    - Novel inserted sequences (NSMFS tag)
    - Microhomology/homeology patterns
    - Repeat element proximity

#### 2. Additional Outputs:

- Raw junction calls for downstream analysis
- Quality metrics for each detected variant
- Intermediate files for debugging and method development

### Technical Innovations

Key algorithmic advancements include:

- Treatment of multiple split-reads as single sequences for accurate two-boundary detection
- Optional non-merging of similar LGRs to preserve mechanistic signatures
- Memory-efficient chromosomal processing
- Specialized handling of Nanopore-specific artifacts

The tool demonstrates particular strength in identifying complex rearrangements that most existing algorithms miss, including multiple tandem duplications, non-reciprocal translocations, and inversions with two defined boundaries. Its performance has been validated across multiple datasets, showing superior accuracy compared to existing tools like Sniffles2, NanoSV, and SVIM, particularly for variants in repetitive regions and those with low read support.

# Installation

## Prerequisites:

    Python 3.7+

    minimap2 (v2.24+)

    vcfanno (v0.3.3+)

    samtools (v1.15+)

    htslib (v1.21+)

    git (for cloning)

> 💡 Pro tip: see full list of program versions in [requirements.txt](/requirements.txt)

## 1. Clone the Repository

    git clone https://github.com/aakechin/eLaRodON.git

## 2. Install Dependencies

    pip install -r requirements.txt

## 3. Verify Installation

    python3 ./eLaRodON-main/main.py -h

> No errors = successful installation

# Quick Start

### Basic Pipeline Execution
```python
python3 main.py \
    -dir ./results_elarodon \
    -bam sample.bam \
    -ref hg38.fa \
    -vcfanno path/to/vcfanno \
    -bed annotations.bed \
    -th 4
```

## Detailed Usage

**Core Arguments**

| Parameter    | Required     | Description          |
| :---        |    :----:   |          ---: |
| `-bam, --bam_file`         | Yes         | BAM file            |
| `-dir, --workdir` |	Yes |	Output directory |
| `-ref, --ref-genome` |	Yes |	Reference genome FASTA|
|`-bed, --bed-file`|	Yes	| Annotation BED file |
| `-vcfanno, --vcf-anno`	|Yes	|vcfanno executable path |

**Processing Parameters**

|Parameter	| Default	 | Required | Description |
| :---        |    :----:   |  :----:   |        ---: |
|`-div, --divide-chroms`|	False | No |	Divide chromosome to analyze |
|`-dvlen, --div-length`|	None | No |	Length of regions for division chromosome to analyze |
|`-len, --minimal-length`|	50| No |Min variant length (bp)|
|`-clip`, `--minimal-clipped-length`|	100	| No | Min clipped length (bp)|
|`-dist`,`-dist-to-join-trl`|	1000	| No | Min clipped length (bp)|
|`-join`,`-maximal-distance-join`|	30	| No | Max distance for fusion joining (bp)|
|`-th, --threads` |	4	| No | CPU threads|
|`-cont, --continue` |	all	| No | Name of stage for start: bam, join, def|

**Special file names**

| Parameter   | Default	 | Required     | Description          |
| :---       |    :----: |    :----:   |          ---: |
| `-in, --input-files`        | auto  | No         | Regular expression for CSV files          |
| `-lrs, --output-lrs` | auto |	No |	CSV output file for LGRs |
| `-ins, --output-ins` | auto |	No |	CSV output file for INS |
| `-sam, --sam_file` | auto |	No | SAM file with INS alignment |

**Output Control**
| Parameter	| Default |	Description |
| :---        |    :----:   |          ---: |
|`-out, --out-vcf` |	auto|	VCF output filename|
|`-nrt_ins, --not-remove-trash-align`|	False|	Keep temp alignment files|
|`-nrt_anno, --not-remove-trash-anno` |	False	|Keep temp annotation files|

*To view all parameters and their descriptions, you can use* 

```python
python3 main.py -h
```

# Output Files

## Main Outputs

    *.junction_stat.LRs_join100.csv - Merged large rearrangements

    *.junction_stat.INS_join100.csv - Insertion calls

    *_all_LGRS.vcf - Final annotated variants
    
## Example
Input Preparation
```bash
samtools index sample.bam
```

Full Analysis
```python
python3 main.py \
    -bam sample.bam \
    -dir lr_results \
    -ref hg38.fa \
    -vcfanno ~/tools/vcfanno \
    -bed repeats.bed \
    -th 8 \
```

Expected Output
vcf

|#CHROM  |POS   |  ID  |    REF | ALT   |   QUAL | FILTER | INFO  | FORMAT	| SAMPLE_NAME |
|    :----:   |    :----: |    :----: |    :----: |    :----: |    :----: |    :----: |    :----: |    :----: |    :----:   
|chr12   |3456789| LR1|       N   | DEL   | 0.78  |  PASS   | SVTYPE=DEL;SVLEN=1200;CHROM2=chr12...| GT:DP:AD:VF|	0/1:330:1,329:0.0|
|chr13   |4123456 | INS1 |     N  |  INS | 0.65 |   PASS |   SVTYPE=INS;SVLEN=350...| GT:DP:AD:VF	| 0/1:244:1,243:0.0 |

# Troubleshooting

## Common Issues

### Missing dependencies:
```bash
Error: minimap2 not found in PATH
```

**Solution**: Install via bioconda or add to PATH

### Memory errors:
```bash
Killed (process exited)
```
**Solution**: Reduce thread count or increase memory

```bash
Too many open files
```
**Solution**: 

Run before eLaRodN: 
```bash
ulimit -n 4096
```

### BAM index missing:
```bash
    [E::idx_find_and_load] Could not retrieve index file for 'sample.bam'
```

**Solution**: Run samtools index sample.bam

# Citation

Please cite:

    eLaRodON: identification of large genomic rearrangements in Oxford Nanopore sequencing data

# License

[MIT License](/LICENSE)
