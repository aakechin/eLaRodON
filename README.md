# eLaRodON

A comprehensive pipeline for detection and analysis of large genomic rearrangements from Oxford Nanopore Technologies (ONT) sequencing data.

![Python](https://img.shields.io/badge/Python-3.7%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
  - [Input Requirements](#input-requirements)
- [Algorithm By Steps](#algorithm-by-steps)
- [Usage](#usage)
  - [Minimal Working Example](#minimal-working-example)
  - [Arguments](#arguments)
  - [Full Analysis](#full-analysis)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

# Overview

The ONTLRcaller pipeline performs detection and characterization of large genomic rearrangements through four integrated modules:

[![](https://mermaid.ink/img/pako:eNpVkFFrgzAUhf9KuM-26NQa8zBodYXC6KAbe5gWCXptHZpIjNu62v--aDfK7lNOzndykpwhlwUCg4Pi7ZG8xKkgZpbJ0_blcZfzuka1J7PZ_bDjn2TU3UBWybusxER0-2tgNTHrqtaosCDGGEiU8Lo6iAaFzjbb53_kRnSodCXFyF2NaDKWQkjNNRYDiZMCy0pgpk8tZrlCs5195GV2q41_awWviex12-uBPCSv0ZqUVY17sMy7qgKYVj1a0KBq-CjhPMZT0EdsMAVmlqaK97VOIRUXE2u5eJOy-Usq2R-OwEped0b1bWGuElfcfNoNQVGgimQvNDAnmI4AdoYvYG7gzD0npAubuiH13MC34ATMC-bUufMdx_M8O3ToxYLvqdOe08C3zQSh6_shpYvLD7JJhTs?type=png)](https://mermaid.live/edit#pako:eNpVkFFrgzAUhf9KuM-26NQa8zBodYXC6KAbe5gWCXptHZpIjNu62v--aDfK7lNOzndykpwhlwUCg4Pi7ZG8xKkgZpbJ0_blcZfzuka1J7PZ_bDjn2TU3UBWybusxER0-2tgNTHrqtaosCDGGEiU8Lo6iAaFzjbb53_kRnSodCXFyF2NaDKWQkjNNRYDiZMCy0pgpk8tZrlCs5195GV2q41_awWviex12-uBPCSv0ZqUVY17sMy7qgKYVj1a0KBq-CjhPMZT0EdsMAVmlqaK97VOIRUXE2u5eJOy-Usq2R-OwEped0b1bWGuElfcfNoNQVGgimQvNDAnmI4AdoYvYG7gzD0npAubuiH13MC34ATMC-bUufMdx_M8O3ToxYLvqdOe08C3zQSh6_shpYvLD7JJhTs)

## Detailed description of the algorithm

The eLaRodON algorithm is a specialized computational pipeline designed for comprehensive detection of large genomic rearrangements (LGRs) from Oxford Nanopore sequencing data. Unlike conventional tools developed primarily for germline variants, eLaRodON incorporates several innovative features specifically optimized for identifying somatic LGRs, including those supported by single reads.

**Input Processing**

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

## Dependencies:

**eLaRodON** requires the following external tools to be installed and available in your PATH:

- minimap2: `conda install -c bioconda minimap2` or download from [GitHub](https://github.com/lh3/minimap2)

- samtools & htslib: `conda install -c bioconda samtools`

- vcfanno: Download the pre-built binary from the [release page](https://github.com/brentp/vcfanno/releases).

> 💡 Pro tip: see full list of program versions in [requirements.txt](/requirements.txt)

## Install eLaRodON via pip
```bash
pip install elarodon
```

### Verify Installation

```bash
elarodon -h
```

## Also you can clone the repository and run without installation, but you'll still need to install the dependencies via pip 

```bash
git clone https://github.com/aakechin/eLaRodON.git
cd eLaRodON
pip install -r requirements.txt
```

### Verify Installation

```bash
./elarodon/main.py -h
```

# Input Requirements

## BAM-file

- **Sorted & Indexed**: Must be coordinate-sorted and indexed (e.g., with `samtools sort` and `samtools index`)
- **Consistent Reference**: Must be aligned to the **exact same reference genome** as provided via `-ref`
- **ONT-specific Alignment**: Optimal results with alignments that preserve soft-clipping (recommended aligners: `minimap2`)

**Example**:

```bash
samtools index sample.bam
```

## Reference genome

- **Indexed**: Must be indexed with `samtools faidx`
- **Consistent**: Must match the reference used for BAM alignment

## BED-file

- **Format**: Standard 3-column BED format (chromosome, start, end) **without header**
- **Purpose**: Used to annotate variants with genomic features (repeats, genes, etc.)
- **Example**:
```bed
chr1  10000  10400  ALU
chr1  15000  15200  LINE1
chr2  5000   5200   SVA
```

# Algorithm By Steps

## Step 1: Breakpoint Detection (ONTLRcaller)
- **Clipped Read Analysis**: Scans BAM for reads with soft/hard-clipped segments exceeding `--minimal-clipped-length`.
- **Breakpoint Clustering**: Groups nearby breakpoints within `--dist-to-join-trl` distance into candidate Structural Variants (SVs).
- **Initial Classification**: Separates candidates into Large Rearrangements (LRs) and potential Insertions (INS).

## Step 2: Variant Consolidation (joinONTLRs)
- **Proximity-based Joining**: Merges neighboring SVs within `--maximal-distance-join` bp into single events.
- **Filtering**: Removes artifacts and low-confidence calls based on supporting read counts and mapping quality.
- **Output Separation**: Generates two CSV files: LRs (deletions, duplications, inversions) and INS (insertion candidates).

## Step 3: Insertion Resolution & Typing
- **Local Assembly**: Extracts unaligned sequences from insertion candidates.
- **Mini-assembly**: Performs local assembly of clipped sequences.
- **Re-alignment**: Aligns assembled contigs to reference using `minimap2` to determine precise insertion points and structure.

## Step 4: Variant Annotation & VCF Generation
- **Type Determination**: Classifies variants based on breakpoint orientation, sequence homology, and alignment patterns.
- **Microhomology Analysis**: Calculates microhomology/homeology at breakpoints.
- **Repeat Annotation**: Annotates with nearby genomic features using provided BED file via `vcfanno`.
- **VCF Export**: Generates `VCF 4.2` compliant output with comprehensive INFO fields.

# Usage

## Minimal Working Example
```bash
elarodon \
    -dir ./results_elarodon \
    -bam sample.bam \
    -ref hg38.fa \
    -vcfanno path/to/vcfanno \
    -th 4
```
## Arguments

### Core Arguments

| Parameter    | Required     | Description          |
| :---        |    :----:   |          ---: |
| `-bam, --bam_file`         | Yes         | BAM file            |
| `-dir, --workdir` |	Yes |	Output directory |
| `-ref, --ref-genome` |	Yes |	Reference genome FASTA|
| `-vcfanno, --vcf-anno`	|Yes	|vcfanno executable path |

### Processing Parameters

|Parameter	| Default	 | Required | Description |
| :---        |    :----:   |  :----:   |        ---: |
|`-bed, --bed-file`|	False | No	| Annotation BED file |
|`-div, --divide-chroms`|	False | No |	Divide chromosome to analyze |
|`-dvlen, --div-length`|	None | No |	Length of regions for division chromosome to analyze |
|`-len, --minimal-length`|	50| No |Min variant length (bp)|
|`-clip`, `--minimal-clipped-length`|	100	| No | Min clipped length (bp)|
|`-dist`,`-dist-to-join-trl`|	1000	| No | Min clipped length (bp)|
|`-join`,`-maximal-distance-join`|	30	| No | Max distance for fusion joining (bp)|
|`-th, --threads` |	4	| No | CPU threads|
|`-cont, --continue` |	all	| No | Name of stage for start: bam, join, def|

About `--continue`:

| Value   | Stage Name	 | Description     | Required Input Files for This Stage          |
| :---       |    :----: |    :----:   |          ---: |
| `all (default)`        | 	Full pipeline |  Runs all stages sequentially: `bam` → `join` → `def` | BAM file, reference genome  |
| `bam` | ONTLRcaller |	**Initial detection**: Parses BAM file to identify clipped reads and cluster breakpoint |	BAM file |
| `join` | joinONTLRs |	**Joining & filtering**: Joins nearby breakpoints, separates LRs from INS candidates |	CSV files from `bam` stage (`*junction_stat.csv`) |
| `def` | define_type_create_vcf |	**Alignment & annotation**: Aligns INS sequences, classifies variants, creates final VCF | CSV files from `join` stage (`*LRs_join*.csv`, `*INS_join*.csv`) |

For example:
```bash
# Re-run only the annotation with a different BED file
elarodon -cont def -dir ./results -ref hg38.fa -vcfanno ./vcfanno -bed new_annotation.bed
```

### Special file names

| Parameter   | Default	 | Required     | Description          |
| :---       |    :----: |    :----:   |          ---: |
| `-in, --input-files`        | auto  | No         | Regular expression for CSV files          |
| `-lrs, --output-lrs` | auto |	No |	CSV output file for LGRs |
| `-ins, --output-ins` | auto |	No |	CSV output file for INS |
| `-sam, --sam_file` | auto |	No | SAM file with INS alignment |

### Output Control
| Parameter	| Default |	Description |
| :---        |    :----:   |          ---: |
|`-out, --out-vcf` |	auto|	VCF output filename|
|`-nrt_ins, --not-remove-trash-align`|	False|	Keep temp alignment files|
|`-nrt_anno, --not-remove-trash-anno` |	False	|Keep temp annotation files|

**To view all parameters and their descriptions, you can use**:

```bash
elarodon -h
```

## Detailed Usage

```bash
elarodon \
    -bam sample.bam \
    -dir lr_results \
    -ref hg38.fa \
    -vcfanno ~/tools/vcfanno \
    -bed repeats.bed \
    -th 8 \
```

# Output Files

## Main Outputs

    *.junction_stat.LRs_join100.csv - Merged large rearrangements

    *.junction_stat.INS_join100.csv - Insertion calls

    *_all_LGRS.vcf - Final annotated variants
    
## VCF Output

eLaRodON enriches VCF output with extensive annotation in the **INFO** column:

### 1. Basic Variant Information
| Field   | Type     | Description          |
| :---       |    :----: |    :----:   |  
| `SVTYPE`        |  String         | Type of structural variant (DEL, INS, INV, etc.)          |
| `CHROM2` | String |	Second chromosome involved (for translocations) |
| `END` | Integer |	End position of the variant |
| `SVLEN` | Integer | Length of variant (exact for DEL/INS, approximate for DUP/INV) |
| `TDRN` | Integer | Number of tandem repeats (for duplications) |

### 2. Breakpoint & Read Evidence
| Field   | Type     | Description          |
| :---       |    :----: |    :----:   |   
| `J1, J2`        |  String         | Clipped side relative to breakpoint (L=left, R=right, ND=not defined)          |
| `S1, S2` | String |	Strand orientation of supporting reads |
| `D1, D2` | Integer |	Distance between breakpoints in paired evidence |
| `RL1, RL2` | Integer | Lengths of supporting read fragments |
| `SR` | Integer | Number of supporting reads |
| `MQ` | Float | Median mapping quality of supporting reads |
| `DP` | Integer | Total read depth at breakpoint |
| `VAF` | Float | Variant allele frequency |

### 3. Sequence Features
| Field   | Type     | Description          |
| :---       |    :----: |    :----:   |   
| `MH, HOM`        |  Integer         | Length of microhomology/homeology at breakpoints          |
| `MHS, HOMS` | String |	Microhomology/homeology sequences |
| `MUTM, MUTV` | Float |	Median and variance of mapping errors in read fragments |

### 4. Variant Classification Metrics
| Field   | Type     | Description          |
| :---       |    :----: |    :----:   |   
| `SBLR, SBTRL, SBINV, SBTD`        |  String         | Evidence for alternative variant classifications   |
| `ISN, NSN, SLN` | Integer |	Counts of different read mapping patterns |
| `ISNMFS, NSNMFS` | Integer |	Length of most frequent overlapping/unmapped sequences |
| `CN, IVN` | Integer | Reads with specific CIGAR patterns or inversions |
| `MNF, PFF, PLF` | Integer | Homopolymeric tract pattern in novel sequence |
| `NSP` | Integer | Metrics on rearrangement fragment distribution |

### 5. Genomic Context Annotation

*Applied only when BED file is provided via* `-bed`

| Field   | Type     | Description          |
| :---       |    :----: |    :----:   | 
| `LERN, RERN`        |  String         | Nearest repeat/element to left/right breakpoint (external)s   |
| `LIRN, RIRN` | Integer |	Nearest repeat/element within left/right segment (internal) |
| `LCRN, RCRN` | Integer |	Repeat/element containing the breakpoint |
| `LERD, RERD` | Integer | Distance to nearest external repeat |
| `LIRD, RIRD` | Integer | Distance to nearest internal repeat |

## Variant Type-Specific Tags

The `ALT` field in VCF header defines all possible variant types:

- `<DEL>`: Deletion (confirmed by multiple reads)
- `<BND_DEL>`: Deletion (single-read evidence)
- `<INS>`: Novel sequence insertion
- `<TD>`: Tandem duplication
- `<BND_TD>`: Possible tandem duplication
- `<INV>`, `<BND_INV>`: Inversion (confirmed/unconfirmed)
- `<INVTD>`, `<BND_INVTD>`: Inverted tandem duplication
- `<TRL>`, `<BND_TRL>`: Translocation

## Quality Filtering

Variants are tagged in the `FILTER` column:

- `PASS`: Confidently typed rearrangement
- `FAIL`: Could not determine rearrangement type (or SV may be *not real*)

## Expected VCF Example

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
**Solution №1**: Reduce thread count or increase memory

```bash
Too many open files
```
**Solution №2**: 

Run before eLaRodON: 
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
