# SeqMatcher Pipeline

A Nextflow DSL2 pipeline for sequence motif discovery, PWM-based scanning, and genomic track generation.

```
input CSVs ──► FETCH_SEQUENCES ──► MATCH_AND_MOTIFS ──► INDEX_BIGWIG
                   │                      │
              fetched TSVs        matches TSV, BED, BigWig,
                                  PSSM, logos, MEME files
```

---

## Table of Contents

1. [What the pipeline does](#what-the-pipeline-does)
2. [Requirements](#requirements)
3. [Installation](#installation)
4. [Input format](#input-format)
5. [Running the pipeline](#running-the-pipeline)
6. [Parameters](#parameters)
7. [Outputs](#outputs)
8. [Execution profiles](#execution-profiles)
9. [Building the Docker image](#building-the-docker-image)
10. [Directory structure](#directory-structure)
11. [JASPAR source-of-truth design](#jaspar-source-of-truth-design)

---

## What the pipeline does

| Step | Process | Script | Description |
|------|---------|--------|-------------|
| 1 | `FETCH_SEQUENCES` | `fetch_sequence.py` | Extracts genomic sequences for each peak interval using pyfaidx. Handles both strands. |
| 2 | `MATCH_AND_MOTIFS` | `sequence_matcher_with_motifs.py` | Single-pass scan of the forward peak sequence with three modes: **exact** (IUPAC regex), **mismatch** (Hamming ≤ N), or **pwm** (log-odds + Monte-Carlo p-value). In all modes, when `--jaspar` is supplied the canonical consensus is derived from the counts matrix via `consensus_from_counts()` — the same function used in `find_motif_both_strands_and_chromsome_wise.py` — so every peak hit is guaranteed to be a strict subset of the genome-wide scan. Minus strand hits are found by scoring `rc(window)` at each position. Coordinates are always reported in + strand BED space (`peak_start + local_offset`). Emits a merged TSV, 6-column BED, BigWig coverage track, per-strand logos, PSSM tables, and MEME motif files. |
| 3 | `INDEX_BIGWIG` | UCSC `bigWigInfo` | Validates the BigWig and writes a summary info file so genome browsers can load it. |

---

## Requirements

| Tool | Minimum version | Notes |
|------|----------------|-------|
| Nextflow | 23.04 | `curl -s get.nextflow.io \| bash` |
| Docker **or** Conda | any | See [Execution profiles](#execution-profiles) |
| Python | 3.11 | Inside container/env only |

Python package dependencies (managed automatically via Docker or Conda):

```
biopython   pandas   numpy   pyfaidx   pyBigWig   logomaker   matplotlib
```

UCSC binary required for `INDEX_BIGWIG`:

```
bigWigInfo   (bioconda: ucsc-bigwiginfo)
```

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/your-lab/seqmatcher-pipeline.git
cd seqmatcher-pipeline

# 2. (Docker users) Build the image once
docker build -t seqmatcher:latest .

# 3. (Conda users) Create the environments
conda env create -f envs/seqmatcher.yml
conda env create -f envs/ucsc.yml
```

---

## Input format

Place one or more CSV files inside `input/` (or the path given by `--input_dir`).

**Required columns** (header names are matched case-insensitively):

| Column | Aliases accepted | Description |
|--------|-----------------|-------------|
| `GeneID` | — | Unique identifier for the peak |
| `Chromosome` | `chrom`, `chr` | Chromosome / contig name matching the FASTA |
| `Peak_Start` | `start`, `chromStart` | 0-based start coordinate |
| `Peak_End` | `end`, `chromEnd` | 0-based exclusive end coordinate |
| `Strand` | — | `+` or `-` |

**Example** (`input/example_peaks.csv`):

```
GeneID,Chromosome,Peak_Start,Peak_End,Strand
gene_001,chr1,1000,1200,+
gene_002,chr1,5000,5250,-
```

Place the genome FASTA at `input/genome/genome.fa` (configurable via `--genome`).  
A `.fai` index is created automatically by pyfaidx if not present.

---

## Running the pipeline

### Recommended: mismatch mode driven by JASPAR file

This is the preferred way to run the pipeline. Providing `--jaspar` means the
canonical consensus is derived from the counts matrix — identical logic to
`find_motif_both_strands_and_chromsome_wise.py` — guaranteeing that every peak
hit is a strict subset of any genome-wide scan run with the same JASPAR file.

```bash
nextflow run main.nf \
    --jaspar     path/to/AGGGGATTTCCC.jaspar \
    --mode       mismatch \
    --mismatches 2 \
    --genome     input/genome/genome.fa \
    -profile     docker
```

> **Note:** `--query` is optional when `--jaspar` is provided. If both are
> given, the JASPAR consensus takes precedence and a warning is printed if they
> differ.

### Mismatch mode — fallback without JASPAR

Only use this when you do not have a JASPAR file. The `--query` string drives
the consensus directly, so cross-tool subset consistency requires you to ensure
`--query` matches the JASPAR consensus used in the genome-wide script manually.

```bash
nextflow run main.nf \
    --query      "AGGGGATTTCCC" \
    --mode       mismatch \
    --mismatches 2 \
    --genome     input/genome/genome.fa \
    -profile     docker
```

### PWM mode with JASPAR file (required)

```bash
nextflow run main.nf \
    --jaspar path/to/AGGGGATTTCCC.jaspar \
    --mode   pwm \
    --pvalue 1e-4 \
    --bg_gc  0.42 \
    -profile docker
```

`--jaspar` is mandatory for PWM mode — the pipeline fast-fails with an
explanatory error if it is omitted.

### Exact mode

```bash
nextflow run main.nf \
    --jaspar path/to/AGGGGATTTCCC.jaspar \
    --mode   exact \
    -profile docker
```

### Forward strand only

```bash
nextflow run main.nf \
    --jaspar   path/to/AGGGGATTTCCC.jaspar \
    --mode     mismatch \
    --fwd_only true \
    -profile   docker
```

### Conda (no Docker)

```bash
nextflow run main.nf \
    --jaspar     path/to/AGGGGATTTCCC.jaspar \
    --mismatches 2 \
    -profile     conda
```

### SLURM cluster

```bash
nextflow run main.nf \
    --jaspar     path/to/AGGGGATTTCCC.jaspar \
    --mismatches 2 \
    --outdir     /scratch/myproject/results \
    -profile     slurm
```

### Resume after failure

```bash
nextflow run main.nf [options] -resume
```

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input_dir` | `input` | Directory containing input CSV files |
| `--genome` | `input/genome/genome.fa` | Reference genome FASTA |
| `--jaspar` | `null` | **Recommended for all modes.** JASPAR-format counts file. The canonical consensus is derived from this file via `consensus_from_counts()`, identical to the genome-wide script, guaranteeing subset consistency. Required for `--mode pwm`. |
| `--query` | `AGGGGATTTCCC` | Fallback consensus string used when `--jaspar` is absent. When both are provided, the JASPAR consensus takes precedence and a warning is printed if they differ. |
| `--mode` | `mismatch` | Search mode: `exact` (IUPAC regex), `mismatch` (Hamming), or `pwm` (log-odds + p-value) |
| `--mismatches` | `2` | Max Hamming mismatches (`mismatch` mode) or hard post-filter cap (`pwm` mode) |
| `--pvalue` | `1e-4` | P-value threshold for PWM hit calling (PWM mode only) |
| `--pseudocount` | `0.1` | Pseudocount added to each PWM cell before log-odds conversion |
| `--bg_gc` | `0.5` | Background GC fraction for the Monte-Carlo null distribution (e.g. `0.42` for human) |
| `--fwd_only` | `false` | Scan forward strand only; default scans both strands |
| `--upstream` | `150` | Flanking bases to extract upstream of each hit |
| `--downstream` | `150` | Flanking bases to extract downstream of each hit |
| `--outdir` | `results` | Output directory |
| `--use_docker` | `true` | Use Docker container (set `false` when running with Conda profile) |

---

## Outputs

```
results/
├── 01_sequences/
│   └── <sample>_with_seq.tsv          # input CSV + Fetched_sequence column
│
├── 02_matches/
│   └── <sample>_matches.tsv           # all hits, one row per hit
│                                      # columns: Genomic_Start, Genomic_End,
│                                      # Strand, Matched_Sequence, Mismatches,
│                                      # Score, Upstream, Downstream, Mode
│                                      # (+ Mismatch_Positions in pwm mode)
│
├── 03_motifs/
│   ├── hit_motif_logo.png             # IC logo for matched motif windows
│   │     (or jaspar_motif_logo.png when --jaspar is supplied in pwm mode)
│   ├── hit_motif_pssm.tsv             # frequency PWM (tab-separated)
│   ├── hit_motif.meme                 # MEME-format motif file
│   ├── +_up_logo.png                  # logo for upstream flanks (plus strand)
│   ├── +_down_logo.png                # logo for downstream flanks (plus strand)
│   ├── -_up_logo.png
│   ├── -_down_logo.png
│   └── *.meme / *_pssm.tsv           # matching MEME and PSSM files
│
├── 04_bed/
│   └── <sample>_hits.bed              # 6-column BED (chrom, start, end, name,
│                                      # score, strand)
│
├── 05_bigwig/
│   ├── <sample>_hits.bw               # per-base score coverage BigWig
│   ├── <sample>_hits_chrom_sizes.txt  # chromosome sizes (for genome browsers)
│   └── <sample>_hits_info.txt         # bigWigInfo validation summary
│
└── pipeline_info/
    ├── timeline.html
    ├── report.html
    ├── trace.txt
    └── dag.svg
```

### TSV output columns

| Column | Description |
|--------|-------------|
| *(input columns)* | All original CSV columns |
| `Strand` | `+` or `-` |
| `Matched_Sequence` | Matched window, 5'→3' on the reported strand |
| `Mismatches` | Hamming distance from consensus |
| `Genomic_Start` | True genomic start = `Peak_Start` + local offset |
| `Genomic_End` | True genomic end = `Genomic_Start` + motif length |
| `Score` | BED score (exact/mismatch) or raw log-odds score (PWM) |
| `Upstream` | Upstream flank, oriented 5'→3' on reported strand |
| `Downstream` | Downstream flank, oriented 5'→3' on reported strand |
| `Mode` | `exact`, `mismatch`, or `pwm` |
| `Mismatch_Positions` | 0-based mismatch positions within window (PWM mode only; `.` if none) |

---

## Execution profiles

| Profile | Executor | Container |
|---------|----------|-----------|
| `docker` | local | Docker (`seqmatcher:latest`) — ARM64 (Apple Silicon) |
| `docker_amd64` | local | Docker (`seqmatcher:latest`) — x86_64 (Intel/AMD) |
| `conda` | local | Conda (`envs/seqmatcher.yml`, `envs/ucsc.yml`) |
| `slurm` | SLURM | Docker |
| `slurm_conda` | SLURM | Conda |

Switch profiles with `-profile <name>` on the CLI.

---

## Building the Docker image

```bash
# From the repository root (auto-detects host architecture)
./build.sh

# Verify
docker run --rm seqmatcher:latest python -c "import pyBigWig, pyfaidx, logomaker; print('OK')"
```

To push to a registry:

```bash
./build.sh --push --registry=ghcr.io/your-org
```

Update `container 'seqmatcher:latest'` in `main.nf` to point to the registry URL.

---

## Directory structure

```
seqmatcher-pipeline/
├── main.nf                         # Nextflow workflow definition
├── nextflow.config                 # Process resources & execution profiles
├── Dockerfile                      # Container definition
├── build.sh                        # Docker build helper (auto-detects arch)
├── README.md
│
├── scripts/
│   ├── fetch_sequence.py           # Step 1 – genomic sequence fetcher
│   └── sequence_matcher_with_motifs.py  # Step 2 – motif search + outputs (v3)
│
├── envs/
│   ├── requirements.txt            # pip requirements (used in Dockerfile)
│   ├── seqmatcher.yml              # Conda env for Steps 1 & 2
│   └── ucsc.yml                    # Conda env for Step 3 (bigWigInfo)
│
└── input/
    ├── example_peaks.csv           # Example input
    └── genome/
        └── genome.fa               # Place reference FASTA here
```

---

## JASPAR source-of-truth design

### Why the JASPAR file must be supplied to both scripts

The pipeline's `MATCH_AND_MOTIFS` process and the standalone genome-wide scanner
(`find_motif_both_strands_and_chromsome_wise.py`) both derive their canonical
motif string from the same function:

```python
def consensus_from_counts(counts: dict) -> str:
    n_pos = len(counts['A'])
    return ''.join(max(BASES, key=lambda b: counts[b][i]) for i in range(n_pos))
```

When both tools are given the **same `.jaspar` file**, `consensus_from_counts()`
produces the **bit-identical** motif string in both. This means:

- Hamming distance comparisons (mismatch mode) use the same reference string.
- Every `(chrom, start, end, strand, matched_seq)` tuple produced by the
  pipeline is guaranteed to also appear in the genome-wide scan output.
- Peak hits are a **strict subset** of genome-wide hits — the fundamental
  correctness property of a ChIP-seq peak enrichment analysis.

### What happens without `--jaspar`

When `--jaspar` is absent, both scripts fall back to the `--query` / `--motif`
CLI string. If those strings happen to differ by even one character, the subset
property silently breaks — peaks hits may reference a motif position that the
genome scanner never evaluated. This is why `--jaspar` is strongly recommended
for mismatch and exact modes, and mandatory for PWM mode.

### Cross-validation warning

If `--jaspar` and `--query` are both supplied and they disagree, the script
prints:

```
[WARNING] --query "AGGGAATTTCCC" differs from JASPAR consensus "AGGGGATTTCCC".
          Using JASPAR consensus as the authoritative motif.
```

The JASPAR consensus always wins. The old default query string `AGGGAATTTCCC`
(with a missing `G` at position 5) was a typo — the corrected default is
`AGGGGATTTCCC`, which matches the attached `.jaspar` file exactly.
