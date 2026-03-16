# Sequence Motif Pipeline

A Nextflow DSL2 pipeline for genome-wide sequence motif discovery and clustering. Given a set of genomic coordinates (e.g., ChIP-seq or ATAC-seq peaks), the pipeline fetches DNA sequences from the mouse mm9 reference genome, searches them for a query motif with configurable mismatch tolerance, generates positional weight matrices (PWMs), sequence logos, and MEME-format motif files, and optionally clusters the flanking sequences using alignment-based distances, UMAP embedding, and K-Means.

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Workflow](#pipeline-workflow)
- [Requirements](#requirements)
- [Installation](#installation)
- [Directory Structure](#directory-structure)
- [Input Files](#input-files)
- [Usage](#usage)
- [Parameters](#parameters)
- [Output Files](#output-files)
- [Scripts](#scripts)
  - [fetch_sequence.py](#fetch_sequencepy)
  - [sequence_matcher_with_motifs.py](#sequence_matcher_with_motifspy)
  - [cluster_sequences.py](#cluster_sequencespy)
- [Docker Image](#docker-image)
- [Troubleshooting](#troubleshooting)
- [Dependencies](#dependencies)

---

## Overview

This pipeline automates three steps:

1. **Sequence Retrieval** — extracts the nucleotide sequence for each genomic region in your input CSV from the mm9 reference genome FASTA, handling strand orientation via reverse complementing.
2. **Motif Matching & Characterization** — searches each retrieved sequence for a user-defined query motif (with allowed mismatches), extracts flanking context around every hit, and produces PWMs, sequence logo PNGs, PSSM tables, and MEME motif files.
3. **Sequence Clustering** *(run independently, post-pipeline)* — clusters the upstream and downstream flanking sequences from the matches TSV using pairwise alignment distances, UMAP dimensionality reduction, and K-Means, producing a labelled results TSV and a 3×3 panel of dendrogram, heatmap, and UMAP scatter plots.

---

## Pipeline Workflow

```
Input CSVs  ──┐
              ├──▶  FETCH_SEQUENCES  ──▶  MATCH_AND_MOTIFS  ──▶  results/
mm9 FASTA   ──┘    (fetch_sequence.py)    (sequence_matcher_with_motifs.py)
                                                   │
                                                   ▼
                                         cluster_sequences.py   ──▶  clustering_results.tsv
                                         (run independently)          clustering_plots.png
```

**FETCH_SEQUENCES** (one task per sample CSV)
- Reads genomic coordinates from the CSV
- Uses `pyfaidx` to extract sequences from the mm9 FASTA by random-access index lookup
- Reverse-complements sequences on the `–` strand
- Outputs a TSV with a new `Fetched_sequence` column

**MATCH_AND_MOTIFS** (one task per FETCH_SEQUENCES output)
- Slides the query motif across every fetched sequence on both strands
- Accepts hits within the configured Hamming distance (mismatch) threshold
- Extracts upstream and downstream flanking sequences for every hit
- Builds a positional weight matrix (PWM) from the pooled flanks
- Writes sequence logo PNGs, PSSM TSV tables, and MEME-format motif files

**cluster_sequences.py** (run manually on any `*_matches.tsv`)
- Computes pairwise alignment distances across all unique upstream and downstream flank sequences
- Reduces dimensionality with UMAP; assigns cluster labels with K-Means
- Generates a 3×3 figure panel: dendrogram, distance heatmap, and UMAP scatter (one row per mode: upstream / downstream / combined)
- Writes a TSV with cluster labels merged back to every row in the matches file

---

## Requirements

| Requirement | Version | Notes |
|---|---|---|
| Nextflow | ≥ 22.10 | `nextflow -version` to check |
| Docker | ≥ 20.10 | Must be running; pipeline uses `seqmatcher:latest` |
| Java | ≥ 11 | Required by Nextflow |
| mm9 genome FASTA | — | See [Installation](#installation) |

> **Apple Silicon (M1/M2/M3) users:** the Docker image is built for `linux/amd64`. Docker Desktop on Apple Silicon runs it transparently via Rosetta emulation. You will see a platform mismatch warning — this is expected and does not affect results. Build the image with `--platform linux/amd64` as shown below.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/bioinfoanalyticsai/sequence-motif-pipeline.git
cd sequence-motif-pipeline
```

### 2. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
mv nextflow /usr/local/bin/   # or any directory on your PATH
```

### 3. Build the Docker image

```bash
docker build --platform linux/amd64 -t seqmatcher:latest .
```

Verify the scripts are baked in correctly:

```bash
docker run --rm seqmatcher:latest \
    grep -A 5 "def pwm_from_seqs" /usr/local/bin/scripts/sequence_matcher_with_motifs.py
```

### 4. Download the mm9 reference genome

```bash
mkdir -p input/genome
cd input/genome

# Download from UCSC (full genome, ~850 MB compressed)
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.fa.gz
gunzip mm9.fa.gz
mv mm9.fa mm9_genome.fa

# The .fai index is generated automatically on first run by pyfaidx.
# To generate it manually:
# pip install pyfaidx
# faidx mm9_genome.fa
```

---

## Directory Structure

```
sequence-motif-pipeline/
├── main.nf                              # Nextflow pipeline definition
├── Dockerfile                           # Docker image build instructions
├── nextflow.config                      # (optional) executor / resource config
├── input/
│   ├── genome/
│   │   └── mm9_genome.fa                # Reference genome (add manually — see above)
│   └── *.csv                            # Your sample input files (one per sample)
├── scripts/
│   ├── fetch_sequence.py                # Step 1: sequence retrieval
│   ├── sequence_matcher_with_motifs.py  # Step 2: motif matching + motif building
│   └── cluster_sequences.py            # Step 3: sequence clustering (run independently)
├── results/                             # Created automatically on pipeline run
└── README.md
```

---

## Input Files

### Sample CSV files

Place one CSV file per sample in `input/` (or the directory set by `--input_dir`). Each file must contain the following columns:

| Column | Type | Description |
|---|---|---|
| `GeneID` | string | Identifier for the region (gene name, peak ID, etc.) |
| `Chromosome` | string | Chromosome name in mm9 format, e.g. `chr1`, `chrX` |
| `Peak_Start` | integer | 0-based start coordinate (inclusive) |
| `Peak_End` | integer | 0-based end coordinate (exclusive) |
| `Strand` | string | `+` or `-` |

> **Coordinate system:** coordinates are 0-based half-open, consistent with BED format. If your peak caller outputs 1-based coordinates, subtract 1 from the start position.

**Example (`sample1.csv`):**

```
GeneID,Chromosome,Peak_Start,Peak_End,Strand
Myc_peak_001,chr8,128746210,128746410,+
Sox2_peak_042,chr3,34891002,34891202,-
Nanog_peak_017,chr6,122334510,122334710,+
```

### Reference genome

The pipeline expects `input/genome/mm9_genome.fa` (configurable via `--genome`). See [Installation](#installation) for download instructions. The `.fai` index is created automatically on first use.

---

## Usage

### Basic run

```bash
nextflow run main.nf
```

### With custom parameters

```bash
nextflow run main.nf \
    --input_dir  my_samples/ \
    --genome     /data/genomes/mm9/mm9_genome.fa \
    --query      TGACTCAGC \
    --mismatches 1 \
    --upstream   100 \
    --downstream 100 \
    --outdir     results_run1/
```

### Resume a failed run (uses Nextflow cache)

```bash
nextflow run main.nf -resume
```

### Clean run (discard all cached work)

```bash
rm -rf work .nextflow .nextflow.log*
nextflow run main.nf
```

> Always do a clean run after rebuilding the Docker image. Nextflow caches task executions by content hash; if the image changes but the input files do not, Nextflow may replay the old (cached) execution from inside the old container.

---

## Parameters

All parameters have defaults set in `main.nf` and can be overridden on the command line with `--param_name value`.

| Parameter | Default | Description |
|---|---|---|
| `--input_dir` | `input` | Directory containing sample CSV files. All `*.csv` files in this directory are processed. |
| `--genome` | `input/genome/mm9_genome.fa` | Path to the mm9 reference genome FASTA file. |
| `--query` | `ACTGATCGATCG` | Query motif sequence (A/T/C/G string). RNA sequences (U bases) are automatically converted to DNA. |
| `--mismatches` | `2` | Maximum number of base mismatches (Hamming distance) allowed when matching the query. `0` = exact matches only. |
| `--upstream` | `50` | Number of bases to extract upstream (5') of each motif hit for flank analysis. |
| `--downstream` | `50` | Number of bases to extract downstream (3') of each motif hit for flank analysis. |
| `--outdir` | `results` | Directory where all output files are published. Created if it does not exist. |
| `--use_docker` | `true` | Whether to run processes inside the Docker container. Set to `false` only if all Python dependencies are installed in the host environment. |

---

## Output Files

All outputs are written to `--outdir` (default: `results/`). For a sample CSV named `rel1h.csv`, Step 1 produces `rel1h_with_seq.tsv`, and Step 2 produces all remaining files using that TSV as the base name.

```
results/
├── rel1h_with_seq.tsv           # Sequences added; input to Step 2
├── rel1h_with_seq_matches.tsv   # All motif hits with full annotation
├── +_up_logo.png                # Sequence logo: upstream context, forward-strand hits
├── +_down_logo.png              # Sequence logo: downstream context, forward-strand hits
├── -_up_logo.png                # Sequence logo: upstream context, reverse-strand hits
├── -_down_logo.png              # Sequence logo: downstream context, reverse-strand hits
├── +_up_pssm.tsv                # PWM table: forward-strand upstream
├── +_down_pssm.tsv              # PWM table: forward-strand downstream
├── -_up_pssm.tsv                # PWM table: reverse-strand upstream
├── -_down_pssm.tsv              # PWM table: reverse-strand downstream
├── +_up.meme                    # MEME-format motif: forward-strand upstream
├── +_down.meme                  # MEME-format motif: forward-strand downstream
├── -_up.meme                    # MEME-format motif: reverse-strand upstream
└── -_down.meme                  # MEME-format motif: reverse-strand downstream
```

### Matches TSV columns

| Column | Description |
|---|---|
| *(all input columns)* | Preserved from the input CSV |
| `Strand` | Strand on which the hit was found (`+` or `-`) |
| `Matched_Sequence` | The actual sequence matched (may differ from query by up to `--mismatches` bases) |
| `Mismatches` | Number of base differences from the query (Hamming distance) |
| `Position` | 0-based position of the hit within the fetched sequence |
| `Alignment_Score` | Biopython PairwiseAligner score for the match |
| `Upstream` | The `--upstream` bases immediately 5' of the hit (padded with `N` if near sequence edge) |
| `Downstream` | The `--downstream` bases immediately 3' of the hit (padded with `N` if near sequence edge) |

### MEME files

The `.meme` files follow the MEME version 4 format and can be used directly with any tool in the [MEME Suite](https://meme-suite.org):

- **TOMTOM** — compare your flanking motif against databases of known transcription factor binding sites
- **FIMO** — scan additional sequences for occurrences of the flanking motif
- **AME** — test for enrichment of the flanking motif in a set of sequences

### Clustering outputs

Running `cluster_sequences.py` on a `*_matches.tsv` file produces two files (prefix is set by `--output-prefix`, default: `clustering`):

| File | Description |
|---|---|
| `clustering_results.tsv` | The full matches TSV with three new columns appended: `kmeans_upstream`, `kmeans_downstream`, and `kmeans_combined` (integer cluster labels). Also includes six UMAP coordinate columns (`umap_upstream_1/2`, `umap_downstream_1/2`, `umap_combined_1/2`) for the unique sequence pairs. |
| `clustering_plots.png` | 3×3 figure panel at 150 dpi. Rows = upstream / downstream / combined modes. Columns = dendrogram (coloured by K-Means cluster) / pairwise distance heatmap / UMAP scatter plot. |

---

## Scripts

### `fetch_sequence.py`

**Purpose:** Retrieve DNA sequences from the mm9 reference genome for each row in the input CSV.

**Location in container:** `/usr/local/bin/scripts/fetch_sequence.py`

**Usage:**
```bash
python fetch_sequence.py \
    -g <genome.fa> \
    -i <input.csv> \
    -o <output.tsv>
```

**Arguments:**

| Flag | Description |
|---|---|
| `-g` / `--genome` | Path to the mm9 FASTA file (must have a `.fai` index alongside it) |
| `-i` / `--input` | Input CSV with `GeneID`, `Chromosome`, `Peak_Start`, `Peak_End`, `Strand` columns |
| `-o` / `--output` | Output TSV path |

**Key logic:**
- Uses `pyfaidx.Fasta` for indexed random-access lookups: `genome[chrom][start:end].seq`
- Applies `reverse_complement()` for `Strand == '-'` rows using `str.maketrans('ATCGatcg', 'TAGCtagc')`
- Writes output as tab-separated with `pandas.DataFrame.to_csv(sep='\t')`

---

### `sequence_matcher_with_motifs.py`

**Purpose:** Search fetched sequences for the query motif and generate PWMs, logos, and MEME files.

**Location in container:** `/usr/local/bin/scripts/sequence_matcher_with_motifs.py`

**Usage:**
```bash
python sequence_matcher_with_motifs.py \
    -i <input_with_seq.tsv> \
    -o <output_matches.tsv> \
    -q <QUERY_MOTIF> \
    -m <max_mismatches> \
    --upstream <N> \
    --downstream <N>
```

**Arguments:**

| Flag | Description |
|---|---|
| `-i` / `--input` | TSV file produced by `fetch_sequence.py` (last column = fetched sequence) |
| `-o` / `--output` | Output TSV for all hits |
| `-q` / `--query` | Query motif string (DNA; U → T conversion applied automatically) |
| `-m` / `--mismatches` | Maximum Hamming distance to accept a hit (default: `0`) |
| `--upstream` | Bases to extract upstream of each hit (default: `50`) |
| `--downstream` | Bases to extract downstream of each hit (default: `50`) |

**Key logic:**

- `normalize(seq)` — uppercases and converts U→T
- `rc(seq)` — reverse complement via `Bio.Seq.reverse_complement()`
- `hamming(a, b)` — counts mismatches between equal-length strings
- `align_match(seq, query, max_mm)` — sliding window: scores each window with `Bio.Align.PairwiseAligner` (global mode; match=+1, mismatch=–1, gap open=–2, gap extend=–1); retains windows where Hamming distance ≤ `max_mm`
- `extract_flanks(seq, pos, k, up, down)` — extracts flanks; pads with `N` when hits are near sequence edges to guarantee consistent length for matrix construction
- `pwm_from_seqs(seqs)` — counts base frequencies per position across all flanks; filters to the most common length to handle any residual length variation; normalizes to probabilities
- `info_content(pwm)` — computes `pwm * log2(pwm / 0.25)` (bits)
- `logo(pwm, outfile, title)` — renders logo at 300 dpi via `logomaker`
- `write_meme(pwm, name, outfile)` — writes MEME v4 letter-probability matrix format

---

### `cluster_sequences.py`

**Purpose:** Cluster the upstream and downstream flanking sequences from the motif matches TSV using alignment-based distances, UMAP embedding, and K-Means. Run this script independently after the Nextflow pipeline completes.

**Location:** `scripts/cluster_sequences.py` (run on the host or inside any environment with the required dependencies — see [Dependencies](#dependencies))

**Usage:**
```bash
python scripts/cluster_sequences.py \
    results/<sample>_matches.tsv \
    --output-prefix  clustering \
    --n-clusters     3 \
    --alignment-mode global \
    --upstream       50 \
    --downstream     50
```

**Positional argument:**

| Argument | Description |
|---|---|
| `input` | Path to the `*_matches.tsv` file produced by `sequence_matcher_with_motifs.py` |

**Optional arguments:**

| Flag | Default | Description |
|---|---|---|
| `--output-prefix` | `clustering` | Prefix for all output files (`<prefix>_results.tsv`, `<prefix>_plots.png`) |
| `--n-clusters` | `3` | Number of K-Means clusters *k*. Adjust based on the expected biological groupings in your data. |
| `--kmeans-seed` | `42` | Random seed for K-Means reproducibility |
| `--linkage` | `average` | Linkage method for the hierarchical dendrogram: `single`, `complete`, `average`, or `ward` |
| `--alignment-mode` | `global` | `global` = Needleman-Wunsch (end-to-end alignment); `local` = Smith-Waterman (best local match) |
| `--w-upstream` | `0.5` | Weight given to upstream distances in combined mode (must sum to 1.0 with `--w-downstream`) |
| `--w-downstream` | `0.5` | Weight given to downstream distances in combined mode |
| `--umap-n-neighbors` | `5` | UMAP neighbourhood size — lower values capture more local structure |
| `--umap-min-dist` | `0.3` | UMAP minimum distance between embedded points — smaller = tighter clusters |
| `--umap-seed` | `42` | Random seed for UMAP reproducibility |
| `--no-plots` | off | Pass this flag to skip figure generation (results TSV is still written) |

**Key logic:**

- `load_input()` — reads the matches TSV, resolves column names case-insensitively, strips and uppercases sequences, drops rows with missing flanks, and assigns a `pair_id` to each unique `(Upstream, Downstream)` combination. Multiple GeneIDs that share identical flank pairs are clustered together and the label is propagated back to all rows in the output.
- `pairwise_distance_matrix()` — computes an N×N distance matrix for a list of sequences using `Bio.Align.PairwiseAligner`. Distance is defined as `1 − score(i,j) / sqrt(score(i,i) × score(j,j))`, which normalises raw alignment scores to the [0, 1] range.
- `combined_distance()` — produces a weighted average of the upstream and downstream distance matrices controlled by `--w-upstream` and `--w-downstream`.
- `run_umap()` — embeds the precomputed distance matrix into 2D using `umap-learn` with `metric='precomputed'`. `n_neighbors` is automatically capped at N−1 when the dataset is small.
- `kmeans_cluster()` — runs `sklearn.cluster.KMeans` on the standardised UMAP coordinates and returns integer cluster labels.
- `hierarchical_linkage()` — converts the distance matrix to a condensed form with `scipy.spatial.distance.squareform` and computes a linkage matrix with `scipy.cluster.hierarchy.linkage` for dendrogram plotting only.
- Three clustering modes are always computed and plotted in parallel: **upstream only**, **downstream only**, and **combined** (weighted average).

**Output figure layout (3 rows × 3 columns):**

| | Column 0: Dendrogram | Column 1: Distance Heatmap | Column 2: UMAP Scatter |
|---|---|---|---|
| **Row 0** | Upstream sequences | Upstream sequences | Upstream sequences |
| **Row 1** | Downstream sequences | Downstream sequences | Downstream sequences |
| **Row 2** | Combined | Combined | Combined |

Leaves/cells/points are coloured by K-Means cluster label. Heatmap diagonal squares are also coloured by cluster for quick cross-reference.

The `Dockerfile` builds an `linux/amd64` image that installs all Python dependencies and copies both scripts to `/usr/local/bin/scripts/`.

**Build:**
```bash
docker build --platform linux/amd64 -t seqmatcher:latest .
```

**Verify scripts inside container:**
```bash
# Check all three scripts are present
docker run --rm seqmatcher:latest ls /usr/local/bin/scripts/

# Check Python dependencies
docker run --rm seqmatcher:latest python -c \
    "import pyfaidx, Bio, logomaker, pandas, numpy; print('All dependencies OK')"
```

> `cluster_sequences.py` uses additional libraries (`umap-learn`, `scikit-learn`, `scipy`, `seaborn`) that may not be in the Docker image. Install them in your host environment or add them to the Dockerfile if you want to run clustering inside the container.

**After any script change, always rebuild before rerunning the pipeline:**
```bash
docker build --platform linux/amd64 -t seqmatcher:latest .
rm -rf work .nextflow .nextflow.log*
nextflow run main.nf
```

---

## Troubleshooting

### `FastaNotFoundError: Cannot read FASTA from file input/genome/mm9_genome.fa`

The genome path is resolved inside the Nextflow work directory (an isolated subdirectory of `work/`), where relative paths do not exist. The genome must be declared as a `path` input and passed via `Channel.fromPath()` so Nextflow stages it. In `main.nf`:

```groovy
genome_ch = Channel.fromPath(params.genome)
FETCH_SEQUENCES(samples_ch.combine(genome_ch))
```

And the process input must be:
```groovy
input:
tuple path(sample_file), path(genome_file)
```

---

### `IndexError: index N is out of bounds for axis 0 with size N`

Occurs in `pwm_from_seqs()` when flanking sequences have inconsistent lengths because a motif hit was found near the start or end of a fetched sequence. Fixed by:

1. Padding short flanks with `N` in `extract_flanks()` (left-pad upstream, right-pad downstream)
2. Filtering to the most common sequence length in `pwm_from_seqs()` before building the matrix

Both fixes are present in the current version of `sequence_matcher_with_motifs.py`. If this error appears, the Docker image likely contains an older version of the script — rebuild the image.

---

### Error persists after fixing a script

Nextflow caches task executions. If you fix a script and rebuild the image but do not clear the cache, Nextflow may re-execute the task from a cached work directory that still contains the old container run. Always run:

```bash
rm -rf work .nextflow .nextflow.log*
nextflow run main.nf
```

---

### `WARNING: The requested image's platform (linux/amd64) does not match the detected host platform (linux/arm64/v8)`

This warning appears on Apple Silicon Macs (M1/M2/M3) running an `amd64` Docker image. Docker Desktop handles the emulation automatically via Rosetta. The warning is informational only and does **not** affect the correctness of results. Build with `--platform linux/amd64` to make this explicit:

```bash
docker build --platform linux/amd64 -t seqmatcher:latest .
```

---

### No hits in output / empty matches TSV

- Verify the query sequence matches the alphabet of your input sequences (DNA, not RNA).
- Check that `--mismatches` is not set too low. Try increasing to `2` or `3`.
- Confirm the `Chromosome` column values in your CSV match the chromosome naming in the mm9 FASTA (e.g., `chr1` not `1`).
- Confirm `Peak_Start` and `Peak_End` are 0-based and within the chromosome length.

---

## Dependencies

### Python libraries — pipeline (inside Docker container)

| Library | Purpose |
|---|---|
| `pyfaidx` | Indexed random-access sequence retrieval from FASTA files |
| `biopython` | Reverse complement (`Bio.Seq`); pairwise alignment scoring (`Bio.Align.PairwiseAligner`) |
| `pandas` | Reading/writing CSV and TSV files; DataFrame operations |
| `numpy` | Position frequency matrix construction and normalization |
| `logomaker` | Rendering sequence logo PNG figures |
| `matplotlib` | Figure backend for logomaker; axis labeling and PNG export |

### Python libraries — clustering (`cluster_sequences.py`, host environment)

Install with:
```bash
pip install biopython scikit-learn scipy matplotlib seaborn pandas numpy umap-learn
```

| Library | Purpose |
|---|---|
| `biopython` | Pairwise sequence alignment (`Bio.Align.PairwiseAligner`) for distance matrix computation |
| `numpy` | Distance matrix operations and numerical arrays |
| `pandas` | TSV I/O; pair-ID assignment; cluster label merging |
| `scipy` | Condensed distance matrix conversion (`squareform`); hierarchical linkage and dendrogram |
| `scikit-learn` | K-Means clustering (`sklearn.cluster.KMeans`); feature scaling (`StandardScaler`) |
| `umap-learn` | UMAP 2D embedding from a precomputed distance matrix |
| `matplotlib` | Figure creation; subplot layout; PNG export |
| `seaborn` | Annotated heatmap rendering (`sns.heatmap`) |

### Infrastructure

| Tool | Purpose |
|---|---|
| Nextflow DSL2 | Workflow orchestration; parallel task execution; data channel management |
| Docker | Containerized, reproducible execution environment |
| Java ≥ 11 | Required runtime for Nextflow |
