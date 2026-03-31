#!/usr/bin/env python3
"""
Biological Sequence Clustering Script
======================================
Designed for TSV files with the following header structure:
  GeneID | Chromosome | Peak_Start | Peak_End | Strand | ... |
  Fetched_sequence | Matched_Sequence | Mismatches | Position |
  Alignment_Score | Upstream | Downstream

The two sequence columns used for clustering are:
  • Upstream    – sequence upstream of the matched motif   (col 1)
  • Downstream  – sequence downstream of the matched motif (col 2)

Multiple GeneIDs can share the same (Upstream, Downstream) pair
(e.g. overlapping genomic peaks). The script clusters on UNIQUE
sequence pairs and maps labels back to every GeneID row in the output.

Clustering strategy
-------------------
  Step 1 – Pairwise alignment distances (Needleman-Wunsch or Smith-Waterman)
  Step 2 – Hierarchical/agglomerative clustering → Dendrogram
  Step 3 – K-Means on UMAP-embedded coordinates → final cluster labels
  Step 4 – Per-cluster sequence logos (upstream, downstream, combined)

  Three modes:
    1. Upstream sequences only
    2. Downstream sequences only
    3. Both columns combined (weighted-average distance)

Output plots (3 rows × 3 columns)
-----------------------------------
  Col 0 – Dendrogram  (hierarchical, coloured by K-Means cluster)
  Col 1 – Pairwise distance heatmap (annotated with K-Means cluster)
  Col 2 – UMAP 2-D scatter (coloured by K-Means cluster)

Sequence logo outputs
---------------------
  For every clustering mode × every cluster, four logo PNGs are written:
    <prefix>_logo_<mode>_cluster<N>_upstream.png
    <prefix>_logo_<mode>_cluster<N>_downstream.png

  A summary mosaic PNG is also written:
    <prefix>_logos_<mode>.png
  showing all clusters side-by-side (upstream on top row, downstream below).

Dependencies
------------
  pip install biopython scikit-learn scipy matplotlib seaborn pandas numpy
              umap-learn logomaker
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")           # non-interactive backend, safe for headless servers
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import logomaker                # sequence-logo rendering
from itertools import combinations
from Bio import Align
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import umap


# ── Alphabet used for PWM / logo construction ──────────────────────────────────
BASES = ["A", "C", "G", "T"]

# ── Default column names (case-insensitive matching applied at load time) ──────
_ID_COL   = "GeneID"
_SEQ1_COL = "Upstream"
_SEQ2_COL = "Downstream"


# ══════════════════════════════════════════════════════════════════════════════
# 1.  INPUT LOADING & VALIDATION
# ══════════════════════════════════════════════════════════════════════════════

def _resolve_col(df: pd.DataFrame, name: str) -> str:
    """
    Return the actual column name in *df* that matches *name*.
    Tries exact match first, then falls back to case-insensitive lookup.
    Raises ValueError if the column is absent entirely.
    """
    if name in df.columns:
        return name
    col_map = {c.lower(): c for c in df.columns}
    if name.lower() in col_map:
        return col_map[name.lower()]
    raise ValueError(
        f"Required column '{name}' not found.\n"
        f"Available columns: {list(df.columns)}"
    )


def load_input(path: str):
    """
    Load the TSV produced by sequence_matcher_with_motifs.py, resolve
    column names, strip/uppercase sequences, and assign pair IDs.

    Why pair IDs?
    -------------
    Multiple GeneID rows can share an identical (Upstream, Downstream)
    pair when genomic peaks overlap.  We cluster on UNIQUE pairs so that
    each sequence pattern is represented once, then map the cluster label
    back to every matching GeneID row in the output TSV.

    Returns
    -------
    full_df  : original DataFrame enriched with a 'pair_id' column
    uniq_df  : one row per unique (Upstream, Downstream) pair
    id_col   : resolved GeneID column name
    seq1_col : resolved Upstream column name
    seq2_col : resolved Downstream column name
    gene_map : dict  pair_id -> list[GeneID]  (for axis labels)
    """
    # Read everything as strings to avoid dtype surprises in sequence columns
    raw = pd.read_csv(path, sep="\t", dtype=str)
    raw.columns = raw.columns.str.strip()   # remove any accidental whitespace

    # Resolve the three required columns with case-insensitive fallback
    id_col   = _resolve_col(raw, _ID_COL)
    seq1_col = _resolve_col(raw, _SEQ1_COL)
    seq2_col = _resolve_col(raw, _SEQ2_COL)

    # ── Clean: drop rows where either flank sequence is missing or blank ──────
    before = len(raw)
    raw = raw.dropna(subset=[seq1_col, seq2_col])
    raw = raw[raw[seq1_col].str.strip() != ""]
    raw = raw[raw[seq2_col].str.strip() != ""]
    dropped = before - len(raw)
    if dropped:
        print(f"  [warn] Dropped {dropped} rows with missing/empty sequences.")

    # Standardise sequences: strip whitespace, force uppercase DNA
    raw[seq1_col] = raw[seq1_col].str.strip().str.upper()
    raw[seq2_col] = raw[seq2_col].str.strip().str.upper()

    print(f"Loaded  : {len(raw)} rows  |  {raw[id_col].nunique()} unique GeneIDs")

    # ── Assign pair IDs – one label per unique (Upstream, Downstream) pair ───
    # Concatenate the two sequences with a separator that cannot appear in DNA
    raw       = raw.copy()
    pair_key  = raw[seq1_col] + "||" + raw[seq2_col]
    uniq_pairs = pair_key.unique()
    # Zero-padded numeric suffix for stable, sortable labels (pair_001, pair_002 …)
    pair_id_map = {p: f"pair_{i+1:03d}" for i, p in enumerate(uniq_pairs)}
    raw["pair_id"] = pair_key.map(pair_id_map)

    # Build a reverse look-up: pair_id → [GeneID, GeneID, …]
    # Used to annotate axes with gene names rather than opaque pair IDs
    gene_map = raw.groupby("pair_id")[id_col].apply(list).to_dict()

    # Deduplicate to one row per unique pair for distance/clustering calculations
    uniq_df = (
        raw[[seq1_col, seq2_col, "pair_id"]]
        .drop_duplicates(subset="pair_id")
        .reset_index(drop=True)
    )

    print(f"Unique (Upstream, Downstream) pairs: {len(uniq_df)}")

    # Need at least 2 sequences to compute a pairwise distance
    if len(uniq_df) < 2:
        sys.exit("Error: need at least 2 unique sequence pairs to cluster.")

    return raw, uniq_df, id_col, seq1_col, seq2_col, gene_map


# ══════════════════════════════════════════════════════════════════════════════
# 2.  PAIRWISE ALIGNMENT DISTANCE
# ══════════════════════════════════════════════════════════════════════════════

def pairwise_distance_matrix(sequences: list, mode: str = "global") -> np.ndarray:
    """
    Build an N×N pairwise distance matrix using Biopython PairwiseAligner.

    Distance formula
    ----------------
    Raw alignment scores are not directly comparable across sequence pairs
    of different lengths or compositions.  We normalise to a [0, 1] range:

        similarity(i, j) = score(i, j) / sqrt( score(i,i) × score(j,j) )
        distance(i, j)   = 1 − similarity(i, j)

    This makes distance(i, i) = 0 and distance(i, j) = 1 for completely
    dissimilar pairs, irrespective of sequence length.

    Parameters
    ----------
    mode : 'global'  → Needleman-Wunsch (end-to-end alignment, good for
                        flanks of the same fixed length)
           'local'   → Smith-Waterman (best local sub-alignment, better
                        when flanks differ substantially in length)
    """
    aligner      = Align.PairwiseAligner()
    aligner.mode = mode   # 'global' or 'local'

    n            = len(sequences)
    dist_matrix  = np.zeros((n, n))   # initialise as all-zeros (diagonal stays 0)

    # Self-alignment scores used as the normalisation denominator
    self_scores = np.array([aligner.score(s, s) for s in sequences])

    # Only compute the upper triangle; mirror to lower triangle (matrix is symmetric)
    for i, j in combinations(range(n), 2):
        raw   = aligner.score(sequences[i], sequences[j])
        denom = np.sqrt(self_scores[i] * self_scores[j])

        # Guard against zero-length or all-gap sequences producing denom = 0
        sim = float(raw) / denom if denom > 0 else 0.0
        sim = min(max(sim, 0.0), 1.0)   # clamp to [0, 1] for floating-point safety
        d   = 1.0 - sim

        dist_matrix[i, j] = d   # upper triangle
        dist_matrix[j, i] = d   # lower triangle (symmetric)

    return dist_matrix


def combined_distance(dist1: np.ndarray, dist2: np.ndarray,
                      w1: float = 0.5, w2: float = 0.5) -> np.ndarray:
    """
    Weighted average of two distance matrices (upstream + downstream).

    Normalising by (w1 + w2) keeps the result in [0, 1] regardless of
    the chosen weights, so the combined matrix is comparable to the
    individual upstream and downstream matrices.
    """
    return (w1 * dist1 + w2 * dist2) / (w1 + w2)


# ══════════════════════════════════════════════════════════════════════════════
# 3.  HIERARCHICAL CLUSTERING  (used only for the dendrogram visualisation)
# ══════════════════════════════════════════════════════════════════════════════

def hierarchical_linkage(dist_matrix: np.ndarray,
                         linkage_method: str = "average") -> np.ndarray:
    """
    Compute a scipy linkage matrix Z from a precomputed square distance matrix.

    Why hierarchical clustering here?
    ----------------------------------
    Hierarchical / agglomerative clustering produces the dendrogram that
    shows how sequences merge into groups step by step.  It is used purely
    for visualisation.  The actual cluster labels used in the output TSV
    come from K-Means (Step 5), which is more stable for non-spherical
    cluster shapes revealed by UMAP.

    linkage_method choices
    ----------------------
    'single'   – minimum pairwise distance between any two points in the clusters
    'complete' – maximum pairwise distance (more compact, sensitive to outliers)
    'average'  – mean pairwise distance (UPGMA, balanced)
    'ward'     – minimises total within-cluster variance (good for equal-sized clusters)
    """
    # Ensure diagonal is exactly zero before converting to condensed form
    np.fill_diagonal(dist_matrix, 0.0)

    # squareform converts the NxN matrix to the N*(N-1)/2 condensed upper-triangle
    # vector that scipy's linkage() expects
    condensed = squareform(dist_matrix, checks=False)
    return linkage(condensed, method=linkage_method)


# ══════════════════════════════════════════════════════════════════════════════
# 4.  UMAP EMBEDDING
# ══════════════════════════════════════════════════════════════════════════════

def run_umap(dist_matrix: np.ndarray,
             n_neighbors: int = 5,
             min_dist: float = 0.3,
             random_state: int = 42) -> np.ndarray:
    """
    Reduce the N×N distance matrix to 2-D coordinates using UMAP.

    Why UMAP?
    ---------
    The pairwise distance matrix lives in high-dimensional space.  UMAP
    (Uniform Manifold Approximation and Projection) learns the underlying
    manifold structure and projects it onto 2 axes that preserve both local
    neighbourhood relationships and broader topology.  This 2-D embedding is
    then used as the input space for K-Means, which works far better in
    low-dimensional space than on raw distance matrices.

    Parameters
    ----------
    n_neighbors : int
        Controls the balance between local and global structure.
        Small values (2–5) reveal fine local clusters; larger values
        (15–50) capture more global topology.
        Automatically capped at n_samples − 1 for small datasets.
    min_dist : float
        How tightly UMAP packs points together in the embedding.
        0.0 = very tight clusters; 0.5+ = more spread-out layout.
    random_state : int
        Seed for reproducibility (UMAP has stochastic elements).
    """
    n  = dist_matrix.shape[0]
    nn = min(n_neighbors, n - 1)   # cap so UMAP doesn't crash on small datasets
    if nn != n_neighbors:
        print(f"  [info] UMAP n_neighbors capped at {nn} (only {n} unique pairs).")

    reducer = umap.UMAP(
        n_components=2,           # project to 2 dimensions for scatter plotting
        metric="precomputed",     # tell UMAP we're supplying a distance matrix
        n_neighbors=nn,
        min_dist=min_dist,
        random_state=random_state,
    )

    # Ensure the diagonal is exactly zero (required by UMAP's precomputed mode)
    np.fill_diagonal(dist_matrix, 0.0)
    return reducer.fit_transform(dist_matrix)   # returns (n, 2) array


# ══════════════════════════════════════════════════════════════════════════════
# 5.  K-MEANS CLUSTERING  (on UMAP coordinates)
# ══════════════════════════════════════════════════════════════════════════════

def kmeans_cluster(umap_coords: np.ndarray,
                   n_clusters: int,
                   random_state: int = 42) -> np.ndarray:
    """
    Assign K-Means cluster labels to the UMAP-embedded points.

    Why K-Means on UMAP rather than directly on the distance matrix?
    ----------------------------------------------------------------
    K-Means assumes Euclidean distances and spherical clusters.  Running it
    directly on a high-dimensional pairwise distance matrix often gives poor
    results.  Instead, UMAP first maps the data to a 2-D Euclidean space
    that respects the underlying sequence similarity structure, and K-Means
    then partitions that space cleanly.

    Returns integer labels that are 1-based (cluster 1, 2, … k) for
    consistency with scipy's fcluster output used in the dendrogram.
    """
    n = umap_coords.shape[0]
    k = min(n_clusters, n)   # can't have more clusters than data points
    if k != n_clusters:
        print(f"  [info] K-Means k capped at {k} (only {n} unique pairs).")

    # n_init=20: run K-Means 20 times with different centroid seeds and keep
    # the best result (minimises sensitivity to random initialisation)
    km     = KMeans(n_clusters=k, random_state=random_state, n_init=20)
    labels = km.fit_predict(umap_coords)
    return labels + 1   # shift from 0-based to 1-based cluster numbering


# ══════════════════════════════════════════════════════════════════════════════
# 6.  PWM / SEQUENCE LOGO CONSTRUCTION
# ══════════════════════════════════════════════════════════════════════════════

def _filter_seqs_by_length(seqs: list) -> list:
    """
    Return only sequences that match the most common length in *seqs*.

    Short sequences arise when a motif hit was found near the edge of a
    fetched region and the flank window was partially out of bounds.  They
    would cause an index error when building the PWM matrix.  Filtering to
    the modal length keeps the matrix rectangular without discarding most
    of the data.
    """
    if not seqs:
        return []
    lengths     = [len(s) for s in seqs]
    modal_len   = max(set(lengths), key=lengths.count)
    # Only keep full-length sequences; partial flanks are silently dropped
    return [s for s in seqs if len(s) == modal_len]


def build_pwm(seqs: list) -> pd.DataFrame:
    """
    Build a position probability matrix (PPM / PWM) from a list of
    equal-length DNA sequences.

    The PWM is a (L × 4) DataFrame where:
      - rows  = positions 0 … L-1 along the sequence
      - columns = A, C, G, T
      - values = fraction of sequences that carry each base at each position

    N-padded bases (used to fill short edge flanks) are ignored because
    'N' is not in BASES, so they do not contribute to any column count.
    A pseudocount of 0 is used here; information content calculations
    add a small epsilon to avoid log(0).
    """
    seqs = _filter_seqs_by_length(seqs)
    if not seqs:
        # Return an empty DataFrame so callers can check `.empty` before plotting
        return pd.DataFrame(columns=BASES)

    L   = len(seqs[0])
    mat = np.zeros((L, 4))   # raw count matrix, shape (positions, bases)

    for s in seqs:
        for i, b in enumerate(s):
            if b in BASES:
                # Increment the count for this base at this position
                mat[i, BASES.index(b)] += 1

    # Normalise each row so counts → frequencies (sum to 1.0 per position)
    row_sums = mat.sum(axis=1, keepdims=True)
    # Avoid division by zero at positions that are all-N (pure padding)
    row_sums[row_sums == 0] = 1.0
    mat /= row_sums

    return pd.DataFrame(mat, columns=BASES)


def pwm_to_information_content(pwm: pd.DataFrame) -> pd.DataFrame:
    """
    Convert a position probability matrix to an information-content matrix.

    Information content (IC) at each position is:
        IC(i, b) = pwm(i, b) × log2( pwm(i, b) / 0.25 )

    The background probability for each base is 0.25 (uniform, i.e. all
    four bases equally likely by chance).  A position where one base
    dominates completely contributes up to 2 bits; a fully random position
    contributes 0 bits.  logomaker uses the IC matrix to size each letter.
    """
    # Small epsilon prevents log2(0); has negligible effect on non-zero values
    return pwm * np.log2((pwm + 1e-9) / 0.25)


def plot_logo_to_ax(ic_df: pd.DataFrame, ax, title: str = "",
                    ylabel: str = "Bits") -> None:
    """
    Draw a sequence logo onto an existing matplotlib Axes object.

    Parameters
    ----------
    ic_df  : information-content DataFrame from pwm_to_information_content()
    ax     : matplotlib Axes to draw into
    title  : subplot title string
    ylabel : y-axis label (default 'Bits')
    """
    if ic_df.empty:
        # No sequences passed quality filters – show a blank placeholder
        ax.text(0.5, 0.5, "No sequences\nin this cluster",
                ha="center", va="center", transform=ax.transAxes, fontsize=9,
                color="grey")
        ax.set_title(title, fontsize=9)
        ax.axis("off")
        return

    # logomaker.Logo draws letter stacks proportional to IC at each position
    logomaker.Logo(ic_df, ax=ax, color_scheme="classic")
    ax.set_title(title, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_xlabel("Position", fontsize=8)
    ax.tick_params(labelsize=7)


def generate_cluster_logos(uniq_df: pd.DataFrame,
                           km_labels_dict: dict,
                           seq1_col: str,
                           seq2_col: str,
                           output_prefix: str,
                           n_clusters: int) -> None:
    """
    Generate per-cluster sequence logos for upstream and downstream flanks.

    For each clustering mode (upstream / downstream / combined) and each
    cluster label, the sequences belonging to that cluster are collected,
    a PWM is built, converted to information content, and rendered as a
    sequence logo.

    Two output types are written:
    1. Individual PNG per cluster per flank type (for easy downstream use):
         <prefix>_logo_<mode>_cluster<N>_upstream.png
         <prefix>_logo_<mode>_cluster<N>_downstream.png

    2. Summary mosaic PNG with all clusters in one figure (easier to compare):
         <prefix>_logos_<mode>.png
       Layout: rows = upstream / downstream; columns = cluster 1, 2 … k

    Parameters
    ----------
    uniq_df        : DataFrame of unique sequence pairs (one row per pair_id)
    km_labels_dict : dict mode -> list of integer cluster labels (1-based)
    seq1_col       : resolved Upstream column name
    seq2_col       : resolved Downstream column name
    output_prefix  : file path prefix for all output files
    n_clusters     : number of K-Means clusters (used for figure sizing)
    """
    seqs_up   = uniq_df[seq1_col].tolist()
    seqs_down = uniq_df[seq2_col].tolist()

    # Iterate over the three clustering modes
    for mode, labels in km_labels_dict.items():

        labels_arr   = np.array(labels)          # convert to numpy for boolean masking
        unique_clust = sorted(set(labels_arr))   # sorted cluster IDs for consistent ordering

        print(f"\n  [{mode:10s}]  Generating per-cluster sequence logos ...")

        # ── 1. Summary mosaic figure (all clusters side-by-side) ─────────────
        # Layout: 2 rows (upstream, downstream) × k columns (one per cluster)
        n_cols  = len(unique_clust)
        fig_logo, axes_logo = plt.subplots(
            2, n_cols,
            figsize=(max(6, 4 * n_cols), 7),   # width scales with cluster count
            squeeze=False                        # always return 2-D axes array
        )

        fig_logo.suptitle(
            f"Per-Cluster Sequence Logos  |  Mode: {mode}",
            fontsize=12, fontweight="bold"
        )

        for col_idx, clust_id in enumerate(unique_clust):

            # Boolean mask selecting all sequence pairs assigned to this cluster
            mask = labels_arr == clust_id

            # Collect the upstream sequences that belong to this cluster
            cluster_seqs_up   = [seqs_up[i]   for i in range(len(seqs_up))   if mask[i]]
            # Collect the downstream sequences that belong to this cluster
            cluster_seqs_down = [seqs_down[i] for i in range(len(seqs_down)) if mask[i]]

            n_seqs = sum(mask)   # number of unique pairs in this cluster
            print(f"    Cluster {clust_id}: {n_seqs} sequences")

            # ── Build PWMs and convert to information content ─────────────────
            pwm_up   = build_pwm(cluster_seqs_up)
            pwm_down = build_pwm(cluster_seqs_down)
            ic_up    = pwm_to_information_content(pwm_up)
            ic_down  = pwm_to_information_content(pwm_down)

            # ── Draw into the mosaic grid ─────────────────────────────────────
            # Row 0 = upstream logo for this cluster
            plot_logo_to_ax(
                ic_up,
                axes_logo[0, col_idx],
                title=f"Cluster {clust_id} — Upstream\n(n={n_seqs})"
            )
            # Row 1 = downstream logo for this cluster
            plot_logo_to_ax(
                ic_down,
                axes_logo[1, col_idx],
                title=f"Cluster {clust_id} — Downstream\n(n={n_seqs})"
            )

            # ── Write individual per-cluster PNG files ────────────────────────
            # Upstream logo
            fig_ind, ax_ind = plt.subplots(figsize=(max(6, len(pwm_up) / 5), 3))
            plot_logo_to_ax(
                ic_up, ax_ind,
                title=f"Upstream Logo | Mode: {mode} | Cluster {clust_id} (n={n_seqs})"
            )
            plt.tight_layout()
            fname_up = f"{output_prefix}_logo_{mode}_cluster{clust_id}_upstream.png"
            fig_ind.savefig(fname_up, dpi=200, bbox_inches="tight")
            plt.close(fig_ind)
            print(f"      Saved: {fname_up}")

            # Downstream logo
            fig_ind, ax_ind = plt.subplots(figsize=(max(6, len(pwm_down) / 5), 3))
            plot_logo_to_ax(
                ic_down, ax_ind,
                title=f"Downstream Logo | Mode: {mode} | Cluster {clust_id} (n={n_seqs})"
            )
            plt.tight_layout()
            fname_down = f"{output_prefix}_logo_{mode}_cluster{clust_id}_downstream.png"
            fig_ind.savefig(fname_down, dpi=200, bbox_inches="tight")
            plt.close(fig_ind)
            print(f"      Saved: {fname_down}")

        # ── Finalise and save the mosaic for this mode ────────────────────────
        # Add row labels on the left side of the first column
        axes_logo[0, 0].set_ylabel("Upstream\n(Bits)", fontsize=9, labelpad=10)
        axes_logo[1, 0].set_ylabel("Downstream\n(Bits)", fontsize=9, labelpad=10)

        plt.tight_layout()
        mosaic_out = f"{output_prefix}_logos_{mode}.png"
        fig_logo.savefig(mosaic_out, dpi=150, bbox_inches="tight")
        plt.close(fig_logo)
        print(f"    Mosaic saved: {mosaic_out}")


# ══════════════════════════════════════════════════════════════════════════════
# 7.  VISUALISATION HELPERS  (clustering plots)
# ══════════════════════════════════════════════════════════════════════════════

# Shared colour palette – tab10 gives 10 visually distinct colours.
# Cluster labels are 1-based, so we subtract 1 for the 0-based palette index.
_CMAP = plt.get_cmap("tab10")

def _cluster_color(lbl: int) -> tuple:
    """Return the RGBA colour tuple for cluster label *lbl* (1-based)."""
    return _CMAP((lbl - 1) % 10)   # modulo 10 handles k > 10 gracefully


def _tick_labels(pair_ids: list, gene_map: dict, max_genes: int = 3) -> list:
    """
    Build readable axis tick labels of the form:
        pair_001
        (Gene1, Gene2, …)
    Truncates to the first *max_genes* gene names to keep labels short.
    """
    out = []
    for pid in pair_ids:
        genes    = gene_map.get(pid, [])
        gene_str = ", ".join(genes[:max_genes])
        if len(genes) > max_genes:
            gene_str += "..."
        out.append(f"{pid}\n({gene_str})")
    return out


def _legend_patches(unique_labels: list) -> list:
    """Build a list of coloured Patch objects for a matplotlib legend."""
    return [
        mpatches.Patch(color=_cluster_color(l), label=f"Cluster {l}")
        for l in sorted(unique_labels)
    ]


# ── Plot A: Dendrogram coloured by K-Means cluster ───────────────────────────

def plot_dendrogram(Z, tick_labels, km_labels, title, ax):
    """
    Draw a hierarchical clustering dendrogram with leaf labels coloured
    by their K-Means cluster assignment.

    A horizontal dashed red line marks the cut height that approximately
    separates the tree into k clusters, providing a visual reference for
    how well the hierarchical tree agrees with K-Means partitioning.
    """
    n = len(tick_labels)

    # Map each tick label string to its cluster colour for leaf colouring
    label_colors = {
        tick_labels[i]: _cluster_color(km_labels[i])
        for i in range(n)
    }

    # Draw the dendrogram with all branches in grey.
    # no_labels=True because we manually re-draw labels with per-cluster colours.
    ddata = dendrogram(
        Z,
        labels=tick_labels,
        ax=ax,
        leaf_rotation=90,
        leaf_font_size=7,
        color_threshold=0,            # disable scipy's own branch colouring
        above_threshold_color="grey",
        no_labels=True,
    )

    # Re-draw x-axis tick labels with cluster colours ─────────────────────────
    # scipy places leaf centres at 5, 15, 25 … (step=10) in dendrogram units
    ax.set_xticks(range(5, n * 10 + 5, 10))
    leaf_order = ddata["ivl"]   # left-to-right leaf label order from scipy
    for tick, lbl in zip(ax.get_xticks(), leaf_order):
        color = label_colors.get(lbl, "black")
        ax.text(tick, -0.02, lbl,
                ha="right", va="top",
                rotation=90, fontsize=6,
                color=color,
                transform=ax.get_xaxis_transform())

    # Draw a dashed horizontal cut line that separates the tree into k clusters
    # The cut height is the midpoint between the (k-1)th and kth merge heights
    k = len(set(km_labels))
    if k < len(Z):
        cut_height = (Z[-(k - 1), 2] + Z[-(k), 2]) / 2
        ax.axhline(cut_height, color="red", linestyle="--",
                   linewidth=1.0, label=f"k={k} cut")

    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Sequence pair  (GeneIDs)", fontsize=8)
    ax.set_ylabel("Distance", fontsize=8)

    # Legend showing cluster → colour mapping
    patches = _legend_patches(sorted(set(km_labels)))
    ax.legend(handles=patches, fontsize=7, loc="upper right", framealpha=0.7)


# ── Plot B: Distance heatmap annotated with cluster bar ─────────────────────

def plot_heatmap(dist_matrix, tick_labels, km_labels, title, ax):
    """
    Seaborn heatmap of the pairwise distance matrix.

    Each diagonal cell is overlaid with a semi-transparent coloured square
    indicating the K-Means cluster of that sequence pair, allowing quick
    visual cross-reference between the heatmap pattern and cluster membership.
    """
    df_hm = pd.DataFrame(dist_matrix,
                         index=tick_labels, columns=tick_labels)

    # viridis_r: low distance = bright yellow (similar), high = dark purple (dissimilar)
    sns.heatmap(df_hm, ax=ax, cmap="viridis_r", square=True,
                linewidths=0.2, cbar_kws={"label": "Distance"},
                xticklabels=True, yticklabels=True)

    # Overlay coloured squares on the diagonal cells to show cluster membership
    n = len(tick_labels)
    for i in range(n):
        c = _cluster_color(km_labels[i])
        ax.add_patch(plt.Rectangle(
            (i, i), 1, 1,           # cell at (row=i, col=i) in heatmap coordinates
            fill=True, color=c,
            alpha=0.45,             # semi-transparent so the distance value shows through
            linewidth=0
        ))

    ax.set_title(title, fontsize=10)
    ax.tick_params(axis="both", labelsize=6)

    patches = _legend_patches(sorted(set(km_labels)))
    ax.legend(handles=patches, fontsize=6, loc="lower right",
              framealpha=0.8, bbox_to_anchor=(1.0, 0.0))


# ── Plot C: UMAP scatter coloured by K-Means cluster ─────────────────────────

def plot_umap(umap_coords, km_labels, pair_ids, gene_map, title, ax):
    """
    2-D scatter plot of the UMAP embedding.

    Each point represents one unique (Upstream, Downstream) sequence pair.
    Points are coloured by K-Means cluster and annotated with their pair_id,
    making it easy to identify which genes belong to each cluster.
    """
    # Draw one scatter series per cluster so each gets its own legend entry
    for lbl in sorted(set(km_labels)):
        mask = np.array(km_labels) == lbl
        ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                   color=_cluster_color(lbl), s=100,
                   edgecolors="k", linewidths=0.5,
                   label=f"Cluster {lbl}", zorder=3)

    # Annotate each point with its pair_id just above the marker
    for idx, pid in enumerate(pair_ids):
        ax.annotate(pid,
                    (umap_coords[idx, 0], umap_coords[idx, 1]),
                    fontsize=6, ha="center", va="bottom")

    ax.set_title(title, fontsize=10)
    ax.set_xlabel("UMAP-1", fontsize=8)
    ax.set_ylabel("UMAP-2", fontsize=8)
    ax.legend(fontsize=7, loc="best", framealpha=0.7)


# ══════════════════════════════════════════════════════════════════════════════
# 8.  MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Cluster Upstream/Downstream sequences from a genomic-peak TSV file.\n"
            "Uses hierarchical clustering for the dendrogram, "
            "K-Means (on UMAP embedding) for final cluster labels, and "
            "logomaker for per-cluster sequence logos."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # ── I/O ───────────────────────────────────────────────────────────────────
    parser.add_argument("input",
                        help="Path to the tab-separated matches file from "
                             "sequence_matcher_with_motifs.py.")
    parser.add_argument("--output-prefix", default="clustering",
                        help="Prefix for all output files (TSV, plots, logos).")
    parser.add_argument("--no-plots", action="store_true",
                        help="Skip generating clustering summary plots.")
    parser.add_argument("--no-logos", action="store_true",
                        help="Skip generating per-cluster sequence logos.")

    # ── K-Means ───────────────────────────────────────────────────────────────
    parser.add_argument("--n-clusters", type=int, default=3,
                        help="Number of clusters k for K-Means.")
    parser.add_argument("--kmeans-seed", type=int, default=42,
                        help="Random seed for K-Means reproducibility.")

    # ── Hierarchical linkage (dendrogram only) ────────────────────────────────
    parser.add_argument("--linkage",
                        choices=["single", "complete", "average", "ward"],
                        default="average",
                        help="Linkage method for the hierarchical dendrogram.")

    # ── UMAP ──────────────────────────────────────────────────────────────────
    parser.add_argument("--umap-n-neighbors", type=int, default=5,
                        help="UMAP n_neighbors — lower = more local structure.")
    parser.add_argument("--umap-min-dist", type=float, default=0.3,
                        help="UMAP min_dist — smaller = tighter clusters.")
    parser.add_argument("--umap-seed", type=int, default=42,
                        help="Random seed for UMAP reproducibility.")

    # ── Alignment ─────────────────────────────────────────────────────────────
    parser.add_argument("--alignment-mode",
                        choices=["global", "local"], default="global",
                        help="global = Needleman-Wunsch,  local = Smith-Waterman.")

    # ── Combined-mode weights ─────────────────────────────────────────────────
    parser.add_argument("--w-upstream", type=float, default=0.5,
                        help="Weight for Upstream distance in combined mode.")
    parser.add_argument("--w-downstream", type=float, default=0.5,
                        help="Weight for Downstream distance in combined mode.")

    args = parser.parse_args()

    # ── Step 1: Load and validate input ───────────────────────────────────────
    # full_df  = all rows from the input TSV, with 'pair_id' column added
    # uniq_df  = one row per unique (Upstream, Downstream) sequence pair
    # gene_map = pair_id → [GeneID, …] for axis annotation
    full_df, uniq_df, id_col, seq1_col, seq2_col, gene_map = load_input(args.input)

    pair_ids  = uniq_df["pair_id"].tolist()
    seqs_up   = uniq_df[seq1_col].tolist()
    seqs_down = uniq_df[seq2_col].tolist()

    # Build axis tick labels of the form "pair_001\n(GeneA, GeneB)"
    tick_labs = _tick_labels(pair_ids, gene_map)

    # ── Step 2: Compute pairwise alignment distance matrices ──────────────────
    # Three matrices: upstream-only, downstream-only, combined (weighted avg)
    print(f"\nComputing Upstream   distances "
          f"({len(seqs_up)} seqs, {args.alignment_mode} alignment) ...")
    dist_up   = pairwise_distance_matrix(seqs_up,   mode=args.alignment_mode)

    print(f"Computing Downstream distances "
          f"({len(seqs_down)} seqs, {args.alignment_mode} alignment) ...")
    dist_down = pairwise_distance_matrix(seqs_down, mode=args.alignment_mode)

    print("Computing combined   distances ...")
    dist_comb = combined_distance(dist_up, dist_down,
                                  w1=args.w_upstream,
                                  w2=args.w_downstream)

    # Bundle matrices with human-readable labels for the plotting loop
    dist_matrices = {
        "upstream":   (dist_up,   "Upstream sequences"),
        "downstream": (dist_down, "Downstream sequences"),
        "combined":   (dist_comb, "Upstream + Downstream (combined)"),
    }

    # ── Step 3: UMAP dimensionality reduction ─────────────────────────────────
    # Each distance matrix is independently embedded into 2-D
    print()
    umap_coords = {}
    for key, (dist, label) in dist_matrices.items():
        print(f"  [{key:10s}]  UMAP embedding ...")
        umap_coords[key] = run_umap(
            dist.copy(),             # copy to prevent in-place diagonal modification
            n_neighbors=args.umap_n_neighbors,
            min_dist=args.umap_min_dist,
            random_state=args.umap_seed,
        )

    # ── Step 4: K-Means clustering on UMAP coordinates ───────────────────────
    print()
    km_labels = {}
    for key in dist_matrices:
        lbls = kmeans_cluster(umap_coords[key],
                              n_clusters=args.n_clusters,
                              random_state=args.kmeans_seed)
        print(f"  [{key:10s}]  K-Means (k={args.n_clusters}) -> "
              f"clusters: {sorted(set(lbls))}")
        km_labels[key] = lbls

    # ── Step 5: Hierarchical linkage for dendrogram visualisation ─────────────
    # Computed after K-Means because the dendrogram leaves are coloured by
    # the K-Means labels, not by the hierarchical cut
    print()
    linkage_Z = {}
    for key, (dist, label) in dist_matrices.items():
        print(f"  [{key:10s}]  Hierarchical linkage ({args.linkage}) ...")
        linkage_Z[key] = hierarchical_linkage(dist.copy(),
                                              linkage_method=args.linkage)

    # ── Step 6: Save cluster-labelled results TSV ─────────────────────────────
    # Start from the unique-pair summary, then attach K-Means labels and
    # UMAP coordinates for all three modes
    pair_summary = uniq_df[[seq1_col, seq2_col, "pair_id"]].copy()
    pair_summary["kmeans_upstream"]   = km_labels["upstream"]
    pair_summary["kmeans_downstream"] = km_labels["downstream"]
    pair_summary["kmeans_combined"]   = km_labels["combined"]

    # Attach UMAP coordinates so downstream tools can reproduce the scatter plots
    for key in dist_matrices:
        pair_summary[f"umap_{key}_1"] = umap_coords[key][:, 0]
        pair_summary[f"umap_{key}_2"] = umap_coords[key][:, 1]

    # Merge cluster labels back to every row in the full input DataFrame
    # (multiple GeneID rows that share a pair_id all receive the same label)
    cluster_cols = ["kmeans_upstream", "kmeans_downstream", "kmeans_combined"]
    out_df = full_df.merge(
        pair_summary[["pair_id"] + cluster_cols],
        on="pair_id", how="left"
    )

    # Re-order columns: GeneID | cluster labels | remaining original columns
    other_cols = [c for c in out_df.columns if c not in cluster_cols + ["pair_id"]]
    out_df     = out_df[other_cols[:1] + cluster_cols + other_cols[1:]]

    tsv_out = f"{args.output_prefix}_results.tsv"
    out_df.to_csv(tsv_out, sep="\t", index=False)
    print(f"\nCluster assignments saved -> {tsv_out}")

    # ── Step 7: Per-cluster sequence logos ────────────────────────────────────
    # Generate upstream and downstream logos for each cluster in each mode.
    # Logos reveal conserved sequence patterns within a cluster, showing whether
    # sequences in the same cluster share a common regulatory motif context.
    if not args.no_logos:
        generate_cluster_logos(
            uniq_df=uniq_df,
            km_labels_dict=km_labels,   # {'upstream': [...], 'downstream': [...], 'combined': [...]}
            seq1_col=seq1_col,
            seq2_col=seq2_col,
            output_prefix=args.output_prefix,
            n_clusters=args.n_clusters,
        )

    # ── Step 8: Clustering summary plots (dendrogram / heatmap / UMAP) ────────
    if not args.no_plots:
        fig, axes = plt.subplots(3, 3, figsize=(24, 20))
        fig.suptitle(
            f"Sequence Clustering  |  Upstream vs Downstream\n"
            f"Alignment: {args.alignment_mode}  |  "
            f"K-Means k={args.n_clusters}  |  "
            f"Hierarchical linkage: {args.linkage}  |  "
            f"UMAP: n_neighbors={args.umap_n_neighbors}, min_dist={args.umap_min_dist}",
            fontsize=11, fontweight="bold", y=1.01
        )

        row_titles = {
            "upstream":   "Upstream sequences",
            "downstream": "Downstream sequences",
            "combined":   "Upstream + Downstream (combined)",
        }

        # Each row of the figure corresponds to one clustering mode
        for row_idx, key in enumerate(["upstream", "downstream", "combined"]):
            dist, title = dist_matrices[key]
            lbls        = km_labels[key]
            Z           = linkage_Z[key]
            coords      = umap_coords[key]

            # Col 0 – Dendrogram with K-Means-coloured leaves
            plot_dendrogram(Z, tick_labs, lbls,
                            f"Dendrogram – {row_titles[key]}", axes[row_idx][0])

            # Col 1 – Pairwise distance heatmap with cluster-coloured diagonal
            plot_heatmap(dist, tick_labs, lbls,
                         f"Distance heatmap – {row_titles[key]}", axes[row_idx][1])

            # Col 2 – UMAP 2-D scatter with K-Means cluster colours
            plot_umap(coords, lbls, pair_ids, gene_map,
                      f"UMAP + K-Means – {row_titles[key]}", axes[row_idx][2])

        plt.tight_layout()
        plot_out = f"{args.output_prefix}_plots.png"
        plt.savefig(plot_out, dpi=150, bbox_inches="tight")
        print(f"Plots saved            -> {plot_out}")
        plt.close()

    # ── Step 9: Console cluster summary ───────────────────────────────────────
    # Print a readable table showing each unique pair, its gene members,
    # and its cluster assignment under each of the three modes
    print("\n── Cluster Summary (unique pairs) ──")
    display = pair_summary[["pair_id"] + cluster_cols].copy()
    display.insert(1, "GeneIDs",
                   display["pair_id"].map(
                       lambda p: ", ".join(gene_map.get(p, []))
                   ))
    print(display.to_string(index=False))


if __name__ == "__main__":
    main()
