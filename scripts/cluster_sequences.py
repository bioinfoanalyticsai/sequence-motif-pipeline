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

  Three modes:
    1. Upstream sequences only
    2. Downstream sequences only
    3. Both columns combined (weighted-average distance)

Output plots (3 rows × 3 columns)
-----------------------------------
  Col 0 – Dendrogram  (hierarchical, coloured by K-Means cluster)
  Col 1 – Pairwise distance heatmap (annotated with K-Means cluster)
  Col 2 – UMAP 2-D scatter (coloured by K-Means cluster)

Dependencies
------------
  pip install biopython scikit-learn scipy matplotlib seaborn pandas numpy umap-learn
"""

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from itertools import combinations
from Bio import Align
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import umap


# ── Default column names (case-insensitive matching applied at load time) ──
_ID_COL   = "GeneID"
_SEQ1_COL = "Upstream"
_SEQ2_COL = "Downstream"


# ══════════════════════════════════════════════
# 1.  INPUT LOADING & VALIDATION
# ══════════════════════════════════════════════

def _resolve_col(df: pd.DataFrame, name: str) -> str:
    """Return the actual column name in *df* matching *name*
    (exact match first, then case-insensitive)."""
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
    Load the TSV, resolve column names, clean sequences.

    Returns
    -------
    full_df  : original DataFrame enriched with 'pair_id'
    uniq_df  : one row per unique (Upstream, Downstream) pair
    id_col   : resolved GeneID column name
    seq1_col : resolved Upstream column name
    seq2_col : resolved Downstream column name
    gene_map : dict  pair_id -> list[GeneID]
    """
    raw = pd.read_csv(path, sep="\t", dtype=str)
    raw.columns = raw.columns.str.strip()

    id_col   = _resolve_col(raw, _ID_COL)
    seq1_col = _resolve_col(raw, _SEQ1_COL)
    seq2_col = _resolve_col(raw, _SEQ2_COL)

    # Drop rows with missing / empty sequences
    before = len(raw)
    raw = raw.dropna(subset=[seq1_col, seq2_col])
    raw = raw[raw[seq1_col].str.strip() != ""]
    raw = raw[raw[seq2_col].str.strip() != ""]
    dropped = before - len(raw)
    if dropped:
        print(f"  [warn] Dropped {dropped} rows with missing/empty sequences.")

    raw[seq1_col] = raw[seq1_col].str.strip().str.upper()
    raw[seq2_col] = raw[seq2_col].str.strip().str.upper()

    print(f"Loaded  : {len(raw)} rows  |  {raw[id_col].nunique()} unique GeneIDs")

    # Assign pair IDs – one per unique (Upstream, Downstream) combination
    raw        = raw.copy()
    pair_key   = raw[seq1_col] + "||" + raw[seq2_col]
    uniq_pairs = pair_key.unique()
    pair_id_map = {p: f"pair_{i+1:03d}" for i, p in enumerate(uniq_pairs)}
    raw["pair_id"] = pair_key.map(pair_id_map)

    # pair_id -> list of GeneIDs (used for axis labels)
    gene_map = raw.groupby("pair_id")[id_col].apply(list).to_dict()

    uniq_df = (
        raw[[seq1_col, seq2_col, "pair_id"]]
        .drop_duplicates(subset="pair_id")
        .reset_index(drop=True)
    )

    print(f"Unique (Upstream, Downstream) pairs: {len(uniq_df)}")

    if len(uniq_df) < 2:
        sys.exit("Error: need at least 2 unique sequence pairs to cluster.")

    return raw, uniq_df, id_col, seq1_col, seq2_col, gene_map


# ══════════════════════════════════════════════
# 2.  PAIRWISE ALIGNMENT DISTANCE
# ══════════════════════════════════════════════

def pairwise_distance_matrix(sequences: list, mode: str = "global") -> np.ndarray:
    """
    NxN distance matrix using BioPython PairwiseAligner.

    distance(i,j) = 1 - score(i,j) / sqrt(score(i,i) * score(j,j))

    mode : 'global'  -> Needleman-Wunsch
           'local'   -> Smith-Waterman
    """
    aligner      = Align.PairwiseAligner()
    aligner.mode = mode

    n           = len(sequences)
    dist_matrix = np.zeros((n, n))
    self_scores = np.array([aligner.score(s, s) for s in sequences])

    for i, j in combinations(range(n), 2):
        raw   = aligner.score(sequences[i], sequences[j])
        denom = np.sqrt(self_scores[i] * self_scores[j])
        sim   = float(raw) / denom if denom > 0 else 0.0
        sim   = min(max(sim, 0.0), 1.0)
        d     = 1.0 - sim
        dist_matrix[i, j] = d
        dist_matrix[j, i] = d

    return dist_matrix


def combined_distance(dist1: np.ndarray, dist2: np.ndarray,
                      w1: float = 0.5, w2: float = 0.5) -> np.ndarray:
    """Normalised weighted average of two distance matrices."""
    return (w1 * dist1 + w2 * dist2) / (w1 + w2)


# ══════════════════════════════════════════════
# 3.  HIERARCHICAL CLUSTERING  (for dendrogram)
# ══════════════════════════════════════════════

def hierarchical_linkage(dist_matrix: np.ndarray,
                         linkage_method: str = "average") -> np.ndarray:
    """Compute linkage matrix Z from a precomputed distance matrix."""
    np.fill_diagonal(dist_matrix, 0.0)
    condensed = squareform(dist_matrix, checks=False)
    return linkage(condensed, method=linkage_method)


# ══════════════════════════════════════════════
# 4.  UMAP EMBEDDING
# ══════════════════════════════════════════════

def run_umap(dist_matrix: np.ndarray,
             n_neighbors: int = 5,
             min_dist: float = 0.3,
             random_state: int = 42) -> np.ndarray:
    """
    2-D UMAP embedding from a precomputed distance matrix.

    n_neighbors : local vs global balance (lower = more local).
                  Auto-capped at n_samples - 1.
    min_dist    : tightness of clusters in embedding (0.0 = very tight).
    """
    n  = dist_matrix.shape[0]
    nn = min(n_neighbors, n - 1)
    if nn != n_neighbors:
        print(f"  [info] UMAP n_neighbors capped at {nn} (only {n} unique pairs).")

    reducer = umap.UMAP(
        n_components=2,
        metric="precomputed",
        n_neighbors=nn,
        min_dist=min_dist,
        random_state=random_state,
    )
    np.fill_diagonal(dist_matrix, 0.0)
    return reducer.fit_transform(dist_matrix)


# ══════════════════════════════════════════════
# 5.  K-MEANS CLUSTERING  (on UMAP coordinates)
# ══════════════════════════════════════════════

def kmeans_cluster(umap_coords: np.ndarray,
                   n_clusters: int,
                   random_state: int = 42) -> np.ndarray:
    """
    K-Means clustering on UMAP 2-D coordinates.

    Running K-Means on the low-dimensional UMAP embedding rather than
    the raw high-dimensional distance matrix respects the non-linear
    structure UMAP has learned and produces more meaningful clusters.

    Returns integer label array (1-indexed for consistency with
    hierarchical fcluster output).
    """
    n = umap_coords.shape[0]
    k = min(n_clusters, n)          # can't have more clusters than points
    if k != n_clusters:
        print(f"  [info] K-Means k capped at {k} (only {n} unique pairs).")

    km     = KMeans(n_clusters=k, random_state=random_state, n_init=20)
    labels = km.fit_predict(umap_coords)
    return labels + 1               # shift to 1-based


# ══════════════════════════════════════════════
# 6.  VISUALISATION HELPERS
# ══════════════════════════════════════════════

# Shared colour palette (tab10, 1-based cluster index)
_CMAP = plt.get_cmap("tab10")

def _cluster_color(lbl: int) -> tuple:
    return _CMAP((lbl - 1) % 10)


def _tick_labels(pair_ids: list, gene_map: dict, max_genes: int = 3) -> list:
    """Build 'pair_001\n(Gene1, Gene2, …)' tick labels."""
    out = []
    for pid in pair_ids:
        genes    = gene_map.get(pid, [])
        gene_str = ", ".join(genes[:max_genes])
        if len(genes) > max_genes:
            gene_str += "..."
        out.append(f"{pid}\n({gene_str})")
    return out


def _legend_patches(unique_labels: list) -> list:
    return [
        mpatches.Patch(color=_cluster_color(l), label=f"Cluster {l}")
        for l in sorted(unique_labels)
    ]


# ── Plot A: Dendrogram coloured by K-Means cluster ──────────────────────────

def plot_dendrogram(Z, tick_labels, km_labels, title, ax):
    """
    Draw dendrogram leaves coloured by their K-Means cluster assignment.
    The cut line is drawn at the height that produces k clusters
    (approximated from the linkage matrix).
    """
    n          = len(tick_labels)
    label_colors = {
        tick_labels[i]: _cluster_color(km_labels[i])
        for i in range(n)
    }

    # scipy dendrogram: set leaf label colours via link_color_func workaround
    # We colour leaves after the fact using the returned x-positions.
    ddata = dendrogram(
        Z,
        labels=tick_labels,
        ax=ax,
        leaf_rotation=90,
        leaf_font_size=7,
        color_threshold=0,          # all links same colour (grey)
        above_threshold_color="grey",
        no_labels=True,             # we draw labels manually
    )

    # Re-draw x-axis tick labels with cluster colours
    ax.set_xticks(range(5, n * 10 + 5, 10))
    leaf_order = ddata["ivl"]      # leaf labels in dendrogram left-to-right order
    for tick, lbl in zip(ax.get_xticks(), leaf_order):
        color = label_colors.get(lbl, "black")
        ax.text(tick, -0.02, lbl,
                ha="right", va="top",
                rotation=90, fontsize=6,
                color=color,
                transform=ax.get_xaxis_transform())

    # Draw a horizontal dashed cut line at the merge height for k clusters
    k = len(set(km_labels))
    if k < len(Z):
        cut_height = (Z[-(k - 1), 2] + Z[-(k), 2]) / 2
        ax.axhline(cut_height, color="red", linestyle="--",
                   linewidth=1.0, label=f"k={k} cut")
        ax.legend(fontsize=7, loc="upper right")

    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Sequence pair  (GeneIDs)", fontsize=8)
    ax.set_ylabel("Distance", fontsize=8)

    # Add colour legend
    patches = _legend_patches(sorted(set(km_labels)))
    ax.legend(handles=patches, fontsize=7, loc="upper right", framealpha=0.7)


# ── Plot B: Distance heatmap annotated with cluster bar ──────────────────────

def plot_heatmap(dist_matrix, tick_labels, km_labels, title, ax):
    """
    Distance heatmap with a coloured cluster-membership strip on the diagonal.
    """
    df_hm = pd.DataFrame(dist_matrix,
                         index=tick_labels, columns=tick_labels)
    sns.heatmap(df_hm, ax=ax, cmap="viridis_r", square=True,
                linewidths=0.2, cbar_kws={"label": "Distance"},
                xticklabels=True, yticklabels=True)

    # Overlay a thin coloured border on each cell diagonal to show cluster
    n = len(tick_labels)
    for i in range(n):
        c = _cluster_color(km_labels[i])
        ax.add_patch(plt.Rectangle(
            (i, i), 1, 1,
            fill=True, color=c, alpha=0.45, linewidth=0
        ))

    ax.set_title(title, fontsize=10)
    ax.tick_params(axis="both", labelsize=6)

    patches = _legend_patches(sorted(set(km_labels)))
    ax.legend(handles=patches, fontsize=6, loc="lower right",
              framealpha=0.8, bbox_to_anchor=(1.0, 0.0))


# ── Plot C: UMAP scatter coloured by K-Means cluster ─────────────────────────

def plot_umap(umap_coords, km_labels, pair_ids, gene_map, title, ax):
    """
    2-D UMAP scatter plot.  Each point = one unique sequence pair.
    Coloured by K-Means cluster.  Annotated with pair_id.
    """
    for lbl in sorted(set(km_labels)):
        mask = np.array(km_labels) == lbl
        ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                   color=_cluster_color(lbl), s=100,
                   edgecolors="k", linewidths=0.5,
                   label=f"Cluster {lbl}", zorder=3)

    for idx, pid in enumerate(pair_ids):
        ax.annotate(pid,
                    (umap_coords[idx, 0], umap_coords[idx, 1]),
                    fontsize=6, ha="center", va="bottom")

    ax.set_title(title, fontsize=10)
    ax.set_xlabel("UMAP-1", fontsize=8)
    ax.set_ylabel("UMAP-2", fontsize=8)
    ax.legend(fontsize=7, loc="best", framealpha=0.7)


# ══════════════════════════════════════════════
# 7.  MAIN
# ══════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Cluster Upstream/Downstream sequences from a genomic-peak TSV file.\n"
            "Uses hierarchical clustering for the dendrogram and "
            "K-Means (on UMAP embedding) for final cluster labels."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # I/O
    parser.add_argument("input",
                        help="Path to the tab-separated input file.")
    parser.add_argument("--output-prefix", default="clustering",
                        help="Prefix for all output files.")
    parser.add_argument("--no-plots", action="store_true",
                        help="Skip generating plots.")

    # K-Means
    parser.add_argument("--n-clusters", type=int, default=3,
                        help="Number of clusters k for K-Means (default: 3).")
    parser.add_argument("--kmeans-seed", type=int, default=42,
                        help="Random seed for K-Means reproducibility.")

    # Hierarchical linkage (dendrogram only)
    parser.add_argument("--linkage",
                        choices=["single", "complete", "average", "ward"],
                        default="average",
                        help="Linkage method for the dendrogram (default: average).")

    # UMAP
    parser.add_argument("--umap-n-neighbors", type=int, default=5,
                        help="UMAP n_neighbors — lower = more local structure.")
    parser.add_argument("--umap-min-dist", type=float, default=0.3,
                        help="UMAP min_dist — smaller = tighter clusters.")
    parser.add_argument("--umap-seed", type=int, default=42,
                        help="Random seed for UMAP reproducibility.")

    # Alignment
    parser.add_argument("--alignment-mode",
                        choices=["global", "local"], default="global",
                        help="global = Needleman-Wunsch,  local = Smith-Waterman.")

    # Combined-mode weights
    parser.add_argument("--w-upstream", type=float, default=0.5,
                        help="Weight for Upstream distance in combined mode.")
    parser.add_argument("--w-downstream", type=float, default=0.5,
                        help="Weight for Downstream distance in combined mode.")

    args = parser.parse_args()

    # ── 1. Load ────────────────────────────────────────────────────────────
    full_df, uniq_df, id_col, seq1_col, seq2_col, gene_map = load_input(args.input)

    pair_ids  = uniq_df["pair_id"].tolist()
    seqs_up   = uniq_df[seq1_col].tolist()
    seqs_down = uniq_df[seq2_col].tolist()
    tick_labs = _tick_labels(pair_ids, gene_map)

    # ── 2. Pairwise distance matrices ──────────────────────────────────────
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

    dist_matrices = {
        "upstream":   (dist_up,   "Upstream sequences"),
        "downstream": (dist_down, "Downstream sequences"),
        "combined":   (dist_comb, "Upstream + Downstream (combined)"),
    }

    # ── 3. UMAP embeddings ─────────────────────────────────────────────────
    print()
    umap_coords = {}
    for key, (dist, label) in dist_matrices.items():
        print(f"  [{key:10s}]  UMAP embedding ...")
        umap_coords[key] = run_umap(
            dist.copy(),
            n_neighbors=args.umap_n_neighbors,
            min_dist=args.umap_min_dist,
            random_state=args.umap_seed,
        )

    # ── 4. K-Means on UMAP coordinates ────────────────────────────────────
    print()
    km_labels = {}
    for key in dist_matrices:
        lbls = kmeans_cluster(umap_coords[key],
                              n_clusters=args.n_clusters,
                              random_state=args.kmeans_seed)
        print(f"  [{key:10s}]  K-Means (k={args.n_clusters}) -> "
              f"clusters: {sorted(set(lbls))}")
        km_labels[key] = lbls

    # ── 5. Hierarchical linkage matrices (for dendrograms) ─────────────────
    print()
    linkage_Z = {}
    for key, (dist, label) in dist_matrices.items():
        print(f"  [{key:10s}]  Hierarchical linkage ({args.linkage}) ...")
        linkage_Z[key] = hierarchical_linkage(dist.copy(),
                                              linkage_method=args.linkage)

    # ── 6. Save TSV results ────────────────────────────────────────────────
    pair_summary = uniq_df[[seq1_col, seq2_col, "pair_id"]].copy()
    pair_summary["kmeans_upstream"]   = km_labels["upstream"]
    pair_summary["kmeans_downstream"] = km_labels["downstream"]
    pair_summary["kmeans_combined"]   = km_labels["combined"]

    # UMAP coordinates saved per mode
    for key in dist_matrices:
        pair_summary[f"umap_{key}_1"] = umap_coords[key][:, 0]
        pair_summary[f"umap_{key}_2"] = umap_coords[key][:, 1]

    cluster_cols = ["kmeans_upstream", "kmeans_downstream", "kmeans_combined"]
    out_df = full_df.merge(
        pair_summary[["pair_id"] + cluster_cols],
        on="pair_id", how="left"
    )
    other_cols = [c for c in out_df.columns if c not in cluster_cols + ["pair_id"]]
    out_df     = out_df[other_cols[:1] + cluster_cols + other_cols[1:]]

    tsv_out = f"{args.output_prefix}_results.tsv"
    out_df.to_csv(tsv_out, sep="\t", index=False)
    print(f"\nCluster assignments saved -> {tsv_out}")

    # ── 7. Plots ───────────────────────────────────────────────────────────
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

        for row_idx, key in enumerate(["upstream", "downstream", "combined"]):
            dist, title = dist_matrices[key]
            lbls        = km_labels[key]
            Z           = linkage_Z[key]
            coords      = umap_coords[key]

            # Col 0 – Dendrogram (coloured by K-Means)
            plot_dendrogram(Z, tick_labs, lbls,
                            f"Dendrogram – {row_titles[key]}", axes[row_idx][0])

            # Col 1 – Distance heatmap (diagonal coloured by K-Means)
            plot_heatmap(dist, tick_labs, lbls,
                         f"Distance heatmap – {row_titles[key]}", axes[row_idx][1])

            # Col 2 – UMAP scatter (coloured by K-Means)
            plot_umap(coords, lbls, pair_ids, gene_map,
                      f"UMAP + K-Means – {row_titles[key]}", axes[row_idx][2])

        plt.tight_layout()
        plot_out = f"{args.output_prefix}_plots.png"
        plt.savefig(plot_out, dpi=150, bbox_inches="tight")
        print(f"Plots saved            -> {plot_out}")
        plt.close()

    # ── 8. Console summary ─────────────────────────────────────────────────
    print("\n── Cluster Summary (unique pairs) ──")
    display = pair_summary[["pair_id"] + cluster_cols].copy()
    display.insert(1, "GeneIDs",
                   display["pair_id"].map(
                       lambda p: ", ".join(gene_map.get(p, []))
                   ))
    print(display.to_string(index=False))


if __name__ == "__main__":
    main()
