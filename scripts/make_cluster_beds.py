#!/usr/bin/env python3
"""
Generate cluster-coloured BED files (and optionally BigWig) from clustering_results.tsv.

Colour palette
--------------
Colours are taken from matplotlib's tab10 colormap — the SAME palette used by
cluster_sequences_updated.py — so BED tracks match the dendrogram / UMAP plots:

  Cluster 1  →  RGB(31,  119, 180)   blue
  Cluster 2  →  RGB(255, 127,  14)   orange
  Cluster 3  →  RGB(44,  160,  44)   green
  Cluster 4  →  RGB(214,  39,  40)   red
  Cluster 5  →  RGB(148, 103, 189)   purple
  Cluster 6  →  RGB(140,  86,  75)   brown
  Cluster 7  →  RGB(227, 119, 194)   pink
  Cluster 8  →  RGB(127, 127, 127)   grey
  Cluster 9  →  RGB(188, 189,  34)   olive
  Cluster 10 →  RGB(23,  190, 207)   cyan

Output files
------------
For each clustering mode (upstream / downstream / combined):

  1. <prefix>_<mode>.bed
       BED9 with itemRgb column → load in UCSC / IGV, colours appear immediately.

  2. <prefix>_<mode>_cluster<N>.bed   (one per cluster)
       Individual BED9 files, one cluster per file, for easy per-cluster loading.

  3. <prefix>_<mode>.bedGraph  (optional, with --bigwig flag)
       bedGraph of cluster IDs (float), converted to BigWig via bedGraphToBigWig.
       Requires bedGraphToBigWig and a chromosome-sizes file (--chrom-sizes).

Usage
-----
  # BED files only (no extra tools needed):
  python make_cluster_beds.py \\
      --input  clustering_results.tsv \\
      --prefix my_clusters

  # Also produce BigWig (needs bedGraphToBigWig + chrom-sizes):
  python make_cluster_beds.py \\
      --input       clustering_results.tsv \\
      --prefix      my_clusters \\
      --bigwig \\
      --chrom-sizes mm9.chrom.sizes

  # Use a specific clustering mode only:
  python make_cluster_beds.py --input clustering_results.tsv --prefix out --mode combined
"""

import argparse
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


# ── Colour palette (mirrors cluster_sequences_updated.py exactly) ─────────────

_CMAP = plt.get_cmap("tab10")

def _cluster_rgb(lbl: int) -> tuple:
    """Return (R, G, B) integer tuple for 1-based cluster label using tab10."""
    rgba = _CMAP((int(lbl) - 1) % 10)
    return (int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255))

def _rgb_str(lbl: int) -> str:
    """Return 'R,G,B' string for BED itemRgb column."""
    return ",".join(str(x) for x in _cluster_rgb(lbl))


# ── BED / bedGraph writers ────────────────────────────────────────────────────

BED_HEADER = (
    'track name="{name}" description="{desc}" '
    'itemRgb="On" visibility=2 useScore=0\n'
)

def _bed9_row(chrom, start, end, name, score, strand, cluster_lbl) -> str:
    """
    Return one BED9 line.

    BED9 columns:
      chrom  chromStart  chromEnd  name  score  strand
      thickStart  thickEnd  itemRgb

    thickStart / thickEnd are set equal to chromStart / chromEnd (solid block).
    itemRgb carries the cluster colour so track viewers render it immediately.
    """
    rgb = _rgb_str(cluster_lbl)
    return (
        f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t"
        f"{start}\t{end}\t{rgb}\n"
    )


def write_combined_bed(df: pd.DataFrame, mode_col: str,
                       output_path: str, track_name: str, desc: str) -> None:
    """Write a single BED9 with all clusters, colour-coded by cluster label."""
    with open(output_path, "w") as fh:
        fh.write(BED_HEADER.format(name=track_name, desc=desc))
        for _, row in df.iterrows():
            fh.write(_bed9_row(
                chrom       = row["Chromosome"],
                start       = int(row["Peak_Start"]),
                end         = int(row["Peak_End"]),
                name        = f"{row['GeneID']}_c{row[mode_col]}",
                score       = 0,
                strand      = row.get("Strand", "."),
                cluster_lbl = int(row[mode_col]),
            ))
    print(f"  [BED]    {output_path}  ({len(df)} regions)")


def write_per_cluster_beds(df: pd.DataFrame, mode_col: str,
                           prefix: str, mode: str) -> None:
    """Write one BED9 file per cluster label."""
    clusters = sorted(df[mode_col].dropna().unique(), key=lambda x: int(x))
    for clust in clusters:
        sub   = df[df[mode_col] == clust]
        r, g, b = _cluster_rgb(int(clust))
        out   = f"{prefix}_{mode}_cluster{clust}.bed"
        desc  = (
            f"Cluster {clust} | mode={mode} | "
            f"n={len(sub)} | RGB({r},{g},{b})"
        )
        with open(out, "w") as fh:
            fh.write(BED_HEADER.format(
                name=f"cluster{clust}_{mode}", desc=desc
            ))
            for _, row in sub.iterrows():
                fh.write(_bed9_row(
                    chrom       = row["Chromosome"],
                    start       = int(row["Peak_Start"]),
                    end         = int(row["Peak_End"]),
                    name        = row["GeneID"],
                    score       = 0,
                    strand      = row.get("Strand", "."),
                    cluster_lbl = int(clust),
                ))
        print(f"  [BED]    {out}  ({len(sub)} regions)")


def write_bedgraph(df: pd.DataFrame, mode_col: str,
                   output_path: str) -> None:
    """
    Write a bedGraph where the score is the integer cluster ID.

    bedGraph format: chrom  chromStart  chromEnd  value
    Sorted by chrom then start (required by bedGraphToBigWig).
    """
    sdf = df.sort_values(["Chromosome", "Peak_Start", "Peak_End"])
    with open(output_path, "w") as fh:
        fh.write('track type=bedGraph\n')
        for _, row in sdf.iterrows():
            fh.write(
                f"{row['Chromosome']}\t"
                f"{int(row['Peak_Start'])}\t"
                f"{int(row['Peak_End'])}\t"
                f"{int(row[mode_col])}\n"
            )
    print(f"  [bedGraph] {output_path}")


def bedgraph_to_bigwig(bedgraph_path: str,
                       chrom_sizes: str,
                       bigwig_path: str) -> bool:
    """
    Convert a bedGraph to BigWig using UCSC's bedGraphToBigWig tool.
    Returns True on success, False on failure.
    """
    cmd = ["bedGraphToBigWig", bedgraph_path, chrom_sizes, bigwig_path]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"  [BigWig]   {bigwig_path}")
            return True
        else:
            print(f"  [ERROR] bedGraphToBigWig failed:\n{result.stderr}",
                  file=sys.stderr)
            return False
    except FileNotFoundError:
        print(
            "  [ERROR] bedGraphToBigWig not found in PATH.\n"
            "          Download from https://hgdownload.soe.ucsc.edu/admin/exe/",
            file=sys.stderr,
        )
        return False


# ── Colour legend printer ─────────────────────────────────────────────────────

def print_color_legend(clusters: list) -> None:
    print("\n── Cluster colour legend (tab10, same as clustering script) ──")
    for c in sorted(clusters, key=int):
        r, g, b = _cluster_rgb(int(c))
        hex_col = "#{:02X}{:02X}{:02X}".format(r, g, b)
        print(f"  Cluster {c:>2}  →  RGB({r:3},{g:3},{b:3})  {hex_col}")
    print()


# ── Main ──────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Create cluster-coloured BED9 (and optionally BigWig) files "
            "from the TSV produced by cluster_sequences_updated.py."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--input",  "-i", required=True,
                   help="clustering_results.tsv from cluster_sequences_updated.py")
    p.add_argument("--prefix", "-p", default="clusters",
                   help="Output file prefix (default: clusters)")
    p.add_argument("--mode",   "-m",
                   choices=["upstream", "downstream", "combined", "all"],
                   default="all",
                   help="Which cluster column(s) to use (default: all three)")
    p.add_argument("--per-cluster", action="store_true",
                   help="Also write one BED file per cluster")
    p.add_argument("--bigwig", action="store_true",
                   help="Convert bedGraphs to BigWig (requires bedGraphToBigWig "
                        "in PATH and --chrom-sizes)")
    p.add_argument("--chrom-sizes", metavar="FILE",
                   help="Chromosome sizes file required for BigWig conversion")
    return p


MODE_COLS = {
    "upstream":   "kmeans_upstream",
    "downstream": "kmeans_downstream",
    "combined":   "kmeans_combined",
}


def main():
    args = build_parser().parse_args()

    # ── Validate inputs ───────────────────────────────────────────────────────
    if not Path(args.input).exists():
        print(f"[ERROR] Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    if args.bigwig and not args.chrom_sizes:
        print("[ERROR] --bigwig requires --chrom-sizes <file>", file=sys.stderr)
        sys.exit(1)

    # ── Load TSV ──────────────────────────────────────────────────────────────
    df = pd.read_csv(args.input, sep="\t", dtype=str)
    df.columns = df.columns.str.strip()

    # Normalise chromosome names (some files use 'chr', others don't)
    # We leave them as-is; user should match their genome assembly.
    required = ["GeneID", "Chromosome", "Peak_Start", "Peak_End",
                "kmeans_upstream", "kmeans_downstream", "kmeans_combined"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"[ERROR] Missing columns: {missing}", file=sys.stderr)
        sys.exit(1)

    # Drop rows without valid coordinates
    df = df.dropna(subset=["Chromosome", "Peak_Start", "Peak_End"])
    df["Peak_Start"] = df["Peak_Start"].astype(int)
    df["Peak_End"]   = df["Peak_End"].astype(int)

    print(f"[INFO] Loaded {len(df):,} regions from {args.input}")

    # Determine which modes to process
    modes = (
        ["upstream", "downstream", "combined"]
        if args.mode == "all"
        else [args.mode]
    )

    # Print colour legend for all unique cluster labels found in the data
    all_clusters = set()
    for mode in modes:
        col = MODE_COLS[mode]
        all_clusters.update(df[col].dropna().unique())
    print_color_legend(list(all_clusters))

    # ── Process each mode ─────────────────────────────────────────────────────
    for mode in modes:
        col = MODE_COLS[mode]
        sub = df.dropna(subset=[col]).copy()

        print(f"── Mode: {mode}  ({len(sub):,} regions, "
              f"{sub[col].nunique()} clusters) ──")

        # 1. Combined BED9 (all clusters in one file)
        bed_out = f"{args.prefix}_{mode}.bed"
        write_combined_bed(
            df        = sub,
            mode_col  = col,
            output_path = bed_out,
            track_name  = f"kmeans_{mode}",
            desc        = f"K-Means clusters ({mode} sequences) – tab10 colours",
        )

        # 2. Per-cluster BED9 files (optional)
        if args.per_cluster:
            write_per_cluster_beds(sub, col, args.prefix, mode)

        # 3. bedGraph + BigWig (optional)
        if args.bigwig:
            bg_out = f"{args.prefix}_{mode}.bedGraph"
            bw_out = f"{args.prefix}_{mode}.bw"
            write_bedgraph(sub, col, bg_out)
            bedgraph_to_bigwig(bg_out, args.chrom_sizes, bw_out)

        print()

    print("[DONE] All output files written.")
    print()
    print("── How to load in genome browsers ──────────────────────────────────")
    print(" IGV  : File → Load from File → select any .bed file")
    print("        Colors appear automatically (itemRgb=On).")
    print(" UCSC : My Data → Custom Tracks → paste/upload .bed")
    print("        Or use a .bigWig for the signal track and a .bed for colours.")
    print("────────────────────────────────────────────────────────────────────")


if __name__ == "__main__":
    main()
