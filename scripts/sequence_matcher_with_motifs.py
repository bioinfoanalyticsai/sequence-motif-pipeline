#!/usr/bin/env python3
"""
sequence_matcher_with_motifs.py
================================
Searches for a query motif inside pre-fetched peak sequences, producing
BED, BigWig, flanking sequence logos, PSSM tables, and MEME motif files.

JASPAR FILE IS THE SOURCE OF TRUTH
------------------------------------
In all three modes, if --jaspar is provided the consensus sequence is derived
directly from the JASPAR counts matrix via consensus_from_counts() — the
SAME function used in find_motif_both_strands_and_chromsome_wise.py.

Providing the same --jaspar file to both scripts guarantees:
  • The motif length is identical (taken from the matrix width).
  • The Hamming reference string for mismatch mode is bit-identical.
  • Every hit found in a peak is guaranteed to also appear in the whole-genome
    scan, i.e. peak hits are a strict subset of genome hits.

The --query / -q string is still accepted for convenience, but is validated
against the JASPAR consensus and a warning is printed when they differ.
The JASPAR consensus always takes precedence when --jaspar is provided.

Search modes (--mode):
  exact     – IUPAC-aware regex search (no mismatches)
  mismatch  – Hamming sliding-window, allows up to --max-mismatches
  pwm       – log-odds PWM with Monte-Carlo p-value threshold

Strand handling — single-pass, consistent with find_motif reference script:
  At every window position i on the FORWARD peak sequence:
    • Forward hit  : score window[i:i+k] against motif / PWM
    • Reverse hit  : score rc(window[i:i+k]) against motif / PWM
  Coordinates are ALWAYS reported as (peak_start + i, peak_start + i + k)
  regardless of strand — standard BED half-open convention on the + strand.
  The matched sequence in the TSV is 5'→3' on the reported strand.

Usage
-----
  # Mismatch mode — JASPAR drives consensus (recommended)
  python sequence_matcher_with_motifs.py \\
      -i input.tsv -o output.tsv \\
      --jaspar AGGGGATTTCCC.jaspar \\
      --mode mismatch --max-mismatches 2

  # PWM mode (JASPAR required)
  python sequence_matcher_with_motifs.py \\
      -i input.tsv -o output.tsv \\
      --jaspar AGGGGATTTCCC.jaspar \\
      --mode pwm --pvalue 1e-4

  # Mismatch mode — fallback without JASPAR (not recommended for cross-tool consistency)
  python sequence_matcher_with_motifs.py \\
      -i input.tsv -o output.tsv \\
      -q AGGGGATTTCCC \\
      --mode mismatch --max-mismatches 2

Dependencies
------------
  biopython, pandas, numpy, logomaker, matplotlib, pyBigWig
"""

import argparse
import bisect
import math
import random
import re
import sys
import pyBigWig
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
import logomaker
import matplotlib.pyplot as plt
from collections import defaultdict
from pathlib import Path

# ── Constants ─────────────────────────────────────────────────────────────────

BASES      = ['A', 'C', 'G', 'T']
BASE_IDX   = {b: i for i, b in enumerate(BASES)}
PSEUDOCOUNT = 0.1      # matches reference script default

IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': '[AG]',  'Y': '[CT]',  'S': '[GC]',
    'W': '[AT]',  'K': '[GT]',  'M': '[AC]',
    'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]',
    'V': '[ACG]', 'N': '[ACGT]',
}

COMPLEMENT = str.maketrans(
    'ACGTRYSWKMBDHVNacgtryswkmbdhvn',
    'TGCAYRSWMKVHDBNtgcayrswmkvhdbn',
)

# ── Column name aliases (lower-cased for matching) ────────────────────────────
_COL_ALIASES = {
    'chromosome': ['chromosome', 'chrom', 'chr'],
    'peak_start': ['peak_start', 'start', 'chromstart'],
    'peak_end':   ['peak_end',   'end',   'chromend'],
}

# ─────────────────────────── Utilities ────────────────────────────────────────

def normalize(seq: str) -> str:
    return seq.upper().replace('U', 'T')

def reverse_complement(seq: str) -> str:
    """Reverse complement handling full IUPAC alphabet."""
    return seq.upper().translate(COMPLEMENT)[::-1]

def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))

def has_n(seq: str) -> bool:
    return 'N' in seq

def iupac_to_regex(motif: str) -> str:
    return ''.join(IUPAC.get(b.upper(), b.upper()) for b in motif)

def _resolve_columns(df: pd.DataFrame) -> dict:
    lower_map = {c.lower(): c for c in df.columns}
    resolved  = {}
    for canonical, aliases in _COL_ALIASES.items():
        match = next((lower_map[a] for a in aliases if a in lower_map), None)
        if match is None:
            raise ValueError(
                f"Required column '{canonical}' not found. "
                f"Expected one of: {aliases}. Got: {list(df.columns)}"
            )
        resolved[canonical] = match
    return resolved

def load_reference(fasta_path: str) -> dict:
    """Return {seq_id: sequence_str} from a FASTA file."""
    ref = {}
    for record in SeqIO.parse(fasta_path, 'fasta'):
        ref[record.id] = normalize(str(record.seq))
    return ref

# ─────────────────────────── Flank extraction ────────────────────────────────

def extract_flanks(seq: str, pos: int, k: int, up: int, down: int,
                   ref_seq=None, ref_offset: int = 0):
    """
    Extract upstream / downstream flanking bases around a hit.
    seq        : forward peak sequence
    pos        : hit start within seq (forward coordinates)
    k          : motif length
    ref_seq    : full chromosome sequence for padding beyond peak boundaries
    ref_offset : peak_start (where seq starts in the chromosome)
    """
    upstream_start   = pos - up
    upstream_end     = pos
    downstream_start = pos + k
    downstream_end   = pos + k + down

    upstream   = seq[max(0, upstream_start):upstream_end]
    downstream = seq[downstream_start:min(len(seq), downstream_end)]

    if len(upstream) < up:
        deficit = up - len(upstream)
        if ref_seq is not None:
            ref_pos = ref_offset + upstream_start
            extra   = ref_seq[max(0, ref_pos):ref_offset + pos - len(upstream)]
            extra   = extra.rjust(deficit, 'N')
        else:
            extra = 'N' * deficit
        upstream = extra + upstream

    if len(downstream) < down:
        deficit = down - len(downstream)
        if ref_seq is not None:
            ref_pos_start = ref_offset + downstream_start + len(downstream)
            ref_pos_end   = ref_offset + downstream_end
            extra = ref_seq[ref_pos_start:min(len(ref_seq), ref_pos_end)]
            extra = extra.ljust(deficit, 'N')
        else:
            extra = 'N' * deficit
        downstream = downstream + extra

    return upstream, downstream

# ─────────────────────────── Search modes ────────────────────────────────────
# All three functions below use the SAME single-pass strategy as the reference:
#   • Iterate once over the FORWARD sequence.
#   • Test rc(window) for the minus strand at the SAME position i.
#   • Yield (local_start, local_end, strand, bed_score, n_mismatches, matched_seq)
#   • Caller adds peak_start to local_start/local_end → genomic BED coordinates.
#
# The motif / consensus parameter in search_exact and search_mismatch MUST be
# derived from the JASPAR file via consensus_from_counts() so that Hamming
# comparisons are bit-identical to find_motif_both_strands_and_chromsome_wise.py.


def search_exact(seq: str, motif: str, both_strands: bool):
    """
    IUPAC regex search.
    Minus strand: search for rc(motif) in the forward sequence.
    Coordinates are always forward-strand positions.

    Yields: (local_start, local_end, strand, bed_score, n_mismatches, matched_seq)
    """
    fwd_pat = re.compile(iupac_to_regex(motif))
    rev_pat = re.compile(iupac_to_regex(reverse_complement(motif)))

    for m in fwd_pat.finditer(seq):
        window = seq[m.start():m.end()]
        yield m.start(), m.end(), '+', 1000, 0, window

    if both_strands:
        for m in rev_pat.finditer(seq):
            # rc(window) is the 5'→3' sequence on the minus strand
            matched = reverse_complement(seq[m.start():m.end()])
            yield m.start(), m.end(), '-', 1000, 0, matched


def search_mismatch(seq: str, motif: str, max_mismatches: int,
                    both_strands: bool):
    """
    Hamming sliding-window scan — mirrors find_motif_both_strands_and_chromsome_wise.py
    exactly.

    At each position i:
      • Compare window to motif            → forward hit if mm ≤ max_mismatches
      • Compare rc(window) to motif        → reverse hit if mm ≤ max_mismatches

    The motif string MUST be derived from the JASPAR file via
    consensus_from_counts() so that Hamming comparisons produce identical
    results in both tools — the necessary condition for peak hits being a
    strict subset of genome hits.

    Score = 1000 × (1 – mismatches / motif_length)
    Coordinates are always forward-strand positions.

    Yields: (local_start, local_end, strand, bed_score, n_mismatches, matched_seq)
    """
    k = len(motif)
    for i in range(len(seq) - k + 1):
        window = seq[i:i + k]
        if has_n(window):
            continue

        # ── Forward strand ────────────────────────────────────────────────────
        mm_fwd = hamming(window, motif)
        if mm_fwd <= max_mismatches:
            score = int(1000 * (1 - mm_fwd / k))
            yield i, i + k, '+', score, mm_fwd, window

        # ── Reverse strand ────────────────────────────────────────────────────
        # rc(window) is scored against the FORWARD motif.
        # Identical strategy to find_motif_both_strands_and_chromsome_wise.py.
        if both_strands:
            rc_window = reverse_complement(window)
            mm_rev    = hamming(rc_window, motif)
            if mm_rev <= max_mismatches:
                score = int(1000 * (1 - mm_rev / k))
                yield i, i + k, '-', score, mm_rev, rc_window


def search_pwm(seq: str, pwm_lom: list, pvalue_threshold: float,
               null_dist: list, both_strands: bool,
               consensus: str = None):
    """
    PWM log-odds scan with Monte-Carlo p-value threshold.
    Mirrors search_pwm_with_mismatch_report() in the reference script.

    At each position i:
      • Score window against forward PWM    → forward hit if p ≤ threshold
      • Score rc(window) against forward PWM → reverse hit if p ≤ threshold

    The consensus parameter MUST be derived from the JASPAR file via
    consensus_from_counts() so mismatch annotations match the genome script.

    Coordinates are always forward-strand positions.

    Yields (local_start, local_end, strand, bed_score, n_mismatches,
            mismatch_positions, raw_score, matched_seq)
    """
    k = len(pwm_lom)
    for i in range(len(seq) - k + 1):
        window = seq[i:i + k]
        if has_n(window):
            continue

        # Forward strand
        fwd_score = _score_sequence(window, pwm_lom)
        if _score_to_pvalue(fwd_score, null_dist) <= pvalue_threshold:
            bed_score    = _score_to_bed_score(fwd_score, pwm_lom)
            n_mm, mm_pos = (_annotate_mm(window, consensus)
                            if consensus else (0, []))
            yield i, i + k, '+', bed_score, n_mm, mm_pos, fwd_score, window

        # Reverse strand (rc of window scored against FORWARD PWM)
        if both_strands:
            rc_window = reverse_complement(window)
            rev_score = _score_sequence(rc_window, pwm_lom)
            if _score_to_pvalue(rev_score, null_dist) <= pvalue_threshold:
                bed_score    = _score_to_bed_score(rev_score, pwm_lom)
                n_mm, mm_pos = (_annotate_mm(rc_window, consensus)
                                if consensus else (0, []))
                yield i, i + k, '-', bed_score, n_mm, mm_pos, rev_score, rc_window


# ─────────────────────────── PWM helpers ─────────────────────────────────────

def parse_jaspar(path: str) -> dict:
    """
    Parse a JASPAR-format counts file.
    Returns {'A': [c0, c1, ...], 'C': [...], 'G': [...], 'T': [...]}
    """
    counts = {b: [] for b in BASES}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            base = line[0].upper()
            if base in BASES:
                nums = [float(x) for x in re.findall(r'[\d.]+', line)]
                counts[base] = nums
    if not all(counts[b] for b in BASES):
        raise ValueError(f"Could not parse all four base rows from: {path}")
    lengths = {len(counts[b]) for b in BASES}
    if len(lengths) != 1:
        raise ValueError(f"Unequal row lengths in PWM: {lengths}")
    return counts


def consensus_from_counts(counts: dict) -> str:
    """
    Derive the consensus sequence from a JASPAR counts dict.

    At each position the base with the highest count wins.
    Ties are broken by BASES order (A > C > G > T).

    This is the SINGLE authoritative way both scripts derive the motif string
    from a JASPAR file, guaranteeing bit-exact agreement between:
      • sequence_matcher_with_motifs.py               (peak-level scan)
      • find_motif_both_strands_and_chromsome_wise.py  (whole-genome scan)
    so that every peak hit is a strict subset of every genome hit.
    """
    n_pos = len(counts['A'])
    return ''.join(max(BASES, key=lambda b: counts[b][i]) for i in range(n_pos))


def counts_to_log_odds(counts: dict, pseudocount: float = PSEUDOCOUNT,
                        bg: dict = None) -> list:
    """
    Convert raw counts → list of per-position log2-odds dicts.
    Identical to the reference script's counts_to_log_odds().
    """
    if bg is None:
        bg = {b: 0.25 for b in BASES}
    pwm = []
    for pos in range(len(counts['A'])):
        col     = {b: counts[b][pos] + pseudocount for b in BASES}
        total   = sum(col.values())
        freq    = {b: col[b] / total for b in BASES}
        logodds = {b: math.log2(freq[b] / bg[b]) for b in BASES}
        pwm.append(logodds)
    return pwm


def _score_sequence(seq: str, pwm: list) -> float:
    return sum(pwm[i].get(base, -10.0) for i, base in enumerate(seq))


def _compute_null_distribution(pwm: list, n_samples: int = 100_000,
                                bg: dict = None) -> list:
    """Monte-Carlo null distribution — identical to reference script."""
    if bg is None:
        bg = {b: 0.25 for b in BASES}
    bases   = list(bg.keys())
    weights = list(bg.values())
    k       = len(pwm)
    scores  = [
        _score_sequence(
            ''.join(random.choices(bases, weights=weights, k=k)), pwm
        )
        for _ in range(n_samples)
    ]
    scores.sort()
    return scores


def _score_to_pvalue(score: float, null_dist: list) -> float:
    return 1.0 - bisect.bisect_left(null_dist, score) / len(null_dist)


def _score_to_bed_score(score: float, pwm: list) -> int:
    min_s = sum(min(pos.values()) for pos in pwm)
    max_s = sum(max(pos.values()) for pos in pwm)
    if max_s == min_s:
        return 500
    return max(0, min(1000, int((score - min_s) / (max_s - min_s) * 1000)))


def _annotate_mm(window: str, consensus: str):
    mm_pos = [i for i, (a, b) in enumerate(zip(window, consensus)) if a != b]
    return len(mm_pos), mm_pos


def _pwm_lom_to_df(pwm_lom: list) -> pd.DataFrame:
    """Convert list-of-dicts PWM back to a DataFrame for logo/MEME output."""
    rows = [{b: d[b] for b in BASES} for d in pwm_lom]
    df   = pd.DataFrame(rows, columns=BASES)
    # Convert log-odds back to frequencies for logo/MEME (2^lo × bg)
    freq = df.apply(lambda col: 2 ** col * 0.25)
    freq = freq.div(freq.sum(axis=1), axis=0)
    return freq


# ─────────────────────────── Motif output helpers ─────────────────────────────

def pwm_from_seqs(seqs: list, pseudocount: float = 0.8) -> pd.DataFrame:
    """Build a frequency PWM DataFrame from equal-length sequences."""
    if not seqs:
        return pd.DataFrame(columns=BASES)
    lengths    = [len(s) for s in seqs]
    target_len = max(set(lengths), key=lengths.count)
    seqs       = [s for s in seqs if len(s) == target_len]
    if not seqs:
        return pd.DataFrame(columns=BASES)
    mat = np.full((target_len, 4), pseudocount)
    for s in seqs:
        for i, b in enumerate(s):
            if b in BASE_IDX:
                mat[i, BASE_IDX[b]] += 1
    mat /= mat.sum(axis=1, keepdims=True)
    return pd.DataFrame(mat, columns=BASES)


def info_content(pwm_df: pd.DataFrame, background: float = 0.25) -> pd.DataFrame:
    pwm = np.clip(pwm_df.values.copy(), 1e-9, 1.0)
    ic  = pwm * np.log2(pwm / background)
    return pd.DataFrame(ic, columns=BASES)


def write_logo(ic_df: pd.DataFrame, outfile: str, title: str):
    plt.figure(figsize=(max(6, len(ic_df) / 2), 3))
    logomaker.Logo(ic_df)
    plt.title(title)
    plt.ylabel('Bits')
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()


def write_meme(pwm_df: pd.DataFrame, name: str, outfile: str):
    with open(outfile, 'w') as f:
        f.write('MEME version 4\n\nALPHABET= ACGT\n\n')
        f.write(f'MOTIF {name}\n')
        f.write(f'letter-probability matrix: alength=4 w={len(pwm_df)}\n')
        for _, row in pwm_df.iterrows():
            f.write(' '.join(f'{v:.4f}' for v in row) + '\n')


# ─────────────────────────── BED / BigWig output ─────────────────────────────

def write_bed(hits: list, outfile: str):
    """
    Write a 6-column BED file.
    chromStart / chromEnd must be TRUE genomic coordinates.
    """
    with open(outfile, 'w') as f:
        for h in hits:
            f.write(
                f"{h['chrom']}\t{h['chromStart']}\t{h['chromEnd']}\t"
                f"{h['name']}\t{h['score']}\t{h['strand']}\n"
            )


def write_bigwig(hits: list, chrom_sizes: dict, outfile: str):
    """
    Write per-base score coverage to a BigWig.
    All coordinates in hits are TRUE genomic coordinates.
    """
    if not hits:
        print('  [BigWig] No hits – skipping BigWig generation.')
        return

    bw     = pyBigWig.open(outfile, 'w')
    header = sorted(chrom_sizes.items())
    bw.addHeader(header)

    cov = defaultdict(lambda: defaultdict(float))
    for h in hits:
        chrom  = h['chrom']
        chrlen = chrom_sizes.get(chrom, 0)
        start  = max(0, h['chromStart'])
        end    = min(chrlen, h['chromEnd'])
        for bp in range(start, end):
            cov[chrom][bp] += float(h['score'])

    for chrom, _ in header:
        if chrom not in cov:
            continue
        positions = sorted(cov[chrom].keys())
        chrlen    = chrom_sizes[chrom]
        starts, ends, values = [], [], []
        i = 0
        while i < len(positions):
            s = positions[i]
            v = cov[chrom][s]
            e = s + 1
            while (i + 1 < len(positions)
                   and positions[i + 1] == e
                   and cov[chrom][positions[i + 1]] == v):
                e += 1
                i += 1
            if s < chrlen:
                starts.append(s)
                ends.append(min(e, chrlen))
                values.append(v)
            i += 1
        if starts:
            bw.addEntries(
                [chrom] * len(starts), starts,
                ends=ends, values=[float(v) for v in values],
            )

    bw.close()
    print(f'  [BigWig] Written → {outfile}')


def write_chrom_sizes(chrom_sizes: dict, outfile: str):
    with open(outfile, 'w') as f:
        for chrom, size in sorted(chrom_sizes.items()):
            f.write(f'{chrom}\t{size}\n')


# ─────────────────────────── Main ─────────────────────────────────────────────

def main(args):
    df      = pd.read_csv(args.input, sep='\t')
    seq_col = df.columns[-1]
    mode    = args.mode

    # Resolve coordinate columns (case-insensitive)
    coord_cols = _resolve_columns(df)
    chrom_col  = coord_cols['chromosome']
    start_col  = coord_cols['peak_start']
    end_col    = coord_cols['peak_end']

    # Load reference FASTA (chromosome name → full sequence)
    reference = {}
    if args.reference:
        reference = load_reference(args.reference)
        if not reference:
            raise ValueError(f'No sequences in reference FASTA: {args.reference}')
        print(f'Loaded {len(reference)} chromosome(s) from reference FASTA.')

    # Build chromosome sizes from reference for BigWig header
    chrom_sizes = {chrom: len(seq) for chrom, seq in reference.items()}

    # ── Resolve canonical motif from JASPAR (all modes) ───────────────────────
    # consensus_from_counts() is the SINGLE authoritative function used by BOTH
    # scripts to derive the consensus from a JASPAR file.  This guarantees that
    # every Hamming comparison below is against the IDENTICAL reference string
    # as find_motif_both_strands_and_chromsome_wise.py, making peak hits a
    # strict subset of genome hits.
    gc = args.bg_gc / 2
    at = (1.0 - args.bg_gc) / 2
    bg = {'A': at, 'C': gc, 'G': gc, 'T': at}

    jaspar_counts   = None
    jaspar_consensus = None

    if args.jaspar:
        jaspar_counts    = parse_jaspar(args.jaspar)
        jaspar_consensus = consensus_from_counts(jaspar_counts)
        print(f'  [JASPAR] Consensus from file : {jaspar_consensus}')

        # Cross-check: warn if --query differs from JASPAR consensus
        if args.query:
            cli_norm = normalize(args.query)
            if cli_norm != jaspar_consensus:
                print(
                    f'  [WARNING] --query "{cli_norm}" differs from JASPAR '
                    f'consensus "{jaspar_consensus}". '
                    f'Using JASPAR consensus as the authoritative motif.',
                    file=sys.stderr,
                )

    # Canonical motif: JASPAR consensus takes precedence; fallback to --query
    if jaspar_consensus:
        motif = jaspar_consensus
    elif args.query:
        motif = normalize(args.query)
    else:
        raise ValueError('Either --jaspar or --query / -q must be provided.')

    k         = len(motif)
    consensus = motif      # used for mismatch annotation in PWM mode

    # ── Build search engine based on mode ─────────────────────────────────────
    both_strands = not args.fwd_only
    pwm_lom      = None
    null_dist    = None
    pvalue_thr   = args.pvalue

    if mode == 'pwm':
        if not args.jaspar:
            raise ValueError('--jaspar is required for --mode pwm')

        print('  [PWM] Building log-odds PWM from JASPAR file …')
        # jaspar_counts already loaded above; reuse it
        pwm_lom = counts_to_log_odds(jaspar_counts, pseudocount=args.pseudocount, bg=bg)
        k       = len(pwm_lom)

        # Verify motif length matches PWM width (must be identical for subset property)
        if k != len(motif):
            print(
                f'  [WARNING] JASPAR PWM width ({k}) differs from consensus '
                f'length ({len(motif)}). Using PWM width.',
                file=sys.stderr,
            )
            motif = jaspar_consensus  # already derived from same JASPAR file

        print(f'  [PWM] Motif width={k}. '
              f'Sampling null distribution ({args.n_samples:,} sequences)…')
        null_dist    = _compute_null_distribution(pwm_lom, args.n_samples, bg)
        cutoff_idx   = int((1 - pvalue_thr) * len(null_dist))
        score_cutoff = null_dist[min(cutoff_idx, len(null_dist) - 1)]
        print(f'  [PWM] p-value threshold={pvalue_thr}, '
              f'log-odds cutoff={score_cutoff:.4f}')

    print(f'Mode: {mode}  |  Motif (consensus): {motif}  |  '
          f'Both strands: {both_strands}  |  '
          f'Max mismatches: {args.max_mismatches}')

    # ── Accumulate results ────────────────────────────────────────────────────
    flanks      = defaultdict(list)   # strand_direction → [sequence strings]
    out_rows    = []
    bed_hits    = []
    hit_windows = []                  # for PWM-from-hits logo (mismatch/exact modes)

    for _, r in df.iterrows():
        seq        = normalize(str(r[seq_col]))
        chrom      = str(r[chrom_col])
        peak_start = int(r[start_col])
        peak_end   = int(r[end_col])
        ref_seq    = reference.get(chrom)

        if chrom not in chrom_sizes:
            chrom_sizes[chrom] = peak_end   # fallback if no reference FASTA

        # ── Dispatch to the correct search function ────────────────────────────
        if mode == 'exact':
            hits_iter = (
                (ls, le, strand, score, n_mm, None, None, matched)
                for ls, le, strand, score, n_mm, matched
                in search_exact(seq, motif, both_strands)
            )
        elif mode == 'mismatch':
            hits_iter = (
                (ls, le, strand, score, n_mm, None, None, matched)
                for ls, le, strand, score, n_mm, matched
                in search_mismatch(seq, motif, args.max_mismatches, both_strands)
            )
        else:  # pwm
            hits_iter = search_pwm(
                seq, pwm_lom, pvalue_thr, null_dist,
                both_strands, consensus
            )

        for ls, le, strand, bed_score, n_mm, mm_pos, raw_score, matched in hits_iter:

            # Apply optional hard mismatch cap in PWM mode
            if mode == 'pwm' and args.max_mismatches is not None:
                if n_mm > args.max_mismatches:
                    continue

            # ── Genomic coordinates (BED half-open, + strand) ─────────────────
            # Single formula for both strands — coordinates are ALWAYS
            # forward-strand positions, exactly as the reference script does:
            #   g_start = peak_start + local_i
            #   g_end   = peak_start + local_i + motif_len
            genomic_start = peak_start + ls
            genomic_end   = peak_start + le

            # ── Flanking sequences (on the forward sequence, then orient) ──────
            up_fwd, dn_fwd = extract_flanks(
                seq, ls, le - ls,
                args.upstream, args.downstream,
                ref_seq=ref_seq, ref_offset=peak_start,
            )
            if strand == '-':
                up_seq = reverse_complement(dn_fwd)
                dn_seq = reverse_complement(up_fwd)
            else:
                up_seq = up_fwd
                dn_seq = dn_fwd

            flanks[f'{strand}_up'].append(up_seq)
            flanks[f'{strand}_down'].append(dn_seq)
            hit_windows.append(matched)

            name = (f"mm{n_mm}_{chrom}:{genomic_start}-{genomic_end}"
                    f"({strand})")

            bed_hits.append({
                'chrom':      chrom,
                'chromStart': genomic_start,
                'chromEnd':   genomic_end,
                'name':       name,
                'score':      bed_score,
                'strand':     strand,
            })

            row = r.to_dict()
            row.update({
                'Strand':           strand,
                'Matched_Sequence': matched,
                'Mismatches':       n_mm,
                'Genomic_Start':    genomic_start,
                'Genomic_End':      genomic_end,
                'Score':            raw_score if raw_score is not None else bed_score,
                'Upstream':         up_seq,
                'Downstream':       dn_seq,
                'Mode':             mode,
            })
            if mode == 'pwm' and mm_pos is not None:
                row['Mismatch_Positions'] = (
                    ','.join(str(p) for p in mm_pos) if mm_pos else '.'
                )
            out_rows.append(row)

    # ── Write TSV ─────────────────────────────────────────────────────────────
    pd.DataFrame(out_rows).to_csv(args.output, sep='\t', index=False)
    print(f'  [TSV] {len(out_rows)} rows → {args.output}')

    # ── Motif logos / PSSMs / MEME (flanking sequences) ───────────────────────
    for key, seqs in flanks.items():
        if not seqs:
            continue
        pwm_df = pwm_from_seqs(seqs)
        if pwm_df.empty:
            continue
        ic = info_content(pwm_df)
        write_logo(ic, f'{key}_logo.png', key)
        pwm_df.to_csv(f'{key}_pssm.tsv', sep='\t', index=False)
        write_meme(pwm_df, key, f'{key}.meme')

    # Emit hit-motif PWM:
    #   PWM mode  → convert lom back to freq PWM for logo/MEME
    #   other modes → build from the matched windows
    if mode == 'pwm' and pwm_lom:
        motif_pwm_df = _pwm_lom_to_df(pwm_lom)
        label = 'jaspar_motif'
    else:
        motif_pwm_df = pwm_from_seqs(hit_windows)
        label = 'hit_motif'

    if not motif_pwm_df.empty:
        ic_hit = info_content(motif_pwm_df)
        write_logo(ic_hit, f'{label}_logo.png',
                   label.replace('_', ' ').title())
        motif_pwm_df.to_csv(f'{label}_pssm.tsv', sep='\t', index=False)
        write_meme(motif_pwm_df, label, f'{label}.meme')
        print(f'  [PWM] PSSM, logo, and MEME written for {label}.')

    # ── BED output ────────────────────────────────────────────────────────────
    bed_path = args.bed or args.output.replace('.tsv', '.bed')
    write_bed(bed_hits, bed_path)
    print(f'  [BED] {len(bed_hits)} intervals → {bed_path}')

    # Summary
    fwd_n = sum(1 for h in bed_hits if h['strand'] == '+')
    rev_n = sum(1 for h in bed_hits if h['strand'] == '-')
    print(f'  [Summary] + strand: {fwd_n:,}  - strand: {rev_n:,}  '
          f'total: {len(bed_hits):,}')

    # ── BigWig output ─────────────────────────────────────────────────────────
    if chrom_sizes:
        bw_path = args.bigwig or args.output.replace('.tsv', '.bw')
        cs_path = bw_path.replace('.bw', '_chrom_sizes.txt')
        write_chrom_sizes(chrom_sizes, cs_path)
        print(f'  [BigWig] Chrom sizes → {cs_path}')
        write_bigwig(bed_hits, chrom_sizes, bw_path)
    else:
        print('  [BigWig] No chromosome sizes available – skipping.')


# ─────────────────────────── CLI ──────────────────────────────────────────────

if __name__ == '__main__':
    p = argparse.ArgumentParser(
        description=(
            'Motif search within pre-fetched peak sequences.\n'
            '\n'
            'JASPAR FILE IS THE SOURCE OF TRUTH\n'
            '  Provide --jaspar in any mode to derive the canonical consensus\n'
            '  from the counts matrix.  Using the same --jaspar file as\n'
            '  find_motif_both_strands_and_chromsome_wise.py guarantees that\n'
            '  every peak hit is a strict subset of genome hits.\n'
            '\n'
            'Uses the same single-pass strand logic as find_motif_both_strands:\n'
            '  forward hit: score window against motif\n'
            '  reverse hit: score rc(window) against motif\n'
            'Coordinates are always reported in + strand BED space.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    p.add_argument('-i', '--input',      required=True,
                   help='Tab-separated input file; last column must contain sequences.')
    p.add_argument('-o', '--output',     required=True,
                   help='Tab-separated output file.')
    p.add_argument('-q', '--query',      default=None,
                   help='Consensus motif string (used for mismatch/exact search '
                        'and as mismatch-annotation reference in PWM mode). '
                        'When --jaspar is also given, this is validated against '
                        'the JASPAR consensus; the JASPAR consensus takes '
                        'precedence if they differ. '
                        'Required only when --jaspar is absent.')
    p.add_argument('-g', '--reference',  default=None,
                   help='Reference FASTA. Used to (1) pad truncated flanks with '
                        'real bases and (2) supply chromosome sizes for BigWig. '
                        'Chromosome column of input must match FASTA record IDs.')

    # ── Mode ─────────────────────────────────────────────────────────────────
    p.add_argument('--mode', choices=['exact', 'mismatch', 'pwm'],
                   default='mismatch',
                   help='Search mode: exact (IUPAC regex), mismatch (Hamming), '
                        'or pwm (log-odds + p-value). Default: mismatch.')
    p.add_argument('--max-mismatches', type=int, default=2,
                   help='Max Hamming mismatches for mismatch mode, or hard '
                        'post-filter cap in PWM mode (default: 2).')

    # ── PWM / JASPAR ─────────────────────────────────────────────────────────
    p.add_argument('--jaspar',      default=None,
                   help='JASPAR-format PWM counts file. '
                        'Required for --mode pwm. '
                        'Strongly recommended for exact/mismatch modes: '
                        'the consensus is derived via consensus_from_counts(), '
                        'the same function used in '
                        'find_motif_both_strands_and_chromsome_wise.py, '
                        'guaranteeing that peak hits are a strict subset of '
                        'genome hits.')
    p.add_argument('--pvalue',      type=float, default=1e-4,
                   help='P-value threshold for PWM hit calling (default: 1e-4).')
    p.add_argument('--pseudocount', type=float, default=0.1,
                   help='Pseudocount added to each PWM cell (default: 0.1).')
    p.add_argument('--n-samples',   type=int,   default=100_000,
                   help='Null-distribution sample count for PWM p-value '
                        'calibration (default: 100000).')
    p.add_argument('--bg-gc',       type=float, default=0.5,
                   help='Background GC fraction for PWM null distribution '
                        '(default: 0.5 = uniform). E.g. 0.42 for human.')

    # ── Strand ───────────────────────────────────────────────────────────────
    p.add_argument('--fwd-only', action='store_true',
                   help='Scan forward (+) strand only (default: both strands).')

    # ── Flanks ───────────────────────────────────────────────────────────────
    p.add_argument('--upstream',   type=int, default=150,
                   help='Bases to extract upstream of each hit (default: 150).')
    p.add_argument('--downstream', type=int, default=150,
                   help='Bases to extract downstream of each hit (default: 150).')

    # ── Output paths ─────────────────────────────────────────────────────────
    p.add_argument('--bed',    default=None,
                   help='Output BED file path (default: <o>.bed).')
    p.add_argument('--bigwig', default=None,
                   help='Output BigWig file path (default: <o>.bw).')

    args = p.parse_args()

    # ── Validate argument combinations ────────────────────────────────────────
    if args.mode == 'pwm' and not args.jaspar:
        p.error('--jaspar is required when --mode pwm')
    if not args.jaspar and not args.query:
        p.error(
            '--jaspar or --query / -q must be provided.\n'
            'Prefer --jaspar: the consensus is derived from the counts matrix,\n'
            'guaranteeing subset consistency with '
            'find_motif_both_strands_and_chromsome_wise.py.'
        )

    main(args)
