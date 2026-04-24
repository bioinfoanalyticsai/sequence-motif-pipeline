#!/usr/bin/env python3
"""
find_motif_both_strands_and_chromsome_wise.py — Genomic motif search with three modes:
  1. exact    — regex with IUPAC ambiguity codes
  2. mismatch — Hamming-distance sliding window, allows k mismatches
  3. pwm      — log-odds Position Weight Matrix with p-value threshold
                + mismatch annotation and optional hard mismatch cap

JASPAR FILE IS THE SOURCE OF TRUTH
------------------------------------
In all three modes, if --jaspar is provided the consensus sequence is derived
directly from the JASPAR counts matrix (argmax per column).  This guarantees
bit-exact agreement with sequence_matcher_with_motifs.py when both tools are
given the same --jaspar file, so peak hits are a strict subset of genome hits.

  --motif / -m  is still accepted for convenience but is validated against the
                JASPAR consensus and a warning is printed if they differ.
                When --jaspar is absent, --motif is the fallback (exact/mismatch).

Chromosome / region targeting:
  --chroms chr1 chr2 chrX   scan only the listed chromosomes
  --region chr1:1000-5000   scan a single genomic region (1-based, inclusive)
  (neither flag)            scan the whole genome (default)

Strand handling (all modes):
  Forward (+) strand : scan sequence directly.
  Reverse (-) strand : at each window position, take reverse_complement(window)
                       and score it against the forward PWM/motif.
                       Coordinates reported are always 0-based half-open on the
                       + strand (BED convention). The matched sequence in the
                       TSV is shown 5'->3' on the reported strand.

Usage examples:
  # Whole genome, mismatch mode driven entirely by JASPAR file
  python find_motif_both_strands_and_chromsome_wise.py \\
      -g genome.fa --jaspar AGGGGATTTCCC.jaspar \\
      -o out.bed --mode mismatch --max-mismatches 2

  # Exact mode with JASPAR (consensus validated against --motif)
  python find_motif_both_strands_and_chromsome_wise.py \\
      -g genome.fa --jaspar AGGGGATTTCCC.jaspar -m AGGGGATTTCCC \\
      -o out.bed --mode exact

  # PWM mode, permissive threshold, hard mismatch cap, single chromosome
  python find_motif_both_strands_and_chromsome_wise.py \\
      -g genome.fa --jaspar pwm.jaspar -m AGGGGATTTCCC \\
      -o out.bed --mode pwm --pvalue 1e-3 --max-mismatches 2 --chroms chr1

  # Mismatch mode, whole genome, both strands
  python find_motif_both_strands_and_chromsome_wise.py \\
      -g genome.fa --jaspar AGGGGATTTCCC.jaspar \\
      -o out.bed --mode mismatch --max-mismatches 2

Dependencies: none beyond the standard library
"""

import re
import sys
import math
import bisect
import random
import argparse
from pathlib import Path


# ── Constants ─────────────────────────────────────────────────────────────────

IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': '[AG]',  'Y': '[CT]',  'S': '[GC]',
    'W': '[AT]',  'K': '[GT]',  'M': '[AC]',
    'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]',
    'V': '[ACG]', 'N': '[ACGT]'
}

COMPLEMENT = str.maketrans(
    'ACGTRYSWKMBDHVNacgtryswkmbdhvn',
    'TGCAYRSWMKVHDBNtgcayrswmkvhdbn'
)

BASES = ['A', 'C', 'G', 'T']


# ── Utilities ─────────────────────────────────────────────────────────────────

def normalize(seq: str) -> str:
    """Upper-case and convert U → T."""
    return seq.upper().replace('U', 'T')

def iupac_to_regex(motif: str) -> str:
    """Convert an IUPAC motif string to a Python regex pattern."""
    return ''.join(IUPAC.get(b.upper(), b.upper()) for b in motif)

def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.upper().translate(COMPLEMENT)[::-1]

def _has_n(seq: str) -> bool:
    """Return True if the sequence contains any ambiguous base."""
    return 'N' in seq

def hamming_distance(a: str, b: str) -> int:
    """Count positions that differ between two equal-length strings."""
    return sum(x != y for x, y in zip(a, b))


# ── JASPAR parsing and consensus extraction ───────────────────────────────────

def parse_jaspar(path: str) -> dict:
    """
    Parse a JASPAR-format counts file.

    Expected layout (brackets and spacing are flexible):
      >MA0001.1  NAME
      A [  10  20  30  15 ]
      C [   5  10   5  20 ]
      G [  30  10  15  25 ]
      T [   5  10  10  10 ]

    Returns {'A': [f0, f1, ...], 'C': [...], 'G': [...], 'T': [...]}
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
      • find_motif_both_strands_and_chromsome_wise.py  (whole-genome scan)
      • sequence_matcher_with_motifs.py               (peak-level scan)
    so that every peak hit is a strict subset of every genome hit.
    """
    n_pos = len(counts['A'])
    return ''.join(max(BASES, key=lambda b: counts[b][i]) for i in range(n_pos))


def resolve_motif(jaspar_path: str | None,
                  cli_motif: str | None,
                  pseudocount: float = 0.1,
                  bg: dict = None) -> tuple:
    """
    Single entry-point that resolves the canonical motif string and (optionally)
    the log-odds PWM from the JASPAR file.

    Priority:
      1. If --jaspar provided → parse counts, derive consensus, build PWM.
         Cross-check against --motif / -m and warn on mismatch.
      2. If only --motif provided → use it directly; pwm_counts = None.

    Returns
    -------
    consensus  : str   — the canonical motif string for Hamming comparison
    pwm_counts : dict | None — raw JASPAR counts (None if no JASPAR file)
    """
    if jaspar_path:
        counts    = parse_jaspar(jaspar_path)
        consensus = consensus_from_counts(counts)
        if cli_motif:
            cli_norm = normalize(cli_motif)
            if cli_norm != consensus:
                print(
                    f"  [WARNING] --motif '{cli_norm}' differs from JASPAR "
                    f"consensus '{consensus}'.  Using JASPAR consensus as the "
                    f"authoritative motif for Hamming comparison.",
                    file=sys.stderr,
                )
        print(f"  [JASPAR] Consensus from file : {consensus}", file=sys.stderr)
        return consensus, counts

    if cli_motif:
        return normalize(cli_motif), None

    raise ValueError("Either --jaspar or --motif / -m must be provided.")


# ── Chromosome / region selection ─────────────────────────────────────────────

def parse_region(region_str: str):
    """
    Parse a region string of the form  chrom:start-end  (1-based, inclusive).

    Returns (chrom, start_0based, end_exclusive) — Python slice coordinates.

    Examples:
      'chr1:1000-2000'  ->  ('chr1', 999, 2000)
      'chrX:1-500'      ->  ('chrX', 0,   500)
    """
    m = re.fullmatch(r'(\S+):(\d+)-(\d+)', region_str.strip())
    if not m:
        raise ValueError(
            f"Region '{region_str}' does not match expected format "
            f"CHROM:START-END (1-based inclusive). "
            f"Example: chr1:1000000-2000000"
        )
    chrom = m.group(1)
    start = int(m.group(2)) - 1   # convert to 0-based
    end   = int(m.group(3))       # already exclusive in Python slice
    if start < 0:
        raise ValueError(f"Region start must be >= 1, got {int(m.group(2))}")
    if end <= start:
        raise ValueError(f"Region end must be > start: {region_str}")
    return chrom, start, end

def build_target_description(chroms: set, region: tuple) -> str:
    """Return a human-readable string describing the scan target."""
    if region:
        chrom, start, end = region
        return (f"region {chrom}:{start+1}-{end} "
                f"(1-based, {end - start:,} bp)")
    if chroms:
        return f"chromosomes: {', '.join(sorted(chroms))}"
    return "whole genome"

def parse_fasta_filtered(fasta_path: str,
                          chroms: set = None,
                          region: tuple = None):
    """
    Yield (chrom, genomic_offset, sequence) from a FASTA file,
    applying chromosome or region filtering.

    Parameters
    ----------
    fasta_path : path to FASTA file
    chroms     : set of chromosome names to include; None means all
    region     : (chrom, start_0based, end_exclusive) tuple; None means
                 whole chromosome

    The returned genomic_offset is the 0-based position of sequence[0] in
    the full chromosome. For whole-chromosome scans this is 0; for region
    scans it equals region_start.  Hit coordinates are reported as
    (genomic_offset + local_position) so they always reflect BED-space
    coordinates regardless of whether a region was sliced.

    Raises ValueError if --region specifies a chromosome absent from the FASTA,
    or if the region extends beyond the chromosome length.
    """
    region_chrom = region[0] if region else None
    region_start = region[1] if region else None
    region_end   = region[2] if region else None
    region_found = False

    current_chrom = None
    parts         = []

    def _emit(chrom, seq_parts):
        nonlocal region_found
        full_seq = ''.join(seq_parts).upper()

        if region:
            # Region mode: only emit the target chromosome
            if chrom != region_chrom:
                return
            region_found = True
            sliced = full_seq[region_start:region_end]
            if not sliced:
                raise ValueError(
                    f"Region {region_chrom}:{region_start+1}-{region_end} "
                    f"is outside chromosome length {len(full_seq):,} bp."
                )
            yield chrom, region_start, sliced

        elif chroms:
            # Chromosome filter mode
            if chrom not in chroms:
                return
            yield chrom, 0, full_seq

        else:
            # Whole-genome mode
            yield chrom, 0, full_seq

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if current_chrom is not None:
                    yield from _emit(current_chrom, parts)
                current_chrom = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)

    if current_chrom is not None:
        yield from _emit(current_chrom, parts)

    if region and not region_found:
        raise ValueError(
            f"Chromosome '{region_chrom}' not found in FASTA: {fasta_path}"
        )


# ── Mode 1: Exact / IUPAC regex ──────────────────────────────────────────────

def search_exact(seq: str, motif: str, both_strands: bool):
    """
    Search by regex (supports IUPAC codes).

    Negative strand strategy:
      Build a separate regex for reverse_complement(motif) and run finditer()
      on the forward sequence. Yields local coordinates; caller adds offset.

    The motif string MUST be derived from the JASPAR file via
    consensus_from_counts() so it is bit-identical to the string used in
    sequence_matcher_with_motifs.py.

    Yields: (local_start, local_end, strand, bed_score, n_mismatches, matched_seq)
    """
    fwd_pat = re.compile(iupac_to_regex(motif))
    rev_pat = re.compile(iupac_to_regex(reverse_complement(motif)))

    for m in fwd_pat.finditer(seq):
        yield m.start(), m.end(), '+', 1000, 0, seq[m.start():m.end()]

    if both_strands:
        for m in rev_pat.finditer(seq):
            matched = reverse_complement(seq[m.start():m.end()])
            yield m.start(), m.end(), '-', 1000, 0, matched


# ── Mode 2: Hamming / mismatch sliding window ─────────────────────────────────

def search_mismatch(seq: str, motif: str, max_mismatches: int,
                    both_strands: bool):
    """
    Hamming-distance sliding window scanner.

    Negative strand strategy:
      At each position i, compute reverse_complement(window) and compare it
      against the forward motif — equivalent to asking whether the minus strand
      sequence at this position matches the motif within k substitutions.

    The motif string MUST be derived from the JASPAR file via
    consensus_from_counts() to guarantee that Hamming comparisons are
    bit-identical to those in sequence_matcher_with_motifs.py, ensuring
    peak hits are a strict subset of genome hits.

    Score = 1000 * (1 - mismatches / motif_length).

    Yields: (local_start, local_end, strand, bed_score, n_mismatches, matched_seq)

    NOTE: n_mismatches is yielded directly (not recalculated by the caller).
    This matches the yield signature of search_mismatch in
    sequence_matcher_with_motifs.py for consistency.
    """
    motif_len = len(motif)

    for i in range(len(seq) - motif_len + 1):
        window = seq[i:i + motif_len]
        if _has_n(window):
            continue

        # ── Forward strand ────────────────────────────────────────────────────
        mm_fwd = hamming_distance(window, motif)
        if mm_fwd <= max_mismatches:
            score = int(1000 * (1 - mm_fwd / motif_len))
            yield i, i + motif_len, '+', score, mm_fwd, window

        # ── Reverse strand ────────────────────────────────────────────────────
        # rc(window) is scored against the FORWARD motif.
        # This is equivalent to asking: does the minus-strand sequence at this
        # position (read 5'→3') match the motif within max_mismatches?
        if both_strands:
            rc_window = reverse_complement(window)
            mm_rev    = hamming_distance(rc_window, motif)
            if mm_rev <= max_mismatches:
                score = int(1000 * (1 - mm_rev / motif_len))
                yield i, i + motif_len, '-', score, mm_rev, rc_window


# ── Mode 3: PWM log-odds scoring ──────────────────────────────────────────────

def counts_to_log_odds(counts: dict, pseudocount: float = 0.1,
                        bg: dict = None) -> list:
    """
    Convert raw counts -> log2 odds weight matrix.
    Returns a list of per-position dicts [{A: score, C: score, ...}, ...]
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

def score_sequence(seq: str, pwm: list) -> float:
    """Sum log-odds scores for a fixed-length sequence against a PWM."""
    return sum(pwm[i].get(base, -10.0) for i, base in enumerate(seq))

def compute_score_distribution(pwm: list, n_samples: int = 100_000,
                                bg: dict = None) -> list:
    """
    Monte-Carlo null distribution under background frequencies.
    Returns a sorted list for bisect-based p-value lookup.
    """
    if bg is None:
        bg = {b: 0.25 for b in BASES}
    bases     = list(bg.keys())
    weights   = list(bg.values())
    motif_len = len(pwm)
    scores    = [
        score_sequence(
            ''.join(random.choices(bases, weights=weights, k=motif_len)), pwm
        )
        for _ in range(n_samples)
    ]
    scores.sort()
    return scores

def score_to_pvalue(score: float, null_dist: list) -> float:
    """Right-tail p-value: fraction of null scores >= observed score."""
    return 1.0 - bisect.bisect_left(null_dist, score) / len(null_dist)

def score_to_bed_score(score: float, pwm: list) -> int:
    """Rescale raw log-odds score to BED integer range [0, 1000]."""
    min_s = sum(min(pos.values()) for pos in pwm)
    max_s = sum(max(pos.values()) for pos in pwm)
    if max_s == min_s:
        return 500
    return max(0, min(1000, int((score - min_s) / (max_s - min_s) * 1000)))

def _annotate_mismatches(window: str, consensus: str):
    """Return (n_mismatches, list_of_0based_mismatch_positions)."""
    mm_pos = [i for i, (a, b) in enumerate(zip(window, consensus)) if a != b]
    return len(mm_pos), mm_pos


# ── PWM search: both strands + mismatch annotation ───────────────────────────

def search_pwm_with_mismatch_report(seq: str, pwm: list,
                                     pvalue_threshold: float,
                                     null_dist: list,
                                     both_strands: bool,
                                     consensus: str = None):
    """
    PWM scanner with explicit positive AND negative strand scoring.

    Negative strand strategy:
      At each position i, compute reverse_complement(window) and score it
      against the FORWARD PWM — mathematically equivalent to scoring the
      window against a reverse-complement PWM, but keeps a single PWM object.
      Mismatch annotation compares rc_window against the forward consensus
      since both are then in the same 5'->3' orientation relative to the motif.

    The consensus parameter MUST be derived from the same JASPAR file via
    consensus_from_counts() so mismatch annotations are consistent with
    sequence_matcher_with_motifs.py.

    Yields:
      (local_start, local_end, strand, bed_score,
       n_mismatches, mismatch_positions, raw_score, matched_seq)

    Coordinates are local to the passed sequence slice; callers add the
    genomic offset before writing BED output.
    """
    motif_len = len(pwm)

    for i in range(len(seq) - motif_len + 1):
        window = seq[i:i + motif_len]
        if _has_n(window):
            continue

        # ── Positive strand ──────────────────────────────────────────────────
        fwd_score = score_sequence(window, pwm)
        if score_to_pvalue(fwd_score, null_dist) <= pvalue_threshold:
            bed_score    = score_to_bed_score(fwd_score, pwm)
            n_mm, mm_pos = (_annotate_mismatches(window, consensus)
                            if consensus else (0, []))
            yield (i, i + motif_len, '+',
                   bed_score, n_mm, mm_pos, fwd_score, window)

        # ── Negative strand ──────────────────────────────────────────────────
        if both_strands:
            rc_window = reverse_complement(window)
            rev_score = score_sequence(rc_window, pwm)
            if score_to_pvalue(rev_score, null_dist) <= pvalue_threshold:
                bed_score    = score_to_bed_score(rev_score, pwm)
                # rc_window vs forward consensus: both 5'->3' on the motif
                n_mm, mm_pos = (_annotate_mismatches(rc_window, consensus)
                                if consensus else (0, []))
                yield (i, i + motif_len, '-',
                       bed_score, n_mm, mm_pos, rev_score, rc_window)


# ── Shared output helpers ─────────────────────────────────────────────────────

def _write_bed_line(fh, chrom, start, end, name, score, strand):
    fh.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")

def _write_tsv_pwm(fh, chrom, start, end, strand,
                   raw_score, bed_score, n_mm, mm_pos, matched_seq):
    mm_str = ','.join(str(p) for p in mm_pos) if mm_pos else '.'
    fh.write(
        f"{chrom}\t{start}\t{end}\t{strand}\t"
        f"{raw_score:.4f}\t{bed_score}\t{n_mm}\t{mm_str}\t{matched_seq}\n"
    )

def _write_tsv_basic(fh, chrom, start, end, strand, bed_score, n_mm, matched_seq):
    fh.write(
        f"{chrom}\t{start}\t{end}\t{strand}\t"
        f"{bed_score}\t{n_mm}\t{matched_seq}\n"
    )

def _print_summary(total, fwd, rev, filtered_mm=None, max_mm=None,
                   bed_path=None, tsv_path=None):
    print("=" * 60, file=sys.stderr)
    print(f"  Total hits         : {total:,}", file=sys.stderr)
    print(f"  Positive strand (+): {fwd:,}", file=sys.stderr)
    print(f"  Negative strand (-): {rev:,}", file=sys.stderr)
    if filtered_mm is not None:
        print(f"  Filtered (mm > {max_mm}) : {filtered_mm:,}", file=sys.stderr)
    if bed_path:
        print(f"  BED output         : {bed_path}", file=sys.stderr)
    if tsv_path:
        print(f"  Detail TSV         : {tsv_path}", file=sys.stderr)
    print("=" * 60, file=sys.stderr)


# ── PWM genome / chromosome / region scanner ─────────────────────────────────

def scan_genome_pwm_mismatch(fasta_path: str, output_bed: str,
                              pwm_path: str,
                              pvalue_threshold: float,
                              consensus: str = None,
                              max_mismatches: int = None,
                              both_strands: bool = True,
                              bg_freqs: dict = None,
                              pseudocount: float = 0.1,
                              n_samples: int = 100_000,
                              chroms: set = None,
                              region: tuple = None):
    """
    Complete PWM scanner supporting whole-genome, per-chromosome, and
    per-region modes.

    The consensus parameter MUST be derived via consensus_from_counts() from
    the same JASPAR file used to build the PWM, guaranteeing that mismatch
    annotations use the identical reference string as sequence_matcher_with_motifs.py.

    Target selection priority (mutually exclusive):
      region  -> scan only that genomic slice (coordinates preserved in BED)
      chroms  -> scan only those chromosomes
      (none)  -> scan the whole genome

    Writes:
      <output_bed>                  -- 6-column BED
      <output_bed stem>_details.tsv -- extended TSV with raw PWM score,
                                       mismatch annotation, and matched sequence
    """
    print("=" * 60, file=sys.stderr)
    print("Building PWM...", file=sys.stderr)
    counts    = parse_jaspar(pwm_path)
    # If caller did not supply consensus, derive it from the JASPAR counts
    if consensus is None:
        consensus = consensus_from_counts(counts)
        print(f"  [JASPAR] Consensus from file : {consensus}", file=sys.stderr)
    pwm       = counts_to_log_odds(counts, pseudocount=pseudocount, bg=bg_freqs)
    motif_len = len(pwm)

    print(f"  Motif length       : {motif_len} bp", file=sys.stderr)
    print(f"  Sampling null dist  ({n_samples:,} sequences)...", file=sys.stderr)
    null_dist = compute_score_distribution(pwm, n_samples=n_samples, bg=bg_freqs)

    cutoff_idx   = int((1 - pvalue_threshold) * len(null_dist))
    score_cutoff = null_dist[min(cutoff_idx, len(null_dist) - 1)]

    print(f"  p-value threshold  : {pvalue_threshold}", file=sys.stderr)
    print(f"  Log-odds cutoff    : {score_cutoff:.4f}", file=sys.stderr)
    if consensus:
        print(f"  Consensus          : {consensus}", file=sys.stderr)
    if max_mismatches is not None:
        print(f"  Hard mismatch cap  : <= {max_mismatches}", file=sys.stderr)
    print(f"  Both strands       : {both_strands}", file=sys.stderr)
    print(f"  Target             : {build_target_description(chroms, region)}",
          file=sys.stderr)
    print("=" * 60, file=sys.stderr)

    tsv_path       = str(Path(output_bed).with_suffix('')) + '_details.tsv'
    total_hits     = 0
    hits_fwd       = 0
    hits_rev       = 0
    filtered_by_mm = 0

    with open(output_bed, 'w') as bed_out, open(tsv_path, 'w') as tsv_out:
        tsv_out.write(
            "chrom\tstart\tend\tstrand\t"
            "pwm_score_raw\tbed_score\t"
            "n_mismatches\tmismatch_positions\tsequence\n"
        )

        for chrom, offset, seq in parse_fasta_filtered(
                fasta_path, chroms=chroms, region=region):

            print(f"  Scanning {chrom}  offset={offset:,}  "
                  f"len={len(seq):,} bp...", file=sys.stderr)
            chrom_fwd = 0
            chrom_rev = 0

            hits = search_pwm_with_mismatch_report(
                seq, pwm, pvalue_threshold, null_dist,
                both_strands, consensus
            )

            for (loc_s, loc_e, strand,
                 bed_score, n_mm, mm_pos,
                 raw_score, matched_seq) in hits:

                if max_mismatches is not None and n_mm > max_mismatches:
                    filtered_by_mm += 1
                    continue

                g_start = offset + loc_s
                g_end   = offset + loc_e
                name    = f"mm{n_mm}_{chrom}:{g_start}-{g_end}({strand})"

                _write_bed_line(bed_out, chrom, g_start, g_end,
                                name, bed_score, strand)
                _write_tsv_pwm(tsv_out, chrom, g_start, g_end, strand,
                               raw_score, bed_score, n_mm, mm_pos, matched_seq)
                total_hits += 1
                if strand == '+':
                    chrom_fwd += 1
                else:
                    chrom_rev += 1

            hits_fwd += chrom_fwd
            hits_rev += chrom_rev
            print(f"    {chrom}: {chrom_fwd} (+) / {chrom_rev} (-)",
                  file=sys.stderr)

    _print_summary(total_hits, hits_fwd, hits_rev,
                   filtered_by_mm, max_mismatches,
                   output_bed, tsv_path)


# ── Exact / mismatch genome / chromosome / region scanner ────────────────────

def scan_genome(fasta_path: str, output_bed: str, mode: str,
                motif: str,
                max_mismatches: int = 0,
                both_strands: bool = True,
                chroms: set = None,
                region: tuple = None):
    """
    Scanner for exact and mismatch modes with chromosome / region filtering.
    Writes a 6-column BED and a companion detail TSV.

    The motif parameter MUST be the consensus derived from the JASPAR file via
    consensus_from_counts() — identical to the string used in
    sequence_matcher_with_motifs.py — so that n_mismatches values and hit
    coordinates are bit-for-bit reproducible across both tools.

    search_mismatch now yields n_mismatches directly (no post-hoc
    recalculation), keeping the logic identical to the peak-search script.
    """
    tsv_path   = str(Path(output_bed).with_suffix('')) + '_details.tsv'
    total_hits = 0
    hits_fwd   = 0
    hits_rev   = 0
    label      = motif if mode == 'exact' else f"{motif}_mm{max_mismatches}"

    print("=" * 60, file=sys.stderr)
    print(f"  Mode               : {mode}", file=sys.stderr)
    print(f"  Motif (consensus)  : {motif}", file=sys.stderr)
    if mode == 'mismatch':
        print(f"  Max mismatches     : {max_mismatches}", file=sys.stderr)
    print(f"  Both strands       : {both_strands}", file=sys.stderr)
    print(f"  Target             : {build_target_description(chroms, region)}",
          file=sys.stderr)
    print("=" * 60, file=sys.stderr)

    with open(output_bed, 'w') as bed_out, open(tsv_path, 'w') as tsv_out:
        tsv_out.write(
            "chrom\tstart\tend\tstrand\t"
            "bed_score\tn_mismatches\tsequence\n"
        )

        for chrom, offset, seq in parse_fasta_filtered(
                fasta_path, chroms=chroms, region=region):

            print(f"  Scanning {chrom}  offset={offset:,}  "
                  f"len={len(seq):,} bp...", file=sys.stderr)
            chrom_fwd = 0
            chrom_rev = 0

            hits = (search_exact(seq, motif, both_strands)
                    if mode == 'exact'
                    else search_mismatch(seq, motif, max_mismatches, both_strands))

            # Both search_exact and search_mismatch now yield:
            #   (local_start, local_end, strand, bed_score, n_mismatches, matched_seq)
            # n_mismatches is yielded directly — no post-hoc recalculation needed.
            for loc_s, loc_e, strand, score, n_mm, matched_seq in hits:
                g_start = offset + loc_s
                g_end   = offset + loc_e
                name    = f"mm{n_mm}_{label}_{chrom}:{g_start}-{g_end}({strand})"

                _write_bed_line(bed_out, chrom, g_start, g_end,
                                name, score, strand)
                _write_tsv_basic(tsv_out, chrom, g_start, g_end,
                                 strand, score, n_mm, matched_seq)
                total_hits += 1
                if strand == '+':
                    chrom_fwd += 1
                else:
                    chrom_rev += 1

            hits_fwd += chrom_fwd
            hits_rev += chrom_rev
            print(f"    {chrom}: {chrom_fwd} (+) / {chrom_rev} (-)",
                  file=sys.stderr)

    _print_summary(total_hits, hits_fwd, hits_rev,
                   bed_path=output_bed, tsv_path=tsv_path)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Motif scanner: exact regex, Hamming mismatch, or PWM mode.\n"
            "\n"
            "JASPAR FILE IS THE SOURCE OF TRUTH\n"
            "  Provide --jaspar in any mode to derive the canonical consensus\n"
            "  from the counts matrix.  This guarantees bit-exact agreement\n"
            "  with sequence_matcher_with_motifs.py so that peak hits are a\n"
            "  strict subset of genome hits.\n"
            "\n"
            "Target selection (mutually exclusive):\n"
            "  --chroms chr1 chr2   scan only the listed chromosomes\n"
            "  --region chr1:S-E    scan one region (1-based, inclusive)\n"
            "  (neither)            scan the whole genome  [default]\n"
            "\n"
            "All modes scan both strands by default.\n"
            "Negative strand hits use + strand coordinates (BED convention).\n"
            "The matched sequence in the TSV is shown 5'->3' on the hit strand."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-g', '--genome', required=True,
                        help='Input genome FASTA file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output BED file path')
    parser.add_argument('--mode', choices=['exact', 'mismatch', 'pwm'],
                        default='exact',
                        help='Search mode (default: exact)')

    # ── Motif / JASPAR ───────────────────────────────────────────────────────
    parser.add_argument('-m', '--motif',
                        help='Consensus motif string, e.g. AGGGGATTTCCC. '
                             'When --jaspar is also given, this is validated '
                             'against the JASPAR consensus and a warning is '
                             'printed if they differ; the JASPAR consensus takes '
                             'precedence. Required only when --jaspar is absent.')
    # --jaspar is the primary JASPAR flag (consistent with sequence_matcher_with_motifs.py).
    # --pwm is kept as a backward-compatible alias.
    parser.add_argument('--jaspar', '--pwm', dest='jaspar',
                        help='Path to JASPAR-format counts file. '
                             'Required for --mode pwm. '
                             'Also accepted for exact/mismatch to derive '
                             'the canonical consensus string and guarantee '
                             'subset consistency with '
                             'sequence_matcher_with_motifs.py.')
    parser.add_argument('--consensus',
                        help='Override consensus for mismatch annotation in '
                             'pwm mode. Falls back to JASPAR consensus if omitted.')

    # ── Mismatch ─────────────────────────────────────────────────────────────
    parser.add_argument('--max-mismatches', type=int, default=None,
                        help='Mismatch mode: Hamming cutoff (default 1 if unset). '
                             'PWM mode: hard post-filter after p-value threshold '
                             '(None = no cap).')

    # ── PWM calibration ───────────────────────────────────────────────────────
    parser.add_argument('--pvalue',      type=float, default=1e-4,
                        help='p-value threshold for PWM hits (default: 1e-4)')
    parser.add_argument('--pseudocount', type=float, default=0.1,
                        help='Pseudocount added to each PWM cell (default: 0.1)')
    parser.add_argument('--n-samples',   type=int,   default=100_000,
                        help='Null distribution sample count (default: 100000)')
    parser.add_argument('--bg-gc',       type=float, default=0.5,
                        help='Background GC fraction e.g. 0.42 for human '
                             '(default: 0.5 = uniform)')

    # ── Strand ────────────────────────────────────────────────────────────────
    parser.add_argument('--fwd-only', action='store_true',
                        help='Scan forward (+) strand only '
                             '(default: both strands)')

    # ── Target selection (mutually exclusive) ─────────────────────────────────
    target = parser.add_mutually_exclusive_group()
    target.add_argument(
        '--chroms', nargs='+', metavar='CHROM',
        help='One or more chromosome names to scan, e.g. --chroms chr1 chr2 chrX. '
             'Names must match FASTA header IDs exactly. '
             'Mutually exclusive with --region.'
    )
    target.add_argument(
        '--region', metavar='CHROM:START-END',
        help='Single genomic region in 1-based inclusive coordinates, '
             'e.g. --region chr1:1000000-2000000. '
             'Mutually exclusive with --chroms.'
    )

    args = parser.parse_args()

    # ── Validate argument combinations ────────────────────────────────────────
    if args.mode == 'pwm' and not args.jaspar:
        parser.error('--jaspar is required for pwm mode')
    if args.mode in ('exact', 'mismatch') and not args.motif and not args.jaspar:
        parser.error(
            '--jaspar or --motif / -m is required for exact and mismatch modes.\n'
            'Prefer --jaspar: the consensus is derived from the counts matrix,\n'
            'guaranteeing subset consistency with sequence_matcher_with_motifs.py.'
        )
    if args.mode == 'mismatch' and args.max_mismatches is None:
        args.max_mismatches = 1
        print("Note: --max-mismatches not set; defaulting to 1 for mismatch mode",
              file=sys.stderr)

    # ── Parse target ──────────────────────────────────────────────────────────
    chroms = set(args.chroms) if args.chroms else None
    region = None
    if args.region:
        try:
            region = parse_region(args.region)
        except ValueError as exc:
            parser.error(str(exc))

    # ── Background frequencies ────────────────────────────────────────────────
    gc = args.bg_gc / 2
    at = (1.0 - args.bg_gc) / 2
    bg = {'A': at, 'C': gc, 'G': gc, 'T': at}

    both_strands = not args.fwd_only

    # ── Resolve canonical motif / consensus from JASPAR ───────────────────────
    # consensus_from_counts() is the single authoritative function for deriving
    # the consensus from a JASPAR file.  Using the same function and the same
    # JASPAR file as sequence_matcher_with_motifs.py guarantees that every
    # Hamming comparison in mismatch mode uses the IDENTICAL reference string,
    # which is the necessary condition for peak hits being a strict subset of
    # genome hits.
    motif, jaspar_counts = resolve_motif(
        jaspar_path = args.jaspar,
        cli_motif   = args.motif,
        pseudocount = args.pseudocount,
        bg          = bg,
    )

    # Consensus used for mismatch annotation in PWM mode:
    # --consensus arg overrides; otherwise use the JASPAR-derived motif.
    consensus = normalize(args.consensus) if args.consensus else motif

    # ── Dispatch ──────────────────────────────────────────────────────────────
    if args.mode == 'pwm':
        scan_genome_pwm_mismatch(
            fasta_path       = args.genome,
            output_bed       = args.output,
            pwm_path         = args.jaspar,
            pvalue_threshold = args.pvalue,
            consensus        = consensus,
            max_mismatches   = args.max_mismatches,
            both_strands     = both_strands,
            bg_freqs         = bg,
            pseudocount      = args.pseudocount,
            n_samples        = args.n_samples,
            chroms           = chroms,
            region           = region,
        )
    else:
        scan_genome(
            fasta_path     = args.genome,
            output_bed     = args.output,
            mode           = args.mode,
            motif          = motif,          # JASPAR-derived consensus
            max_mismatches = args.max_mismatches if args.max_mismatches is not None else 0,
            both_strands   = both_strands,
            chroms         = chroms,
            region         = region,
        )


if __name__ == '__main__':
    main()
