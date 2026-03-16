#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
import logomaker
import matplotlib.pyplot as plt
from collections import defaultdict

BASES = ["A", "C", "G", "T"]

# ---------------- Utilities ----------------

def normalize(seq):
    return seq.upper().replace("U", "T")

def rc(seq):
    return str(Seq(seq).reverse_complement())

def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))

def extract_flanks(seq, pos, k, up, down):
    return (
        seq[max(0, pos - up):pos],
        seq[pos + k:pos + k + down]
    )

# ---------------- Matching ----------------

def align_match(seq, query, max_mm):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1

    hits = []
    k = len(query)

    for i in range(len(seq) - k + 1):
        window = seq[i:i+k]
        aln = aligner.align(query, window)
        mm = hamming(query, window)
        if aln and mm <= max_mm:
            hits.append((window, mm, i, aln[0].score))
    return hits

# ---------------- Motif functions ----------------

def pwm_from_seqs(seqs):
    mat = np.zeros((len(seqs[0]), 4))
    for s in seqs:
        for i, b in enumerate(s):
            if b in BASES:
                mat[i, BASES.index(b)] += 1
    mat /= mat.sum(axis=1, keepdims=True)
    return pd.DataFrame(mat, columns=BASES)

def info_content(pwm):
    return pwm * np.log2((pwm + 1e-9) / 0.25)

def logo(pwm, outfile, title):
    plt.figure(figsize=(max(6, len(pwm)/2), 3))
    logomaker.Logo(pwm)
    plt.title(title)
    plt.ylabel("Bits")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()

def write_meme(pwm, name, outfile):
    with open(outfile, "w") as f:
        f.write("MEME version 4\n\nALPHABET= ACGT\n\n")
        f.write(f"MOTIF {name}\n")
        f.write(f"letter-probability matrix: alength=4 w={len(pwm)}\n")
        for _, row in pwm.iterrows():
            f.write(" ".join(f"{v:.4f}" for v in row) + "\n")

# ---------------- Main ----------------

def main(args):
    df = pd.read_csv(args.input, sep="\t")
    seq_col = df.columns[-1]
    query = normalize(args.query)

    flanks = defaultdict(list)
    out = []

    for _, r in df.iterrows():
        seq = normalize(str(r[seq_col]))

        for strand, s in [("+", seq), ("-", rc(seq))]:
            hits = align_match(s, query, args.mismatches)

            for w, mm, pos, score in hits:
                up, down = extract_flanks(s, pos, len(query),
                                          args.upstream, args.downstream)

                if strand == "-":
                    w, up, down = rc(w), rc(down), rc(up)

                flanks[f"{strand}_up"].append(up)
                flanks[f"{strand}_down"].append(down)

                row = r.to_dict()
                row.update({
                    "Strand": strand,
                    "Matched_Sequence": w,
                    "Mismatches": mm,
                    "Position": pos,
                    "Alignment_Score": score,
                    "Upstream": up,
                    "Downstream": down
                })
                out.append(row)

    pd.DataFrame(out).to_csv(args.output, sep="\t", index=False)

    # ---- Motifs ----
    for key, seqs in flanks.items():
        if not seqs:
            continue
        pwm = pwm_from_seqs(seqs)
        ic = info_content(pwm)

        logo(ic, f"{key}_logo.png", key)
        pwm.to_csv(f"{key}_pssm.tsv", sep="\t", index=False)
        write_meme(pwm, key, f"{key}.meme")

# ---------------- CLI ----------------

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input", required=True)
    p.add_argument("-o", "--output", required=True)
    p.add_argument("-q", "--query", required=True)
    p.add_argument("-m", "--mismatches", type=int, default=0)
    p.add_argument("--upstream", type=int, default=50)
    p.add_argument("--downstream", type=int, default=50)
    args = p.parse_args()
    main(args)
