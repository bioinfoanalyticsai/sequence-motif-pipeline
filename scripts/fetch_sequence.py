#!/usr/bin/env python3

import pandas as pd
from pyfaidx import Fasta
import argparse

def fetch_sequences(genome_fasta, input_csv, output_file):
    # Load genome fasta
    genome = Fasta(genome_fasta)

    # Load input CSV
    df = pd.read_csv(input_csv)

    sequences = []
    for idx, row in df.iterrows():
        chrom = str(row['Chromosome'])
        start = int(row['Peak_Start'])
        end = int(row['Peak_End'])
        strand = row['Strand']

        # pyfaidx uses 0-based indexing internally
        seq = genome[chrom][start:end].seq

        # Reverse complement if strand is '-'
        if strand == '-':
            seq = reverse_complement(seq)

        sequences.append(seq)

    # Add sequences to DataFrame
    df['Fetched_sequence'] = sequences

    # Write output as tab-separated file
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Output written to {output_file}")

def reverse_complement(seq):
    complement = str.maketrans('ATCGatcg', 'TAGCtagc')
    return seq.translate(complement)[::-1]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch sequences from genome FASTA based on peaks")
    parser.add_argument('-g', '--genome', required=True, help='Genome FASTA file')
    parser.add_argument('-i', '--input', required=True, help='Input CSV file with GeneID,Chromosome,Peak_Start,Peak_End,Strand')
    parser.add_argument('-o', '--output', required=True, help='Output tab-separated file')
    args = parser.parse_args()

    fetch_sequences(args.genome, args.input, args.output)
