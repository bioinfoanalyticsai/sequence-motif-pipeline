nextflow run main.nf \                
    --jaspar     ../../matching_logic/matching_logic_scripts/AGGGGATTTCCC.jaspar \
    --mode       mismatch \
    --mismatches 2 \
    --genome     input/genome/mm9_genome.fa \
    -profile     docker


python sequence_matcher_with_motifs.py \
    -i ../../pipeline_fixed/results/01_sequences/p651hcombined_with_seq.tsv -o p651hcombined_peak_hits.tsv \    --jaspar AGGGGATTTCCC.jaspar \
    --mode mismatch --max-mismatches 2


python cluster_sequences_updated.py p651hcombined_peak_hits.tsv --n-clusters 7 --output-prefix p651hcombined_cluster_7_re

python make_cluster_beds.py --input p651hcombined_cluster_7_re_results.tsv --prefix p651hcombined_cluster_7_re_beds --mode combined

Command used to get the motifs from mm9 fasta genome file:
python find_motif_both_strands_and_chromsome_wise.py \
    -g genome.fa --jaspar AGGGGATTTCCC.jaspar \
    -o genome_hits.bed --mode mismatch --max-mismatches 2
