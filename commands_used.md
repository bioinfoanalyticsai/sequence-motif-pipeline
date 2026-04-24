# Command used in testing and running SeqMatcher Pipeline

### be sure to use `./build.sh` to build the docker container


## What run main script to fetch the sequence and find the motif with the pipeline

```bash
nextflow run main.nf \                
    --jaspar     ../../matching_logic/matching_logic_scripts/AGGGGATTTCCC.jaspar \
    --mode       mismatch \
    --mismatches 2 \
    --genome     input/genome/mm9_genome.fa \
    -profile     docker
```

### Manual run with fetched sequence file only

```bash
python sequence_matcher_with_motifs.py \
    -i ../../pipeline_fixed/results/01_sequences/p651hcombined_with_seq.tsv -o p651hcombined_peak_hits.tsv \    --jaspar AGGGGATTTCCC.jaspar \
    --mode mismatch --max-mismatches 2
```

## Cluster upstream 150 and downstream 150 base pairs

```bash
python cluster_sequences_updated.py p651hcombined_peak_hits.tsv --n-clusters 7 --output-prefix p651hcombined_cluster_7_re
```
## make cluster beds

```bash
python make_cluster_beds.py --input p651hcombined_cluster_7_re_results.tsv --prefix p651hcombined_cluster_7_re_beds --mode combined
```

## Command used to get the motifs from mm9 fasta genome file:
```bash
python find_motif_both_strands_and_chromsome_wise.py \
    -g genome.fa --jaspar AGGGGATTTCCC.jaspar \
    -o genome_hits.bed --mode mismatch --max-mismatches 2
```
