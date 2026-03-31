nextflow.enable.dsl=2

params.input_dir   = "input"
params.query       = "AGGGAATTTCCC"
params.mismatches  = 2
params.upstream    = 150
params.downstream  = 150
params.outdir      = "results"
params.use_docker  = true
params.genome      = "input/genome/mm9_genome.fa"

workflow {

    genome_ch  = Channel.fromPath(params.genome)
    samples_ch = Channel.fromPath("${params.input_dir}/*.csv")

    // Pair every sample with the genome file
    FETCH_SEQUENCES(samples_ch.combine(genome_ch))

    // Pass FETCH_SEQUENCES output into MATCH_AND_MOTIFS
    MATCH_AND_MOTIFS(FETCH_SEQUENCES.out)
}

process FETCH_SEQUENCES {

    tag { sample_file.baseName }
    container params.use_docker ? 'seqmatcher:latest' : null

    input:
    tuple path(sample_file), path(genome_file)   // ← combined tuple input

    output:
    path "${sample_file.baseName}_with_seq.tsv"

    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    python /usr/local/bin/scripts/fetch_sequence.py \
        -i ${sample_file} \
        -g ${genome_file} \
        -o ${sample_file.baseName}_with_seq.tsv
    """
}

process MATCH_AND_MOTIFS {

    tag { sample_file.baseName }
    container params.use_docker ? 'seqmatcher:latest' : null

    input:
    path sample_file

    output:
    path "*.tsv"
    path "*.png"
    path "*.meme"

    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    python /usr/local/bin/scripts/sequence_matcher_with_motifs.py \
        -i ${sample_file} \
        -o ${sample_file.baseName}_matches.tsv \
        -q ${params.query} \
        -m ${params.mismatches} \
        --upstream ${params.upstream} \
        --downstream ${params.downstream}
    """
}
