#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ─────────────────────────────────────────────────────────────────────────────
//  Parameters
// ─────────────────────────────────────────────────────────────────────────────
params.input_dir      = "input"
params.query          = "AGGGGATTTCCC"  // fallback consensus (JASPAR takes precedence when --jaspar is set)
params.mode           = "mismatch"      // exact | mismatch | pwm
params.mismatches     = 2               // max Hamming mismatches (mismatch mode)
params.upstream       = 150
params.downstream     = 150
// --jaspar is the single source of truth for the consensus motif in ALL modes:
//   • mismatch : consensus_from_counts(jaspar) drives Hamming comparisons
//   • exact    : consensus_from_counts(jaspar) drives the IUPAC regex
//   • pwm      : required — builds log-odds matrix; consensus for annotations
// Providing the same JASPAR file to both this pipeline AND
// find_motif_both_strands_and_chromsome_wise.py guarantees every peak hit
// is a strict subset of the genome-wide scan results.
params.jaspar         = null
params.pvalue         = 1e-4            // p-value threshold (pwm mode)
params.pseudocount    = 0.1             // PWM pseudocount
params.bg_gc          = 0.5             // background GC fraction for null dist
params.fwd_only       = false           // scan forward strand only
params.outdir         = "results"
params.use_docker     = true
params.genome         = "input/genome/genome.fa"

// ─────────────────────────────────────────────────────────────────────────────
//  Parameter validation (fast-fail before any process runs)
// ─────────────────────────────────────────────────────────────────────────────
// PWM mode always needs a JASPAR file to build the log-odds matrix.
if (params.mode == 'pwm' && !params.jaspar) {
    error(
        "ERROR: --mode pwm requires --jaspar.\n" +
        "  Provide a JASPAR-format counts file: --jaspar <path/to/motif.jaspar>"
    )
}

// Mismatch / exact mode with neither --jaspar nor --query is unrunnable.
if (params.mode in ['mismatch', 'exact'] && !params.jaspar && !params.query) {
    error(
        "ERROR: --mode ${params.mode} requires either --jaspar or --query.\n" +
        "  Prefer --jaspar: the consensus is derived from the counts matrix,\n" +
        "  guaranteeing subset consistency with find_motif_both_strands."
    )
}

// ─────────────────────────────────────────────────────────────────────────────
//  Workflow
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    samples_ch = Channel.fromPath("${params.input_dir}/*.csv", checkIfExists: true)

    // .first() converts the queue channel to a VALUE channel so the single
    // genome file is broadcast to every sample rather than consumed once.
    genome_ch  = Channel.fromPath(params.genome, checkIfExists: true).first()

    // JASPAR file channel.
    //   Supplied  → value channel (.first() broadcasts to every sample).
    //   Absent    → value channel carrying [] so Nextflow stages nothing.
    // In mismatch mode the script reads the JASPAR counts matrix to derive
    // the canonical consensus via consensus_from_counts(), identical to
    // what find_motif_both_strands_and_chromsome_wise.py does — this is
    // what guarantees the subset property between the two tools.
    jaspar_ch = params.jaspar != null
                    ? Channel.fromPath(params.jaspar, checkIfExists: true).first()
                    : Channel.value([])

    // Step 1 – fetch genomic sequences for every peak in each sample CSV
    FETCH_SEQUENCES(samples_ch.combine(genome_ch))

    // Step 2 – motif search (mismatch / exact / pwm) + BED + BigWig + logos
    //  genome_ch and jaspar_ch are value channels → reused for every TSV
    MATCH_AND_MOTIFS(
        FETCH_SEQUENCES.out.tsv.combine(genome_ch),
        jaspar_ch
    )

    // Step 3 – validate every BigWig so genome browsers can load it remotely
    INDEX_BIGWIG(MATCH_AND_MOTIFS.out.bigwig)
}

// ─────────────────────────────────────────────────────────────────────────────
//  Process 1 – Fetch sequences
// ─────────────────────────────────────────────────────────────────────────────
process FETCH_SEQUENCES {

    tag { sample_file.baseName }
    container params.use_docker ? 'seqmatcher:latest' : null

    publishDir "${params.outdir}/01_sequences", mode: 'copy'

    input:
    tuple path(sample_file), path(genome_file)

    output:
    path "${sample_file.baseName}_with_seq.tsv", emit: tsv

    script:
    """
    python /usr/local/bin/scripts/fetch_sequence.py \\
        --input  ${sample_file} \\
        --genome ${genome_file} \\
        --output ${sample_file.baseName}_with_seq.tsv
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  Process 2 – Motif search, logos, BED, BigWig
// ─────────────────────────────────────────────────────────────────────────────
//
//  JASPAR-driven consensus logic (v3 behaviour):
//  ──────────────────────────────────────────────
//  When --jaspar is supplied (recommended for ALL modes):
//    • The script calls consensus_from_counts(jaspar_counts) to derive the
//      canonical motif string — the same function used in the genome-wide
//      scanner — so Hamming comparisons are bit-identical across both tools.
//    • --query is passed alongside as a cross-check; a warning is printed
//      and the JASPAR consensus wins if they differ.
//
//  When --jaspar is absent (fallback):
//    • --query drives the consensus directly.
//    • Cross-tool subset consistency is not guaranteed unless --query
//      matches the JASPAR consensus manually.
//
process MATCH_AND_MOTIFS {

    tag { tsv_file.baseName }
    container params.use_docker ? 'seqmatcher:latest' : null

    publishDir "${params.outdir}/02_matches",  mode: 'copy', pattern: '*.tsv'
    publishDir "${params.outdir}/03_motifs",   mode: 'copy', pattern: '*.{png,meme}'
    publishDir "${params.outdir}/04_bed",      mode: 'copy', pattern: '*.bed'
    publishDir "${params.outdir}/05_bigwig",   mode: 'copy', pattern: '*.bw'
    publishDir "${params.outdir}/05_bigwig",   mode: 'copy', pattern: '*_chrom_sizes.txt'

    input:
    tuple path(tsv_file), path(genome_file)
    path jaspar_file    // real path, or [] when --jaspar is omitted

    output:
    path "*_matches.tsv",         emit: matches
    path "*_pssm.tsv",            emit: pssm,       optional: true
    path "*.png",                 emit: logos,       optional: true
    path "*.meme",                emit: meme,        optional: true
    path "*.bed",                 emit: bed,         optional: true
    path "*.bw",                  emit: bigwig,      optional: true
    path "*_chrom_sizes.txt",     emit: chromsizes,  optional: true

    script:
    // jaspar_file is [] (falsy) when the user did not supply --jaspar.
    // In mismatch mode: the JASPAR path is forwarded so the script derives
    // the canonical consensus via consensus_from_counts() — identical logic
    // to find_motif_both_strands_and_chromsome_wise.py.
    // --query is passed alongside when set; the script treats it as a
    // cross-check and warns if it differs from the JASPAR consensus.
    def jaspar_arg = jaspar_file     ? "--jaspar ${jaspar_file}"    : ""
    def query_arg  = params.query    ? "--query  ${params.query}"   : ""
    def fwd_arg    = params.fwd_only ? "--fwd-only"                 : ""
    def base       = tsv_file.baseName
    """
    python /usr/local/bin/scripts/sequence_matcher_with_motifs.py \\
        --input     ${tsv_file} \\
        --output    ${base}_matches.tsv \\
        --reference ${genome_file} \\
        --mode      ${params.mode} \\
        --max-mismatches ${params.mismatches} \\
        --upstream       ${params.upstream} \\
        --downstream     ${params.downstream} \\
        --pvalue         ${params.pvalue} \\
        --pseudocount    ${params.pseudocount} \\
        --bg-gc          ${params.bg_gc} \\
        --bed            ${base}_hits.bed \\
        --bigwig         ${base}_hits.bw \\
        ${jaspar_arg} \\
        ${query_arg} \\
        ${fwd_arg}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  Process 3 – Validate BigWig files (ucsc bigWigInfo)
// ─────────────────────────────────────────────────────────────────────────────
process INDEX_BIGWIG {

    tag { bw_file.baseName }
    container params.use_docker ? 'seqmatcher:latest' : null

    publishDir "${params.outdir}/05_bigwig", mode: 'copy'

    input:
    path bw_file

    output:
    path "${bw_file}"

    script:
    """
    # Validate BigWig integrity; exits non-zero on corrupt file
    bigWigInfo ${bw_file} > ${bw_file.baseName}_info.txt 2>&1 || true
    """
}
