process FASTP {
    tag "${meta.id}"
    container 'quay.io/biocontainers/fastp@sha256:99eb308e4c4f1a6467beb775019bce0f88414303395fa8f27ad6da9a014bf14b'
    cpus 16
    memory { "${16 * task.attempt} GB" }

    input:
    tuple val(meta), path(reads, stageAs: 'input/*')
    val trim

    output:
    tuple val(meta), path("${meta.id}.fastp.json"), emit: json
    tuple val(meta), path("${meta.id}.fastp.html"), emit: html
    tuple val(meta), path("*.trimmed.fastq.gz"), emit: reads, optional: true

    script:
    def args = fastpArgs(meta, reads, [
        trim: trim, clip_r1: params.clip_r1, clip_r2: params.clip_r2,
        three_prime_clip_r1: params.three_prime_clip_r1, three_prime_clip_r2: params.three_prime_clip_r2,
        poly_g: params.poly_g_trimming, length_required: params.length_required, threads: task.cpus
    ])
    """
    set -euo pipefail
    fastp ${args}
    """

    stub:
    def n = meta.single_end ? 1 : 2
    def outs = fastpWritesReads(trim, params.poly_g_trimming) ? (1..n).collect { i -> "${meta.id}_${i}.trimmed.fastq.gz" } : []
    """
    echo '{"summary":{"before_filtering":{"total_reads":0},"after_filtering":{"total_reads":0}}}' > ${meta.id}.fastp.json
    touch ${meta.id}.fastp.html ${outs.join(' ')}
    """
}

def fastpArgs(Map meta, Object reads, Map opts) {
    def r = (reads instanceof List ? reads : [reads]).collect { it.toString() }
    def a = ["-i ${r[0]}"]
    if (!meta.single_end) a << "-I ${r[1]}"
    if (fastpWritesReads(opts.trim, opts.poly_g)) {
        a << "-o ${meta.id}_1.trimmed.fastq.gz"
        if (!meta.single_end) a << "-O ${meta.id}_2.trimmed.fastq.gz"
        a << '--disable_quality_filtering'
        if (opts.poly_g == 'on') a << '--trim_poly_g'
        if (opts.poly_g == 'off') a << '--disable_trim_poly_g'
        if (opts.trim) {
            if (!meta.single_end) a << '--detect_adapter_for_pe'
            if (opts.clip_r1) a << "--trim_front1 ${opts.clip_r1}"
            if (opts.three_prime_clip_r1) a << "--trim_tail1 ${opts.three_prime_clip_r1}"
            if (!meta.single_end && opts.clip_r2) a << "--trim_front2 ${opts.clip_r2}"
            if (!meta.single_end && opts.three_prime_clip_r2) a << "--trim_tail2 ${opts.three_prime_clip_r2}"
        } else {
            a << '--disable_adapter_trimming'
        }
        a << "--length_required ${opts.length_required}"
        a << '-z 1'
    } else {
        a += ['--disable_quality_filtering', '--disable_adapter_trimming', '--disable_trim_poly_g', '--disable_length_filtering']
    }
    a << "--thread ${opts.threads}"
    a << "-j ${meta.id}.fastp.json -h ${meta.id}.fastp.html"
    a.join(' ')
}

def fastpWritesReads(Object trim, Object polyG) {
    trim.toString() == 'true' || polyG.toString() != 'off'
}
