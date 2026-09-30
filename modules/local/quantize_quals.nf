process QUANTIZE_QUALS {
    tag "${meta.sample}"
    container params.quantize_quals_container
    cpus 8
    memory '8 GB'

    input:
    tuple val(meta), path(bam, stageAs: 'input/*')
    tuple path(ref_fasta), path(ref_fasta_fai)
    val output_fmt

    output:
    tuple val(meta), path("${meta.sample}.recal.${output_fmt}"), path("${meta.sample}.recal.${output_fmt}.{bai,crai}")

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail
    quantize_quals \\
        --in ${bam} \\
        --out ${meta.sample}.recal.${output_fmt} \\
        --ref ${ref_fasta} \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    def idx = output_fmt == 'cram' ? 'crai' : 'bai'
    """
    touch ${meta.sample}.recal.${output_fmt} ${meta.sample}.recal.${output_fmt}.${idx}
    """
}
