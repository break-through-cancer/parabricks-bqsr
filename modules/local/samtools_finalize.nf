process SAMTOOLS_FINALIZE {
    tag "${meta.sample}"
    container 'community.wave.seqera.io/library/htslib_samtools@sha256:a55ddea590e567a91df592300a960aa534cfc1bd16e7623e3938ec21f4f3df15'
    cpus 8
    memory '8 GB'

    input:
    tuple val(meta), path(bam)
    tuple path(ref_fasta), path(ref_fasta_fai)
    val output_fmt

    output:
    tuple val(meta), path("${meta.sample}.recal.${output_fmt}"), path("${meta.sample}.recal.${output_fmt}.{bai,crai}"), emit: alignment

    script:
    if (output_fmt == 'bam') {
        """
        set -euo pipefail
        if [ "${bam}" != "${meta.sample}.recal.bam" ]; then
            ln -s "${bam}" "${meta.sample}.recal.bam"
        fi
        samtools index -@ ${task.cpus} ${meta.sample}.recal.bam
        """
    } else {
        """
        set -euo pipefail
        samtools view -@ ${task.cpus} -C -T ${ref_fasta} --write-index \\
            -o ${meta.sample}.recal.cram ${bam}
        """
    }

    stub:
    def idx = output_fmt == 'cram' ? 'crai' : 'bai'
    """
    touch ${meta.sample}.recal.${output_fmt} ${meta.sample}.recal.${output_fmt}.${idx}
    """
}
