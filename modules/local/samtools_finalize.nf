process SAMTOOLS_FINALIZE {
    tag "${meta.sample}"
    container 'community.wave.seqera.io/library/htslib_samtools@sha256:a55ddea590e567a91df592300a960aa534cfc1bd16e7623e3938ec21f4f3df15'
    cpus 8
    memory '8 GB'

    input:
    tuple val(meta), path(alignment, stageAs: 'input/*')
    tuple path(ref_fasta), path(ref_fasta_fai)
    val output_fmt

    output:
    tuple val(meta), path("${meta.sample}.${meta.suffix}.${output_fmt}"), path("${meta.sample}.${meta.suffix}.${output_fmt}.${output_fmt == 'cram' ? 'crai' : 'bai'}"), emit: alignment

    script:
    def out = "${meta.sample}.${meta.suffix}.${output_fmt}"
    def in_fmt = alignment.name.endsWith('.cram') ? 'cram' : 'bam'
    if (in_fmt == output_fmt) {
        """
        set -euo pipefail
        ln -s ${alignment} ${out}
        samtools index -@ ${task.cpus} ${out}
        """
    } else {
        """
        set -euo pipefail
        samtools view -@ ${task.cpus} ${output_fmt == 'cram' ? '-C' : '-b'} -T ${ref_fasta} --write-index -o ${out}##idx##${out}.${output_fmt == 'cram' ? 'crai' : 'bai'} ${alignment}
        """
    }

    stub:
    def idx = output_fmt == 'cram' ? 'crai' : 'bai'
    """
    touch ${meta.sample}.${meta.suffix}.${output_fmt} ${meta.sample}.${meta.suffix}.${output_fmt}.${idx}
    """
}
