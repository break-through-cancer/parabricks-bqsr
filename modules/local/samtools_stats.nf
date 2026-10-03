process SAMTOOLS_STATS {
    tag "${meta.sample}"
    container 'community.wave.seqera.io/library/htslib_samtools@sha256:a55ddea590e567a91df592300a960aa534cfc1bd16e7623e3938ec21f4f3df15'
    cpus 4
    memory '4 GB'

    input:
    tuple val(meta), path(alignment), path(index)
    tuple path(ref_fasta), path(ref_fasta_fai)

    output:
    tuple val(meta), path("${meta.sample}.stats"), emit: stats

    script:
    """
    samtools stats -@ ${task.cpus} --reference ${ref_fasta} ${alignment} > ${meta.sample}.stats
    """

    stub:
    """
    printf 'SN\\traw total sequences:\\t0\\nSN\\tnon-primary alignments:\\t0\\nSN\\tsupplementary alignments:\\t0\\n' > ${meta.sample}.stats
    """
}
