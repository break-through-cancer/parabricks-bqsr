process SAMTOOLS_FLAGSTAT {
    tag "${meta.sample}"
    container 'community.wave.seqera.io/library/htslib_samtools@sha256:a55ddea590e567a91df592300a960aa534cfc1bd16e7623e3938ec21f4f3df15'
    cpus 2
    memory '2 GB'

    input:
    tuple val(meta), path(alignment)
    tuple path(ref_fasta), path(ref_fasta_fai)

    output:
    tuple val(meta), path("${meta.sample}.flagstat"), emit: flagstat

    script:
    """
    samtools flagstat -@ ${task.cpus} --input-fmt-option reference=${ref_fasta} -O tsv ${alignment} > ${meta.sample}.flagstat
    """

    stub:
    """
    printf '0\\t0\\tprimary\\n0\\t0\\tsecondary\\n0\\t0\\tsupplementary\\n' > ${meta.sample}.flagstat
    """
}
