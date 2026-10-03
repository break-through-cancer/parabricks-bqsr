process MOSDEPTH {
    tag "${meta.sample}"
    container 'quay.io/biocontainers/mosdepth@sha256:94b54a6a5610da01030307377ae75963820f661df8c251668c4156b0fdd0b76c'
    cpus 4
    memory '4 GB'

    input:
    tuple val(meta), path(alignment), path(index)
    tuple path(ref_fasta), path(ref_fasta_fai)

    output:
    tuple val(meta), path("${meta.sample}.mosdepth.{summary,global.dist}.txt"), emit: reports

    script:
    """
    mosdepth --fast-mode --no-per-base -t ${task.cpus} --fasta ${ref_fasta} ${meta.sample} ${alignment}
    """

    stub:
    """
    touch ${meta.sample}.mosdepth.summary.txt ${meta.sample}.mosdepth.global.dist.txt
    """
}
