process PARABRICKS_APPLYBQSR {
    tag "${meta.sample}"
    stageInMode 'copy'
    container 'nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1'
    cpus 12
    memory '72 GB'

    input:
    tuple val(meta), path(alignment, stageAs: 'input/*'), path(alignment_index, stageAs: 'input/*'), path(recal_table)
    tuple path(ref_fasta), path(ref_fasta_fai)

    output:
    tuple val(meta), path("${meta.sample}.recal.bam"), emit: bam

    script:
    def args     = task.ext.args ?: ''
    def num_gpus = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ''
    def expected_index = "${alignment}.${alignment.name.endsWith('.cram') ? 'crai' : 'bai'}"
    """
    set -euo pipefail
    if [ "${alignment_index}" != "${expected_index}" ]; then
        ln -sf "\$(basename ${alignment_index})" "${expected_index}"
    fi

    pbrun \\
        applybqsr \\
        --ref ${ref_fasta} \\
        --in-bam ${alignment} \\
        --in-recal-file ${recal_table} \\
        --out-bam ${meta.sample}.recal.bam \\
        --num-threads ${task.cpus} \\
        ${num_gpus} \\
        ${args}
    """

    stub:
    """
    touch ${meta.sample}.recal.bam
    """
}
