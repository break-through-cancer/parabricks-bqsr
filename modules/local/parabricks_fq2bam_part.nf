include { pbrunFunction } from './pbrun_wrapper'
include { fq2bamArgs; fq2bamResources } from './parabricks_fq2bam'

process PARABRICKS_FQ2BAM_PART {
    tag "${meta.sample}:${meta.part}"
    stageInMode 'copy'
    container 'nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1'
    cpus { fq2bamResources(params.fq2bam_gpus).cpus }
    memory { "${fq2bamResources(params.fq2bam_gpus).memory_gb * task.attempt} GB" }

    input:
    tuple val(meta), path(reads, stageAs: 'reads/?/*')
    tuple path(fasta), path(fai)
    path bwa_index

    output:
    tuple val(meta), path("${meta.sample}.${meta.part}.bam"), emit: bam

    script:
    def args = fq2bamArgs(meta, reads, [], [], [
        cpus: task.cpus, memory_gb: task.memory.toGiga(), num_gpus: task.accelerator ? task.accelerator.request : 1,
        low_memory: params.fq2bam_low_memory, gpuwrite: params.fq2bam_gpuwrite,
        align_only: true, out_bam: "${meta.sample}.${meta.part}.bam"
    ])
    """
    set -euo pipefail
    ${pbrunFunction()}
    if command -v nvidia-smi >/dev/null 2>&1; then
        echo "GPU: \$(nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader | paste -sd ';' -)" >&2
    else
        echo "GPU: nvidia-smi not available" >&2
    fi
    INDEX=\$(find -L ${bwa_index}/ -name '*.amb' | sed 's/\\.amb\$//')
    cp -L ${fasta} "\$INDEX"
    cp -L ${fai} "\$INDEX.fai"
    pbrun fq2bam --ref "\$INDEX" ${args}
    """

    stub:
    """
    touch ${meta.sample}.${meta.part}.bam
    """
}
