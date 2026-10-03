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
    tuple val(meta), path("${meta.sample}.input.idxstats"), emit: idxstats

    script:
    def expected_index = "${alignment}.${alignment.name.endsWith('.cram') ? 'crai' : 'bai'}"
    def args = applybqsrArgs(alignment, recal_table, "${meta.sample}.recal.bam",
        [cpus: task.cpus, num_gpus: task.accelerator ? task.accelerator.request : 1])
    """
    set -euo pipefail
    if [ "${alignment_index}" != "${expected_index}" ]; then
        ln -sf "\$(basename ${alignment_index})" "${expected_index}"
    fi

    samtools view -H ${alignment} | awk -F'\\t' '\$1=="@SQ"{sn="";ln="";for(i=2;i<=NF;i++){if(\$i~/^SN:/)sn=substr(\$i,4);if(\$i~/^LN:/)ln=substr(\$i,4)};print sn"\\t"ln}' > alignment.contigs
    cut -f1,2 ${ref_fasta_fai} > reference.contigs
    if ! cmp -s alignment.contigs reference.contigs; then
        echo "Alignment ${alignment.name} was not aligned to ${ref_fasta.name}: header contigs (names, lengths, order) differ from the reference .fai" >&2
        diff alignment.contigs reference.contigs | head -20 >&2 || true
        exit 1
    fi
    samtools idxstats ${alignment} > ${meta.sample}.input.idxstats

    pbrun applybqsr --ref ${ref_fasta} ${args}
    """

    stub:
    """
    touch ${meta.sample}.recal.bam
    printf 'chrT\\t500\\t0\\t0\\n*\\t0\\t0\\t0\\n' > ${meta.sample}.input.idxstats
    """
}

def applybqsrArgs(Object alignment, Object table, String out, Map opts) {
    "--in-bam ${alignment} --in-recal-file ${table} --out-bam ${out} --num-threads ${opts.cpus} --num-gpus ${opts.num_gpus}".toString()
}
