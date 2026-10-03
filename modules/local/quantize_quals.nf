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
    tuple val(meta), path("${meta.sample}.${meta.suffix}.${output_fmt}"), path("${meta.sample}.${meta.suffix}.${output_fmt}.{bai,crai}"), emit: alignment
    path "${meta.sample}.quantize.log", emit: log
    path "${meta.sample}.quantize_mqc.tsv", emit: mqc

    script:
    def args = task.ext.args ?: ''
    def out = "${meta.sample}.${meta.suffix}.${output_fmt}"
    """
    set -euo pipefail
    quantize_quals \\
        --in ${bam} \\
        --out ${out} \\
        --ref ${ref_fasta} \\
        --threads ${task.cpus} \\
        ${args} 2> ${meta.sample}.quantize.log || { cat ${meta.sample}.quantize.log >&2; exit 1; }
    cat ${meta.sample}.quantize.log >&2
    {
        echo "# id: 'quantize_quals'"
        echo "# section_name: 'Quality-score quantization'"
        echo "# plot_type: 'table'"
        printf 'Sample\\tRecords\\tQualities\\tChanged\\tChanged_pct\\n'
        sed -n 's/^quantize_quals: done: \\([0-9]*\\) records ([0-9]* without qualities), \\([0-9]*\\) qualities, \\([0-9]*\\) changed (\\([0-9.]*\\)%).*/${meta.sample}\\t\\1\\t\\2\\t\\3\\t\\4/p' ${meta.sample}.quantize.log
    } > ${meta.sample}.quantize_mqc.tsv
    """

    stub:
    def idx = output_fmt == 'cram' ? 'crai' : 'bai'
    """
    touch ${meta.sample}.${meta.suffix}.${output_fmt} ${meta.sample}.${meta.suffix}.${output_fmt}.${idx} ${meta.sample}.quantize.log
    printf 'Sample\\tRecords\\n${meta.sample}\\t0\\n' > ${meta.sample}.quantize_mqc.tsv
    """
}
