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
    def idx = output_fmt == 'cram' ? 'crai' : 'bai'
    def convert = "samtools view -@ ${task.cpus} ${output_fmt == 'cram' ? "-C --output-fmt-option version=${params.cram_version}" : '-b'} -T ${ref_fasta} --write-index -o ${out}##idx##${out}.${idx} ${alignment}"
    def link = "ln -s ${alignment} ${out} && samtools index -@ ${task.cpus} ${out}"
    if (in_fmt != output_fmt) {
        """
        set -euo pipefail
        ${convert}
        """
    } else if (output_fmt == 'bam') {
        """
        set -euo pipefail
        ${link}
        """
    } else {
        """
        set -euo pipefail
        version=\$(head -c 6 ${alignment} | tail -c 2 | od -An -tu1 | tr -s ' ' | sed 's/^ //; s/ \$//; s/ /./')
        if [ "\$version" = "${params.cram_version}" ]; then
            ${link}
        else
            ${convert}
        fi
        """
    }

    stub:
    def idx = output_fmt == 'cram' ? 'crai' : 'bai'
    """
    touch ${meta.sample}.${meta.suffix}.${output_fmt} ${meta.sample}.${meta.suffix}.${output_fmt}.${idx}
    """
}
