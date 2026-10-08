process MOSDEPTH {
    tag "${meta.sample}"
    container 'quay.io/biocontainers/mosdepth@sha256:94b54a6a5610da01030307377ae75963820f661df8c251668c4156b0fdd0b76c'
    cpus 4
    memory { "${8 * task.attempt} GB" }

    input:
    tuple val(meta), path(alignment), path(index)
    tuple path(ref_fasta), path(ref_fasta_fai)

    output:
    tuple val(meta), path("${meta.sample}.mosdepth.{summary,global.dist}.txt"), emit: reports

    script:
    """
    set -euo pipefail
    mosdepth --fast-mode --no-per-base -t ${task.cpus} --fasta ${ref_fasta} ${meta.sample} ${alignment} 2> mosdepth.err \
        || { rc=\$?; cat mosdepth.err >&2; exit \$rc; }
    cat mosdepth.err >&2
    covered=\$(awk '\$1 == "total" { print \$2 }' ${meta.sample}.mosdepth.summary.txt)
    if grep -qE '\\[E::|hts-nim\\] error' mosdepth.err || [ "\${covered:-0}" -eq 0 ]; then
        echo "mosdepth could not read ${alignment}: decode errors or no coverage (an unsupported CRAM version?)" >&2
        exit 1
    fi
    """

    stub:
    """
    touch ${meta.sample}.mosdepth.summary.txt ${meta.sample}.mosdepth.global.dist.txt
    """
}
