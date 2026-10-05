include { pbrunFunction } from './pbrun_wrapper'

process PARABRICKS_MARKDUP {
    tag "${meta.sample}"
    stageInMode 'copy'
    container 'nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1'
    cpus 16
    memory { "${64 * task.attempt} GB" }

    input:
    tuple val(meta), path(parts, stageAs: 'parts/*')
    tuple path(fasta), path(fai)
    tuple path(vcfs, stageAs: 'known_sites/?/*'), path(tbis, stageAs: 'known_sites/?/*')
    path intervals

    output:
    tuple val(meta), path("${meta.sample}.md.${params.fq2bam_intermediate_fmt}"), path("${meta.sample}.md.${params.fq2bam_intermediate_fmt}.${params.fq2bam_intermediate_fmt == 'bam' ? 'bai' : 'crai'}"), emit: cram
    tuple val(meta), path("${meta.sample}.table"), emit: table, optional: true
    tuple val(meta), path("${meta.sample}.fq2bam.idxstats"), emit: idxstats
    tuple val(meta), path("${meta.sample}.duplicate-metrics.txt"), emit: duplicate_metrics
    tuple val(meta), path("${meta.sample}_qc_metrics"), emit: qc_metrics

    script:
    def s = meta.sample
    def mem = Math.max(1, (task.memory.toGiga() as long).intdiv(2))
    def gpuwrite = params.fq2bam_gpuwrite ? '--gpuwrite' : ''
    def markdup = markdupArgs(meta, [
        optical_distance: params.optical_duplicate_pixel_distance, markdups_se_mode: params.markdups_se_mode,
        cpus: task.cpus, memory_gb: task.memory.toGiga(), gpuwrite: params.fq2bam_gpuwrite, intermediate_fmt: params.fq2bam_intermediate_fmt
    ])
    def out = "${s}.md.${params.fq2bam_intermediate_fmt}"
    def index = "${out}.${params.fq2bam_intermediate_fmt == 'bam' ? 'bai' : 'crai'}"
    def bqsr = bqsrArgs(meta, vcfs, intervals, out)
    """
    set -euo pipefail
    ${pbrunFunction()}
    if command -v nvidia-smi >/dev/null 2>&1; then
        echo "GPU: \$(nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader | paste -sd ';' -)" >&2
    else
        echo "GPU: nvidia-smi not available" >&2
    fi
    samtools merge -@ ${task.cpus} -o ${s}.merged.bam parts/*
    pbrun bamsort --ref ${fasta} --in-bam ${s}.merged.bam --out-bam ${s}.qsorted.bam \\
        --sort-order queryname --sort-compatibility picard --gpusort ${gpuwrite} --mem-limit ${mem} --tmp-dir .
    rm ${s}.merged.bam
    pbrun markdup --ref ${fasta} ${markdup}
    rm ${s}.qsorted.bam
    HD=\$(samtools view -H ${out} | grep '^@HD' || true)
    if [[ "\$HD" != *SO:coordinate* ]]; then
        echo "markdup output is not coordinate-sorted; sorting it" >&2
        mv ${out} unsorted.${out}
        rm -f ${index}
        pbrun bamsort --ref ${fasta} --in-bam unsorted.${out} --out-bam ${out} \\
            --sort-order coordinate --gpusort ${gpuwrite} --mem-limit ${mem} --tmp-dir .
        rm unsorted.${out}
    fi
    [ -f ${index} ] || samtools index -@ ${task.cpus} ${out}
    ${bqsr ? "pbrun bqsr --ref ${fasta} ${bqsr} --tmp-dir ." : ''}
    pbrun collectmultiplemetrics --ref ${fasta} --bam ${out} --out-qc-metrics-dir ${s}_qc_metrics --gen-all-metrics --tmp-dir .
    samtools idxstats -@ ${task.cpus} ${out} > ${s}.fq2bam.idxstats
    """

    stub:
    def out = "${meta.sample}.md.${params.fq2bam_intermediate_fmt}"
    def table = bqsrArgs(meta, vcfs, intervals, out) ? "touch ${meta.sample}.table" : ''
    """
    touch ${out} ${out}.${params.fq2bam_intermediate_fmt == 'bam' ? 'bai' : 'crai'} ${meta.sample}.duplicate-metrics.txt
    printf 'chrT\\t500\\t0\\t0\\n*\\t0\\t0\\t0\\n' > ${meta.sample}.fq2bam.idxstats
    mkdir ${meta.sample}_qc_metrics
    ${table}
    """
}

def markdupArgs(Map meta, Map opts) {
    def a = []
    a << "--in-bam ${meta.sample}.qsorted.bam"
    a << "--out-bam ${meta.sample}.md.${opts.intermediate_fmt ?: 'cram'}"
    a << "--out-duplicate-metrics ${meta.sample}.duplicate-metrics.txt"
    a << "--optical-duplicate-pixel-distance ${opts.optical_distance}"
    if (asList(meta.lane_single_end).any() && opts.markdups_se_mode == 'start-end') a << '--markdups-single-ended-start-end'
    a << "--num-worker-threads ${opts.cpus}"
    a << "--num-zip-threads ${opts.cpus}"
    a << "--mem-limit ${Math.max(1, (opts.memory_gb as long).intdiv(2))}"
    a << '--gpusort'
    if (opts.gpuwrite) a << '--gpuwrite'
    a << '--tmp-dir .'
    a.join(' ')
}

def bqsrArgs(Map meta, Object vcfs, Object intervals, Object alignment) {
    def sites = asList(vcfs).collect { it.toString() }.findAll { it.endsWith('.vcf.gz') }
    if (!sites) return ''
    def a = ["--in-bam ${alignment}"]
    sites.each { a << "--knownSites ${it}" }
    def iv = asList(intervals)
    if (iv) a << "--interval-file ${iv[0]}"
    a << "--out-recal-file ${meta.sample}.table"
    a.join(' ')
}

def asList(Object x) {
    x == null ? [] : (x instanceof List ? x : [x]).findAll { it }
}
