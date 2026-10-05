include { pbrunFunction } from './pbrun_wrapper'

process PARABRICKS_FQ2BAM {
    tag "${meta.sample}"
    stageInMode 'copy'
    container 'nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1'
    cpus { fq2bamResources(params.fq2bam_gpus).cpus }
    memory { "${fq2bamResources(params.fq2bam_gpus).memory_gb * task.attempt} GB" }

    input:
    tuple val(meta), path(reads, stageAs: 'reads/?/*')
    tuple path(fasta), path(fai)
    path bwa_index
    tuple path(vcfs, stageAs: 'known_sites/?/*'), path(tbis, stageAs: 'known_sites/?/*')
    path intervals

    output:
    tuple val(meta), path("${meta.sample}.md.cram"), path("${meta.sample}.md.cram.crai"), emit: cram
    tuple val(meta), path("${meta.sample}.table"), emit: table, optional: true
    tuple val(meta), path("${meta.sample}.fq2bam.idxstats"), emit: idxstats
    tuple val(meta), path("${meta.sample}.duplicate-metrics.txt"), emit: duplicate_metrics
    tuple val(meta), path("${meta.sample}_qc_metrics"), emit: qc_metrics

    script:
    def args = fq2bamArgs(meta, reads, vcfs, intervals, [
        optical_distance: params.optical_duplicate_pixel_distance, markdups_se_mode: params.markdups_se_mode,
        cpus: task.cpus, memory_gb: task.memory.toGiga(), num_gpus: task.accelerator ? task.accelerator.request : 1,
        low_memory: params.fq2bam_low_memory, gpuwrite: params.fq2bam_gpuwrite
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
    samtools idxstats -@ ${task.cpus} ${meta.sample}.md.cram > ${meta.sample}.fq2bam.idxstats
    """

    stub:
    def table = (vcfs instanceof List ? vcfs : [vcfs]).findAll() ? "touch ${meta.sample}.table" : ''
    """
    touch ${meta.sample}.md.cram ${meta.sample}.md.cram.crai ${meta.sample}.duplicate-metrics.txt
    printf 'chrT\\t500\\t0\\t0\\n*\\t0\\t0\\t0\\n' > ${meta.sample}.fq2bam.idxstats
    mkdir ${meta.sample}_qc_metrics
    ${table}
    """
}

def asList(Object x) {
    x == null ? [] : (x instanceof List ? x : [x]).findAll { it }
}

def fq2bamArgs(Map meta, Object reads, Object vcfs, Object intervals, Map opts) {
    def r = asList(reads).collect { it.toString() }
    def lane_se = meta.lane_single_end ?: meta.read_groups.collect { meta.single_end }
    if (lane_se.any() && !lane_se.every()) {
        error "fq2bam takes either paired or single-end FASTQs per call; ${meta.sample} has both and must be aligned in parts"
    }
    def a = []
    def offset = 0
    meta.read_groups.eachWithIndex { rg, i ->
        def per = lane_se[i] ? 1 : 2
        def files = r.subList(offset, offset + per).join(' ')
        offset += per
        a << (lane_se[i] ? "--in-se-fq ${files} \"${rg}\"" : "--in-fq ${files} \"${rg}\"")
    }
    if (opts.align_only) {
        a << '--no-markdups'
        a << "--out-bam ${opts.out_bam}"
        a.addAll(fq2bamPerformanceArgs(opts, meta.status == 1))
        return a.join(' ')
    }
    def sites = asList(vcfs).collect { it.toString() }.findAll { it.endsWith('.vcf.gz') }
    sites.each { a << "--knownSites ${it}" }
    if (sites) a << "--out-recal-file ${meta.sample}.table"
    def iv = asList(intervals)
    if (iv) a << "--interval-file ${iv[0]}"
    a << "--out-bam ${meta.sample}.md.cram"
    a << "--out-duplicate-metrics ${meta.sample}.duplicate-metrics.txt"
    a << "--out-qc-metrics-dir ${meta.sample}_qc_metrics"
    a << "--optical-duplicate-pixel-distance ${opts.optical_distance}"
    if (lane_se.any() && opts.markdups_se_mode == 'start-end') a << '--markdups-single-ended-start-end'
    a.addAll(fq2bamPerformanceArgs(opts, meta.status == 1))
    a.join(' ')
}

def fq2bamPerformanceArgs(Map opts, Boolean tumor) {
    def a = []
    a << "--bwa-options=\"-K 100000000 -Y${tumor ? ' -B 3' : ''}\""
    a << "--bwa-cpu-thread-pool ${opts.cpus}"
    a << "--memory-limit ${Math.max(1, (opts.memory_gb as long).intdiv(2))}"
    if (opts.gpuwrite) a << '--gpuwrite'
    a << '--gpusort'
    if (opts.low_memory) a << '--low-memory'
    a << '--monitor-usage'
    a << "--num-gpus ${opts.num_gpus}"
    a << '--tmp-dir .'
    a
}

def fq2bamResources(Object gpus) {
    def n = gpus.toString().toInteger()
    [cpus: Math.max(16, 12 * n), memory_gb: Math.max(64, 44 * n)]
}
