process PARABRICKS_FQ2BAM {
    tag "${meta.sample}"
    stageInMode 'copy'
    container 'nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1'
    cpus 16
    memory '64 GB'

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
        cpus: task.cpus, memory_gb: task.memory.toGiga(), num_gpus: task.accelerator ? task.accelerator.request : 1
    ])
    """
    set -euo pipefail
    INDEX=\$(find -L ${bwa_index}/ -name '*.amb' | sed 's/\\.amb\$//')
    cp -L ${fasta} "\$INDEX"
    cp -L ${fai} "\$INDEX.fai"
    pbrun fq2bam --ref "\$INDEX" ${args}
    samtools idxstats ${meta.sample}.md.cram > ${meta.sample}.fq2bam.idxstats
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
    def per = meta.single_end ? 1 : 2
    def a = []
    meta.read_groups.eachWithIndex { rg, i ->
        def files = r.subList(i * per, i * per + per).join(' ')
        a << (meta.single_end ? "--in-se-fq ${files} \"${rg}\"" : "--in-fq ${files} \"${rg}\"")
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
    a << "--bwa-options=\"-K 100000000 -Y${meta.status == 1 ? ' -B 3' : ''}\""
    a << "--bwa-cpu-thread-pool ${opts.cpus}"
    a << "--memory-limit ${Math.max(1, (opts.memory_gb as long).intdiv(2))}"
    a += ['--gpuwrite', '--gpusort', '--low-memory', '--monitor-usage']
    if (meta.single_end && opts.markdups_se_mode == 'start-end') a << '--markdups-single-ended-start-end'
    a << "--num-gpus ${opts.num_gpus}"
    a << '--tmp-dir .'
    a.join(' ')
}
