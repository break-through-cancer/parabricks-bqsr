workflow READ_CHECKS {
    take:
    final_stats     // [ meta, stats ]
    before_idxstats // [ meta, idxstats ]
    fastq_counts    // [ sample, long ] -- empty for the alignment entry

    main:
    finals = final_stats.map { meta, stats -> [meta.sample, samtoolsStatsCounts(stats)] }
    befores = before_idxstats.map { meta, idx -> [meta.sample, idxstatsTotal(idx)] }
    fastqs = fastq_counts.toList().map { rows -> rows.collectEntries { [(it[0]): it[1]] } }

    results = finals
        .join(befores, failOnDuplicate: true, failOnMismatch: true)
        .combine(fastqs.map { [counts: it] })
        .map { sample, fin, before, wrapper ->
            def c = [sample: sample, fastq_reads: wrapper.counts[sample], before_records: before,
                     final_primary: fin.primary, final_records: fin.records]
            def failures = readCheckFailures(c)
            if (failures) error failures.join('\n')
            log.info "Read checks passed for '${sample}': ${c.final_primary} primary reads, ${c.final_records} records"
            "${sample}\t${c.fastq_reads ?: 'NA'}\t${c.before_records}\t${c.final_primary}\t${c.final_records}\tPASS"
        }

    mqc = results.collectFile(
        name: 'read_checks_mqc.tsv', newLine: true, sort: true,
        seed: "# id: 'read_checks'\n# section_name: 'Read-count checks'\n# plot_type: 'table'\nSample\tFASTQ_reads\tRecords_before\tFinal_primary\tFinal_records\tResult"
    )

    emit:
    mqc
}

def fastpReadCount(Object json, boolean trimmed) {
    def s = new groovy.json.JsonSlurper().parseText(toPath(json).text).summary
    (trimmed ? s.after_filtering.total_reads : s.before_filtering.total_reads) as long
}

def samtoolsStatsCounts(Object stats) {
    def v = [:]
    toPath(stats).eachLine { line ->
        def t = line.tokenize('\t')
        if (t.size() >= 3 && t[0] == 'SN') v[t[1]] = t[2] as long
    }
    def raw = v['raw total sequences:'] ?: 0L
    [primary: raw, records: raw + (v['non-primary alignments:'] ?: 0L) + (v['supplementary alignments:'] ?: 0L)]
}

def idxstatsTotal(Object idxstats) {
    def total = 0L
    toPath(idxstats).eachLine { line ->
        def t = line.tokenize('\t')
        if (t.size() >= 4) total += (t[2] as long) + (t[3] as long)
    }
    total
}

def readCheckFailures(Map c) {
    def f = []
    if (c.fastq_reads != null && c.fastq_reads != c.final_primary) {
        f << "Sample '${c.sample}': FASTQ reads (${c.fastq_reads}) != final primary reads (${c.final_primary}). Likely causes: reads longer than fq2bam --max-read-length (480 bp) or zero-length reads."
    }
    if (c.before_records != c.final_records) {
        f << "Sample '${c.sample}': records before BQSR/finalize (${c.before_records}) != final records (${c.final_records})."
    }
    f
}

def toPath(Object p) {
    p instanceof CharSequence ? file(p.toString()) : p
}
