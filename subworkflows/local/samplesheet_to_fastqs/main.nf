include { resolveRelativeToProjectDir } from '../samplesheet_to_samples/main'

workflow SAMPLESHEET_TO_FASTQS {
    take:
    samplesheet

    main:
    rows = Channel.fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> validateFastqRow(row) }

    checked_ch = rows.toList().map { all -> checkFastqSamples(all); all }

    lanes = checked_ch
        .flatMap { all ->
            def nLanes = all.countBy { it[0].sample }
            all.collect { meta, reads -> [meta + [n_lanes: nLanes[meta.sample]], reads] }
        }
        .ifEmpty { error "Samplesheet produced zero FASTQ rows -- check that it has at least one data row below the header" }

    emit:
    lanes // [ meta, reads ]
}

def samplesheetEntry(Object path) {
    def header = toPath(path).withReader { it.readLine() }?.tokenize(',')*.trim() ?: []
    def fastq = 'fastq_1' in header
    def alignment = 'alignment' in header
    if (fastq && alignment) error "Samplesheet ${path} has both FASTQ and alignment columns; use one entry point per run"
    if (!fastq && !alignment) error "Samplesheet ${path} has neither a 'fastq_1' nor an 'alignment' column"
    fastq ? 'fastq' : 'alignment'
}

def flowcellFromFastq(Object path) {
    def f = toPath(path)
    def line = f.withInputStream { is ->
        def stream = f.name.endsWith('.gz') ? new java.util.zip.GZIPInputStream(is) : is
        new BufferedReader(new InputStreamReader(stream)).readLine()
    }
    def fields = line?.startsWith('@') ? line.substring(1).tokenize(' ')[0].tokenize(':') : []
    fields.size() >= 7 ? fields[2] : 'unknown'
}

def validateFastqRow(Map row) {
    if (!row.sample) error "FASTQ samplesheet row missing 'sample': ${row}"
    if (!row.lane) error "FASTQ samplesheet row '${row.sample}' missing 'lane'"
    if (!row.fastq_1) error "FASTQ samplesheet row '${row.sample}' lane '${row.lane}' missing 'fastq_1'"
    def status = row.status ?: '0'
    if (!(status in ['0', '1'])) error "FASTQ samplesheet row '${row.sample}': status must be 0 or 1, got '${status}'"
    def patient = row.patient ?: row.sample
    def reads = [file(resolveRelativeToProjectDir(row.fastq_1), checkIfExists: true)]
    if (row.fastq_2) reads << file(resolveRelativeToProjectDir(row.fastq_2), checkIfExists: true)
    def flowcell = flowcellFromFastq(reads[0])
    def unit = "${flowcell}.${row.sample}.${row.lane}"
    def rg = "@RG\\tID:${unit}\\tPU:${unit}\\tSM:${patient}_${row.sample}\\tLB:${row.sample}\\tPL:${params.seq_platform}".toString()
    def meta = [id: "${row.sample}_${row.lane}".toString(), sample: row.sample, patient: patient, status: status.toInteger(),
                lane: row.lane, single_end: reads.size() == 1, read_group: rg]
    [meta, reads]
}

def checkFastqSamples(List rows) {
    rows.groupBy { it[0].sample }.each { sample, lanes ->
        if (lanes.collect { it[0].patient }.unique().size() > 1) error "Sample '${sample}' has more than one patient"
        if (lanes.collect { it[0].status }.unique().size() > 1) error "Sample '${sample}' has more than one status"
        if (lanes.collect { it[0].single_end }.unique().size() > 1) {
            error "Sample '${sample}' mixes paired-end and single-end lanes; this is not supported until a single fq2bam call is validated with both"
        }
        lanes.countBy { it[0].lane }.findAll { _lane, n -> n > 1 }.each { lane, _n ->
            error "Sample '${sample}' has lane '${lane}' more than once"
        }
    }
}

def toPath(Object p) {
    p instanceof CharSequence ? file(p.toString()) : p
}
