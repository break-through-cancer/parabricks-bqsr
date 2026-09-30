workflow SAMPLESHEET_TO_SAMPLES {
    take:
    samplesheet // path to samplesheet CSV (sample,alignment,alignment_index,recal_table)

    main:
    rows = Channel.fromPath(resolveRelativeToProjectDir(samplesheet), checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> validateSamplesheetRow(row) }

    unique_ids_ch = rows
        .map { meta, _alignment, _index, _table -> meta.sample }
        .toList()
        .map { ids ->
            def dupes = ids.countBy { it }.findAll { _id, n -> n > 1 }.keySet().sort()
            if (dupes) {
                error "Duplicate sample ID(s) in samplesheet: ${dupes.join(', ')}"
            }
            [ ok: true ]
        }

    samples = rows
        .combine(unique_ids_ch)
        .map { meta, alignment, index, table, _ok -> [ meta, alignment, index, table ] }
        .ifEmpty { error "Samplesheet produced zero samples -- check that it has at least one data row below the header" }

    emit:
    samples // [ meta(sample), alignment, alignment_index, recal_table ]
}

def validateSamplesheetRow(row) {
    if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
    ['alignment', 'alignment_index', 'recal_table'].each { col ->
        if (!row[col]) error "Samplesheet row '${row.sample}' missing '${col}'"
    }

    def meta      = [ sample: row.sample ]
    def alignment = file(resolveRelativeToProjectDir(row.alignment), checkIfExists: true)
    def index     = file(resolveRelativeToProjectDir(row.alignment_index), checkIfExists: true)
    def table     = file(resolveRelativeToProjectDir(row.recal_table), checkIfExists: true)
    [ meta, alignment, index, table ]
}

def resolveRelativeToProjectDir(p) {
    def s = p.toString()
    def isAbsolute = s.startsWith('/') || s ==~ /^[a-zA-Z][a-zA-Z0-9+.\-]*:\/\/.*/
    isAbsolute ? file(s) : file("${projectDir}/${s}")
}
