include { TABIX_VCF } from '../../../modules/local/tabix_vcf'

workflow PREPARE_KNOWN_SITES {
    take:
    known_sites // List<String> of .vcf.gz paths

    main:
    indexed = known_sites.withIndex().collect { v, i ->
        def vcf = file(v, checkIfExists: true)
        def tbi = file("${v}.tbi")
        [i, vcf, tbi.exists() ? tbi : null]
    }
    have_ch = Channel.fromList(indexed.findAll { it[2] })
    TABIX_VCF(Channel.fromList(indexed.findAll { !it[2] }.collect { [it[0], it[1]] }))
    made_ch = TABIX_VCF.out

    known_sites_out = have_ch.mix(made_ch)
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { rows -> [rows.collect { it[1] }, rows.collect { it[2] }] }

    emit:
    known_sites = known_sites_out
}
