nextflow.enable.dsl = 2

include { SAMPLESHEET_TO_SAMPLES                    } from './subworkflows/local/samplesheet_to_samples/main'
include { SAMPLESHEET_TO_FASTQS; samplesheetEntry   } from './subworkflows/local/samplesheet_to_fastqs/main'
include { resolveReferences; checkReferences; intervalsMessage } from './subworkflows/local/references/main'
include { PREPARE_KNOWN_SITES                       } from './subworkflows/local/prepare_known_sites/main'
include { READ_CHECKS; fastpReadCount               } from './subworkflows/local/read_checks/main'
include { FASTP                                     } from './modules/local/fastp'
include { PARABRICKS_FQ2BAM                         } from './modules/local/parabricks_fq2bam'
include { PARABRICKS_APPLYBQSR                      } from './modules/local/parabricks_applybqsr'
include { QUANTIZE_QUALS                            } from './modules/local/quantize_quals'
include { SAMTOOLS_FINALIZE                         } from './modules/local/samtools_finalize'
include { SAMTOOLS_STATS                            } from './modules/local/samtools_stats'
include { MOSDEPTH                                  } from './modules/local/mosdepth'
include { MULTIQC                                   } from './modules/local/multiqc'

workflow {
    validateParams()
    def entry = samplesheetEntry(params.input)
    def refs = resolveReferences(params)
    def checks = checkReferences(refs, entry == 'fastq')
    checks.warnings.each { log.warn it }
    if (checks.errors) error checks.errors.join('\n')

    ref_ch = Channel.value([file(refs.ref_fasta, checkIfExists: true), file(refs.ref_fasta_fai, checkIfExists: true)])
    def quantize = params.quantize_quals_enabled.toString() == 'true'

    if (entry == 'fastq') {
        def bqsr = refs.known_sites as boolean
        log.info intervalsMessage(refs)
        if (!bqsr) log.warn "No known sites supplied: BQSR is skipped; outputs are published as .md files under preprocessing/markduplicates/"

        SAMPLESHEET_TO_FASTQS(params.input)
        def trim = params.trim_fastq.toString() == 'true'
        FASTP(SAMPLESHEET_TO_FASTQS.out.lanes, trim)

        fastq_counts = FASTP.out.json
            .map { meta, json -> [groupKey(meta.sample, meta.n_lanes), fastpReadCount(json, trim)] }
            .groupTuple()
            .map { sample, counts -> [sample.toString(), counts.sum()] }

        sample_reads = FASTP.out.reads
            .map { meta, reads -> [groupKey(meta.sample, meta.n_lanes), meta, reads instanceof List ? reads : [reads]] }
            .groupTuple()
            .map { sample, metas, reads ->
                def ordered = [metas, reads].transpose().sort { it[0].lane.toString() }
                def m = ordered[0][0]
                [[sample: m.sample, patient: m.patient, status: m.status, single_end: m.single_end,
                  read_groups: ordered.collect { it[0].read_group }], ordered.collectMany { it[1] }]
            }

        PREPARE_KNOWN_SITES(refs.known_sites)
        PARABRICKS_FQ2BAM(
            sample_reads, ref_ch, Channel.value(file(refs.bwa_index, checkIfExists: true)),
            PREPARE_KNOWN_SITES.out.known_sites, Channel.value(refs.intervals ? file(refs.intervals, checkIfExists: true) : [])
        )

        if (bqsr) {
            PARABRICKS_APPLYBQSR(
                PARABRICKS_FQ2BAM.out.cram.join(PARABRICKS_FQ2BAM.out.table, failOnMismatch: true)
                    .map { meta, cram, crai, table -> [[sample: meta.sample], cram, crai, table] },
                ref_ch
            )
            to_finish = PARABRICKS_APPLYBQSR.out.bam.map { meta, bam -> [meta + [suffix: 'recal'], bam] }
        } else {
            to_finish = PARABRICKS_FQ2BAM.out.cram.map { meta, cram, _crai -> [[sample: meta.sample, suffix: 'md'], cram] }
        }
        before_idxstats = PARABRICKS_FQ2BAM.out.idxstats.map { meta, idx -> [[sample: meta.sample], idx] }
        extra_qc = FASTP.out.json.map { it[1] }
            .mix(PARABRICKS_FQ2BAM.out.duplicate_metrics.map { it[1] })
    } else {
        fastq_counts = Channel.empty()
        SAMPLESHEET_TO_SAMPLES(params.input)
        PARABRICKS_APPLYBQSR(SAMPLESHEET_TO_SAMPLES.out.samples, ref_ch)
        to_finish = PARABRICKS_APPLYBQSR.out.bam.map { meta, bam -> [meta + [suffix: 'recal'], bam] }
        before_idxstats = PARABRICKS_APPLYBQSR.out.idxstats
        extra_qc = Channel.empty()
    }

    if (quantize) {
        QUANTIZE_QUALS(to_finish, ref_ch, params.output_fmt)
        final_ch = QUANTIZE_QUALS.out.alignment
        quant_mqc = QUANTIZE_QUALS.out.mqc
    } else {
        SAMTOOLS_FINALIZE(to_finish, ref_ch, params.output_fmt)
        final_ch = SAMTOOLS_FINALIZE.out.alignment
        quant_mqc = Channel.empty()
    }

    SAMTOOLS_STATS(final_ch, ref_ch)
    MOSDEPTH(final_ch, ref_ch)
    READ_CHECKS(
        SAMTOOLS_STATS.out.stats.map { meta, s -> [[sample: meta.sample], s] },
        before_idxstats.map { meta, idx -> [[sample: meta.sample], idx] },
        fastq_counts
    )
    MULTIQC(
        extra_qc
            .mix(SAMTOOLS_STATS.out.stats.map { it[1] })
            .mix(MOSDEPTH.out.reports.flatMap { it[1] })
            .mix(quant_mqc)
            .mix(READ_CHECKS.out.mqc)
            .collect()
    )
}

def validateParams() {
    if (!params.input) error "Missing required parameter: input (path to samplesheet CSV). See README.md and assets/samplesheet.csv."
    if (!(params.output_fmt in ['bam', 'cram'])) error "Invalid output_fmt '${params.output_fmt}': must be 'bam' or 'cram'"
    if (!(params.quantize_quals_enabled.toString() in ['true', 'false'])) error "Invalid quantize_quals_enabled '${params.quantize_quals_enabled}': must be true or false"
    if (!(params.markdups_se_mode in ['5prime', 'start-end'])) error "Invalid markdups_se_mode '${params.markdups_se_mode}': must be '5prime' or 'start-end'"
    if (!params.fq2bam_gpus.toString().isInteger() || params.fq2bam_gpus.toString().toInteger() < 1) error "Invalid fq2bam_gpus '${params.fq2bam_gpus}': must be a positive integer"
    if (params.quantize_quals_enabled.toString() == 'true') {
        def bins = params.static_quantized_quals.toString().tokenize(',')*.trim()
        if (!bins || bins.any { !it.isInteger() }) error "Invalid static_quantized_quals '${params.static_quantized_quals}': must be a comma-separated list of integers, e.g. '10,20,30'"
        if (!params.preserve_qscores_less_than.toString().isInteger()) error "Invalid preserve_qscores_less_than '${params.preserve_qscores_less_than}': must be an integer"
    }
    ['clip_r1', 'clip_r2', 'three_prime_clip_r1', 'three_prime_clip_r2', 'length_required'].each { p ->
        if (!params[p].toString().isInteger() || params[p].toString().toInteger() < 0) error "Invalid ${p} '${params[p]}': must be a non-negative integer"
    }
}
