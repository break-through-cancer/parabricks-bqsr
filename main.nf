nextflow.enable.dsl = 2

include { SAMPLESHEET_TO_SAMPLES                    } from './subworkflows/local/samplesheet_to_samples/main'
include { SAMPLESHEET_TO_FASTQS; samplesheetEntry   } from './subworkflows/local/samplesheet_to_fastqs/main'
include { resolveReferences; checkReferences; intervalsMessage } from './subworkflows/local/references/main'
include { PREPARE_KNOWN_SITES                       } from './subworkflows/local/prepare_known_sites/main'
include { READ_CHECKS; fastpReadCount               } from './subworkflows/local/read_checks/main'
include { FASTP; fastpWritesReads                   } from './modules/local/fastp'
include { PARABRICKS_FQ2BAM                         } from './modules/local/parabricks_fq2bam'
include { PARABRICKS_FQ2BAM_PART                    } from './modules/local/parabricks_fq2bam_part'
include { PARABRICKS_MARKDUP                        } from './modules/local/parabricks_markdup'
include { PARABRICKS_APPLYBQSR                      } from './modules/local/parabricks_applybqsr'
include { QUANTIZE_QUALS; quantizeHistJson          } from './modules/local/quantize_quals'
include { SAMTOOLS_FINALIZE                         } from './modules/local/samtools_finalize'
include { SAMTOOLS_FINALIZE as SAMTOOLS_FINALIZE_RAW } from './modules/local/samtools_finalize'
include { SAMTOOLS_STATS                            } from './modules/local/samtools_stats'
include { SAMTOOLS_FLAGSTAT                         } from './modules/local/samtools_flagstat'
include { MOSDEPTH                                  } from './modules/local/mosdepth'
include { MULTIQC                                   } from './modules/local/multiqc'
include { paramsSummaryLog                          } from 'plugin/nf-schema'

workflow {
    log.info paramsSummaryLog(workflow)
    validateParams()
    def entry = samplesheetEntry(params.input)
    def refs = resolveReferences(params)
    def checks = checkReferences(refs, entry == 'fastq')
    checks.warnings.each { log.warn it }
    if (checks.errors) error checks.errors.join('\n')

    ref_ch = Channel.value([file(refs.ref_fasta, checkIfExists: true), file(refs.ref_fasta_fai, checkIfExists: true)])
    def quantize = params.quantize_quals_enabled.toString() == 'true'
    // Skips SAMTOOLS_FINALIZE when applybqsr's own BAM output already needs no change.
    def applybqsr_direct = params.output_fmt == 'bam' && !quantize
    def use_direct_final = false
    def direct_final = channel.empty()

    if (entry == 'fastq') {
        def bqsr = refs.known_sites as boolean
        def apply_bqsr = bqsr && params.apply_bqsr.toString() != 'false'
        // Skips SAMTOOLS_FINALIZE when fq2bam/markdup's own output already needs no change.
        def md_direct = !apply_bqsr && params.fq2bam_intermediate_fmt == params.output_fmt
        log.info intervalsMessage(refs)
        if (!bqsr) log.warn "No known sites supplied: BQSR is skipped; outputs are published as .md files under preprocessing/markduplicates/"

        SAMPLESHEET_TO_FASTQS(params.input)
        def trim = params.trim_fastq.toString() == 'true'
        def fastp_reads = fastpWritesReads(trim, params.poly_g_trimming)
        FASTP(SAMPLESHEET_TO_FASTQS.out.lanes, trim)

        fastq_counts = FASTP.out.json
            .map { meta, json -> [groupKey(meta.sample, meta.n_lanes), fastpReadCount(json, fastp_reads)] }
            .groupTuple()
            .map { sample, counts -> [sample.toString(), counts.sum()] }

        sample_reads = (fastp_reads ? FASTP.out.reads : SAMPLESHEET_TO_FASTQS.out.lanes)
            .map { meta, reads -> [groupKey(meta.sample, meta.n_lanes), meta, reads instanceof List ? reads : [reads]] }
            .groupTuple()
            .map { sample, metas, reads ->
                def ordered = [metas, reads].transpose().sort { it[0].lane.toString() }
                def m = ordered[0][0]
                [[sample: m.sample, patient: m.patient, status: m.status,
                  single_end: ordered.every { it[0].single_end }, lane_single_end: ordered.collect { it[0].single_end },
                  read_groups: ordered.collect { it[0].read_group }, publish_directly: md_direct], ordered.collectMany { it[1] }]
            }

        PREPARE_KNOWN_SITES(refs.known_sites)
        def bwa_index_ch = channel.value(file(refs.bwa_index, checkIfExists: true))
        def intervals_ch = channel.value(refs.intervals ? file(refs.intervals, checkIfExists: true) : [])
        by_kind = sample_reads.branch { meta, _reads ->
            mixed: meta.lane_single_end.any() && !meta.lane_single_end.every()
            single: true
        }
        PARABRICKS_FQ2BAM(by_kind.single, ref_ch, bwa_index_ch, PREPARE_KNOWN_SITES.out.known_sites, intervals_ch)

        PARABRICKS_FQ2BAM_PART(by_kind.mixed.flatMap { meta, reads -> splitByReadKind(meta, reads) }, ref_ch, bwa_index_ch)
        PARABRICKS_MARKDUP(
            PARABRICKS_FQ2BAM_PART.out.bam
                .map { meta, bam -> [groupKey(meta.sample, 2), meta.parent, bam] }
                .groupTuple()
                .map { _sample, parents, bams -> [parents[0], bams.sort { bam -> bam.name }] },
            ref_ch, PREPARE_KNOWN_SITES.out.known_sites, intervals_ch
        )
        aligned_cram = PARABRICKS_FQ2BAM.out.cram.mix(PARABRICKS_MARKDUP.out.cram)
        aligned_table = PARABRICKS_FQ2BAM.out.table.mix(PARABRICKS_MARKDUP.out.table)

        // The pre-recalibration alignment is otherwise discarded once BQSR runs.
        if (apply_bqsr && params.publish_markduplicates.toString() == 'true') {
            SAMTOOLS_FINALIZE_RAW(
                aligned_cram.map { meta, cram, _crai -> [[sample: meta.sample, suffix: 'md'], cram] },
                ref_ch, params.output_fmt
            )
        }

        if (apply_bqsr) {
            PARABRICKS_APPLYBQSR(
                aligned_cram.join(aligned_table, failOnMismatch: true)
                    .map { meta, cram, crai, table -> [[sample: meta.sample], cram, crai, table] },
                ref_ch
            )
            to_finish = PARABRICKS_APPLYBQSR.out.bam.map { meta, bam -> [meta + [suffix: 'recal'], bam] }
            if (applybqsr_direct) {
                use_direct_final = true
                direct_final = PARABRICKS_APPLYBQSR.out.bam.join(PARABRICKS_APPLYBQSR.out.bai)
                    .map { meta, bam, bai -> [meta + [suffix: 'recal'], bam, bai] }
            }
        } else {
            to_finish = aligned_cram.map { meta, cram, _crai -> [[sample: meta.sample, suffix: 'md'], cram] }
            if (md_direct) {
                use_direct_final = true
                direct_final = aligned_cram.map { meta, cram, crai -> [[sample: meta.sample, suffix: 'md'], cram, crai] }
            }
        }
        // Quantization is deferred along with BQSR when apply_bqsr=false.
        quantize = quantize && (apply_bqsr || !bqsr)
        before_idxstats = PARABRICKS_FQ2BAM.out.idxstats.mix(PARABRICKS_MARKDUP.out.idxstats)
            .map { meta, idx -> [[sample: meta.sample], idx] }
        extra_qc = FASTP.out.json.map { it[1] }
            .mix(PARABRICKS_FQ2BAM.out.duplicate_metrics.mix(PARABRICKS_MARKDUP.out.duplicate_metrics).map { _meta, metrics -> metrics })
    } else {
        fastq_counts = Channel.empty()
        SAMPLESHEET_TO_SAMPLES(params.input)
        PARABRICKS_APPLYBQSR(SAMPLESHEET_TO_SAMPLES.out.samples, ref_ch)
        to_finish = PARABRICKS_APPLYBQSR.out.bam.map { meta, bam -> [meta + [suffix: 'recal'], bam] }
        if (applybqsr_direct) {
            use_direct_final = true
            direct_final = PARABRICKS_APPLYBQSR.out.bam.join(PARABRICKS_APPLYBQSR.out.bai)
                .map { meta, bam, bai -> [meta + [suffix: 'recal'], bam, bai] }
        }
        before_idxstats = PARABRICKS_APPLYBQSR.out.idxstats
        extra_qc = Channel.empty()
    }

    if (quantize) {
        QUANTIZE_QUALS(to_finish, ref_ch, params.output_fmt)
        final_ch = QUANTIZE_QUALS.out.alignment
        quant_mqc = QUANTIZE_QUALS.out.mqc.mix(
            QUANTIZE_QUALS.out.hist.collect().map { files -> quantizeHistJson(files) }
                .collectFile(name: 'quantize_quals_hist_mqc.json', newLine: false)
        )
    } else if (use_direct_final) {
        final_ch = direct_final
        quant_mqc = channel.empty()
    } else {
        SAMTOOLS_FINALIZE(to_finish, ref_ch, params.output_fmt)
        final_ch = SAMTOOLS_FINALIZE.out.alignment
        quant_mqc = Channel.empty()
    }

    // stats reads the pre-quantization alignment (unaffected by BQSR/quantization except
    // quality bytes); flagstat and mosdepth read the true published file.
    SAMTOOLS_STATS(to_finish, ref_ch)
    SAMTOOLS_FLAGSTAT(final_ch.map { meta, aln, _idx -> [meta, aln] }, ref_ch)
    if (params.output_fmt == 'cram' && params.cram_version == '3.1') {
        log.warn "mosdepth cannot read CRAM 3.1: coverage QC is skipped (use --cram_version 3.0 to keep it)"
        mosdepth_reports = Channel.empty()
    } else {
        MOSDEPTH(final_ch, ref_ch)
        mosdepth_reports = MOSDEPTH.out.reports.flatMap { it[1] }
    }
    READ_CHECKS(
        SAMTOOLS_FLAGSTAT.out.flagstat.map { meta, s -> [[sample: meta.sample], s] },
        before_idxstats.map { meta, idx -> [[sample: meta.sample], idx] },
        fastq_counts
    )
    MULTIQC(
        extra_qc
            .mix(SAMTOOLS_STATS.out.stats.map { it[1] })
            .mix(SAMTOOLS_FLAGSTAT.out.flagstat.map { _meta, flagstat -> flagstat })
            .mix(mosdepth_reports)
            .mix(quant_mqc)
            .mix(READ_CHECKS.out.mqc)
            .collect()
    )
}

def validateParams() {
    if (!params.input) error "Missing required parameter: input (path to samplesheet CSV). See README.md and assets/samplesheet.csv."
    if (!(params.output_fmt in ['bam', 'cram'])) error "Invalid output_fmt '${params.output_fmt}': must be 'bam' or 'cram'"
    if (!(params.quantize_quals_enabled.toString() in ['true', 'false'])) error "Invalid quantize_quals_enabled '${params.quantize_quals_enabled}': must be true or false"
    if (!(params.publish_markduplicates.toString() in ['true', 'false'])) error "Invalid publish_markduplicates '${params.publish_markduplicates}': must be true or false"
    if (!(params.apply_bqsr.toString() in ['true', 'false'])) error "Invalid apply_bqsr '${params.apply_bqsr}': must be true or false"
    if (!(params.cram_version.toString() in ['3.0', '3.1'])) error "Invalid cram_version '${params.cram_version}': must be '3.0' or '3.1'"
    if (!(params.poly_g_trimming in ['auto', 'on', 'off'])) error "Invalid poly_g_trimming '${params.poly_g_trimming}': must be 'auto', 'on' or 'off'"
    if (!(params.markdups_se_mode in ['5prime', 'start-end'])) error "Invalid markdups_se_mode '${params.markdups_se_mode}': must be '5prime' or 'start-end'"
    if (params.fq2bam_memory_gb != null && params.fq2bam_memory_gb.toString() != '') {
        def m = params.fq2bam_memory_gb.toString()
        if (!m.isInteger() || !(m.toInteger() in 16..768)) error "Invalid fq2bam_memory_gb '${params.fq2bam_memory_gb}': must be an integer between 16 and 768"
    }
    if (!(params.fq2bam_intermediate_fmt in ['cram', 'bam'])) error "Invalid fq2bam_intermediate_fmt '${params.fq2bam_intermediate_fmt}': must be 'cram' or 'bam'"
    if (!params.fq2bam_gpus.toString().isInteger() || !(params.fq2bam_gpus.toString().toInteger() in 1..4)) error "Invalid fq2bam_gpus '${params.fq2bam_gpus}': must be between 1 and 4"
    if (params.quantize_quals_enabled.toString() == 'true') {
        def bins = params.static_quantized_quals.toString().tokenize(',')*.trim()
        if (!bins || bins.any { !it.isInteger() }) error "Invalid static_quantized_quals '${params.static_quantized_quals}': must be a comma-separated list of integers, e.g. '10,20,30'"
        if (!params.preserve_qscores_less_than.toString().isInteger()) error "Invalid preserve_qscores_less_than '${params.preserve_qscores_less_than}': must be an integer"
    }
    ['clip_r1', 'clip_r2', 'three_prime_clip_r1', 'three_prime_clip_r2', 'length_required'].each { p ->
        if (!params[p].toString().isInteger() || params[p].toString().toInteger() < 0) error "Invalid ${p} '${params[p]}': must be a non-negative integer"
    }
}

def splitByReadKind(Map meta, Object reads) {
    def files = reads instanceof List ? reads : [reads]
    def lanes = []
    def offset = 0
    meta.lane_single_end.eachWithIndex { se, i ->
        def per = se ? 1 : 2
        lanes << [se: se, read_group: meta.read_groups[i], files: files.subList(offset, offset + per)]
        offset += per
    }
    [false, true].collect { se ->
        def own = lanes.findAll { lane -> lane.se == se }
        [meta + [part: se ? 'singleton' : 'paired', single_end: se, lane_single_end: own.collect { lane -> lane.se },
                 read_groups: own.collect { lane -> lane.read_group }, parent: meta],
         own.collectMany { lane -> lane.files }]
    }
}
