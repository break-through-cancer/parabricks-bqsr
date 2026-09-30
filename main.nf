nextflow.enable.dsl = 2

include { SAMPLESHEET_TO_SAMPLES } from './subworkflows/local/samplesheet_to_samples/main'
include { PARABRICKS_APPLYBQSR   } from './modules/local/parabricks_applybqsr'
include { QUANTIZE_QUALS         } from './modules/local/quantize_quals'
include { SAMTOOLS_FINALIZE      } from './modules/local/samtools_finalize'

workflow {
    if (!params.input) {
        error "Missing required parameter: input (path to samplesheet CSV). See README.md and assets/samplesheet.csv."
    }
    if (!(params.output_fmt in ['bam', 'cram'])) {
        error "Invalid output_fmt '${params.output_fmt}': must be 'bam' or 'cram'"
    }
    if (!(params.quantize_quals_enabled.toString() in ['true', 'false'])) {
        error "Invalid quantize_quals_enabled '${params.quantize_quals_enabled}': must be true or false"
    }
    def quantize = params.quantize_quals_enabled.toString() == 'true'
    if (quantize) {
        def bins = params.static_quantized_quals.toString().tokenize(',')*.trim()
        if (!bins || bins.any { !it.isInteger() }) {
            error "Invalid static_quantized_quals '${params.static_quantized_quals}': must be a comma-separated list of integers, e.g. '10,20,30'"
        }
        if (!params.preserve_qscores_less_than.toString().isInteger()) {
            error "Invalid preserve_qscores_less_than '${params.preserve_qscores_less_than}': must be an integer"
        }
    }

    SAMPLESHEET_TO_SAMPLES(params.input)

    ref_ch = Channel.value([
        file(params.ref_fasta, checkIfExists: true),
        file(params.ref_fasta_fai, checkIfExists: true)
    ])

    PARABRICKS_APPLYBQSR(SAMPLESHEET_TO_SAMPLES.out.samples, ref_ch)

    if (quantize) {
        QUANTIZE_QUALS(PARABRICKS_APPLYBQSR.out.bam, ref_ch, params.output_fmt)
    } else {
        SAMTOOLS_FINALIZE(PARABRICKS_APPLYBQSR.out.bam, ref_ch, params.output_fmt)
    }
}
