process MULTIQC {
    container 'community.wave.seqera.io/library/multiqc@sha256:ae074e961979ea85c72c6c67a8ed432fa3368c736b8fe37799f7223901ac84ab'
    cpus 2
    memory '4 GB'

    input:
    path files, stageAs: 'inputs/?/*'

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_report_data', emit: data

    script:
    """
    multiqc --force --filename multiqc_report.html inputs
    """

    stub:
    """
    touch multiqc_report.html
    mkdir multiqc_report_data
    """
}
