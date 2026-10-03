process TABIX_VCF {
    tag "${vcf.name}"
    container 'community.wave.seqera.io/library/htslib_samtools@sha256:a55ddea590e567a91df592300a960aa534cfc1bd16e7623e3938ec21f4f3df15'
    cpus 1
    memory '2 GB'

    input:
    tuple val(i), path(vcf)

    output:
    tuple val(i), path(vcf), path("${vcf}.tbi")

    script:
    """
    set -euo pipefail
    if ! htsfile ${vcf} | grep -q BGZF; then
        echo "Known-sites VCF ${vcf} is not bgzip-compressed; re-compress it with bgzip (e.g. zcat ${vcf} | bgzip > fixed.vcf.gz) and index it with tabix" >&2
        exit 1
    fi
    tabix -p vcf ${vcf}
    """

    stub:
    """
    touch ${vcf}.tbi
    """
}
