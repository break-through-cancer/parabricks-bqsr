process QUANTIZE_QUALS {
    tag "${meta.sample}"
    container params.quantize_quals_container
    cpus 8
    memory '8 GB'

    input:
    tuple val(meta), path(bam, stageAs: 'input/*')
    tuple path(ref_fasta), path(ref_fasta_fai)
    val output_fmt

    output:
    tuple val(meta), path("${meta.sample}.${meta.suffix}.${output_fmt}"), path("${meta.sample}.${meta.suffix}.${output_fmt}.${output_fmt == 'cram' ? 'crai' : 'bai'}"), emit: alignment
    path "${meta.sample}.quantize.log", emit: log
    path "${meta.sample}.quantize_mqc.tsv", emit: mqc
    path "${meta.sample}.quantize_hist_mqc.tsv", emit: hist

    script:
    def args = task.ext.args ?: ''
    def out = "${meta.sample}.${meta.suffix}.${output_fmt}"
    """
    set -euo pipefail
    quantize_quals \\
        --in ${bam} \\
        --out ${out} \\
        --ref ${ref_fasta} \\
        --threads ${task.cpus} \\
        ${args} 2> ${meta.sample}.quantize.log || { cat ${meta.sample}.quantize.log >&2; exit 1; }
    cat ${meta.sample}.quantize.log >&2
    {
        echo "# id: 'quantize_quals'"
        echo "# section_name: 'Quality-score quantization'"
        echo "# plot_type: 'table'"
        printf 'Sample\\tRecords\\tQualities\\tChanged\\tChanged_pct\\n'
        sed -n 's/^quantize_quals: done: \\([0-9]*\\) records ([0-9]* without qualities), \\([0-9]*\\) qualities, \\([0-9]*\\) changed (\\([0-9.]*\\)%).*/${meta.sample}\\t\\1\\t\\2\\t\\3\\t\\4/p' ${meta.sample}.quantize.log
    } > ${meta.sample}.quantize_mqc.tsv
    # Before/after per-quality-bin percentages, already computed and logged by the tool but
    # otherwise discarded. Zero-padded to 2 digits so plain (non-numeric) `sort`/`join` agree
    # with numeric order; `join` needs both sides sorted and compares byte-for-byte.
    awk '
        /^quantize_quals: input qualities:/ {
            line=\$0; sub(/^quantize_quals: input qualities: */, "", line)
            n = split(line, toks, " ")
            for (i=1;i<=n;i++) { split(toks[i], kv, ":"); q=kv[1]; sub(/^Q/,"",q); pct=kv[2]; sub(/%\$/,"",pct); printf "%02d\\t%s\\n", q, pct > "before.tsv" }
        }
        /^quantize_quals: output qualities:/ {
            line=\$0; sub(/^quantize_quals: output qualities: */, "", line)
            n = split(line, toks, " ")
            for (i=1;i<=n;i++) { split(toks[i], kv, ":"); q=kv[1]; sub(/^Q/,"",q); pct=kv[2]; sub(/%\$/,"",pct); printf "%02d\\t%s\\n", q, pct > "after.tsv" }
        }
    ' ${meta.sample}.quantize.log
    sort -t '\t' -k1,1 before.tsv -o before.sorted.tsv
    sort -t '\t' -k1,1 after.tsv -o after.sorted.tsv
    {
        echo "# id: 'quantize_quals_hist'"
        echo "# section_name: 'Quality-score distribution, before vs after quantization'"
        echo "# plot_type: 'linegraph'"
        echo "# pconfig:"
        echo "#     xlab: 'Quality score'"
        echo "#     ylab: '% of bases'"
        printf 'Category\\tbefore\\tafter\\n'
        join -t '\t' -a1 -a2 -e 0 -o 0,1.2,2.2 before.sorted.tsv after.sorted.tsv | awk -F'\\t' '{sub(/^0/,"",\$1); print "Q"\$1"\\t"\$2"\\t"\$3}'
    } > ${meta.sample}.quantize_hist_mqc.tsv
    """

    stub:
    def idx = output_fmt == 'cram' ? 'crai' : 'bai'
    """
    touch ${meta.sample}.${meta.suffix}.${output_fmt} ${meta.sample}.${meta.suffix}.${output_fmt}.${idx} ${meta.sample}.quantize.log
    printf 'Sample\\tRecords\\n${meta.sample}\\t0\\n' > ${meta.sample}.quantize_mqc.tsv
    printf 'Category\\tbefore\\tafter\\n' > ${meta.sample}.quantize_hist_mqc.tsv
    """
}
