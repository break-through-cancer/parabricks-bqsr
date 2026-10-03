def genomes() {
    [
    'GATK.GRCh38': [
        ref_fasta  : 'Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta',
        bwa_index  : 'Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/',
        known_sites: [
            'Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz',
            'Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz',
            'Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/beta/Homo_sapiens_assembly38.known_indels.vcf.gz'
        ],
        intervals  : 'Homo_sapiens/GATK/GRCh38/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed'
    ]
    ]
}

def resolveReferences(Map p) {
    def genome = p.genome in [null, '', 'null', false] ? null : p.genome.toString()
    def defaults = [:]
    if (genome) {
        if (!genomes().containsKey(genome)) {
            error "Unknown genome '${genome}'. Built-in genomes: ${genomes().keySet().join(', ')}. Use --genome null with --ref_fasta and --bwa_index for a custom genome."
        }
        def base = p.igenomes_base.toString().replaceAll('/+$', '')
        def g = genomes()[genome]
        defaults = [
            ref_fasta  : "${base}/${g.ref_fasta}".toString(),
            bwa_index  : "${base}/${g.bwa_index}".toString(),
            known_sites: g.known_sites.collect { "${base}/${it}".toString() },
            intervals  : "${base}/${g.intervals}".toString()
        ]
    }
    if (p.no_intervals && p.intervals) {
        error "--intervals and --no_intervals cannot be combined"
    }
    def refFasta = p.ref_fasta ?: defaults.ref_fasta
    if (!refFasta) {
        error "--ref_fasta is required with a custom genome (--genome null)"
    }
    def knownSites = p.known_sites ? p.known_sites.toString().tokenize(',')*.trim().findAll() : (defaults.known_sites ?: [])
    def intervals = p.no_intervals ? null : (p.intervals ?: defaults.intervals)
    [
        genome          : genome,
        ref_fasta       : refFasta.toString(),
        ref_fasta_fai   : (p.ref_fasta_fai ?: "${refFasta}.fai").toString(),
        bwa_index       : (p.bwa_index ?: defaults.bwa_index)?.toString(),
        known_sites     : knownSites,
        intervals       : intervals?.toString(),
        intervals_source: p.no_intervals ? 'disabled' : p.intervals ? 'explicit' : defaults.intervals ? 'genome' : 'none'
    ]
}

def intervalsMessage(Map refs) {
    if (refs.intervals_source == 'disabled') return "--no_intervals set: building the BQSR table genome-wide by request"
    if (refs.intervals_source == 'none') return "No intervals supplied: building the BQSR table genome-wide"
    "Building the BQSR table on intervals ${refs.intervals}".toString()
}

def readerFor(Object path) {
    def f = toPath(path)
    def stream = f.newInputStream()
    new BufferedReader(new InputStreamReader(f.name.endsWith('.gz') ? new java.util.zip.GZIPInputStream(stream) : stream))
}

def faiContigs(Object path) {
    def contigs = [:]
    toPath(path).eachLine { line ->
        def t = line.tokenize('\t')
        if (t.size() >= 2) contigs[t[0]] = t[1] as Long
    }
    contigs
}

def annContigs(Object path) {
    def lines = toPath(path).readLines()
    def contigs = [:]
    (1..<lines.size()).step(2).findAll { i -> i + 1 < lines.size() }.each { i ->
        contigs[lines[i].tokenize(' ')[1]] = lines[i + 1].tokenize(' ')[1] as Long
    }
    contigs
}

def vcfContigs(Object path) {
    def lines = readerFor(path).withCloseable { reader -> reader.iterator().take(50000).toList() }
    def header = lines.findAll { it.startsWith('##contig=') }
        .collect { line -> def m = line =~ /ID=([^,>]+)/; m.find() ? m.group(1) : null }
        .findAll() as Set
    def sampled = lines.findAll { !it.startsWith('#') }.take(10000).collect { it.tokenize('\t')[0] } as Set
    def lengths = lines.findAll { it.startsWith('##contig=') }
        .collect { line ->
            def id = line =~ /ID=([^,>]+)/
            def len = line =~ /length=(\d+)/
            id.find() && len.find() ? [id.group(1), len.group(1) as Long] : null
        }
        .findAll()
        .collectEntries()
    header ? [contigs: header, complete: true, lengths: lengths] : [contigs: sampled, complete: false, lengths: [:]]
}

def intervalContigs(Object path) {
    def contigs = [] as Set
    readerFor(path).withCloseable { reader ->
        reader.eachLine { line ->
            if (!line || line.startsWith('@') || line.startsWith('#') || line.startsWith('track') || line.startsWith('browser')) return
            def m = line =~ /^([^\t:]+)(?::\d+(?:-\d+)?)?(?:\t|$)/
            if (m.find()) contigs << m.group(1)
        }
    }
    contigs
}

def checkReferences(Map refs, boolean fastqEntry) {
    def errors = []
    def warnings = []
    def fai = faiContigs(refs.ref_fasta_fai)
    if (!fastqEntry) return [errors: errors, warnings: warnings]

    if (!refs.bwa_index) {
        errors << "--bwa_index is required for FASTQ input"
    } else {
        def amb = file("${refs.bwa_index.toString().replaceAll('/+$', '')}/*.amb")
        def ambs = amb instanceof List ? amb : [amb].findAll { it.exists() }
        if (ambs.size() != 1) {
            errors << "BWA index ${refs.bwa_index}: expected exactly one *.amb file, found ${ambs.size()}"
        } else if (annContigs(ambs[0].resolveSibling(ambs[0].name.replaceAll(/\.amb$/, '.ann'))) != fai) {
            errors << "BWA index ${refs.bwa_index} does not match ${refs.ref_fasta}: contig names or lengths differ"
        }
    }
    refs.known_sites.each { vcf ->
        def v = vcfContigs(vcf)
        def shared = v.contigs.intersect(fai.keySet())
        def wrong = shared.find { c -> v.lengths[c] != null && v.lengths[c] != fai[c] }
        if (wrong) {
            errors << "Known-sites VCF ${vcf} does not match ${refs.ref_fasta}: contig '${wrong}' has length ${v.lengths[wrong]} in the VCF but ${fai[wrong]} in the FASTA (a different assembly? use --genome null with matching --known_sites for a custom genome)"
        } else if (!shared) {
            errors << "Known-sites VCF ${vcf} has no contig names in common with ${refs.ref_fasta} (e.g. VCF '${v.contigs.take(3).join("', '")}' vs FASTA '${fai.keySet().take(3).join("', '")}')"
        } else if (v.complete && shared.size() < v.contigs.size()) {
            warnings << "Known-sites VCF ${vcf}: ${v.contigs.size() - shared.size()} of ${v.contigs.size()} contigs are not in ${refs.ref_fasta}"
        }
    }
    if (refs.intervals) {
        def bed = intervalContigs(refs.intervals)
        if (!bed.intersect(fai.keySet())) {
            errors << "Intervals ${refs.intervals} have no contig names in common with ${refs.ref_fasta}"
        }
    }
    [errors: errors, warnings: warnings]
}

def toPath(Object p) {
    p instanceof CharSequence ? file(p.toString()) : p
}
