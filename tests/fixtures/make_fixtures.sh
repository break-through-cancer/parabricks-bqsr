#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
SAMTOOLS=community.wave.seqera.io/library/htslib_samtools@sha256:a55ddea590e567a91df592300a960aa534cfc1bd16e7623e3938ec21f4f3df15
BWA=quay.io/biocontainers/bwa:0.7.18--he4a0461_1
run() { docker run --rm --platform linux/amd64 -v "$PWD":/w -w /w "$@"; }

mkdir -p genome/bwa genome/bwa_mismatch fastq
cp ../../tools/quantize_quals/tests/fixtures/ref.fa genome/reference.fasta
run "$SAMTOOLS" samtools faidx genome/reference.fasta
cp genome/reference.fasta genome/bwa/reference.fasta.64
run "$BWA" bwa index genome/bwa/reference.fasta.64
rm genome/bwa/reference.fasta.64
sed 's/^>chrT/>chrOther/' genome/reference.fasta > genome/bwa_mismatch/other.fasta
run "$BWA" bwa index genome/bwa_mismatch/other.fasta
rm genome/bwa_mismatch/other.fasta

printf '##fileformat=VCFv4.2\n##contig=<ID=chrT,length=500>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrT\t100\trs1\tA\tG\t.\tPASS\t.\n' > genome/known_sites.vcf
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrT\t200\trs2\tC\tT\t.\tPASS\t.\n' > genome/known_sites_nocontig.vcf
printf '##fileformat=VCFv4.2\n##contig=<ID=1,length=500>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\trs1\tA\tG\t.\tPASS\t.\n' > genome/known_sites_wrongcontig.vcf
printf '##fileformat=VCFv4.2\n##contig=<ID=chrT,length=999>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrT\t100\trs1\tA\tG\t.\tPASS\t.\n' > genome/known_sites_wronglength.vcf
for v in known_sites known_sites_nocontig known_sites_wrongcontig known_sites_wronglength; do
    run "$SAMTOOLS" bgzip -f genome/$v.vcf
done
run "$SAMTOOLS" tabix -f -p vcf genome/known_sites.vcf.gz
run "$SAMTOOLS" tabix -f -p vcf genome/known_sites_wrongcontig.vcf.gz
run "$SAMTOOLS" tabix -f -p vcf genome/known_sites_wronglength.vcf.gz
printf '##fileformat=VCFv4.2\n##contig=<ID=chrT,length=500>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrT\t300\trs3\tG\tA\t.\tPASS\t.\n' | gzip -c > genome/known_sites_plaingzip.vcf.gz
printf 'chrT\t0\t400\n' > genome/intervals.bed
printf '1\t0\t400\n' > genome/intervals_wrongcontig.bed

python3 - <<'PY'
import gzip, random
r = random.Random(20261003)
ref = ''.join(l.strip() for l in open('genome/reference.fasta') if not l.startswith('>'))
comp = str.maketrans('ACGT', 'TGCA')
def write(path, recs):
    with gzip.open(path, 'wt') as f:
        for name, seq in recs:
            f.write(f"@{name}\n{seq}\n+\n{'F' * len(seq)}\n")
def pairs(flowcell, lane, n):
    r1, r2 = [], []
    for i in range(n):
        p = r.randint(0, len(ref) - 200)
        frag = ref[p:p + 200]
        h = f"A00123:8:{flowcell}:{lane}:1101:{1000 + i}:2000"
        r1.append((f"{h} 1:N:0:ACGT", frag[:100]))
        r2.append((f"{h} 2:N:0:ACGT", frag[-100:].translate(comp)[::-1]))
    return r1, r2
a1, a2 = pairs('HFLOWCLXX', 1, 40); write('fastq/S1_L1_R1.fastq.gz', a1); write('fastq/S1_L1_R2.fastq.gz', a2)
b1, b2 = pairs('HFLOWCLXX', 2, 30); write('fastq/S1_L2_R1.fastq.gz', b1); write('fastq/S1_L2_R2.fastq.gz', b2)
c1, _ = pairs('HFLOWCLYY', 1, 25); write('fastq/S2_L1.fastq.gz', c1)
write('fastq/S3_sra.fastq.gz', [(f"SRR000001.{i} {i} length=100", ref[i:i + 100]) for i in range(20)])
PY

D=tests/fixtures/fastq
cat > fastq/valid.csv <<CSV
patient,sample,status,lane,fastq_1,fastq_2
P1,S1,0,1,$D/S1_L1_R1.fastq.gz,$D/S1_L1_R2.fastq.gz
P1,S1,0,2,$D/S1_L2_R1.fastq.gz,$D/S1_L2_R2.fastq.gz
P2,S2,1,1,$D/S2_L1.fastq.gz,
CSV
cat > fastq/defaults.csv <<CSV
sample,lane,fastq_1,fastq_2
S3,1,$D/S3_sra.fastq.gz,
CSV
cat > fastq/mixed_pe_se.csv <<CSV
patient,sample,status,lane,fastq_1,fastq_2
P1,S1,0,1,$D/S1_L1_R1.fastq.gz,$D/S1_L1_R2.fastq.gz
P1,S1,0,2,$D/S2_L1.fastq.gz,
CSV
cat > fastq/duplicate_lane.csv <<CSV
patient,sample,status,lane,fastq_1,fastq_2
P1,S1,0,1,$D/S1_L1_R1.fastq.gz,$D/S1_L1_R2.fastq.gz
P1,S1,0,1,$D/S1_L2_R1.fastq.gz,$D/S1_L2_R2.fastq.gz
CSV
cat > fastq/two_patients.csv <<CSV
patient,sample,status,lane,fastq_1,fastq_2
P1,S1,0,1,$D/S1_L1_R1.fastq.gz,$D/S1_L1_R2.fastq.gz
P9,S1,0,2,$D/S1_L2_R1.fastq.gz,$D/S1_L2_R2.fastq.gz
CSV
cat > fastq/missing_fastq.csv <<CSV
patient,sample,status,lane,fastq_1,fastq_2
P1,S1,0,1,$D/does_not_exist.fastq.gz,
CSV
cat > fastq/both_entries.csv <<CSV
sample,lane,fastq_1,fastq_2,alignment,alignment_index,recal_table
S1,1,$D/S1_L1_R1.fastq.gz,,x.bam,x.bam.bai,x.table
CSV
