# parabricks_bqsr

Nextflow pipeline for GPU alignment and base quality score recalibration with NVIDIA
Parabricks. Starting from FASTQs (the main entry point) or from an existing alignment plus
its BQSR table, it applies BQSR genome-wide with `pbrun applybqsr` (never restricted to
intervals, so unmapped and off-target reads are recalibrated too) and can replicate GATK
`ApplyBQSR`'s static quality-score quantization, a step Parabricks does not provide.

```
FASTQ samplesheet                          alignment + table samplesheet
      │                                          │
FASTP (per lane: read counts, read QC,           │
       optional trimming)                        │
      │                                          │
PARABRICKS_FQ2BAM (per sample, GPU:              │
  alignment, duplicate marking, BQSR table)      │
      │                                          │
      ├─ known sites ─▶ PARABRICKS_APPLYBQSR (GPU, genome-wide) ◀─┘
      │                        │
      │                QUANTIZE_QUALS or SAMTOOLS_FINALIZE
      │                        └─▶ preprocessing/recalibrated/<s>/<s>.recal.<cram|bam>
      └─ no known sites ─▶ QUANTIZE_QUALS or SAMTOOLS_FINALIZE
                               └─▶ preprocessing/markduplicates/<s>/<s>.md.<cram|bam>

QC: samtools stats, mosdepth (skipped for CRAM 3.1), read-count checks, MultiQC
```

## Usage

Requires Nextflow 26.04.0 or later.

```bash
nextflow run main.nf --input fastq_samplesheet.csv
nextflow run main.nf --input alignment_samplesheet.csv --quantize_quals_enabled false --output_fmt bam
```

`--input` is the only required parameter. The samplesheet header selects the entry point:
a `fastq_1` column means FASTQ input, an `alignment` column means alignment + table
input; both together is an error. Invalid parameters, missing samplesheet values,
nonexistent paths and mismatched reference files fail at startup with a message naming
the problem.

`--input` may be local or `s3://`; a relative `--input` resolves against the launch
directory. Paths inside samplesheets may be `s3://` and are checked at startup. Use
absolute or `s3://` paths inside samplesheets for AWS Batch runs: relative row paths
resolve against the pipeline directory (`projectDir`), which suits the test fixtures only.

### FASTQ samplesheet

One row per lane. See `assets/samplesheet_fastq.csv`.

| Column | Required | Meaning |
| --- | --- | --- |
| `patient` | no (defaults to `sample`) | Groups samples; part of the read group `SM`. |
| `sample` | yes | Sample ID. |
| `status` | no (defaults to `0`) | `0` normal, `1` tumor (adds BWA `-B 3`). |
| `lane` | yes | Lane identifier, unique within a sample. |
| `fastq_1` | yes | Read 1, or the single-end FASTQ. |
| `fastq_2` | no | Read 2; empty means single-end. |

Read groups follow nf-core/sarek: `ID` and `PU` are `<flowcell>.<sample>.<lane>` (flowcell
from the first Illumina read header, `unknown` otherwise), `SM` is `<patient>_<sample>`,
`LB` is `<sample>`, `PL` is `--seq_platform`. All lanes of a sample go into one `fq2bam`
call, giving one alignment and one BQSR table per sample. A sample cannot mix paired-end
and single-end lanes.

### Alignment + table samplesheet

One row per sample. See `assets/samplesheet.csv`.

| Column | Meaning |
| --- | --- |
| `sample` | Sample ID; unique within the samplesheet. |
| `alignment` | Pre-BQSR BAM or CRAM, e.g. the fq2bam output of `sarek_align`. |
| `alignment_index` | Its `.bai` / `.crai`. |
| `recal_table` | The BQSR table for that alignment. |

The alignment header contigs must match the selected genome; otherwise the task stops
before `applybqsr` runs.

### Genomes

`--genome` selects a built-in genome; explicit parameters override its defaults.

- **`--genome GATK.GRCh38` (default):** FASTA, BWA index, known sites (dbSNP 146, Mills
  and 1000G gold-standard indels, known indels) and WGS calling intervals from iGenomes
  under `--igenomes_base` (default `s3://ngi-igenomes/igenomes/`).
- **Custom genome (`--genome null`):** `--ref_fasta` (and `--ref_fasta_fai`, default
  `<fasta>.fai`), `--bwa_index` (directory holding a classic BWA index of that FASTA),
  optional `--known_sites` (comma-separated `.vcf.gz`), optional `--intervals`.

**Resolution:** every reference takes the explicit parameter if given, else the genome
default. `--no_intervals` turns intervals off, including a genome default; combining it
with an explicit `--intervals` is an error.

**Known sites:** only the `.vcf.gz` files are passed to `fq2bam`. A VCF without a `.tbi`
is indexed automatically (it must be bgzip-compressed). With no known sites, BQSR is
skipped and outputs are published as `.md` files under `preprocessing/markduplicates/`.

**Intervals** restrict BQSR table building only; reads outside them are still aligned,
output and recalibrated. Without intervals the table is built genome-wide, with a warning.

**Startup checks** (headers and index files only): the BWA index matches the FASTA; each
known-sites VCF and the intervals share contig names with the FASTA (zero overlap fails;
partial VCF overlap warns).

### GPU

`PARABRICKS_FQ2BAM` requests `--fq2bam_gpus` GPUs (default 1); `PARABRICKS_APPLYBQSR`
requests 1 (`pbrun applybqsr` accepts only 1 or 2). On AWS Batch, `accelerator` places
the task on a GPU compute environment; with the local executor and Docker, `--gpus all`
is added. Local runs need a Linux host with an NVIDIA GPU.

### Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--outdir` | `./results` | Output directory. |
| `--genome` | `GATK.GRCh38` | Built-in genome, or `null` for a custom genome. |
| `--igenomes_base` | `s3://ngi-igenomes/igenomes/` | Base path for built-in genome files. |
| `--ref_fasta`, `--ref_fasta_fai`, `--bwa_index`, `--known_sites`, `--intervals` | from `--genome` | Explicit reference overrides. |
| `--no_intervals` | `false` | Build the BQSR table genome-wide. |
| `--seq_platform` | `ILLUMINA` | Read group `PL`. |
| `--trim_fastq` | `false` | Enable fastp trimming (fastp always runs for counts and QC). |
| `--clip_r1`, `--clip_r2` | `0` | Bases removed from the 5′ end of read 1 / read 2. |
| `--three_prime_clip_r1`, `--three_prime_clip_r2` | `0` | Bases removed from the 3′ end. |
| `--trim_nextseq` | `false` | Trim poly-G tails (two-colour instruments). |
| `--length_required` | `15` | Minimum read length after trimming. |
| `--save_trimmed` | `false` | Publish trimmed FASTQs. |
| `--markdups_se_mode` | `5prime` | Single-end duplicate marking: `5prime` (standard) or `start-end` (adapter-trimmed short fragments such as cfDNA). |
| `--optical_duplicate_pixel_distance` | `100` | Optical-duplicate metrics only; 2500 is usual for patterned flowcells. |
| `--fq2bam_gpus` | `1` | GPUs for alignment. |
| `--fq2bam_low_memory` | `true` | `--low-memory` for fq2bam (one BWA stream per GPU, fits a 24 GB GPU). Off can be faster on larger GPUs or fail on smaller ones. |
| `--fq2bam_gpuwrite` | `true` | `--gpuwrite` for fq2bam. With one GPU it shares the device with BQSR. |
| `--output_fmt` | `cram` | `bam` or `cram`. |
| `--cram_version` | `3.0` | CRAM version for CRAM output. `3.0` is readable by essentially all tools; `3.1` is smaller, but older HTSlib builds and htsjdk-based tools may not read it, and mosdepth coverage QC is skipped. |
| `--quantize_quals_enabled` | `true` | Quantize quality scores; `false` publishes unquantized output. |
| `--static_quantized_quals` | `10,20,30` | Comma-separated static bins. |
| `--preserve_qscores_less_than` | `6` | Qualities below this value remain unchanged. |
| `--round_down_quantized` | `false` | Round down to a bin instead of the nearest bin in probability space. |
| `--quantize_quals_container` | `ghcr.io/break-through-cancer/parabricks-bqsr:0.1.2` | Quantizer image; see `tools/quantize_quals/README.md`. |

Parameters are declared in `nextflow_schema.json` (nf-schema); the startup log prints the
parameters that differ from their defaults, as nf-core pipelines do. A new parameter must be
added to both `nextflow.config` and the schema; `tests/config` checks they match.

Quality filtering in fastp is always disabled (BQSR handles base qualities); with
trimming on, `--length_required` is the only filter that removes reads.

### Outputs

```
preprocessing/recalibrated/<s>/<s>.recal.<cram|bam> (+ index)   BQSR ran
preprocessing/markduplicates/<s>/<s>.md.<cram|bam> (+ index)    no known sites
preprocessing/recal_table/<s>/<s>.table                          FASTQ entry with known sites
preprocessing/fastp/<s>/                                         --save_trimmed only
reports/fastp/<s>/, reports/markduplicates/<s>/, reports/parabricks_qc/<s>/,
reports/samtools/<s>/, reports/mosdepth/<s>/, reports/quantize/<s>/
multiqc/multiqc_report.html
pipeline_info/
```

The layout matches nf-core/sarek, so sarek-based variant calling (including Cirro's
`sarek_call_variants`) recognises the outputs. The pre-BQSR `fq2bam` CRAM and the
intermediate `applybqsr` BAM are not published.

### Read-count checks

Each run fails if reads go missing:
- **FASTQ → final:** fastp read total (before filtering, or after filtering when trimming)
  equals the final file's primary reads (`samtools stats` raw total).
- **Before → final:** records in the `fq2bam` output (FASTQ entry) or the input alignment
  (alignment entry) equal the final file's records.

The FASTQ check names `fq2bam`'s 480 bp maximum read length and zero-length reads as the
likely causes of a mismatch. Results appear in the log and in MultiQC.

## Testing

```bash
make -C tools/quantize_quals test      # C unit tests + CLI tests (needs HTSlib, samtools)
nf-test test                           # modules, subworkflows, pipeline (needs Docker)
python -m pytest tests/config .cirro   # GPU config checks, Cirro preprocess (needs nextflow, pandas)
```

Parabricks steps run under `-stub` locally; real execution needs an NVIDIA GPU. nf-test
pulls the quantizer image from GHCR; to test local changes to the tool, build and tag it
under the same name first:

```bash
docker build --platform linux/amd64 -t ghcr.io/break-through-cancer/parabricks-bqsr:0.1.2 tools/quantize_quals
```

Output content is decoded inside nf-test with the
[nft-bam](https://github.com/nvnieuwk/nft-bam) plugin. All fixtures are small synthetic
files (`tests/fixtures/make_fixtures.sh` regenerates the genome and FASTQ fixtures).

## Cirro

Two process registrations come from this repository:

| Directory | Entry point | Input dataset |
| --- | --- | --- |
| `.cirro/fastq/` | FASTQ | Paired or single-end FASTQ datasets |
| `.cirro/alignment/` | alignment + table | A `sarek_align` dataset with `preprocessing/parabricks/<sample>/<sample>.{bam,bam.bai,table}` (Parabricks aligner, known sites, `save_mapped` on, `baserecalibrator` skipped) |

Each directory holds `process-form.json`, `process-input.json`, `preprocess.py`,
`process-compute.config` and `process-output.json`. Both Parabricks processes run on the
on-demand GPU queue (`PW_ONDEMAND_JOB_QUEUE`) with retries on resource-related exit codes.

**FASTQ form:** genome source (iGenomes GATK.GRCh38 with an intervals checkbox, or a
custom genome: a BWA index dataset containing `genome.fasta`, optional known-sites VCFs
from the references library under `germline_resource`, optional BED under `genome_bed`),
output format, quantization, trimming, single-end duplicate marking, optical pixel
distance and alignment GPUs.

**Registration settings:** repository `break-through-cancer/parabricks-bqsr`, entry
script `main.nf`, configuration directory `.cirro/fastq` or `.cirro/alignment`, Nextflow
`26.04.0` or later (required; nf-schema 2.8.0 needs it). A registration created
before `.cirro/alignment/` existed must be re-pointed to that directory. Output file
mapping can reuse `sarek_align`'s patterns:
- `preprocessing/(?P<bamType>recalibrated|markduplicates)/(?P<sampleName>[^/]+)/[^/]+\.(?:bam|cram)$`
- `preprocessing/(?P<bamType>recalibrated|markduplicates)/(?P<sampleName>[^/]+)/[^/]+\.(?:bam\.bai|cram\.crai)$`

## Known gaps

1. **The FASTQ entry point has not run on a GPU.** Validation on Cirro, in order:
   - V1: a low-pass human sample from FASTQ (read checks, `fq2bam` writing `.crai`, QC
     directory contents, `--memory-limit`, resource sizing).
   - V2: a single-lane sample compared with `sarek_align` on the same FASTQs (read groups,
     mapping and duplicate rates, identical table, matching pre-BQSR records).
   - V3: a single-end sample with both duplicate-marking modes.
   - V4: one `fq2bam` call with paired and single-end inputs (decides whether the
     mixed-sample rejection can be lifted).
   - V5: a canine custom genome with and without known sites.
   - V6: low-pass tables (observations per read group below 1x).
   - V7: quantization parity with GATK `ApplyBQSR --static-quantized-quals`.
2. **`stageInMode 'copy'`** is inherited from the nf-core Parabricks modules;
   `--preserve-file-symlinks` may make the copy unnecessary.
3. **`SAMTOOLS_FINALIZE` exists because `applybqsr` writes BAM only** (NVIDIA documents
   its `--out-bam` as "Output BAM file"). With quantization on, `QUANTIZE_QUALS` writes
   CRAM directly and this step does not run.
