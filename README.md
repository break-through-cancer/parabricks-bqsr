# parabricks-fq2bam-bqsr

Nextflow pipeline for GPU alignment and base quality score recalibration with NVIDIA
Parabricks. Starting from FASTQs (the main entry point) or from an existing alignment plus
its BQSR table, it applies BQSR genome-wide with `pbrun applybqsr` (never restricted to
intervals, so unmapped and off-target reads are recalibrated too) and can replicate GATK
`ApplyBQSR`'s static quality-score quantization, a step Parabricks does not provide.
Measured on a 105x WGS sample: 7.2h / $6.99, against 3d15h / $115.80 for sarek/CPU on the
same sample.

```
FASTQ samplesheet                          alignment + table samplesheet
      │                                          │
FASTP (per lane: read counts, read QC,           │
       optional trimming)                        │
      │                                          │
PARABRICKS_FQ2BAM (per sample, GPU:              │
  alignment, duplicate marking, BQSR table)      │
      │                                          │
      ├─ known sites, apply_bqsr ─▶ PARABRICKS_APPLYBQSR (GPU, genome-wide) ◀─┘
      │                        │
      │                QUANTIZE_QUALS or SAMTOOLS_FINALIZE
      │                        └─▶ preprocessing/recalibrated/<s>/<s>.recal.<cram|bam>
      └─ no known sites, or apply_bqsr=false ─▶ QUANTIZE_QUALS or SAMTOOLS_FINALIZE
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
call, giving one alignment and one BQSR table per sample. A sample may mix paired lanes and
single-end lanes (e.g. singletons left over from a BAM-to-FASTQ conversion); they share the
sample's `SM` and `LB`, and `--markdups_se_mode` applies to the single-end reads. `fq2bam`
takes one input kind per call, so such a sample is aligned in two `fq2bam` calls (paired,
single-end) without duplicate marking. The two are merged and queryname-sorted (`bamsort`), then
duplicates are marked once across all reads (`markdup`, matching `fq2bam`'s GATK
MarkDuplicates behaviour). The BQSR table (`bqsr`) and QC metrics (`collectmultiplemetrics`)
are built from that file. Outputs match a single `fq2bam` call.

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

`PARABRICKS_FQ2BAM` requests `--fq2bam_gpus` GPUs (default 2); `PARABRICKS_APPLYBQSR`
requests 1 (`pbrun applybqsr` accepts only 1 or 2). On AWS Batch, `accelerator` places
the task on a GPU compute environment; with the local executor and Docker, `--gpus all`
is added. Local runs need a Linux host with an NVIDIA GPU.

### Settings by sequencing depth

`fq2bam` is the costly step. Its host memory peaks in its last phase (duplicate marking,
BQSR table, writing), and that peak grows with depth, not with `--fq2bam_gpus` or
`--fq2bam_low_memory` (those change BWA speed and BWA memory only). A memory kill restarts
`fq2bam` from the beginning; retries use 2× then 3× the memory. For deep samples, set
`--fq2bam_memory_gb` (on Cirro: *fq2bam memory override*) so the first attempt succeeds
instead of paying for one that doesn't.

| Depth (WGS) | GPUs | Low-memory | `--fq2bam_memory_gb` |
|---|---|---|---|
| Low-pass (≤5x) | 1 | off | unset (64 GB) |
| About 20–30x | 2 (default) | off | unset (88 GB) |
| About 60x | 2 | off | 176 (estimate) |
| About 100x and deeper | 4 | off | 352 |

Placement on g5: a g5.12xlarge has 4 GPUs, 48 vCPUs and 192 GB, and holds two default
2-GPU jobs. A 352 GB job needs a whole g5.24xlarge regardless of GPU count, so 4 GPUs is
the better choice at that memory tier: same cost, faster alignment.

### Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--outdir` | `./results` | Output directory. |
| `--genome` | `GATK.GRCh38` | Built-in genome, or `null` for a custom genome. |
| `--igenomes_base` | `s3://ngi-igenomes/igenomes/` | Base path for built-in genome files. |
| `--ref_fasta`, `--ref_fasta_fai`, `--bwa_index`, `--known_sites`, `--intervals` | from `--genome` | Explicit reference overrides. |
| `--no_intervals` | `false` | Build the BQSR table genome-wide. |
| `--seq_platform` | `ILLUMINA` | Read group `PL`. |
| `--trim_fastq` | `false` | Enable fastp adapter trimming and clipping (fastp always runs for counts and QC). |
| `--clip_r1`, `--clip_r2` | `0` | Bases removed from the 5′ end of read 1 / read 2. |
| `--three_prime_clip_r1`, `--three_prime_clip_r2` | `0` | Bases removed from the 3′ end. |
| `--poly_g_trimming` | `auto` | Poly-G tail trimming by fastp, applied whether or not `--trim_fastq` is set: `auto` (fastp detects two-colour instruments such as NovaSeq and NextSeq), `on` or `off`. When not `off`, fq2bam aligns fastp's output instead of the original FASTQs. |
| `--length_required` | `15` | Minimum read length after trimming (including poly-G trimming). |
| `--save_trimmed` | `false` | Publish trimmed FASTQs. |
| `--markdups_se_mode` | `5prime` | Single-end duplicate marking: `5prime` (standard) or `start-end` (adapter-trimmed short fragments such as cfDNA). |
| `--optical_duplicate_pixel_distance` | `100` | Optical-duplicate metrics only; 2500 is usual for patterned flowcells. |
| `--fq2bam_gpus` | `2` | GPUs for alignment, 1–4; set 1 to minimise cost, more for speed. CPUs and memory scale with it: 12 CPUs and 44 GB per GPU, at least 16 CPUs / 64 GB. |
| `--fq2bam_low_memory` | `false` | `--low-memory` for fq2bam (one BWA stream per GPU). Turn on if a smaller GPU runs out of memory. |
| `--fq2bam_gpuwrite` | `true` | `--gpuwrite` for fq2bam. Not shown on the Cirro form. |
| `--fq2bam_intermediate_fmt` | `bam` | Format of the duplicate-marked alignment passed from fq2bam (or markdup) to applybqsr: `bam` or `cram`. Not published; final outputs follow `--output_fmt`. Not shown on the Cirro form. |
| `--fq2bam_memory_gb` | unset | Override `fq2bam`'s first-attempt host memory in GB (16–768); retries multiply it by the attempt number and `--memory-limit` stays at half. Unset: 44 GB per GPU, at least 64. See *Settings by sequencing depth*. |
| `--apply_bqsr` | `true` | FASTQ entry with known sites: run `applybqsr` (and quantization) in this same run. `false` aligns only, publishing the markduplicates alignment and BQSR table as this run's output, to apply BQSR later with a separate alignment-entry run against them. No effect without known sites, or on the alignment entry (which always applies BQSR). |
| `--output_fmt` | `bam` | `bam` or `cram`. |
| `--cram_version` | `3.0` | CRAM version for CRAM output. `3.0` is readable by essentially all tools; `3.1` is smaller, but older HTSlib builds and htsjdk-based tools may not read it, and mosdepth coverage QC is skipped. |
| `--publish_markduplicates` | `false` | FASTQ entry with known sites only: also publish the pre-BQSR, duplicate-marked alignment as an indexed file in `--output_fmt`, alongside the BQSR table that's already published there. Lets a later alignment-entry run apply BQSR without repeating alignment. No effect without known sites. |
| `--quantize_quals_enabled` | `true` | Quantize quality scores; `false` publishes unquantized output. |
| `--static_quantized_quals` | `10,20,30` | Comma-separated static bins. |
| `--preserve_qscores_less_than` | `6` | Qualities below this value remain unchanged. |
| `--round_down_quantized` | `false` | Round down to a bin instead of the nearest bin in probability space. |
| `--quantize_quals_container` | `ghcr.io/break-through-cancer/parabricks-fq2bam-bqsr:0.1.2` | Quantizer image; see `tools/quantize_quals/README.md`. |

Parameters are declared in `nextflow_schema.json` (nf-schema); the startup log prints the
parameters that differ from their defaults, as nf-core pipelines do. A new parameter must be
added to both `nextflow.config` and the schema; `tests/config` checks they match.

Quality filtering in fastp is always disabled: GATK discourages quality trimming, as base
qualities are handled by soft-clipping, BQSR and the variant callers. Poly-G tails are
different: on two-colour instruments "no signal" is called as a high-quality G, producing
tails that alignment does not clip and that materially affect mapping and duplicate rates
untrimmed (see `--poly_g_trimming`). `--length_required` is the only filter that removes
reads (pairs emptied by trimming), and the read checks count fastp's output whenever fq2bam
aligns it.

When fq2bam aligns fastp's output, it waits for fastp to finish (a CPU task; no GPU is
held meanwhile). fastp runs with 16 threads, the most fastp 0.24 uses, and fast output
compression (`-z 1`).

### Outputs

```
preprocessing/recalibrated/<s>/<s>.recal.<cram|bam> (+ index)   BQSR ran
preprocessing/markduplicates/<s>/<s>.md.<cram|bam> (+ index)    no known sites, --apply_bqsr false, or --publish_markduplicates
preprocessing/recal_table/<s>/<s>.table                          FASTQ entry with known sites
preprocessing/fastp/<s>/                                         --save_trimmed only
reports/fastp/<s>/, reports/markduplicates/<s>/, reports/parabricks_qc/<s>/,
reports/samtools/<s>/ (<s>.stats, <s>.flagstat), reports/mosdepth/<s>/, reports/quantize/<s>/
multiqc/multiqc_report.html
pipeline_info/
```

The layout matches nf-core/sarek, so sarek-based variant calling (including Cirro's
`sarek_call_variants`) recognises the outputs. The intermediate `applybqsr` BAM is never
published. With `--apply_bqsr false` (known sites, FASTQ entry), the `fq2bam`/`markdup`
alignment and BQSR table *are* this run's output; otherwise that alignment is published only
via `--publish_markduplicates`, or discarded once BQSR has run.

When the alignment to publish is already in `--output_fmt` with no quantization pending, the
step that produced it indexes and publishes it directly; `SAMTOOLS_FINALIZE` then does not
appear in that run's task list at all. This is expected, not an error.

`reports/samtools/<s>/<s>.stats` reads the pre-quantization alignment rather than the
published file (BQSR and quantization only rewrite quality bytes), while `<s>.flagstat`
reads the true published file and backs the read-count check below. MultiQC also plots the
quality-score distribution before and after quantization (Before/After tabs, one stacked
bar per sample), from the `reports/quantize/<s>/` log.

### Read-count checks

Each run fails if reads go missing:
- **FASTQ → final:** fastp read total (before filtering, or after filtering when trimming)
  equals the final file's primary reads (`samtools flagstat` primary count, on the true
  published file).
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

Each nf-test case pays its own Nextflow engine start-up (several seconds), so the full suite
is slower than its test count suggests. `scripts/test-parallel.sh [N]` splits it across N
parallel `nf-test --shard` processes (default 4) for a full local run; extra arguments pass
through to nf-test, e.g. `scripts/test-parallel.sh 4 --tag quantize`. For fast iteration on a
single change, prefer `nf-test test --changedSince HEAD` or `--relatedTests` over a full run.

Parabricks steps run under `-stub` locally; real execution needs an NVIDIA GPU. nf-test
pulls the quantizer image from GHCR; to test local changes to the tool, build and tag it
under the same name first:

```bash
docker build --platform linux/amd64 -t ghcr.io/break-through-cancer/parabricks-fq2bam-bqsr:0.1.2 tools/quantize_quals
```

Output content is decoded inside nf-test with the
[nft-bam](https://github.com/nvnieuwk/nft-bam) plugin. All fixtures are small synthetic
files (`tests/fixtures/make_fixtures.sh` regenerates the genome and FASTQ fixtures).

## Cirro

Two process registrations come from this repository:

| Directory | Entry point | Input dataset |
| --- | --- | --- |
| `.cirro/align/` | FASTQ | Paired or single-end FASTQ datasets |
| `.cirro/apply_bqsr/` | alignment + table | A `sarek_align` dataset (`preprocessing/parabricks/<sample>/<sample>.{bam,bam.bai,table}`), or this pipeline's own `--apply_bqsr false` output (`preprocessing/markduplicates/` + `preprocessing/recal_table/`) |

Each directory holds `process-form.json`, `process-input.json`, `preprocess.py`,
`process-compute.config` and `process-output.json`. Both Parabricks processes run on the
on-demand GPU queue (`PW_ONDEMAND_JOB_QUEUE`) with retries on resource-related exit codes.

**FASTQ form:** genome source (iGenomes GATK.GRCh38, or a custom genome via a BWA index
dataset plus known-sites/intervals from the references library), output format, `apply_bqsr`
(hides quantization when off), quantization, trimming, single-end duplicate marking, optical
pixel distance, alignment GPUs, low-memory mode and an fq2bam memory override. A
`fastq_singleton` column in the dataset's samplesheet adds that file to the sample as a
single-end lane, resolved next to the sample's read 1.

**Registration settings:** repository `break-through-cancer/parabricks-fq2bam-bqsr`, entry
script `main.nf`, configuration directory `.cirro/align` or `.cirro/apply_bqsr`, Nextflow
`26.04.0` or later. The align process needs `file_mapping_rules` so its output dataset's
files are discoverable by a later apply_bqsr run — Cirro's engine uses .NET-style named
groups (`(?<sample>...)`, not Python's `(?P<sample>...)`), and `is_sample: true` is what
makes a rule's captured group populate per-file sample metadata:
- Aligned reads (`is_sample: true`): `preprocessing/(?<bamType>recalibrated|markduplicates)/(?<sample>[^/]+)/[^/]+\.(?:bam|cram)$` and the matching `.bam.bai|.cram.crai` pattern
- BQSR table (`is_sample: true`): `preprocessing/recal_table/(?<sample>[^/]+)/[^/]+\.table$`

## Known gaps

1. **`stageInMode 'copy'`** is inherited from the nf-core Parabricks modules;
   `--preserve-file-symlinks` may make the copy unnecessary.
2. **`SAMTOOLS_FINALIZE` exists because `applybqsr` writes BAM only** (NVIDIA documents
   its `--out-bam` as "Output BAM file"). With quantization on, `QUANTIZE_QUALS` writes
   CRAM directly and this step does not run.
