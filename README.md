# parabricks_bqsr

Standalone Nextflow pipeline that applies a `sarek_align` BQSR recalibration table
genome-wide with NVIDIA Parabricks `applybqsr`. It can also replicate GATK `ApplyBQSR`'s
static quality-score quantization (`--static-quantized-quals`), a step Parabricks does not
provide.

```
samplesheet ─▶ SAMPLESHEET_TO_SAMPLES ─▶ PARABRICKS_APPLYBQSR (no intervals, BAM)
                                              │
                   quantize_quals_enabled ────┤
                                  true        ▼        false
                          QUANTIZE_QUALS           SAMTOOLS_FINALIZE
                                   └──── <sample>.recal.{bam,cram} + index ────┘
```

## Local usage

```bash
nextflow run main.nf --input samplesheet.csv
nextflow run main.nf --input samplesheet.csv --quantize_quals_enabled false --output_fmt bam
```

`--input` is the only required parameter. Invalid parameter values, missing samplesheet
values and nonexistent paths all fail at startup with a message naming the problem.

### Samplesheet format

CSV, one row per sample. See `assets/samplesheet.csv`.

| Column | Meaning |
| --- | --- |
| `sample` | Sample ID; must be unique within the samplesheet. |
| `alignment` | Aligned BAM or CRAM from `sarek_align`. |
| `alignment_index` | Its `.bai` / `.crai`. |
| `recal_table` | The BQSR recalibration table `sarek_align` publishes for that sample. |

`--input` itself may be local or `s3://`; a relative `--input` resolves against the launch
directory. Paths inside the samplesheet may be `s3://` (or any URI Nextflow supports)
and are checked for existence at startup. Use absolute or `s3://` paths inside the
samplesheet for AWS Batch runs: relative row paths resolve against the pipeline
directory (`projectDir`) on the head node, which suits the committed test fixtures only.

### GPU

`PARABRICKS_APPLYBQSR` requests one GPU (`accelerator = 1` in `conf/modules.config`).
On AWS Batch this places the task on a GPU compute environment and passes
`--num-gpus 1` to `pbrun`. With the local executor and Docker, `--gpus all` is added to
the container options. `applybqsr` uses at most two GPUs; raise `accelerator` to `2` to
trade cost for speed.

### Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--outdir` | `./results` | Output directory. |
| `--ref_fasta` / `--ref_fasta_fai` | hg38 from `s3://broad-references` | Reference used for alignment. |
| `--output_fmt` | `cram` | `bam` or `cram`; independent of input format. |
| `--quantize_quals_enabled` | `true` | Run `QUANTIZE_QUALS` after `applybqsr`; `false` publishes the unquantized recalibrated output instead. |
| `--static_quantized_quals` | `10,20,30` | Comma-separated static bins. |
| `--preserve_qscores_less_than` | `6` | Qualities below this value remain unchanged. |
| `--round_down_quantized` | `false` | Round down to a bin instead of the nearest bin in probability space. |
| `--quantize_quals_container` | `ghcr.io/break-through-cancer/parabricks-bqsr:0.1.1` | Quantizer image; see `tools/quantize_quals/README.md`. |

### Outputs

`<outdir>/preprocessing/recalibrated/<sample>/<sample>.recal.<bam|cram>` plus its index,
the same layout nf-core/sarek uses for recalibrated alignments, so downstream sarek-based
variant calling (including Cirro's `sarek_call_variants`) recognises the output. With
quantization on, this file is the quantized one. The intermediate `applybqsr` BAM is not
published.

## Testing

```bash
make -C tools/quantize_quals test   # C unit tests + CLI end-to-end tests (needs HTSlib, samtools)
nf-test test                        # module, subworkflow and pipeline tests (needs Docker)
```

nf-test pulls the quantizer image from GHCR. To test local changes to the tool,
build and tag it under the same name first:

```bash
docker build --platform linux/amd64 -t ghcr.io/break-through-cancer/parabricks-bqsr:0.1.1 tools/quantize_quals
```

Output content is decoded inside nf-test with the
[nft-bam](https://github.com/nvnieuwk/nft-bam) plugin (`nf-test.config`); nf-test
downloads it on first run, so no host `samtools` is needed. All fixtures are tiny
synthetic files, not patient data.

## Cirro

`.cirro/` holds the Cirro custom-pipeline configuration:

| File | Purpose |
| --- | --- |
| `process-form.json` | Run form: output format and quantization settings. |
| `process-input.json` | Maps form values to pipeline parameters; sets the iGenomes GATK.GRCh38 reference from Cirro's references bucket (the reference `sarek_align` uses). |
| `preprocess.py` | Builds the samplesheet from the input dataset's `preprocessing/parabricks/<sample>/` files: the pre-BQSR fq2bam alignment, its index and its `.table`. Other stages are never used. |
| `process-compute.config` | AWS Batch overrides: `applybqsr` gets 1 GPU on the on-demand queue (`PW_ONDEMAND_JOB_QUEUE`), as in `sarek_align`; retries on resource-related exit codes. |
| `process-output.json` | No post-processing commands. |

Registration settings for the custom pipeline:

- **Repository:** `break-through-cancer/parabricks-bqsr`, entry script `main.nf`, configuration directory `.cirro`.
- **Nextflow version:** `25.10.4` (stub-run verified).
- **Input dataset:** a `sarek_align` run with the Parabricks aligner and known sites supplied, `save_mapped` on and `baserecalibrator` skipped. That combination publishes `preprocessing/parabricks/<sample>/<sample>.{bam,bam.bai,table}`.
- **Output file mapping:** the same patterns `sarek_align` uses for recalibrated alignments:
  - `preprocessing/(?P<bamType>recalibrated)/(?P<sampleName>[^/]+)/[^/]+\.(?:bam|cram)$`
  - `preprocessing/(?P<bamType>recalibrated)/(?P<sampleName>[^/]+)/[^/]+\.(?:bam\.bai|cram\.crai)$`

`preprocess.py` tests run outside Cirro:

```bash
python -m pytest .cirro
```

## Known gaps

1. **`CirroBio/Cirro-pipelines` PR #115 is not merged.** Until it merges, `sarek_align`
   does not publish the recalibration table, so a real `recal_table` input requires a
   manual `sarek_align` run on a branch carrying commit `47a75115`.
2. **No real Parabricks execution has been validated.** The dev machine has no NVIDIA
   GPU, so `PARABRICKS_APPLYBQSR` is tested only through `-stub`. Real validation must
   happen on Cirro against a GPU instance. Items to confirm there: the `accelerator`
   request co-existing with Cirro's own GPU compute config, `stageInMode 'copy'`
   (inherited from the nf-core Parabricks modules; `--preserve-file-symlinks` may remove
   the copy) and CPU/memory sizing.
3. **GATK parity is not yet checked.** The reference path is GATK `ApplyBQSR
   --static-quantized-quals 10 20 30 --preserve-qscores-less-than 6`; the test path is
   Parabricks `applybqsr` followed by `quantize_quals` with the same settings, on the same
   pre-BQSR BAM and table. Records should match exactly. This needs a GATK install and
   real data; the quantizer is so far verified against its documented mapping tables
   only.
4. **`SAMTOOLS_FINALIZE` exists because of a Parabricks limitation.** The pipeline assumes
   `applybqsr` writes BAM only: NVIDIA documents its `--out-bam` as "Output BAM file",
   while `fq2bam` documents "Path of a BAM/CRAM file". With quantization off, this step
   indexes the BAM or converts it to indexed CRAM. With quantization on,
   `QUANTIZE_QUALS` writes CRAM directly and this step does not run.
5. **The Cirro configuration has not run on Cirro yet.** `preprocess.py` is tested
   locally and produces the expected samplesheet from a real `sarek_align` dataset
   listing; the form, input mapping and compute config are untested on the platform.
