# parabricks_bqsr

Standalone Nextflow pipeline that applies a `sarek_align` BQSR recalibration table
genome-wide with NVIDIA Parabricks `applybqsr`. It can also replicate GATK `ApplyBQSR`'s
static quality-score quantization (`--static-quantized-quals`), a step Parabricks does not
provide. See `SPEC.md` for the full design.

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
nextflow run main.nf --input samplesheet.csv --quantize_quals_enabled true --output_fmt cram
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

Relative paths resolve against the pipeline directory (`projectDir`). Absolute paths,
local or `s3://`, are recommended for real runs.

### Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--outdir` | `./results` | Output directory. |
| `--ref_fasta` / `--ref_fasta_fai` | hg38 from `s3://broad-references` | Reference used for alignment. |
| `--output_fmt` | `cram` | `bam` or `cram`; independent of input format. |
| `--quantize_quals_enabled` | `false` | Run `QUANTIZE_QUALS` after `applybqsr`. |
| `--static_quantized_quals` | `10,20,30` | Comma-separated static bins. |
| `--preserve_qscores_less_than` | `6` | Qualities below this value remain unchanged. |
| `--round_down_quantized` | `false` | Round down to a bin instead of the nearest bin in probability space. |
| `--quantize_quals_container` | `quantize-quals:0.1.0` | Quantizer image; see `tools/quantize_quals/README.md`. |

### Outputs

`<outdir>/quantize_quals/` (quantization on) or `<outdir>/samtools_finalize/`
(quantization off) holds `<sample>.recal.<bam|cram>` plus its index. The intermediate
`applybqsr` BAM is not published.

## Testing

```bash
make -C tools/quantize_quals test   # C unit tests + CLI end-to-end tests (needs HTSlib, samtools)
nf-test test                        # module, subworkflow and pipeline tests (needs Docker)
```

Before running nf-test, build the quantizer image once:

```bash
docker build --platform linux/amd64 -t quantize-quals:0.1.0 tools/quantize_quals
```

Output content is decoded inside nf-test with the
[nft-bam](https://github.com/nvnieuwk/nft-bam) plugin (`nf-test.config`); nf-test
downloads it on first run, so no host `samtools` is needed. All fixtures are tiny
synthetic files, not patient data.

## Known gaps

1. **`CirroBio/Cirro-pipelines` PR #115 is not merged.** Until it merges, `sarek_align`
   does not publish the recalibration table, so a real `recal_table` input requires a
   manual `sarek_align` run on a branch carrying commit `47a75115`.
2. **No real Parabricks execution has been validated.** The dev machine has no NVIDIA
   GPU, so `PARABRICKS_APPLYBQSR` is tested only through `-stub`. Real validation must
   happen on Cirro against a GPU instance. Items to confirm there: `stageInMode 'copy'`
   (inherited from the nf-core Parabricks modules; `--preserve-file-symlinks` may remove
   the copy) and CPU/memory sizing.
3. **GATK parity (`SPEC.md` §4.5) is not yet checked.** This check requires a real pre-BQSR
   BAM, its recalibration table and a GATK install. The quantizer is verified only
   against the spec's documented mapping tables.
4. **`SAMTOOLS_FINALIZE` is an addition to `SPEC.md` §2.** The pipeline assumes
   `applybqsr` writes BAM only: NVIDIA documents its `--out-bam` as "Output BAM file",
   while `fq2bam` documents "Path of a BAM/CRAM file". With quantization off, this step
   indexes the BAM or converts it to indexed CRAM. With quantization on,
   `QUANTIZE_QUALS` writes CRAM directly and this step does not run.
5. **Cirro wiring (`.cirro/`) is not added yet** (`SPEC.md` §8).
6. **The quantizer image is local only.** Push `quantize-quals:0.1.0` to a registry
   and set `--quantize_quals_container` before any remote run.
