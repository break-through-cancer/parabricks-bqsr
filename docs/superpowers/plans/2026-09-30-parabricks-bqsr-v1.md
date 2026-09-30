# parabricks_bqsr v1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Standalone Nextflow pipeline that applies a `sarek_align` recalibration table genome-wide with `pbrun applybqsr`, then optionally runs a purpose-built C/HTSlib static quality-score quantizer, publishing an indexed BAM or CRAM per sample.

**Architecture:** `SAMPLESHEET_TO_SAMPLES` (thin CSV validation) → `PARABRICKS_APPLYBQSR` (no intervals, always writes BAM) → either `QUANTIZE_QUALS` (writes final `bam|cram` + index) or `SAMTOOLS_FINALIZE` (indexes BAM, or converts to indexed CRAM). Scaffolding follows ICONICC-NF: flat `cpus`/`memory` per process, digest-pinned containers, `publishDir` default in `nextflow.config`, `projectDir`-relative samplesheet paths, nf-test under `tests/local`, `tests/subworkflows`, `tests/main.nf.test`, committed tiny synthetic fixtures.

**Tech Stack:** Nextflow 26 (DSL2), nf-test 0.9.5, NVIDIA Parabricks 4.7.1-1, HTSlib 1.22.1 (C99), samtools 1.24, Docker.

**Spec:** `SPEC.md` (moved from `docs/superpowers/specs/2026-09-30-parabricks-bqsr-design.md`, matching ICONICC-NF's root-level `SPEC.md`).

## Global Constraints

- No dependency on `nf-core/sarek` or `ICONICC-NF` in either direction.
- `pbrun applybqsr` receives no `--interval-file` / `-L`.
- Static quantization only; no `--quantize-quals`, `--bqsr-recal-file`, or dynamic-mode scaffolding.
- `--preserve-qscores-less-than` default `6`; `--round-down-quantized` default false.
- Nearest-boundary mode compares `Pcorrect(Q) = 1 - 10^(-Q/10)`, not raw Phred.
- Only `QUAL` changes; all other record fields and aux tags preserved; `@PG` header line added.
- Quantizer container built fresh and versioned (`quantize-quals:0.1.0`), never `:latest`.
- `params.output_fmt` ∈ {`bam`, `cram`}, independent of input format.
- Samplesheet errors fail at startup naming the row/path.
- Minimal inline comments (user preference); rationale goes in commit messages.
- No GPU on the dev machine: Parabricks is tested via `-stub` only.

## Deviation from spec §2 (flagged)

Parabricks documents `applybqsr --out-bam` as a BAM output only. CRAM output from
`applybqsr` is unverified, so `PARABRICKS_APPLYBQSR` always writes BAM, and a third
small module, `SAMTOOLS_FINALIZE`, produces the final indexed `bam|cram` whenever
quantization is disabled. With quantization enabled, `QUANTIZE_QUALS` writes the final
format directly and no extra pass occurs.

## Review Focus

- Records with no qualities (`QUAL = *`, stored as `0xff`) must pass through untouched, not be mapped to the top bin.
- Input BAM staged under the same name as the output (for example `S.recal.bam`) must not collide: inputs are staged into `input/`.
- Duplicate `sample` values in the samplesheet must fail at startup, not overwrite each other in `results/`.
- Unsorted or duplicate `--static-quantized-quals` values (`30 10 20 20`) must behave like the sorted, de-duplicated set.
- CRAM input or output without `--ref` must fail with a clear message, not a deep HTSlib error.

---

### Task 1: Scaffolding and `SAMPLESHEET_TO_SAMPLES`

**Files:**
- Create: `main.nf` (param guard + subworkflow call), `nextflow.config`, `assets/samplesheet.csv`, `subworkflows/local/samplesheet_to_samples/main.nf`, `tests/subworkflows/samplesheet_to_samples.nf.test`, `tests/fixtures/samplesheet_to_samples/*.csv` + placeholder files
- Modify: `.gitignore`; `git mv` the spec to `SPEC.md`

**Interfaces:**
- Produces: `SAMPLESHEET_TO_SAMPLES(samplesheet)` → `emit: samples // [ meta(sample), alignment, alignment_index, recal_table ]`

- [ ] Write tests: valid two-row sheet emits two tuples with `meta.sample`; missing `recal_table` value errors naming the sample; non-existent path errors; duplicate sample errors ("Duplicate sample"); header-only sheet errors ("zero samples").
- [ ] Run `nf-test test tests/subworkflows/samplesheet_to_samples.nf.test` → FAIL (script missing).
- [ ] Implement `validateSamplesheetRow` and `resolveRelativeToProjectDir` as in ICONICC-NF's `SAMPLESHEET_TO_PAIRS`; duplicates via `.map{it[0].sample}.toList()` check; `.ifEmpty { error }`.
- [ ] Re-run → PASS. Commit.

### Task 2: `quantize_quals` C tool

**Files:**
- Create: `tools/quantize_quals/src/{quant.h,quant.c,main.c}`, `tools/quantize_quals/tests/{test_quant.c,test_cli.sh,fixtures/*}`, `tools/quantize_quals/Makefile`, `tools/quantize_quals/Dockerfile`, `tools/quantize_quals/README.md`

**Interfaces:**
- `int qq_build_mapping(int preserve, const int *bins, size_t n_bins, int round_down, uint8_t mapping[256], char *err, size_t err_len)` → 0 on success, -1 with message on invalid input.
- CLI: `quantize_quals --in <sam|bam|cram> --out <bam|cram|sam> [--ref fasta] --static-quantized-quals Q [Q...] [--static-quantized-quals Q] [--preserve-qscores-less-than 6] [--round-down-quantized] [--threads N]`. Writes `<out>.bai` / `<out>.crai` for BAM/CRAM output.

- [ ] Write `test_quant.c`: default mode expected table `0-5 id, 6-7→6, 8-12→10, 13-22→20, ≥23→30`; round-down `6-9→6, 10-19→10, 20-29→20, ≥30→30`; unsorted/duplicate bins equal sorted set; `mapping[255]==255`; empty bins rejected; `preserve >= min bin` rejected; negative / >93 bin rejected.
- [ ] `make test` → FAIL (no `quant.c`).
- [ ] Implement `quant.c`; `make test` → PASS.
- [ ] Write `test_cli.sh`: synthetic SAM (mapped, unmapped, `QUAL=*`, aux tags, secondary) → BAM and CRAM; assert `samtools view` of output equals input with only column 11 changed per the table; `.bai`/`.crai` exist; `@PG` present; CLI error cases (no bins, CRAM without `--ref`, bad threshold, missing `--in`).
- [ ] Implement `main.c` (shared `hts_tpool`, `sam_idx_init`, `sam_hdr_add_pg`); `make test` → PASS.
- [ ] Dockerfile: multi-stage `debian:bookworm-slim`, HTSlib 1.22.1 from source, runs `make test` in build stage; runtime stage carries binary, libhts, `procps`. `docker build -t quantize-quals:0.1.0 tools/quantize_quals`. Commit.

### Task 3: `QUANTIZE_QUALS` module

**Files:** Create `modules/local/quantize_quals.nf`, `tests/local/quantize_quals.nf.test`, `tests/fixtures/modules/{reference.fasta,.fai,sample.bam,.bai}`

**Interfaces:**
- Input: `tuple val(meta), path(bam, stageAs: 'input/*')`, `tuple path(ref_fasta), path(ref_fasta_fai)`, `val(output_fmt)`
- Output: `tuple val(meta), path("${meta.sample}.recal.${output_fmt}"), path("${meta.sample}.recal.${output_fmt}.{bai,crai}")`
- Args from `task.ext.args` (default set in `nextflow.config` from params).

- [ ] Tests (real run, local `quantize-quals:0.1.0`): BAM output has quantized quals (`samtools view` in test via `path(...).text` is not possible for BAM, so assert on file names + non-empty and run a CRAM case); stub test for shape. Implement, PASS, commit.

### Task 4: `PARABRICKS_APPLYBQSR` and `SAMTOOLS_FINALIZE`

**Files:** Create `modules/local/parabricks_applybqsr.nf`, `modules/local/samtools_finalize.nf`, tests `tests/local/parabricks_applybqsr.nf.test` (+ `.nf.config` sizing override), `tests/local/samtools_finalize.nf.test`

**Interfaces:**
- `PARABRICKS_APPLYBQSR`: input `tuple val(meta), path(alignment, stageAs: 'input/*'), path(alignment_index, stageAs: 'input/*'), path(recal_table)`, `tuple path(ref_fasta), path(ref_fasta_fai)`; output `tuple val(meta), path("${meta.sample}.recal.bam")`.
- `SAMTOOLS_FINALIZE`: input `tuple val(meta), path(bam)`, `tuple path(ref_fasta), path(ref_fasta_fai)`, `val(output_fmt)`; output same shape as `QUANTIZE_QUALS`.

- [ ] Tests: applybqsr stub emits one BAM and the script contains no `--interval-file`; finalize real run for `bam` (index only) and `cram` (converted + `.crai`). Implement, PASS, commit.

### Task 5: `main.nf` wiring

- [ ] Validate `params.output_fmt` and `params.quantize_quals_enabled` at startup.
- [ ] `tests/main.nf.test` (`-stub`): quantization on → `QUANTIZE_QUALS` succeeds once per sample, `SAMTOOLS_FINALIZE` never runs; off → the reverse; `--output_fmt sam` fails naming the param. PASS, commit.

### Task 6: README

- [ ] Usage, samplesheet table, params table, container build, known gaps (PR #115 unmerged, no GPU validation, GATK-diff validation pending fixture data). Commit.
