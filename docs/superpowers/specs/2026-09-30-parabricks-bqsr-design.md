# parabricks_bqsr: Technical Spec

Status: design/planning notes, no code written yet. Standalone Nextflow pipeline,
deliberately separate from `nf-core/sarek`'s own `sarek_align` and from `ICONICC-NF` —
neither of those repos should gain any dependency on this one, or vice versa.

## 1. Background and scope

Org convention runs alignment + BQSR via `nf-core/sarek`'s `sarek_align` Cirro process,
using the Parabricks aligner (`PARABRICKS_FQ2BAM`) for GPU acceleration. A recent change
(`CirroBio/Cirro-pipelines` PR #115, commit `47a75115`) makes `sarek_align` publish the
BQSR recalibration `.table` file that `PARABRICKS_FQ2BAM` already produces as a side
effect of alignment whenever known-sites resources (`dbsnp`/`known_indels`) are supplied
— built efficiently against a WGS calling-intervals file, matching GATK/Broad's own
best-practice pattern of building the recalibration model on a representative interval
subset for speed.

That table is not itself the end of BQSR: applying it needs to happen genome-wide,
**without** interval restriction, so reads outside the model-building intervals —
including unmapped reads — still get recalibrated quality scores. Parabricks' own
`fq2bam`-integrated application path is interval-scoped, so it cannot be reused for this
step; a separate `pbrun applybqsr` pass is required.

**External dependency, unresolved as of this writing (2026-09-30, re-checked via `gh pr
view 115 --repo CirroBio/Cirro-pipelines`): PR #115 is still open/draft, not merged.**
Real end-to-end testing of this pipeline against an actual `sarek_align` output is
blocked until it merges — until then, a `recal_table` samplesheet input can only come
from a manually-run `sarek_align` on a branch carrying that commit, or a hand-constructed
fixture. Code development and unit-level testing don't need to wait on it.

Additionally, NVIDIA Parabricks currently has no equivalent of GATK `ApplyBQSR`'s
`--static-quantized-quals` (or its dynamic `--quantize-quals`) — GATK's own quality-score
quantization step doesn't exist in Parabricks at all. Replicating it requires a
purpose-built post-processing tool.

**This pipeline's job, precisely:**
1. Given a sample's already-aligned BAM/CRAM + its recalibration table (both from
   `sarek_align`), run `pbrun applybqsr` with no interval restriction, producing a fully
   recalibrated BAM/CRAM.
2. Optionally (toggleable), run a custom quality-score quantization tool over that output,
   replicating GATK `ApplyBQSR`'s static quantization behavior exactly.

Explicitly out of scope: alignment itself (stays in `sarek_align`), recalibration-table
*generation* (already happens inside `PARABRICKS_FQ2BAM`, published by the org's own
already-in-flight PR), and dynamic quantization (`--quantize-quals`, see §5 — deliberately
deferred).

## 2. Architecture

```
samplesheet.csv (sample, alignment, alignment_index, recal_table)
        │
        ▼
  PARABRICKS_APPLYBQSR   -- hand-written local module (no nf-core/sarek module exists
        │                  for this yet -- vendored PARABRICKS_FQ2BAM/HAPLOTYPECALLER
        │                  are the only Parabricks modules nf-core/modules currently
        │                  ships). No --interval-file passed.
        ▼
  [ recalibrated BAM/CRAM, per sample ]
        │
        ▼ (if params.quantize_quals_enabled)
  QUANTIZE_QUALS          -- Nextflow wrapper around a new, purpose-built C/HTSlib
        │                  binary (see §4). Container built and versioned fresh,
        │                  same convention ICONICC-NF used for its R container:
        │                  no reverse-pinning an existing unversioned image.
        ▼
  final BAM/CRAM (params.output_fmt: bam|cram, mirroring sarek's own toggle)
```

Two new local Nextflow modules (`PARABRICKS_APPLYBQSR`, `QUANTIZE_QUALS`), one thin
samplesheet-parsing subworkflow (no pairing logic — this pipeline has no tumor/normal
concept, unlike ICONICC-NF), and one new tool + container (the quantization binary).

**Why a separate post-processing step, not folded into a single custom BQSR+quantize
tool:** considered and rejected. Reusing Parabricks' own tested `applybqsr` for
recalibration and keeping quantization as an independent, focused tool is lower-risk than
re-deriving GATK/Parabricks' recalibration math ourselves — and it keeps the quantization
tool small enough to validate in isolation (feed it any BAM/CRAM with known qualities,
check the output against the documented lookup table, no BQSR logic involved at all). The
cost is one extra full read/write pass over a WGS-scale file, mitigated by multi-threaded
BGZF I/O (§4).

## 3. Samplesheet schema

CSV, one row per sample — no pairing, no patient/tumor/normal concept:

```
sample,alignment,alignment_index,recal_table
SAMPLE_A,SAMPLE_A.cram,SAMPLE_A.cram.crai,SAMPLE_A.table
```

`alignment`/`alignment_index` rather than `cram`/`crai` deliberately: a row may hold
either BAM or CRAM input, since Parabricks' `--in-bam` flag accepts CRAM directly (already
verified during ICONICC-NF's Phase 2 work, same NVIDIA tool family). `params.output_fmt`
(`bam`|`cram`) controls the *output* format independent of what format came in, mirroring
sarek's own toggle — not derived from the input's extension.

## 4. The quantization tool (v1: static mode only)

Full behavior spec supplied directly by the dev (`bqsr_flag_specs.txt`, read and folded
in verbatim below) — this is not something to re-derive or approximate, it is already
GATK-verified.

### 4.1 CLI (v1 scope)

```
quantize_quals --in <bam|cram> --out <bam|cram> --ref <fasta, required for CRAM> \
               --static-quantized-quals <Q...>    (repeatable, e.g. 10 20 30)
               --preserve-qscores-less-than <Q>    (default: 6)
               --round-down-quantized               (boolean flag, default: false)
               --threads <N>                        (BGZF compression threads, ~8-16)
```

**Deliberately not implemented in v1** (see §5): `--quantize-quals`, `--bqsr-recal-file`,
and the mutual-exclusion validation between static and dynamic modes. The dev confirmed:
build only the static-mode flags now: add dynamic mode later, as a genuinely separate
follow-on effort — no scaffolding for it needs to exist in v1's argument parser.

### 4.2 Algorithm

Build a lookup table once at program startup (`mapping[0..MAX_QUAL]`), not per-base
floating-point math:

- For `q < preserve_threshold`: `mapping[q] = q` (unchanged).
- Otherwise, boundaries = `[preserve_threshold] + sorted(static_quantized_quals)`.
- **`--round-down-quantized` true**: `mapping[q]` = the largest boundary `<= q`.
- **`--round-down-quantized` false (default)**: `mapping[q]` = the boundary minimizing
  `abs(Pcorrect(q) - Pcorrect(boundary))`, where `Pcorrect(Q) = 1 - 10^(-Q/10)` (nearest
  boundary in probability space, not in raw Phred-score space).
- Values above the largest static bin map to the largest bin.

Worked example, `--preserve-qscores-less-than 6 --static-quantized-quals 10 20 30`
(boundaries: `6, 10, 20, 30`) — the two modes have genuinely different bucket edges, not
just different target values, so they're quoted separately rather than merged into one
table (a merged grid was tried while writing this section and produced wrong entries):

**Nearest-in-probability-space (default, `--round-down-quantized` not passed):**
`0-5` unchanged, `6-7` → `6`, `8-12` → `10`, `13-22` → `20`, `≥23` → `30`.

**Round-down (`--round-down-quantized` passed):**
`0-5` unchanged, `6-9` → `6`, `10-19` → `10`, `20-29` → `20`, `≥30` → `30`.

### 4.3 Record processing

Per alignment record: replace only the per-base `QUAL` array via the lookup table.
Preserve everything else exactly — `QNAME`, `FLAG`, `RNAME`, `POS`, `MAPQ`, `CIGAR`,
`RNEXT`, `PNEXT`, `TLEN`, `SEQ`, all auxiliary tags. Output stays coordinate-sorted if
input is. No BQSR recalculation, no alignment modification of any kind.

### 4.4 Performance

HTSlib (or equivalent), multithreaded BGZF decompression/compression (`hts_set_threads()`
or equivalent), ~8-16 threads. The transform itself is a trivial array lookup per base;
runtime should be dominated by I/O and (de)compression, not computation. No GPU needed.

### 4.5 Validation

Against real GATK, not synthetic-only: same original pre-BQSR BAM + same recalibration
table, through both paths —

- **Reference**: `GATK ApplyBQSR --static-quantized-quals 10 20 30 --preserve-qscores-less-than 6`
- **Test**: `Parabricks ApplyBQSR` → `quantize_quals --static-quantized-quals 10 20 30 --preserve-qscores-less-than 6`

Compare `QNAME`/`SEQ`/`QUAL` record-by-record; the only expected divergence from the
un-quantized Parabricks output is the `QUAL` quantization itself. This needs a real
recalibration table and a real GATK install to run — an accepted gap until that fixture
data exists, same posture ICONICC-NF already uses for its own CNV fixtures.

## 5. Explicitly deferred: dynamic quantization (`--quantize-quals`)

Full behavior already spec'd by the dev (`future-feature_quantize-quals.txt`) for when
this is picked up later — recorded here so the boundary is explicit, not implemented now:

- `--quantize-quals N`: `N == 0` disables, `N > 0` recomputes an N-level quantization
  mapping from the BQSR report's own quality-score distribution (GATK's
  `QuantizationInfo`/`QualQuantizer` bin-merging algorithm, not an approximation), `N < 0`
  reuses the mapping already encoded in the recalibration report.
- Mutually exclusive with `--static-quantized-quals`/`--round-down-quantized` — the tool
  must reject combining them once dynamic mode exists.
- Requires a second input (`--bqsr-recal-file`) beyond what static mode needs.
- Both modes share `--preserve-qscores-less-than`.

v1 builds none of this — no argument-parser scaffolding, no mode-dispatch groundwork.
When this is picked up, it is its own brainstorm/plan cycle, not a v1 afterthought.

## 6. Error handling

- Samplesheet: missing `alignment`/`alignment_index`/`recal_table` on any row, or a path
  that doesn't exist, fails fast at startup naming the offending row/path — matching
  ICONICC-NF's established convention, not failing deep inside a container.
- Quantization tool argument validation (v1 scope only): reject an empty
  `--static-quantized-quals` list; reject a `--preserve-qscores-less-than` value that
  isn't strictly below the smallest supplied static bin.

## 7. Testing strategy

Two distinct layers, not one:

- **Nextflow layer**: nf-test per module/subworkflow, following the org's established
  "known accepted gap" convention (placeholder fixtures where real data isn't available,
  documented rather than faked, mirroring ICONICC-NF throughout).
- **Quantization tool layer** (separate from nf-test): a C-level test harness —
  synthetic BAMs with known input qualities, verifying the lookup table produces the
  exact documented mappings (§4.2, both modes) — plus §4.5's real-GATK-diff validation,
  which is the accepted gap noted there until real fixture data (a pre-BQSR BAM + recal
  table + GATK install) is available.

## 8. Deployment target

Nextflow pipeline now; Cirro wiring (`.cirro/` process-form.json etc.) added later, once
the core pipeline works, living in this same repo (not split across a separate
Cirro-pipelines-style repo) — confirmed with the dev.

## 9. Repo layout (informs the implementation plan, not a task breakdown itself)

```
main.nf
nextflow.config
assets/samplesheet.csv
subworkflows/local/samplesheet_to_samples/main.nf   -- thin CSV parsing, no pairing
modules/local/parabricks_applybqsr.nf               -- hand-written, no vendor exists
modules/local/quantize_quals.nf                     -- Nextflow wrapper around the tool
tools/quantize_quals/                               -- C/HTSlib source tree
  src/                                              -- implementation
  tests/                                            -- C-level unit tests
  Makefile
  Dockerfile                                        -- fresh, versioned container
tests/                                              -- nf-test suite
SPEC.md
README.md
```

## Rulings / decisions log

- **Deployment**: Cirro process eventually; `.cirro/` files live in this same repo, added
  once the core pipeline works (not upfront).
- **Input contract**: samplesheet, one row per sample — matches sarek/ICONICC-NF
  convention over flat CLI params or Cirro-dataset-driven input.
- **Output format**: configurable (`params.output_fmt`, bam|cram), matching sarek's own
  toggle rather than hardcoding either format.
- **Quantization scope**: in v1, toggleable — but only the *static* mode (§4); dynamic
  mode (§5) is fully deferred, no scaffolding built now.
- **Quantization tool ownership**: built as part of this effort (not wrapping an
  externally-supplied binary).
- **Quantization architecture**: a separate post-processing step after Parabricks'
  `applybqsr` (Approach A), not folded into one custom recalibration+quantization tool —
  reuses Parabricks' tested recalibration logic, keeps the quantization tool
  independently testable, at the cost of one extra I/O pass.
