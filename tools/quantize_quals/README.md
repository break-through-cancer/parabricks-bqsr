# tools/quantize_quals

C/HTSlib tool that replicates GATK `ApplyBQSR --static-quantized-quals` on an already
recalibrated BAM/CRAM (Parabricks `applybqsr` has no quantization step). Only per-base
`QUAL` values change.

```bash
quantize_quals --in in.bam --out out.cram --ref ref.fa \
    --static-quantized-quals 10 20 30 --preserve-qscores-less-than 6 --threads 8
```

Coordinate-sorted BAM/CRAM output is indexed alongside (`out.bam.bai` / `out.cram.crai`).
Records without qualities (`QUAL = *`) pass through unchanged. An `@PG` line is added.

## Algorithm

A 256-entry lookup table is built once at startup; each base quality is then replaced by
`mapping[q]`.

- `q < --preserve-qscores-less-than`: unchanged.
- Otherwise the boundaries are the threshold plus the sorted, de-duplicated
  `--static-quantized-quals` values.
- Default: the boundary nearest to `q` in probability space, comparing
  `Pcorrect(Q) = 1 - 10^(-Q/10)`, not raw Phred values.
- `--round-down-quantized`: the largest boundary `<= q`.
- Values above the largest bin map to the largest bin. `QUAL = *` is left alone.

With `--preserve-qscores-less-than 6 --static-quantized-quals 10 20 30`
(boundaries 6, 10, 20, 30):

| Input Q | Default (nearest) | `--round-down-quantized` |
| --- | --- | --- |
| 0-5 | unchanged | unchanged |
| 6-7 | 6 | 6 |
| 8-9 | 10 | 6 |
| 10-12 | 10 | 10 |
| 13-19 | 20 | 10 |
| 20-22 | 20 | 20 |
| 23-29 | 30 | 20 |
| 30+ | 30 | 30 |

## Log output

Everything is written to stderr (stdout stays empty), so Nextflow captures it in the
task's `.command.log`:

```
quantize_quals 0.1.2 (htslib 1.22.1)
quantize_quals: input: input/S.recal.bam
quantize_quals: mode: nearest bin in probability space; preserve Q<6; bins 10,20,30; threads 8
quantize_quals: mapping: 6-7->6 8-12->10 13-22->20 23+->30
quantize_quals: output: S.recal.cram (CRAM 3.0, indexed, reference Homo_sapiens_assembly38.fasta)
quantize_quals: progress: 10000000 records, at chr1:23905377, 41.2 s, 242718 records/s
...
quantize_quals: done: 812345678 records (1234 without qualities), 121851851700 qualities, 98765432100 changed (81.1%) in 3402.6 s
quantize_quals: input qualities: Q2:0.1% ... Q34:20.3%
quantize_quals: output qualities: Q2:0.1% ... Q30:89.4%
```

The values above are illustrative. `mapping` shows the exact lookup table applied.
`--cram-version 3.0|3.1` selects the CRAM version (default 3.0, the most widely readable).
`--progress-every N` sets the progress interval in records (default 10,000,000; `0`
disables it).

## Build and test locally

Requires HTSlib (found through `pkg-config`) and `samtools` for the CLI tests.

```bash
make test
```

`tests/test_quant.c` checks the lookup table against both documented mappings;
`tests/test_cli.sh` runs the binary on `tests/fixtures/input.sam` (synthetic, one 500 bp
contig) and diffs every record against the expected outputs, which were generated
independently from the spec tables.

## Container

The build stage runs `make test`, so a failing test fails the image build. Always tag
with the version in `src/main.c` (`QQ_VERSION`), never `latest`:

```bash
docker build --platform linux/amd64 -t ghcr.io/break-through-cancer/parabricks-bqsr:<version> tools/quantize_quals
docker push ghcr.io/break-through-cancer/parabricks-bqsr:<version>
```

Then update `params.quantize_quals_container` in `nextflow.config` to the new tag.
