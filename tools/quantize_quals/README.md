# tools/quantize_quals

C/HTSlib tool that replicates GATK `ApplyBQSR --static-quantized-quals` on an already
recalibrated BAM/CRAM (Parabricks `applybqsr` has no quantization step). Only per-base
`QUAL` values change. See `SPEC.md` §4 for the algorithm and worked examples.

```bash
quantize_quals --in in.bam --out out.cram --ref ref.fa \
    --static-quantized-quals 10 20 30 --preserve-qscores-less-than 6 --threads 8
```

Coordinate-sorted BAM/CRAM output is indexed alongside (`out.bam.bai` / `out.cram.crai`).
Records without qualities (`QUAL = *`) pass through unchanged. An `@PG` line is added.

## Log output

Everything is written to stderr (stdout stays empty), so Nextflow captures it in the
task's `.command.log`:

```
quantize_quals 0.1.1 (htslib 1.22.1)
quantize_quals: input: input/S.recal.bam
quantize_quals: mode: nearest bin in probability space; preserve Q<6; bins 10,20,30; threads 8
quantize_quals: mapping: 6-7->6 8-12->10 13-22->20 23+->30
quantize_quals: output: S.recal.cram (CRAM, indexed, reference Homo_sapiens_assembly38.fasta)
quantize_quals: progress: 10000000 records, at chr1:23905377, 41.2 s, 242718 records/s
...
quantize_quals: done: 812345678 records (1234 without qualities), 121851851700 qualities, 98765432100 changed (81.1%) in 3402.6 s
quantize_quals: input qualities: Q2:0.1% ... Q34:20.3%
quantize_quals: output qualities: Q2:0.1% ... Q30:89.4%
```

The values above are illustrative. `mapping` shows the exact lookup table applied.
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
