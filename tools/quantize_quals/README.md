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
docker build --platform linux/amd64 -t quantize-quals:0.1.0 tools/quantize_quals
# then push where the executor can pull it, e.g.:
#   docker tag quantize-quals:0.1.0 ghcr.io/<org>/quantize-quals:0.1.0
#   docker push ghcr.io/<org>/quantize-quals:0.1.0
```

Then point the pipeline at the pushed image with `--quantize_quals_container`.
