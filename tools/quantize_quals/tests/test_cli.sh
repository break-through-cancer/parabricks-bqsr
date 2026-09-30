#!/usr/bin/env bash
set -uo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
bin="${QUANTIZE_QUALS:?set QUANTIZE_QUALS to the quantize_quals binary}"
fx="$here/fixtures"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
failures=0

pass() { echo "ok   $1"; }
fail() { echo "FAIL $1"; failures=$((failures + 1)); }

expect_body() {
    local name="$1" out="$2" expected="$3"
    shift 3
    if diff <(samtools view "$@" "$out") "$expected" > "$tmp/diff.txt"; then
        pass "$name"
    else
        fail "$name"; head -20 "$tmp/diff.txt"
    fi
}

expect_error() {
    local name="$1" pattern="$2"
    shift 2
    if "$bin" "$@" > "$tmp/stdout.txt" 2> "$tmp/stderr.txt"; then
        fail "$name (exited 0)"
    elif grep -q -- "$pattern" "$tmp/stderr.txt"; then
        pass "$name"
    else
        fail "$name (stderr did not match '$pattern')"; cat "$tmp/stderr.txt"
    fi
}

"$bin" --in "$fx/input.sam" --out "$tmp/nearest.bam" --static-quantized-quals 10 20 30 --threads 2 \
    && pass "SAM -> BAM, default mode exits 0" || fail "SAM -> BAM, default mode exits 0"
expect_body "default mode: only QUAL changes, nearest in probability space" \
    "$tmp/nearest.bam" "$fx/expected_nearest.body.sam"
[[ -s "$tmp/nearest.bam.bai" ]] && pass "BAM output is indexed" || fail "BAM output is indexed"

"$bin" --in "$fx/input.sam" --out "$tmp/rdown.bam" \
    --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
    --preserve-qscores-less-than 6 --round-down-quantized \
    && pass "repeated-flag form exits 0" || fail "repeated-flag form exits 0"
expect_body "round-down mode" "$tmp/rdown.bam" "$fx/expected_round_down.body.sam"

"$bin" --in "$fx/input.sam" --out "$tmp/messy.bam" --static-quantized-quals 30 10 20 20 \
    && expect_body "unsorted/duplicate bins equal the sorted set" "$tmp/messy.bam" "$fx/expected_nearest.body.sam" \
    || fail "unsorted/duplicate bins exits 0"

"$bin" --in "$tmp/nearest.bam" --out "$tmp/idempotent.bam" --static-quantized-quals 10 20 30 \
    && expect_body "quantizing already-quantized output is a no-op" "$tmp/idempotent.bam" "$fx/expected_nearest.body.sam" \
    || fail "idempotent run exits 0"

samtools view -b -o "$tmp/input.bam" "$fx/input.sam"
"$bin" --in "$tmp/input.bam" --out "$tmp/nearest.cram" --ref "$fx/ref.fa" --static-quantized-quals 10 20 30 \
    && pass "BAM -> CRAM exits 0" || fail "BAM -> CRAM exits 0"
[[ -s "$tmp/nearest.cram.crai" ]] && pass "CRAM output is indexed" || fail "CRAM output is indexed"
if diff <(samtools view -T "$fx/ref.fa" "$tmp/nearest.cram" | cut -f1-11 | sort) \
        <(cut -f1-11 "$fx/expected_nearest.body.sam" | sort) > "$tmp/diff.txt"; then
    pass "CRAM output: core fields and QUAL match"
else
    fail "CRAM output: core fields and QUAL match"; head -20 "$tmp/diff.txt"
fi

"$bin" --in "$tmp/nearest.cram" --out "$tmp/from_cram.bam" --ref "$fx/ref.fa" --static-quantized-quals 10 20 30 \
    && pass "CRAM -> BAM exits 0" || fail "CRAM -> BAM exits 0"

if samtools view -H "$tmp/nearest.bam" | grep -q $'^@PG\tID:quantize_quals'; then
    pass "@PG header line added"
else
    fail "@PG header line added"
fi
if diff <(samtools view -H "$tmp/nearest.bam" | grep -v '^@PG') <(grep '^@' "$fx/input.sam") > /dev/null; then
    pass "non-@PG header lines preserved"
else
    fail "non-@PG header lines preserved"
fi

grep -v '^@HD' "$fx/input.sam" > "$tmp/unsorted.sam"
"$bin" --in "$tmp/unsorted.sam" --out "$tmp/unsorted.bam" --static-quantized-quals 10 20 30 2> "$tmp/stderr.txt" \
    && [[ ! -e "$tmp/unsorted.bam.bai" ]] && grep -q "not coordinate-sorted" "$tmp/stderr.txt" \
    && pass "input without SO:coordinate is processed but not indexed, with a warning" \
    || fail "input without SO:coordinate is processed but not indexed, with a warning"

expect_error "no --static-quantized-quals"  "--static-quantized-quals" --in "$fx/input.sam" --out "$tmp/x.bam"
expect_error "flag with no values"          "--static-quantized-quals" --in "$fx/input.sam" --out "$tmp/x.bam" --static-quantized-quals --round-down-quantized
expect_error "threshold not below min bin"  "strictly below"           --in "$fx/input.sam" --out "$tmp/x.bam" --static-quantized-quals 10 20 30 --preserve-qscores-less-than 10
expect_error "non-integer bin"              "not an integer"           --in "$fx/input.sam" --out "$tmp/x.bam" --static-quantized-quals 10 2x0
expect_error "missing --in"                 "--in"                     --out "$tmp/x.bam" --static-quantized-quals 10 20 30
expect_error "missing --out"                "--out"                    --in "$fx/input.sam" --static-quantized-quals 10 20 30
expect_error "CRAM output without --ref"    "--ref"                    --in "$fx/input.sam" --out "$tmp/x.cram" --static-quantized-quals 10 20 30
expect_error "CRAM input without --ref"     "--ref"                    --in "$tmp/nearest.cram" --out "$tmp/x.bam" --static-quantized-quals 10 20 30
expect_error "unknown output extension"     "extension"                --in "$fx/input.sam" --out "$tmp/x.txt" --static-quantized-quals 10 20 30
expect_error "unknown option"               "unknown option"           --in "$fx/input.sam" --out "$tmp/x.bam" --static-quantized-quals 10 20 30 --quantize-quals 4
expect_error "nonexistent input"            "cannot open"              --in "$tmp/nope.bam" --out "$tmp/x.bam" --static-quantized-quals 10 20 30
expect_error "bad --threads"                "--threads"                --in "$fx/input.sam" --out "$tmp/x.bam" --static-quantized-quals 10 20 30 --threads 0

if [[ $failures -gt 0 ]]; then
    echo "test_cli: $failures failure(s)"
    exit 1
fi
echo "test_cli: all checks passed"
