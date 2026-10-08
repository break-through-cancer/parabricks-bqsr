#!/usr/bin/env bash
# Runs the nf-test suite across N parallel shards (default 4) and reports a
# combined pass/fail. Extra arguments are passed through to each shard's
# nf-test invocation, e.g.:
#   scripts/test-parallel.sh 4 --tag quantize
#   scripts/test-parallel.sh --changedSince HEAD^   # N defaults to 4
set -euo pipefail

N=4
if [[ "${1:-}" =~ ^[0-9]+$ ]]; then
    N="$1"
    shift
fi

log_dir=$(mktemp -d)
trap 'rm -rf "$log_dir"' EXIT

pids=()
for i in $(seq 1 "$N"); do
    nf-test test --shard "${i}/${N}" "$@" > "${log_dir}/shard-${i}.log" 2>&1 &
    pids+=("$!")
done

status=0
for pid in "${pids[@]}"; do
    wait "$pid" || status=1
done

for i in $(seq 1 "$N"); do
    echo "===== shard ${i}/${N} ====="
    cat "${log_dir}/shard-${i}.log"
done

exit "$status"
