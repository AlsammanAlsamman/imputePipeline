#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <dataset> [snakemake-args...]" >&2
    exit 2
fi

dataset="$1"
shift || true

snakemake \
    --snakefile rules/merge_plink2.smk \
    "results/10_plink_merged/${dataset}_merged.done" \
    "$@"
