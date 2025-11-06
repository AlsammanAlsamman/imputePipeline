#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 11 ]]; then
    echo "Usage: $0 <dataset> <chr> <input_vcf> <output_prefix> <conversion_log> <info_score> <maf> <missing_rate> <threads> <module> <done_file>" >&2
    exit 2
fi

dataset="$1"
chromosome="$2"
input_vcf="$3"
output_prefix="$4"
conversion_log="$5"
info_score="$6"
maf_threshold="$7"
missing_rate="$8"
thread_count="$9"
software_module="${10}"
done_file="${11}"

mkdir -p "$(dirname "$output_prefix")"
mkdir -p "$(dirname "$conversion_log")"

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

load_module() {
    if [[ -n "${MODULESHOME:-}" ]]; then
        # shellcheck disable=SC1090
        source "${MODULESHOME}/init/bash"
    fi

    if command -v module >/dev/null 2>&1; then
        module load "$software_module"
    else
        log "Module command not available; proceeding without loading $software_module"
    fi
}

log "=== VCF to PLINK1 Conversion ==="
log "Dataset: $dataset"
log "Chromosome: $chromosome"
log "Input VCF: $input_vcf"
log "PLINK prefix: $output_prefix"
log "Conversion log: $conversion_log"
log "QC filters - INFO: $info_score, MAF: $maf_threshold, Missing: $missing_rate"
log "==============================="

log "Preparing environment..."
load_module

if ! command -v plink >/dev/null 2>&1; then
    log "ERROR: plink command not found after loading module $software_module"
    exit 1
fi

log "Running PLINK conversion..."
plink \
    --vcf "$input_vcf" \
    --allow-extra-chr \
    --chr "$chromosome" \
    --double-id \
    --make-bed \
    --maf "$maf_threshold" \
    --geno "$missing_rate" \
    --threads "$thread_count" \
    --memory 30000 \
    --out "$output_prefix" \
    >"$conversion_log" 2>&1

log "PLINK conversion finished successfully, writing done marker"

touch "$done_file"

log "Done marker created at $done_file"
