#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 9 ]]; then
    echo "Usage: $0 <dataset> <module> <pfile_prefix> <out_prefix> <log_file> <out_bed> <out_bim> <out_fam> <done_file>" >&2
    exit 2
fi

dataset="$1"
software_module="$2"
pfile_prefix="$3"
out_prefix="$4"
log_file="$5"
out_bed="$6"
out_bim="$7"
out_fam="$8"
done_file="$9"

mkdir -p "$(dirname "$log_file")"
mkdir -p "$(dirname "$out_prefix")"

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$log_file"
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

ensure_command() {
    if ! command -v "$1" >/dev/null 2>&1; then
        log "ERROR: required command '$1' not found"
        exit 1
    fi
}

require_file() {
    if [[ ! -s "$1" ]]; then
        log "ERROR: required file missing or empty: $1"
        exit 1
    fi
}

log "=== PLINK2 to PLINK1 conversion ==="
log "Dataset: $dataset"
log "PFILE prefix: $pfile_prefix"
log "Output prefix: $out_prefix"

require_file "${pfile_prefix}.pgen"
require_file "${pfile_prefix}.pvar"
require_file "${pfile_prefix}.psam"

load_module
ensure_command plink2

log "Running plink2 --make-bed"
plink2 \
    --pfile "$pfile_prefix" \
    --max-alleles 2 \
    --make-bed \
    --memory 16000 \
    --out "$out_prefix" \
    >> "$log_file" 2>&1

log "PLINK2 conversion complete"

if [[ ! -s "${out_prefix}.bed" || ! -s "${out_prefix}.bim" || ! -s "${out_prefix}.fam" ]]; then
    log "ERROR: missing PLINK1 outputs after conversion"
    exit 1
fi

if [[ "${out_prefix}.bed" != "$out_bed" ]]; then
    cp "${out_prefix}.bed" "$out_bed"
fi
if [[ "${out_prefix}.bim" != "$out_bim" ]]; then
    cp "${out_prefix}.bim" "$out_bim"
fi
if [[ "${out_prefix}.fam" != "$out_fam" ]]; then
    cp "${out_prefix}.fam" "$out_fam"
fi

log "Outputs ready: $out_bed, $out_bim, $out_fam"

touch "$done_file"
log "Done marker created at $done_file"
