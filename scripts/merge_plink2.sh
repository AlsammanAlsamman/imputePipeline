#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 10 ]]; then
    echo "Usage: $0 <dataset> <space_chr_list> <input_dir> <output_prefix> <log_file> <module> <out_pgen> <out_pvar> <out_psam> <done_file>" >&2
    exit 2
fi

dataset="$1"
chr_list="$2"
input_dir="$3"
output_prefix="$4"
log_file="$5"
software_module="$6"
expected_pgen="$7"
expected_pvar="$8"
expected_psam="$9"
done_file="${10}"

work_dir="$(dirname "$output_prefix")"
intermediate_dir="${work_dir}/plink2_intermediate"
merge_list="${intermediate_dir}/pmerge_list.txt"
current_prefix="${intermediate_dir}/current_merged"

mkdir -p "$work_dir"
mkdir -p "$(dirname "$log_file")"
mkdir -p "$intermediate_dir"

touch "$log_file"

read -r -a chromosomes <<< "$chr_list"

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

convert_to_pgen() {
    local chr="$1"
    local prefix="${input_dir}/chr${chr}"
    local out_prefix="${intermediate_dir}/chr${chr}"

    for ext in bed bim fam; do
        if [[ ! -s "${prefix}.${ext}" ]]; then
            log "ERROR: missing ${ext} file for chromosome ${chr}: ${prefix}.${ext}"
            exit 1
        fi
    done

    log "Converting chr${chr} to PLINK2 format"
    plink2 \
        --bfile "$prefix" \
        --make-pgen \
        --out "$out_prefix" \
        >> "$log_file" 2>&1
}

attempt_direct_merge() {
    local target_prefix="$1"

    log "Attempting direct PLINK2 merge via --pmerge-list"
    if plink2 \
        --pmerge-list "$merge_list" \
        --make-pgen \
        --multiallelics-already-joined \
        --out "$target_prefix" \
        >> "$log_file" 2>&1; then
        if [[ -s "${target_prefix}.pgen" ]]; then
            log "Direct merge succeeded"
            return 0
        fi
    fi

    log "Direct merge failed"
    return 1
}

sequential_merge() {
    local target_prefix="$1"

    cp "${intermediate_dir}/chr${chromosomes[0]}.pgen" "${current_prefix}.pgen"
    cp "${intermediate_dir}/chr${chromosomes[0]}.pvar" "${current_prefix}.pvar"
    cp "${intermediate_dir}/chr${chromosomes[0]}.psam" "${current_prefix}.psam"

    for (( idx=1; idx<${#chromosomes[@]}; idx++ )); do
        local chr="${chromosomes[idx]}"
        local source_prefix="${intermediate_dir}/chr${chr}"
        local step_prefix="${intermediate_dir}/step_${idx}"

        log "Sequential merge: adding chromosome ${chr}"
        if plink2 \
            --pfile "$current_prefix" \
            --pmerge "$source_prefix" \
            --make-pgen \
            --multiallelics-already-joined \
            --out "$step_prefix" \
            >> "$log_file" 2>&1; then
            mv "${step_prefix}.pgen" "${current_prefix}.pgen"
            mv "${step_prefix}.pvar" "${current_prefix}.pvar"
            mv "${step_prefix}.psam" "${current_prefix}.psam"
            continue
        fi

        log "Sequential merge failed for chr${chr}, retrying with allele limit"
        if plink2 \
            --pfile "$current_prefix" \
            --pmerge "$source_prefix" \
            --make-pgen \
            --multiallelics-already-joined \
            --merge-max-allele-ct 2 \
            --new-id-max-allele-ct 2 \
            --out "$step_prefix" \
            >> "$log_file" 2>&1; then
            mv "${step_prefix}.pgen" "${current_prefix}.pgen"
            mv "${step_prefix}.pvar" "${current_prefix}.pvar"
            mv "${step_prefix}.psam" "${current_prefix}.psam"
            continue
        fi

        log "ERROR: unable to merge chromosome ${chr} even after retries"
        return 1
    done

    mv "${current_prefix}.pgen" "${target_prefix}.pgen"
    mv "${current_prefix}.pvar" "${target_prefix}.pvar"
    mv "${current_prefix}.psam" "${target_prefix}.psam"
    return 0
}

cleanup_intermediate() {
    rm -f "${intermediate_dir}"/chr*.pgen \
          "${intermediate_dir}"/chr*.pvar \
          "${intermediate_dir}"/chr*.psam \
          "$merge_list" \
          "${intermediate_dir}"/step_*.pgen \
          "${intermediate_dir}"/step_*.pvar \
          "${intermediate_dir}"/step_*.psam \
          "${current_prefix}.pgen" \
          "${current_prefix}.pvar" \
          "${current_prefix}.psam"
}

trap cleanup_intermediate EXIT

log "=== PLINK2 Merge Workflow ==="
log "Dataset: $dataset"
log "Chromosomes: $chr_list"
log "Arguments received: $# (expected 10)"
log "Working directory: $work_dir"
log "Intermediate directory: $intermediate_dir"

load_module
ensure_command plink2

log "Converting all chromosomes to PLINK2 intermediate format"
> "$merge_list"
for chr in "${chromosomes[@]}"; do
    convert_to_pgen "$chr"
    printf '%s\n' "${intermediate_dir}/chr${chr}" >> "$merge_list"
done

temp_prefix="${output_prefix}.tmp"
rm -f "${temp_prefix}.pgen" "${temp_prefix}.pvar" "${temp_prefix}.psam"

if ! attempt_direct_merge "$temp_prefix"; then
    log "Falling back to sequential merge"
    if ! sequential_merge "$temp_prefix"; then
        log "ERROR: all merge strategies failed"
        cleanup_intermediate
        exit 1
    fi
fi

log "Preparing final outputs"
mv "${temp_prefix}.pgen" "${output_prefix}.pgen"
mv "${temp_prefix}.pvar" "${output_prefix}.pvar"
mv "${temp_prefix}.psam" "${output_prefix}.psam"

if [[ "${output_prefix}.pgen" != "$expected_pgen" ]]; then
    cp "${output_prefix}.pgen" "$expected_pgen"
fi
if [[ "${output_prefix}.pvar" != "$expected_pvar" ]]; then
    cp "${output_prefix}.pvar" "$expected_pvar"
fi
if [[ "${output_prefix}.psam" != "$expected_psam" ]]; then
    cp "${output_prefix}.psam" "$expected_psam"
fi

log "Output pgen: ${output_prefix}.pgen"
log "Output pvar: ${output_prefix}.pvar"
log "Output psam: ${output_prefix}.psam"

cleanup_intermediate

log "Merge complete"
touch "$done_file"
log "Done marker created at $done_file"
