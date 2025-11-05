"""
Snakemake rule for preparing filtered VCF files for PLINK conversion.

This rule processes filtered Michigan VCF files to remove problematic variants
that could cause PLINK failures (multi-allelic variants, duplicates, etc.).

Author: Genomics Pipeline
Date: November 2025
"""

import os
import sys

# Add utils directory to Python path
utils_path = "utils"
if utils_path not in sys.path:
    sys.path.insert(0, utils_path)

from bioconfigme import get_config_value, get_software_module

# Load configuration values
def get_dataset_info():
    """Extract dataset information from configuration"""
    dataset_name = get_config_value("configs/analysis.yml", "dataset.name")
    chromosomes = get_config_value("configs/analysis.yml", "dataset.chromosomes")
    return dataset_name, chromosomes

# Get configuration values
DATASET_NAME, N_CHROMOSOMES = get_dataset_info()

# Get software modules
BCFTOOLS_MODULE = get_software_module("bcftools")

# Define input and output directories
FILTERED_DIR = "results/07_filtered"
PLINK_READY_DIR = "results/08_plink_ready"

rule prepare_vcf_for_plink_chr:
    """
    Prepare individual chromosome VCF files for PLINK conversion.
    
    This rule removes problematic variants that cause PLINK failures:
    - Multi-allelic variants
    - Duplicate variant IDs
    - Non-SNP variants
    """
    input:
        # Filtered VCF file from previous filtering step
        filtered_vcf = f"{FILTERED_DIR}/chr{{chr}}.filtered.vcf.gz",
        filtered_index = f"{FILTERED_DIR}/chr{{chr}}.filtered.vcf.gz.csi"
    output:
        # PLINK-ready VCF file
        plink_ready_vcf = f"{PLINK_READY_DIR}/chr{{chr}}.plink_ready.vcf.gz",
        plink_ready_index = f"{PLINK_READY_DIR}/chr{{chr}}.plink_ready.vcf.gz.csi",
        # Brief summary of changes
        prep_summary = f"{PLINK_READY_DIR}/chr{{chr}}.prep_summary.txt"
    params:
        # Chromosome number
        chromosome = "{chr}",
        # Dataset name for logging
        dataset = DATASET_NAME,
        # Software module to load
        bcftools_module = BCFTOOLS_MODULE
    log:
        f"results/log/prepare_plink_chr{{chr}}_{DATASET_NAME}.log"
    threads: 2
    resources:
        mem_mb = 32000,
        time = "00:30:00"
    shell:
        """
        # Load required software module
        module load {params.bcftools_module}
        
        # Create output directory
        mkdir -p {PLINK_READY_DIR}
        mkdir -p results/log
        
        # Start logging
        echo "Starting VCF preparation for PLINK - chromosome {params.chromosome}" > {log}
        echo "Dataset: {params.dataset}" >> {log}
        echo "Loaded module: {params.bcftools_module}" >> {log}
        echo "Started at: $(date)" >> {log}
        echo "" >> {log}
        
        # Use helper script for VCF preparation
        bash scripts/prepare_vcf_for_plink.sh \\
            {input.filtered_vcf} \\
            {output.plink_ready_vcf} \\
            {output.prep_summary} \\
            {params.chromosome} \\
            {threads} 2>&1 | tee -a {log}
        
        echo "" >> {log}
        echo "VCF preparation completed at: $(date)" >> {log}
        """

rule prepare_vcf_for_plink_all:
    """
    Prepare all chromosomes for PLINK conversion in parallel.
    
    This rule coordinates the preparation of all autosomes (chromosomes 1-22)
    and creates a completion marker with overall summary.
    """
    input:
        # Wait for all individual chromosome preparations to complete
        expand(f"{PLINK_READY_DIR}/chr{{chr}}.plink_ready.vcf.gz", chr=range(1, N_CHROMOSOMES + 1)),
        expand(f"{PLINK_READY_DIR}/chr{{chr}}.prep_summary.txt", chr=range(1, N_CHROMOSOMES + 1))
    output:
        # Completion marker and brief summary report
        done_marker = f"{PLINK_READY_DIR}/vcf_preparation_complete.done"
    params:
        dataset = DATASET_NAME,
        n_chr = N_CHROMOSOMES
    log:
        f"results/log/prepare_plink_all_{DATASET_NAME}.log"
    resources:
        mem_mb = 32000,
        time = "00:30:00"
    shell:
        """
        echo "Generating VCF preparation summary for dataset {params.dataset}" > {log}
        echo "Processed {params.n_chr} chromosomes" >> {log}
        echo "Started at: $(date)" >> {log}
        
        # Create brief summary report inside done marker
        echo "VCF Preparation for PLINK Complete" > {output.done_marker}
        echo "===================================" >> {output.done_marker}
        echo "Dataset: {params.dataset}" >> {output.done_marker}
        echo "Date: $(date)" >> {output.done_marker}
        echo "Chromosomes processed: {params.n_chr}" >> {output.done_marker}
        echo "" >> {output.done_marker}
        
        # Brief summary statistics
        echo "Summary Statistics:" >> {output.done_marker}
        TOTAL_INPUT=0
        TOTAL_OUTPUT=0
        
        for chr in {{1..{params.n_chr}}}; do
            if [[ -f {PLINK_READY_DIR}/chr${{chr}}.prep_summary.txt ]]; then
                # Extract key numbers from individual summaries
                INPUT_VAR=$(grep "Input variants:" {PLINK_READY_DIR}/chr${{chr}}.prep_summary.txt | awk '{{print $NF}}' || echo "0")
                OUTPUT_VAR=$(grep "Output variants:" {PLINK_READY_DIR}/chr${{chr}}.prep_summary.txt | awk '{{print $NF}}' || echo "0")
                
                if [[ "$INPUT_VAR" =~ ^[0-9]+$ ]]; then
                    TOTAL_INPUT=$((TOTAL_INPUT + INPUT_VAR))
                fi
                if [[ "$OUTPUT_VAR" =~ ^[0-9]+$ ]]; then
                    TOTAL_OUTPUT=$((TOTAL_OUTPUT + OUTPUT_VAR))
                fi
            fi
        done
        
        echo "Total input variants: $TOTAL_INPUT" >> {output.done_marker}
        echo "Total output variants: $TOTAL_OUTPUT" >> {output.done_marker}
        echo "Total variants removed: $((TOTAL_INPUT - TOTAL_OUTPUT))" >> {output.done_marker}
        if [[ $TOTAL_INPUT -gt 0 ]]; then
            echo "Retention rate: $(echo "scale=1; $TOTAL_OUTPUT * 100 / $TOTAL_INPUT" | bc)%" >> {output.done_marker}
        fi
        echo "" >> {output.done_marker}
        echo "Output directory: {PLINK_READY_DIR}" >> {output.done_marker}
        echo "Files ready for PLINK conversion." >> {output.done_marker}
        
        echo "VCF preparation summary completed at: $(date)" >> {log}
        echo "Completion marker: {output.done_marker}" >> {log}
        """