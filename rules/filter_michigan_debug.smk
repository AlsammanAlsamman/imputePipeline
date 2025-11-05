"""
Simple debug version of Michigan filtering rule.

This version has minimal dependencies and better error handling.
"""

import os
import yaml

# Load configuration directly
def load_config():
    with open("configs/analysis.yml", 'r') as f:
        return yaml.safe_load(f)

CONFIG = load_config()
DATASET_NAME = CONFIG['dataset']['name']
N_CHROMOSOMES = CONFIG['dataset']['chromosomes']
FILTER_PARAMS = CONFIG['qc_filters']['post_imputation']

# Define directories
EXTRACTED_DIR = "results/06_extracted_data"
FILTERED_DIR = "results/07_filtered"

rule filter_michigan_chr_debug:
    """
    Debug version of Michigan filtering rule.
    """
    input:
        dose_vcf = f"{EXTRACTED_DIR}/chr{{chr}}/chr{{chr}}.dose.vcf.gz",
        info_file = f"{EXTRACTED_DIR}/chr{{chr}}/chr{{chr}}.info.gz"
    output:
        filtered_vcf = f"{FILTERED_DIR}/chr{{chr}}.filtered.vcf.gz",
        filtered_index = f"{FILTERED_DIR}/chr{{chr}}.filtered.vcf.gz.csi",
        filter_summary = f"{FILTERED_DIR}/chr{{chr}}.filter_summary.txt"
    params:
        info_threshold = FILTER_PARAMS['info_score'],
        maf_threshold = FILTER_PARAMS['maf'],
        missing_threshold = FILTER_PARAMS['missing_rate'],
        chromosome = "{chr}",
        dataset = DATASET_NAME
    log:
        f"logs/filter_michigan_debug_chr{{chr}}_{DATASET_NAME}.log"
    threads: 2
    resources:
        mem_mb = 8000,
        time = "01:00:00"
    shell:
        """
        # Create output directory and logs
        mkdir -p {FILTERED_DIR}
        mkdir -p logs
        
        # Start logging
        echo "=== Michigan Filtering Debug ===" > {log}
        echo "Chromosome: {params.chromosome}" >> {log}
        echo "Dataset: {params.dataset}" >> {log}
        echo "Started at: $(date)" >> {log}
        echo "" >> {log}
        
        # Check input files
        echo "Checking input files..." >> {log}
        if [[ -f {input.dose_vcf} ]]; then
            echo "✓ VCF file found: {input.dose_vcf}" >> {log}
            echo "  Size: $(stat -c%s {input.dose_vcf}) bytes" >> {log}
        else
            echo "✗ VCF file missing: {input.dose_vcf}" >> {log}
            exit 1
        fi
        
        if [[ -f {input.info_file} ]]; then
            echo "✓ INFO file found: {input.info_file}" >> {log}
            echo "  Size: $(stat -c%s {input.info_file}) bytes" >> {log}
        else
            echo "✗ INFO file missing: {input.info_file}" >> {log}
            exit 1
        fi
        
        # Check dependencies
        echo "" >> {log}
        echo "Checking dependencies..." >> {log}
        
        if command -v bcftools &> /dev/null; then
            echo "✓ bcftools found: $(which bcftools)" >> {log}
            echo "  Version: $(bcftools --version | head -1)" >> {log}
        else
            echo "✗ bcftools not found" >> {log}
            echo "Loading bcftools module..." >> {log}
            module load bcftools 2>/dev/null || echo "Failed to load bcftools module" >> {log}
        fi
        
        if command -v awk &> /dev/null; then
            echo "✓ awk found: $(which awk)" >> {log}
        else
            echo "✗ awk not found" >> {log}
            exit 1
        fi
        
        # Check first few lines of INFO file
        echo "" >> {log}
        echo "Checking INFO file format..." >> {log}
        echo "First 5 lines of INFO file:" >> {log}
        zcat {input.info_file} | head -5 >> {log}
        
        # Count total variants
        TOTAL_VARIANTS=$(zcat {input.info_file} | wc -l)
        TOTAL_VARIANTS=$((TOTAL_VARIANTS - 1))  # Subtract header
        echo "" >> {log}
        echo "Total variants in INFO file: $TOTAL_VARIANTS" >> {log}
        
        # Create filtered variant list
        echo "" >> {log}
        echo "Creating filtered variant list..." >> {log}
        echo "Filter parameters:" >> {log}
        echo "  INFO/R² threshold: {params.info_threshold}" >> {log}
        echo "  MAF threshold: {params.maf_threshold}" >> {log}
        echo "  Missing rate threshold: {params.missing_threshold}" >> {log}
        
        # Extract high-quality variants
        zcat {input.info_file} | awk '
        BEGIN {{
            OFS="\\t"
            total = 0
            passed = 0
        }}
        NR == 1 {{ next }}  # Skip header
        {{
            total++
            snp_id = $1
            maf = $5
            avg_call = $6
            rsq = $7
            missing_rate = 1 - avg_call
            
            if (rsq >= {params.info_threshold} && 
                maf >= {params.maf_threshold} && 
                missing_rate <= {params.missing_threshold}) {{
                print snp_id
                passed++
            }}
        }}
        END {{
            print "Processed " total " variants, " passed " passed filters" > "/dev/stderr"
        }}
        ' > {FILTERED_DIR}/chr{params.chromosome}.high_quality_variants.txt 2>> {log}
        
        # Count filtered variants
        PASSED_VARIANTS=$(wc -l < {FILTERED_DIR}/chr{params.chromosome}.high_quality_variants.txt)
        echo "Variants passing filters: $PASSED_VARIANTS" >> {log}
        
        if [[ $PASSED_VARIANTS -eq 0 ]]; then
            echo "ERROR: No variants passed filters!" >> {log}
            echo "This suggests the filter parameters may be too strict" >> {log}
            exit 1
        fi
        
        # Filter VCF file
        echo "" >> {log}
        echo "Filtering VCF file..." >> {log}
        
        bcftools view \\
            --include "ID=@{FILTERED_DIR}/chr{params.chromosome}.high_quality_variants.txt" \\
            --output-type z \\
            --threads {threads} \\
            {input.dose_vcf} > {output.filtered_vcf}
        
        if [[ $? -ne 0 ]]; then
            echo "ERROR: bcftools filtering failed" >> {log}
            exit 1
        fi
        
        # Index filtered VCF
        echo "Indexing filtered VCF..." >> {log}
        bcftools index {output.filtered_vcf}
        
        if [[ $? -ne 0 ]]; then
            echo "ERROR: bcftools indexing failed" >> {log}
            exit 1
        fi
        
        # Count final variants
        FINAL_VARIANTS=$(bcftools view -H {output.filtered_vcf} | wc -l)
        echo "Final variants in filtered VCF: $FINAL_VARIANTS" >> {log}
        
        # Create summary
        {{
            echo "Michigan Filtering Summary - Chromosome {params.chromosome}"
            echo "=========================================="
            echo "Dataset: {params.dataset}"
            echo "Processing date: $(date)"
            echo ""
            echo "Input files:"
            echo "  VCF: {input.dose_vcf}"
            echo "  INFO: {input.info_file}"
            echo ""
            echo "Filter parameters:"
            echo "  INFO/R² threshold: {params.info_threshold}"
            echo "  MAF threshold: {params.maf_threshold}"
            echo "  Missing rate threshold: {params.missing_threshold}"
            echo ""
            echo "Results:"
            echo "  Total variants: $TOTAL_VARIANTS"
            echo "  Variants passing filters: $PASSED_VARIANTS"
            echo "  Final variants in VCF: $FINAL_VARIANTS"
            echo "  Filter rate: $(awk "BEGIN {{if($TOTAL_VARIANTS>0) printf \\"%.2f\\", 100.0*$PASSED_VARIANTS/$TOTAL_VARIANTS; else print \\"0\\"}}")%"
        }} > {output.filter_summary}
        
        # Clean up temporary files
        rm -f {FILTERED_DIR}/chr{params.chromosome}.high_quality_variants.txt
        
        echo "" >> {log}
        echo "Filtering completed successfully at: $(date)" >> {log}
        """

rule filter_michigan_all_debug:
    """
    Debug version - filter all chromosomes.
    """
    input:
        expand(f"{FILTERED_DIR}/chr{{chr}}.filtered.vcf.gz", chr=range(1, N_CHROMOSOMES + 1)),
        expand(f"{FILTERED_DIR}/chr{{chr}}.filter_summary.txt", chr=range(1, N_CHROMOSOMES + 1))
    output:
        done_marker = f"{FILTERED_DIR}/michigan_filtering_debug_complete.done"
    shell:
        """
        echo "Michigan filtering debug completed successfully!" > {output.done_marker}
        echo "Completion time: $(date)" >> {output.done_marker}
        """