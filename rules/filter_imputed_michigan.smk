"""
Snakemake rule for filtering Michigan Imputation Server results.

This rule filters imputed VCF files based on quality metrics from .info.gz files.
Filtering criteria are extracted from analysis.yml configuration.

Author: Genomics Pipeline
Date: November 2025
"""

"""
Snakemake rule for filtering Michigan Imputation Server results.

This rule filters imputed VCF files based on quality metrics from .info.gz files.
Filtering criteria are extracted from analysis.yml configuration.

Author: Genomics Pipeline
Date: November 2025
"""

import os
import yaml
import sys

# Add utils directory to Python path
utils_path = "utils"
if utils_path not in sys.path:
    sys.path.insert(0, utils_path)

from bioconfigme import get_config_value, get_software_module

# Load configuration directly with YAML
def load_config():
    """Load configuration from analysis.yml"""
    config_file = "configs/analysis.yml"
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config

# Load configuration
CONFIG = load_config()

def get_filter_params():
    """Extract filtering parameters from configuration"""
    params = {}
    params['info_score'] = CONFIG['qc_filters']['post_imputation']['info_score']
    params['maf'] = CONFIG['qc_filters']['post_imputation']['maf']
    params['missing_rate'] = CONFIG['qc_filters']['post_imputation']['missing_rate']
    return params

def get_dataset_info():
    """Extract dataset information from configuration"""
    dataset_name = CONFIG['dataset']['name']
    chromosomes = CONFIG['dataset']['chromosomes']
    return dataset_name, chromosomes

# Get configuration values
DATASET_NAME, N_CHROMOSOMES = get_dataset_info()
FILTER_PARAMS = get_filter_params()

# Get software modules
BCFTOOLS_MODULE = get_software_module("bcftools")

# Define input and output directories
EXTRACTED_DIR = "results/06_extracted_data"
FILTERED_DIR = "results/07_filtered"

rule filter_imputed_michigan_chr:
    """
    Filter individual chromosome VCF files from Michigan Imputation Server.
    
    This rule:
    1. Reads quality metrics from .info.gz files
    2. Filters variants based on R², MAF, and INFO score thresholds
    3. Creates filtered VCF files with only high-quality variants
    4. Indexes the filtered VCF files for downstream analysis
    """
    input:
        # Michigan server produces these files per chromosome in subdirectories
        dose_vcf = f"{EXTRACTED_DIR}/chr{{chr}}/chr{{chr}}.dose.vcf.gz",
        info_file = f"{EXTRACTED_DIR}/chr{{chr}}/chr{{chr}}.info.gz"
    output:
        # Filtered VCF file
        filtered_vcf = f"{FILTERED_DIR}/chr{{chr}}.filtered.vcf.gz",
        # Index file
        filtered_index = f"{FILTERED_DIR}/chr{{chr}}.filtered.vcf.gz.csi",
        # Filtering summary
        filter_summary = f"{FILTERED_DIR}/chr{{chr}}.filter_summary.txt",
        # Detailed SNP statistics table
        snp_stats = f"{FILTERED_DIR}/chr{{chr}}.snp_stats.tsv"
    params:
        # Extract filtering thresholds from config
        info_threshold = FILTER_PARAMS['info_score'],
        maf_threshold = FILTER_PARAMS['maf'],
        missing_threshold = FILTER_PARAMS['missing_rate'],
        # Chromosome number
        chromosome = "{chr}",
        # Dataset name for logging
        dataset = DATASET_NAME,
        # Software module to load
        bcftools_module = BCFTOOLS_MODULE
    log:
        f"logs/filter_michigan_chr{{chr}}_{DATASET_NAME}.log"
    threads: 2
    resources:
        mem_mb = 8000,
        time = "01:00:00"
    shell:
        """
        # Load required software modules
        module load {params.bcftools_module}
        
        # Create output directory
        mkdir -p {FILTERED_DIR}
        
        # Log start of filtering
        echo "Starting Michigan variant filtering for chromosome {params.chromosome}" > {log}
        echo "Dataset: {params.dataset}" >> {log}
        echo "Loaded module: {params.bcftools_module}" >> {log}
        echo "Started at: $(date)" >> {log}
        echo "" >> {log}
        
        # Log filtering parameters
        echo "Filter parameters:" >> {log}
        echo "  INFO/R² threshold: {params.info_threshold}" >> {log}
        echo "  MAF threshold: {params.maf_threshold}" >> {log}
        echo "  Missing rate threshold: {params.missing_threshold}" >> {log}
        echo "" >> {log}
        
        # Check input files
        echo "Input files:" >> {log}
        echo "  VCF: {input.dose_vcf}" >> {log}
        echo "  INFO: {input.info_file}" >> {log}
        
        # Step 1: Extract high-quality variants and create SNP stats table
        echo "Step 1: Creating variant filter list and SNP statistics table..." >> {log}
        
        # Create detailed SNP statistics table with header
        echo -e "SNP_ID\\tCHROM\\tPOS\\tREF\\tALT\\tAF\\tMAF\\tAVG_CS\\tR2\\tMISSING_RATE\\tIMPUTED\\tSTATUS" > {output.snp_stats}
        
        zcat {input.info_file} | awk -v info_thresh="{params.info_threshold}" \\
                                    -v maf_thresh="{params.maf_threshold}" \\
                                    -v miss_thresh="{params.missing_threshold}" \\
                                    -v stats_file="{output.snp_stats}" '
        BEGIN {{
            total_variants = 0
            passed_variants = 0
            
            # Initialize counters for detailed statistics
            low_r2 = 0
            low_maf = 0
            high_missing = 0
            multiple_reasons = 0
        }}
        
        # Skip VCF header lines
        /^#/ {{ next }}
        
        {{
            total_variants++
            
            # Extract fields: CHROM POS ID REF ALT QUAL FILTER INFO
            chrom = $1
            pos = $2
            snp_id = $3
            ref = $4
            alt = $5
            info_field = $8
            
            # Parse INFO field to extract quality metrics
            # Example: IMPUTED;AF=0.341481;MAF=0.341481;AVG_CS=0.664044;R2=0.0388952
            
            af = 0
            maf = 0
            avg_cs = 0
            r2 = 0
            imputed = "FALSE"
            
            # Extract AF (Allele Frequency)
            if (match(info_field, /AF=([0-9.e-]+)/, arr)) {{
                af = arr[1]
            }}
            
            # Extract MAF
            if (match(info_field, /MAF=([0-9.e-]+)/, arr)) {{
                maf = arr[1]
            }}
            
            # Extract AVG_CS (call score)
            if (match(info_field, /AVG_CS=([0-9.e-]+)/, arr)) {{
                avg_cs = arr[1]
            }}
            
            # Extract R2 (imputation quality)
            if (match(info_field, /R2=([0-9.e-]+)/, arr)) {{
                r2 = arr[1]
            }}
            
            # Check if imputed
            if (match(info_field, /IMPUTED/)) {{
                imputed = "TRUE"
            }}
            
            # Calculate missing rate
            missing_rate = 1 - avg_cs
            
            # Determine filter status and reasons
            status = "PASS"
            reasons = ""
            reason_count = 0
            
            if (r2 < info_thresh) {{
                status = "FAIL"
                reasons = reasons "LOW_R2;"
                low_r2++
                reason_count++
            }}
            
            if (maf < maf_thresh) {{
                status = "FAIL"
                reasons = reasons "LOW_MAF;"
                low_maf++
                reason_count++
            }}
            
            if (missing_rate > miss_thresh) {{
                status = "FAIL"
                reasons = reasons "HIGH_MISSING;"
                high_missing++
                reason_count++
            }}
            
            if (reason_count > 1) {{
                multiple_reasons++
            }}
            
            # Remove trailing semicolon
            gsub(/;$/, "", reasons)
            if (reasons == "") reasons = "PASS"
            
            # Write to SNP stats table
            printf "%s\\t%s\\t%s\\t%s\\t%s\\t%.6f\\t%.6f\\t%.6f\\t%.6f\\t%.6f\\t%s\\t%s\\n", \\
                snp_id, chrom, pos, ref, alt, af, maf, avg_cs, r2, missing_rate, imputed, reasons >> stats_file
            
            # Add to high-quality variant list if passed
            if (status == "PASS") {{
                print snp_id
                passed_variants++
            }}
        }}
        
        END {{
            print "Total variants: " total_variants > "/dev/stderr"
            print "Variants passing filters: " passed_variants > "/dev/stderr"
            print "Variants failing R2 filter: " low_r2 > "/dev/stderr"
            print "Variants failing MAF filter: " low_maf > "/dev/stderr"
            print "Variants failing missing rate filter: " high_missing > "/dev/stderr"
            print "Variants failing multiple filters: " multiple_reasons > "/dev/stderr"
            if (total_variants > 0) {{
                print "Filter rate: " (100.0 * passed_variants / total_variants) "%" > "/dev/stderr"
            }}
        }}
        ' > {FILTERED_DIR}/chr{params.chromosome}.high_quality_variants.txt 2>> {log}
        
        # Check if any variants passed
        PASSED_VARIANTS=$(wc -l < {FILTERED_DIR}/chr{params.chromosome}.high_quality_variants.txt)
        echo "Variants passing quality filters: $PASSED_VARIANTS" >> {log}
        
        if [[ $PASSED_VARIANTS -eq 0 ]]; then
            echo "ERROR: No variants passed the quality filters!" >> {log}
            echo "This may indicate filter parameters are too strict" >> {log}
            exit 1
        fi
        
        # Step 2: Filter VCF file using bcftools
        echo "Step 2: Filtering VCF file with bcftools..." >> {log}
        
        bcftools view \\
            --include "ID=@{FILTERED_DIR}/chr{params.chromosome}.high_quality_variants.txt" \\
            --output-type z \\
            --threads {threads} \\
            {input.dose_vcf} > {output.filtered_vcf}
        
        # Step 3: Index filtered VCF
        echo "Step 3: Indexing filtered VCF..." >> {log}
        bcftools index {output.filtered_vcf}
        
        # Step 4: Generate comprehensive summary statistics
        FINAL_VARIANTS=$(bcftools view -H {output.filtered_vcf} | wc -l)
        TOTAL_SNPS=$(tail -n +2 {output.snp_stats} | wc -l)
        PASS_SNPS=$(tail -n +2 {output.snp_stats} | awk '$12=="PASS"' | wc -l)
        LOW_R2_SNPS=$(tail -n +2 {output.snp_stats} | awk '$12~"LOW_R2"' | wc -l)
        LOW_MAF_SNPS=$(tail -n +2 {output.snp_stats} | awk '$12~"LOW_MAF"' | wc -l)
        HIGH_MISS_SNPS=$(tail -n +2 {output.snp_stats} | awk '$12~"HIGH_MISSING"' | wc -l)
        IMPUTED_SNPS=$(tail -n +2 {output.snp_stats} | awk '$11=="TRUE"' | wc -l)
        
        # Calculate statistics
        AVG_R2=$(tail -n +2 {output.snp_stats} | awk '$12=="PASS" {{sum+=$9; count++}} END {{if(count>0) printf "%.4f", sum/count; else print "0"}}')
        AVG_MAF=$(tail -n +2 {output.snp_stats} | awk '$12=="PASS" {{sum+=$7; count++}} END {{if(count>0) printf "%.4f", sum/count; else print "0"}}')
        AVG_CALL_RATE=$(tail -n +2 {output.snp_stats} | awk '$12=="PASS" {{sum+=$8; count++}} END {{if(count>0) printf "%.4f", sum/count; else print "0"}}')
        
        {{
            echo "Michigan Filtering Summary - Chromosome {params.chromosome}"
            echo "=========================================="
            echo "Dataset: {params.dataset}"
            echo "Processing date: $(date)"
            echo ""
            echo "Filter parameters:"
            echo "  INFO/R² threshold: {params.info_threshold}"
            echo "  MAF threshold: {params.maf_threshold}"
            echo "  Missing rate threshold: {params.missing_threshold}"
            echo ""
            echo "Variant counts:"
            echo "  Total variants processed: $TOTAL_SNPS"
            echo "  Variants passing all filters: $PASS_SNPS"
            echo "  Variants in final VCF: $FINAL_VARIANTS"
            echo "  Variants failing R² filter: $LOW_R2_SNPS"
            echo "  Variants failing MAF filter: $LOW_MAF_SNPS"
            echo "  Variants failing missing rate filter: $HIGH_MISS_SNPS"
            echo "  Imputed variants: $IMPUTED_SNPS"
            echo ""
            echo "Quality statistics (passing variants only):"
            echo "  Average R²: $AVG_R2"
            echo "  Average MAF: $AVG_MAF"
            echo "  Average call rate: $AVG_CALL_RATE"
            echo ""
            echo "Output files:"
            echo "  Filtered VCF: {output.filtered_vcf}"
            echo "  VCF Index: {output.filtered_index}"
            echo "  SNP Statistics Table: {output.snp_stats}"
            echo "  Summary Report: {output.filter_summary}"
        }} > {output.filter_summary}
        
        # Clean up temporary files
        rm -f {FILTERED_DIR}/chr{params.chromosome}.high_quality_variants.txt
        
        echo "Filtering completed successfully at: $(date)" >> {log}
        """

rule filter_imputed_michigan_all:
    """
    Filter all chromosomes from Michigan Imputation Server in parallel.
    
    This rule coordinates the filtering of all autosomes (chromosomes 1-22)
    and creates a completion marker when all filtering is done.
    """
    input:
        # Wait for all individual chromosome filtering to complete
        expand(f"{FILTERED_DIR}/chr{{chr}}.filtered.vcf.gz", chr=range(1, N_CHROMOSOMES + 1)),
        expand(f"{FILTERED_DIR}/chr{{chr}}.filtered.vcf.gz.csi", chr=range(1, N_CHROMOSOMES + 1)),
        expand(f"{FILTERED_DIR}/chr{{chr}}.filter_summary.txt", chr=range(1, N_CHROMOSOMES + 1)),
        expand(f"{FILTERED_DIR}/chr{{chr}}.snp_stats.tsv", chr=range(1, N_CHROMOSOMES + 1))
    output:
        # Completion marker and summary report
        done_marker = f"{FILTERED_DIR}/michigan_filtering_complete.done",
        summary_report = f"{FILTERED_DIR}/filtering_summary_report.txt",
        # Combined SNP statistics table for all chromosomes
        combined_snp_stats = f"{FILTERED_DIR}/all_chromosomes_snp_stats.tsv"
    params:
        dataset = DATASET_NAME,
        n_chr = N_CHROMOSOMES,
        filter_params = FILTER_PARAMS
    log:
        f"logs/filter_michigan_all_{DATASET_NAME}.log"
    shell:
        """
        echo "Generating filtering summary report for dataset {params.dataset}" > {log}
        echo "Processed {params.n_chr} chromosomes" >> {log}
        echo "Started at: $(date)" >> {log}
        
        # Create combined SNP statistics table
        echo "Creating combined SNP statistics table..." >> {log}
        
        # Start with header from first chromosome
        head -1 {FILTERED_DIR}/chr1.snp_stats.tsv > {output.combined_snp_stats}
        
        # Append data from all chromosomes (skip headers)
        for chr in {{1..{params.n_chr}}}; do
            if [[ -f {FILTERED_DIR}/chr${{chr}}.snp_stats.tsv ]]; then
                tail -n +2 {FILTERED_DIR}/chr${{chr}}.snp_stats.tsv >> {output.combined_snp_stats}
            else
                echo "WARNING: Missing SNP stats for chromosome ${{chr}}" >> {log}
            fi
        done
        
        # Calculate overall statistics from combined table
        TOTAL_SNPS=$(tail -n +2 {output.combined_snp_stats} | wc -l)
        PASS_SNPS=$(tail -n +2 {output.combined_snp_stats} | awk '$12=="PASS"' | wc -l)
        LOW_R2_SNPS=$(tail -n +2 {output.combined_snp_stats} | awk '$12~"LOW_R2"' | wc -l)
        LOW_MAF_SNPS=$(tail -n +2 {output.combined_snp_stats} | awk '$12~"LOW_MAF"' | wc -l)
        HIGH_MISS_SNPS=$(tail -n +2 {output.combined_snp_stats} | awk '$12~"HIGH_MISSING"' | wc -l)
        IMPUTED_SNPS=$(tail -n +2 {output.combined_snp_stats} | awk '$11=="TRUE"' | wc -l)
        
        # Calculate genome-wide averages
        GENOME_AVG_R2=$(tail -n +2 {output.combined_snp_stats} | awk '$12=="PASS" {{sum+=$9; count++}} END {{if(count>0) printf "%.4f", sum/count; else print "0"}}')
        GENOME_AVG_MAF=$(tail -n +2 {output.combined_snp_stats} | awk '$12=="PASS" {{sum+=$7; count++}} END {{if(count>0) printf "%.4f", sum/count; else print "0"}}')
        GENOME_AVG_CALL_RATE=$(tail -n +2 {output.combined_snp_stats} | awk '$12=="PASS" {{sum+=$8; count++}} END {{if(count>0) printf "%.4f", sum/count; else print "0"}}')
        
        echo "Combined SNP statistics:" >> {log}
        echo "  Total SNPs across all chromosomes: $TOTAL_SNPS" >> {log}
        echo "  SNPs passing filters: $PASS_SNPS" >> {log}
        echo "  Genome-wide average R²: $GENOME_AVG_R2" >> {log}
        echo "  Genome-wide average MAF: $GENOME_AVG_MAF" >> {log}
        
        # Create comprehensive summary report
        echo "Michigan Imputation Server - Filtering Summary Report" > {output.summary_report}
        echo "====================================================" >> {output.summary_report}
        echo "Dataset: {params.dataset}" >> {output.summary_report}
        echo "Date: $(date)" >> {output.summary_report}
        echo "Chromosomes processed: {params.n_chr}" >> {output.summary_report}
        echo "" >> {output.summary_report}
        
        echo "Filtering Parameters:" >> {output.summary_report}
        echo "  INFO/R² threshold: {params.filter_params[info_score]}" >> {output.summary_report}
        echo "  MAF threshold: {params.filter_params[maf]}" >> {output.summary_report}
        echo "  Missing rate threshold: {params.filter_params[missing_rate]}" >> {output.summary_report}
        echo "" >> {output.summary_report}
        
        echo "Per-Chromosome Results:" >> {output.summary_report}
        echo "======================" >> {output.summary_report}
        
        # Aggregate statistics from all chromosomes
        TOTAL_VARIANTS=0
        TOTAL_PASSED=0
        
        for chr in {{1..{params.n_chr}}}; do
            if [[ -f {FILTERED_DIR}/chr${{chr}}.filter_stats.txt ]]; then
                echo "Chromosome ${{chr}}:" >> {output.summary_report}
                cat {FILTERED_DIR}/chr${{chr}}.filter_stats.txt | grep -E "(Total variants|Variants passing|Filter rate)" | sed 's/^/  /' >> {output.summary_report}
                echo "" >> {output.summary_report}
                
                # Extract numbers for overall statistics
                CHR_TOTAL=$(grep "Total variants:" {FILTERED_DIR}/chr${{chr}}.filter_stats.txt | awk '{{print $3}}')
                CHR_PASSED=$(grep "Variants passing filters:" {FILTERED_DIR}/chr${{chr}}.filter_stats.txt | awk '{{print $4}}')
                
                TOTAL_VARIANTS=$((TOTAL_VARIANTS + CHR_TOTAL))
                TOTAL_PASSED=$((TOTAL_PASSED + CHR_PASSED))
            else
                echo "WARNING: Missing filter stats for chromosome ${{chr}}" >> {log}
                echo "Chromosome ${{chr}}: ERROR - Missing filter statistics" >> {output.summary_report}
            fi
        done
        
        # Overall summary statistics
        echo "Overall Summary:" >> {output.summary_report}
        echo "===============" >> {output.summary_report}
        echo "Total variants across all chromosomes: $TOTAL_VARIANTS" >> {output.summary_report}
        echo "Total variants passing filters: $TOTAL_PASSED" >> {output.summary_report}
        if [[ $TOTAL_VARIANTS -gt 0 ]]; then
            OVERALL_RATE=$(awk "BEGIN {{printf \\"%.2f\\", 100.0 * $TOTAL_PASSED / $TOTAL_VARIANTS}}")
            echo "Overall filter rate: ${{OVERALL_RATE}}%" >> {output.summary_report}
        fi
        echo "" >> {output.summary_report}
        
        echo "Output Files:" >> {output.summary_report}
        echo "============" >> {output.summary_report}
        for chr in {{1..{params.n_chr}}}; do
            echo "Chromosome ${{chr}}:" >> {output.summary_report}
            echo "  Filtered VCF: {FILTERED_DIR}/chr${{chr}}.filtered.vcf.gz" >> {output.summary_report}
            echo "  Index: {FILTERED_DIR}/chr${{chr}}.filtered.vcf.gz.csi" >> {output.summary_report}
            echo "  Statistics: {FILTERED_DIR}/chr${{chr}}.filter_stats.txt" >> {output.summary_report}
        done
        
        # Create completion marker
        echo "Michigan Imputation Server filtering completed successfully" > {output.done_marker}
        echo "Dataset: {params.dataset}" >> {output.done_marker}
        echo "Chromosomes: {params.n_chr}" >> {output.done_marker}
        echo "Completion time: $(date)" >> {output.done_marker}
        echo "Total variants: $TOTAL_VARIANTS" >> {output.done_marker}
        echo "Filtered variants: $TOTAL_PASSED" >> {output.done_marker}
        
        echo "Filtering summary completed at: $(date)" >> {log}
        echo "Summary report: {output.summary_report}" >> {log}
        echo "Completion marker: {output.done_marker}" >> {log}
        """