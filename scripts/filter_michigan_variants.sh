#!/bin/bash

# Filter Michigan Imputation Server variants based on quality metrics
# Usage: filter_michigan_variants.sh <info_file> <dose_vcf> <output_vcf> <info_threshold> <maf_threshold> <missing_threshold> <chr>

set -e  # Exit on any error

# Parse command line arguments
INFO_FILE="$1"
DOSE_VCF="$2"
OUTPUT_VCF="$3"
INFO_THRESHOLD="$4"
MAF_THRESHOLD="$5"
MISSING_THRESHOLD="$6"
CHROMOSOME="$7"
THREADS="${8:-2}"

# Validate input arguments
if [[ $# -lt 7 ]]; then
    echo "Error: Insufficient arguments provided"
    echo "Usage: $0 <info_file> <dose_vcf> <output_vcf> <info_threshold> <maf_threshold> <missing_threshold> <chr> [threads]"
    exit 1
fi

# Check if input files exist
if [[ ! -f "$INFO_FILE" ]]; then
    echo "Error: INFO file not found: $INFO_FILE"
    exit 1
fi

if [[ ! -f "$DOSE_VCF" ]]; then
    echo "Error: VCF file not found: $DOSE_VCF"
    exit 1
fi

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_VCF")
mkdir -p "$OUTPUT_DIR"

# Temporary files
TEMP_DIR=$(mktemp -d)
HIGH_QUALITY_VARIANTS="${TEMP_DIR}/high_quality_variants.txt"
FILTER_STATS="${OUTPUT_DIR}/chr${CHROMOSOME}.filter_stats.txt"

# Cleanup function
cleanup() {
    rm -rf "$TEMP_DIR"
}
trap cleanup EXIT

echo "=========================================="
echo "Michigan Imputation Server Variant Filter"
echo "=========================================="
echo "Input INFO file: $INFO_FILE"
echo "Input VCF file: $DOSE_VCF"
echo "Output VCF file: $OUTPUT_VCF"
echo "Chromosome: $CHROMOSOME"
echo "Filter thresholds:"
echo "  INFO/R² ≥ $INFO_THRESHOLD"
echo "  MAF ≥ $MAF_THRESHOLD"
echo "  Missing rate ≤ $MISSING_THRESHOLD"
echo "Threads: $THREADS"
echo "=========================================="

# Step 1: Extract high-quality variants from INFO VCF file
echo "Step 1: Analyzing quality metrics from VCF INFO file and creating variant filter list..."

# Michigan .info.gz file is actually a VCF file with INFO fields:
# Format: CHROM POS ID REF ALT QUAL FILTER INFO
# INFO field contains: IMPUTED;AF=0.341481;MAF=0.341481;AVG_CS=0.664044;R2=0.0388952

zcat "$INFO_FILE" | awk -v info_thresh="$INFO_THRESHOLD" \
                        -v maf_thresh="$MAF_THRESHOLD" \
                        -v miss_thresh="$MISSING_THRESHOLD" '
BEGIN {
    OFS="\t"
    total_variants = 0
    passed_variants = 0
    print "Processing Michigan VCF INFO file..." > "/dev/stderr"
}

# Skip VCF header lines
/^#/ { 
    if (NR <= 5) print "Header: " $0 > "/dev/stderr"
    next 
}

{
    total_variants++
    
    # Extract fields from VCF
    chrom = $1
    pos = $2
    snp_id = $3
    ref = $4
    alt = $5
    info_field = $8
    
    # Parse INFO field to extract quality metrics
    # Example: IMPUTED;AF=0.341481;MAF=0.341481;AVG_CS=0.664044;R2=0.0388952
    
    # Initialize variables
    maf = 0
    avg_cs = 0
    r2 = 0
    
    # Extract MAF
    if (match(info_field, /MAF=([0-9.e-]+)/, arr)) {
        maf = arr[1]
    }
    
    # Extract AVG_CS (call score - higher is better, 1.0 = 100% called)
    if (match(info_field, /AVG_CS=([0-9.e-]+)/, arr)) {
        avg_cs = arr[1]
    }
    
    # Extract R2 (imputation quality)
    if (match(info_field, /R2=([0-9.e-]+)/, arr)) {
        r2 = arr[1]
    }
    
    # Calculate missing rate (1 - call rate)
    missing_rate = 1 - avg_cs
    
    # Debug: Print first few variants
    if (NR <= 10) {
        printf "Variant %s: MAF=%.4f, R²=%.4f, AVG_CS=%.4f, Miss=%.4f\n", snp_id, maf, r2, avg_cs, missing_rate > "/dev/stderr"
    }
    
    # Apply quality filters
    if (r2 >= info_thresh && maf >= maf_thresh && missing_rate <= miss_thresh) {
        print snp_id
        passed_variants++
    }
}

END {
    printf "Total variants: %d\n", total_variants > "/dev/stderr"
    printf "Variants passing filters: %d\n", passed_variants > "/dev/stderr"
    if (total_variants > 0) {
        printf "Filter rate: %.2f%%\n", (100.0 * passed_variants / total_variants) > "/dev/stderr"
    }
    
    # Save statistics
    printf "variant_filtering_stats\n" > "'$FILTER_STATS'"
    printf "total_variants\t%d\n", total_variants >> "'$FILTER_STATS'"
    printf "passed_variants\t%d\n", passed_variants >> "'$FILTER_STATS'"
    if (total_variants > 0) {
        printf "filter_rate_percent\t%.2f\n", (100.0 * passed_variants / total_variants) >> "'$FILTER_STATS'"
    }
    printf "info_threshold\t%.3f\n", info_thresh >> "'$FILTER_STATS'"
    printf "maf_threshold\t%.3f\n", maf_thresh >> "'$FILTER_STATS'"
    printf "missing_threshold\t%.3f\n", miss_thresh >> "'$FILTER_STATS'"
}' > "$HIGH_QUALITY_VARIANTS"

# Check if any variants passed the filters
VARIANT_COUNT=$(wc -l < "$HIGH_QUALITY_VARIANTS")
echo "Variants passing quality filters: $VARIANT_COUNT"

if [[ $VARIANT_COUNT -eq 0 ]]; then
    echo "Warning: No variants passed the quality filters!"
    echo "This might indicate:"
    echo "  - Thresholds are too stringent"
    echo "  - Input files have quality issues"
    echo "  - File format issues"
    
    # Create empty output file
    touch "$OUTPUT_VCF"
    bcftools index "$OUTPUT_VCF" 2>/dev/null || true
    exit 0
fi

# Step 2: Filter VCF file using high-quality variant list
echo "Step 2: Filtering VCF file with bcftools..."

# Use bcftools to filter VCF based on variant IDs
bcftools view \
    --include "ID=@${HIGH_QUALITY_VARIANTS}" \
    --output-type z \
    --threads "$THREADS" \
    "$DOSE_VCF" > "$OUTPUT_VCF"

# Check if output file was created successfully
if [[ ! -f "$OUTPUT_VCF" ]] || [[ ! -s "$OUTPUT_VCF" ]]; then
    echo "Error: Failed to create filtered VCF file"
    exit 1
fi

# Step 3: Index the filtered VCF file
echo "Step 3: Indexing filtered VCF file..."
bcftools index "$OUTPUT_VCF"

# Step 4: Generate detailed statistics
echo "Step 4: Generating filtering statistics..."

# Count variants in output VCF
OUTPUT_VARIANT_COUNT=$(bcftools view -H "$OUTPUT_VCF" | wc -l)

# Create comprehensive statistics file
{
    echo "Michigan Imputation Server - Filtering Statistics"
    echo "================================================="
    echo "Chromosome: $CHROMOSOME"
    echo "Processing date: $(date)"
    echo ""
    echo "Input files:"
    echo "  INFO file: $INFO_FILE"
    echo "  VCF file: $DOSE_VCF"
    echo ""
    echo "Output files:"
    echo "  Filtered VCF: $OUTPUT_VCF"
    echo "  Index: ${OUTPUT_VCF}.csi"
    echo ""
    echo "Filtering parameters:"
    echo "  INFO/R² threshold: ≥ $INFO_THRESHOLD"
    echo "  MAF threshold: ≥ $MAF_THRESHOLD"  
    echo "  Missing rate threshold: ≤ $MISSING_THRESHOLD"
    echo ""
    echo "Results:"
    if [[ -f "$FILTER_STATS" ]]; then
        TOTAL=$(grep "total_variants" "$FILTER_STATS" | cut -f2)
        PASSED=$(grep "passed_variants" "$FILTER_STATS" | cut -f2)
        RATE=$(grep "filter_rate_percent" "$FILTER_STATS" | cut -f2)
        
        echo "  Total variants in INFO file: $TOTAL"
        echo "  Variants passing filters: $PASSED"
        echo "  Variants in output VCF: $OUTPUT_VARIANT_COUNT"
        echo "  Filter rate: ${RATE}%"
        
        if [[ $PASSED -ne $OUTPUT_VARIANT_COUNT ]]; then
            echo "  Warning: Mismatch between filtered variants ($PASSED) and output VCF ($OUTPUT_VARIANT_COUNT)"
        fi
    else
        echo "  Error: Could not read filtering statistics"
    fi
} > "${OUTPUT_DIR}/chr${CHROMOSOME}.filter_summary.txt"

echo "=========================================="
echo "Filtering completed successfully!"
echo "Output VCF: $OUTPUT_VCF"
echo "Index: ${OUTPUT_VCF}.csi"
echo "Statistics: ${OUTPUT_DIR}/chr${CHROMOSOME}.filter_summary.txt"
echo "Variants retained: $OUTPUT_VARIANT_COUNT"
echo "=========================================="