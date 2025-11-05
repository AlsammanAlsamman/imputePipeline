#!/bin/bash

# Create harmonization summary after all chromosomes are processed
# Usage: harmonize_summary.sh <dataset> <build> <stats_file> <done_file>

set -e
set -u

# Parse command line arguments
DATASET="$1"
BUILD="$2"
STATS_FILE="$3"
DONE_FILE="$4"

# Setup directories
HARMONIZE_DIR=$(dirname "$STATS_FILE")

echo "=== Generating harmonization summary ==="
echo "Dataset: $DATASET"
echo "Build: $BUILD"
echo "Directory: $HARMONIZE_DIR"

# Count total variants across all chromosomes
TOTAL_INPUT_VARIANTS=0
TOTAL_HARMONIZED_VARIANTS=0
PROCESSED_CHROMOSOMES=()
FAILED_CHROMOSOMES=()

for CHR in {1..22}; do
    HARMONIZED_VCF="${HARMONIZE_DIR}/chr${CHR}.vcf.gz"
    
    if [[ -f "$HARMONIZED_VCF" ]]; then
        HARMONIZED_VARIANTS=$(bcftools view -H "$HARMONIZED_VCF" | wc -l)
        TOTAL_HARMONIZED_VARIANTS=$((TOTAL_HARMONIZED_VARIANTS + HARMONIZED_VARIANTS))
        PROCESSED_CHROMOSOMES+=("$CHR")
        echo "Chromosome $CHR: $HARMONIZED_VARIANTS variants"
    else
        FAILED_CHROMOSOMES+=("$CHR")
        echo "WARNING: Missing harmonized VCF for chromosome $CHR"
    fi
done

# Count samples from first harmonized VCF
FIRST_HARMONIZED="${HARMONIZE_DIR}/chr1.vcf.gz"
SAMPLES=0
if [[ -f "$FIRST_HARMONIZED" ]]; then
    SAMPLES=$(bcftools query -l "$FIRST_HARMONIZED" | wc -l)
fi

# Generate harmonization statistics
{
    echo "=== VCF Harmonization Statistics ==="
    echo "Dataset: $DATASET"
    echo "Genome Build: $BUILD"
    echo "Date: $(date)"
    echo ""
    
    echo "Summary:"
    echo "  - Total harmonized variants: $TOTAL_HARMONIZED_VARIANTS"
    echo "  - Samples in harmonized VCFs: $SAMPLES"
    echo "  - Successfully processed chromosomes: ${#PROCESSED_CHROMOSOMES[@]}"
    
    if [[ ${#FAILED_CHROMOSOMES[@]} -gt 0 ]]; then
        echo "  - Failed chromosomes: ${FAILED_CHROMOSOMES[*]}"
    else
        echo "  - All chromosomes processed successfully"
    fi
    
    echo ""
    echo "Chromosome details:"
    for CHR in "${PROCESSED_CHROMOSOMES[@]}"; do
        HARMONIZED_VCF="${HARMONIZE_DIR}/chr${CHR}.vcf.gz"
        if [[ -f "$HARMONIZED_VCF" ]]; then
            VARIANTS=$(bcftools view -H "$HARMONIZED_VCF" | wc -l)
            echo "  - Chromosome $CHR: $VARIANTS variants"
        fi
    done
    
    echo ""
    echo "Harmonization steps performed:"
    echo "  - Indel removal (SNPs only retained)"
    echo "  - Multiallelic variant splitting"
    echo "  - Reference genome normalization"
    echo "  - Chromosome naming for $BUILD build"
    if [[ "$BUILD" == "hg19" || "$BUILD" == "GRCh37" ]]; then
        echo "    * Chromosomes without prefix (e.g., '1', '2', '22')"
    elif [[ "$BUILD" == "hg38" || "$BUILD" == "GRCh38" ]]; then
        echo "    * Chromosomes with 'chr' prefix (e.g., 'chr1', 'chr2', 'chr22')"
    fi
    
    echo ""
    echo "Output files:"
    for CHR in "${PROCESSED_CHROMOSOMES[@]}"; do
        echo "  - ${HARMONIZE_DIR}/chr${CHR}.vcf.gz"
        echo "  - ${HARMONIZE_DIR}/chr${CHR}.vcf.gz.tbi"
    done
} > "$STATS_FILE"

# Create done file
{
    echo "VCF harmonization completed successfully"
    echo "Dataset: $DATASET"
    echo "Build: $BUILD"
    echo "Chromosomes processed: ${PROCESSED_CHROMOSOMES[*]}"
    echo "Total harmonized variants: $TOTAL_HARMONIZED_VARIANTS"
    echo "Samples: $SAMPLES"
    echo "Completed at: $(date)"
} > "$DONE_FILE"

echo "=== Harmonization summary completed ==="
echo "Processed chromosomes: ${PROCESSED_CHROMOSOMES[*]}"
echo "Total harmonized variants: $TOTAL_HARMONIZED_VARIANTS"
echo "Statistics saved to: $STATS_FILE"
echo "Done marker created: $DONE_FILE"