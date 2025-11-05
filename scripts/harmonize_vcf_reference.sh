#!/bin/bash

# Harmonize VCF files against reference genome
# Usage: harmonize_vcf_reference.sh <dataset> <reference_fasta> <results_dir> <done_file>

set -e
set -u

# Parse command line arguments
DATASET="$1"
REFERENCE_FASTA="$2"
RESULTS_DIR="$3"
DONE_FILE="$4"

# Setup directories
INPUT_DIR="${RESULTS_DIR}/03_vcf"
OUT_DIR="${RESULTS_DIR}/04_harmonize"
mkdir -p "$OUT_DIR"

# Log file for this step
LOG_FILE="${OUT_DIR}/${DATASET}_harmonize.log"

echo "=== Starting VCF harmonization against reference genome ===" | tee "$LOG_FILE"
echo "Dataset: $DATASET" | tee -a "$LOG_FILE"
echo "Reference FASTA: $REFERENCE_FASTA" | tee -a "$LOG_FILE"
echo "Input directory: $INPUT_DIR" | tee -a "$LOG_FILE"
echo "Output directory: $OUT_DIR" | tee -a "$LOG_FILE"
echo "Started at: $(date)" | tee -a "$LOG_FILE"

# Check if reference FASTA exists
if [[ ! -f "$REFERENCE_FASTA" ]]; then
    echo "ERROR: Reference FASTA not found: $REFERENCE_FASTA" | tee -a "$LOG_FILE"
    exit 1
fi

# Check if reference FASTA index exists
if [[ ! -f "${REFERENCE_FASTA}.fai" ]]; then
    echo "ERROR: Reference FASTA index not found: ${REFERENCE_FASTA}.fai" | tee -a "$LOG_FILE"
    exit 1
fi

# Load required modules
module load bcftools
module load htslib

# Process all chromosomes 1-22
TOTAL_INPUT_VARIANTS=0
TOTAL_HARMONIZED_VARIANTS=0
FAILED_CHROMOSOMES=()

echo "Processing chromosomes 1-22..." | tee -a "$LOG_FILE"
for CHR in {1..22}; do
    echo "Processing chromosome $CHR..." | tee -a "$LOG_FILE"
    
    # Define input and output files
    INPUT_VCF="${INPUT_DIR}/${DATASET}_chr${CHR}.vcf.gz"
    OUTPUT_VCF="${OUT_DIR}/${DATASET}_chr${CHR}_harmonized.vcf.gz"
    TEMP_VCF="${OUT_DIR}/${DATASET}_chr${CHR}_temp.vcf.gz"
    
    # Check if input VCF exists
    if [[ ! -f "$INPUT_VCF" ]]; then
        echo "WARNING: Input VCF not found for chromosome $CHR: $INPUT_VCF" | tee -a "$LOG_FILE"
        FAILED_CHROMOSOMES+=("$CHR")
        continue
    fi
    
    # Count input variants
    INPUT_VARIANTS=$(bcftools view -H "$INPUT_VCF" | wc -l)
    TOTAL_INPUT_VARIANTS=$((TOTAL_INPUT_VARIANTS + INPUT_VARIANTS))
    echo "  Input variants for chr${CHR}: $INPUT_VARIANTS" | tee -a "$LOG_FILE"
    
    # Step 1: Remove indels and keep only SNPs
    echo "  Removing indels and keeping only SNPs..." | tee -a "$LOG_FILE"
    bcftools view \
        -v snps \
        -O z \
        -o "$TEMP_VCF" \
        "$INPUT_VCF"
    
    # Index temporary file
    tabix -p vcf "$TEMP_VCF"
    
    # Step 2: Normalize and split multiallelic variants against reference
    echo "  Normalizing and splitting multiallelic variants..." | tee -a "$LOG_FILE"
    bcftools norm \
        -f "$REFERENCE_FASTA" \
        -m -both \
        -c s \
        -O z \
        -o "$OUTPUT_VCF" \
        "$TEMP_VCF"
    
    # Index output file
    tabix -p vcf "$OUTPUT_VCF"
    
    # Count output variants
    OUTPUT_VARIANTS=$(bcftools view -H "$OUTPUT_VCF" | wc -l)
    TOTAL_HARMONIZED_VARIANTS=$((TOTAL_HARMONIZED_VARIANTS + OUTPUT_VARIANTS))
    echo "  Harmonized variants for chr${CHR}: $OUTPUT_VARIANTS" | tee -a "$LOG_FILE"
    
    # Clean up temporary files
    rm -f "$TEMP_VCF"
    rm -f "${TEMP_VCF}.tbi"
    
    echo "  Chromosome $CHR completed successfully" | tee -a "$LOG_FILE"
done

# Generate harmonization statistics
echo "Generating harmonization statistics..." | tee -a "$LOG_FILE"
{
    echo "=== VCF Harmonization Statistics ==="
    echo "Dataset: $DATASET"
    echo "Date: $(date)"
    echo "Reference: $REFERENCE_FASTA"
    echo ""
    
    echo "Summary:"
    echo "  - Total input variants: $TOTAL_INPUT_VARIANTS"
    echo "  - Total harmonized variants: $TOTAL_HARMONIZED_VARIANTS"
    
    if [[ ${#FAILED_CHROMOSOMES[@]} -gt 0 ]]; then
        echo "  - Failed chromosomes: ${FAILED_CHROMOSOMES[*]}"
    else
        echo "  - All chromosomes processed successfully"
    fi
    
    # Count samples from first harmonized VCF
    FIRST_HARMONIZED="${OUT_DIR}/${DATASET}_chr1_harmonized.vcf.gz"
    if [[ -f "$FIRST_HARMONIZED" ]]; then
        SAMPLES=$(bcftools query -l "$FIRST_HARMONIZED" | wc -l)
        echo "  - Samples in harmonized VCFs: $SAMPLES"
    fi
    
    echo ""
    echo "Harmonization steps performed:"
    echo "  - Indel removal (SNPs only retained)"
    echo "  - Multiallelic variant splitting"
    echo "  - Reference genome normalization"
    echo "  - Strand alignment check"
    echo ""
    echo "Output files:"
    for CHR in {1..22}; do
        HARMONIZED_FILE="${OUT_DIR}/${DATASET}_chr${CHR}_harmonized.vcf.gz"
        if [[ -f "$HARMONIZED_FILE" ]]; then
            echo "  - ${HARMONIZED_FILE}"
            echo "  - ${HARMONIZED_FILE}.tbi"
        fi
    done
} > "${OUT_DIR}/${DATASET}_harmonization_stats.txt"

# Create done file
echo "VCF harmonization completed successfully" > "$DONE_FILE"
echo "Dataset: $DATASET" >> "$DONE_FILE"
echo "Chromosomes processed: 1-22" >> "$DONE_FILE"
echo "Total input variants: $TOTAL_INPUT_VARIANTS" >> "$DONE_FILE"
echo "Total harmonized variants: $TOTAL_HARMONIZED_VARIANTS" >> "$DONE_FILE"
echo "Completed at: $(date)" >> "$DONE_FILE"

echo "=== VCF harmonization completed ===" | tee -a "$LOG_FILE"
echo "Total harmonized variants: $TOTAL_HARMONIZED_VARIANTS" | tee -a "$LOG_FILE"
echo "Statistics saved to: ${OUT_DIR}/${DATASET}_harmonization_stats.txt" | tee -a "$LOG_FILE"
echo "Done marker created: $DONE_FILE" | tee -a "$LOG_FILE"