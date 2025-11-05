#!/bin/bash
#
# VCF conversion and harmonization script for all chromosomes
# Usage: vcf_convert_harmonize.sh <reference_fasta> <input_bed> <input_bim> <input_fam> <dataset> <done_file> <stats_file>
#

set -e

# Parse arguments
REFERENCE_FASTA="$1"
INPUT_BED="$2"
INPUT_BIM="$3"
INPUT_FAM="$4"
DATASET="$5"
DONE_FILE="$6"
STATS_FILE="$7"

# Create output directory
OUT_DIR="results/03_vcf"
mkdir -p "$OUT_DIR"
mkdir -p "results/logs"

echo "=== VCF Conversion and Harmonization for All Chromosomes ==="
echo "Dataset: $DATASET"
echo "Input prefix: ${INPUT_BED%.bed}"
echo "Output directory: $OUT_DIR"
echo "Reference: $REFERENCE_FASTA"

# Validate reference genome exists
if [[ ! -f "$REFERENCE_FASTA" ]]; then
    echo "ERROR: Reference FASTA not found: $REFERENCE_FASTA"
    exit 1
fi

# Get input prefix (remove .bed extension)
INPUT_PREFIX="${INPUT_BED%.bed}"

# Load required modules
module load plink2/1.90b3w
module load bcftools
module load htslib

# Process all chromosomes 1-22
TOTAL_VARIANTS=0
TOTAL_HARMONIZED=0

echo "Processing chromosomes 1-22..."
for CHR in {1..22}; do
    echo "Processing chromosome $CHR..."
    
    TEMP_VCF="${OUT_DIR}/${DATASET}_chr${CHR}_temp.vcf"
    OUT_VCF="${OUT_DIR}/${DATASET}_chr${CHR}.vcf.gz"
    
    # Step 1: Extract chromosome with PLINK
    plink --bfile "$INPUT_PREFIX" \
          --chr "$CHR" \
          --recode vcf \
          --out "${OUT_DIR}/${DATASET}_chr${CHR}_temp" \
          --memory 32000 \
          --threads 4 \
          --allow-no-sex
    
    # Step 2: Harmonize VCF against reference genome
    bcftools norm \
        -f "$REFERENCE_FASTA" \
        -c s \
        -m - \
        -O z \
        -o "$OUT_VCF" \
        "$TEMP_VCF"
    
    # Step 3: Create index file
    tabix -p vcf "$OUT_VCF"
    
    # Count variants for this chromosome
    if [[ -f "$TEMP_VCF" ]]; then
        CHR_VARIANTS=$(grep -v "^#" "$TEMP_VCF" | wc -l)
        TOTAL_VARIANTS=$((TOTAL_VARIANTS + CHR_VARIANTS))
    fi
    
    if [[ -f "$OUT_VCF" ]]; then
        CHR_HARMONIZED=$(bcftools view -H "$OUT_VCF" | wc -l)
        TOTAL_HARMONIZED=$((TOTAL_HARMONIZED + CHR_HARMONIZED))
    fi
    
    # Clean up temporary files
    rm -f "$TEMP_VCF"
    rm -f "${OUT_DIR}/${DATASET}_chr${CHR}_temp.log"
done

# Generate overall VCF conversion statistics
echo "Generating overall VCF conversion statistics..."
{
    echo "=== VCF Conversion and Harmonization Statistics ==="
    echo "Dataset: $DATASET"
    echo "Date: $(date)"
    echo "Chromosomes processed: 1-22"
    echo ""
    
    echo "Summary:"
    echo "  - Total variants from PLINK: $TOTAL_VARIANTS"
    echo "  - Total variants after harmonization: $TOTAL_HARMONIZED"
    
    # Count samples from first VCF
    FIRST_VCF="${OUT_DIR}/${DATASET}_chr1.vcf.gz"
    if [[ -f "$FIRST_VCF" ]]; then
        SAMPLES=$(bcftools query -l "$FIRST_VCF" | wc -l)
        echo "  - Samples in VCFs: $SAMPLES"
    fi
    
    echo ""
    echo "Harmonization steps performed:"
    echo "  - Chromosome extraction (1-22)"
    echo "  - PLINK to VCF conversion"
    echo "  - bcftools normalization against reference"
    echo "  - Multiallelic variant splitting"
    echo "  - VCF compression and indexing"
    echo ""
    echo "Output files:"
    for CHR in {1..22}; do
        echo "  - ${OUT_DIR}/${DATASET}_chr${CHR}.vcf.gz"
        echo "  - ${OUT_DIR}/${DATASET}_chr${CHR}.vcf.gz.tbi"
    done
} > "$STATS_FILE"

# Create done file
echo "VCF conversion and harmonization completed successfully for all chromosomes" > "$DONE_FILE"
echo "Dataset: $DATASET" >> "$DONE_FILE"
echo "Chromosomes: 1-22" >> "$DONE_FILE"
echo "Total variants: $TOTAL_HARMONIZED" >> "$DONE_FILE"
echo "Completed at: $(date)" >> "$DONE_FILE"

echo "=== VCF conversion and harmonization completed for all chromosomes ==="
echo "Total harmonized variants: $TOTAL_HARMONIZED"
echo "Statistics saved to: $STATS_FILE"
echo "Done marker created: $DONE_FILE"