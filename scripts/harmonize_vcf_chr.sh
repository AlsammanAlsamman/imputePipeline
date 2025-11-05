#!/bin/bash

# Harmonize individual chromosome VCF against reference genome
# Usage: harmonize_vcf_chr.sh <dataset> <chr> <reference_fasta> <build> <input_vcf> <output_vcf>

set -e
set -u

# Parse command line arguments
DATASET="$1"
CHR="$2"
REFERENCE_FASTA="$3"
BUILD="$4"
INPUT_VCF="$5"
OUTPUT_VCF="$6"

# Setup directories
OUT_DIR=$(dirname "$OUTPUT_VCF")
mkdir -p "$OUT_DIR"

echo "=== Starting VCF harmonization for chromosome $CHR ==="
echo "Dataset: $DATASET"
echo "Chromosome: $CHR"
echo "Build: $BUILD"
echo "Reference FASTA: $REFERENCE_FASTA"
echo "Input VCF: $INPUT_VCF"
echo "Output VCF: $OUTPUT_VCF"
echo "Started at: $(date)"

# Check if reference FASTA exists
if [[ ! -f "$REFERENCE_FASTA" ]]; then
    echo "ERROR: Reference FASTA not found: $REFERENCE_FASTA"
    exit 1
fi

# Check if reference FASTA index exists
if [[ ! -f "${REFERENCE_FASTA}.fai" ]]; then
    echo "ERROR: Reference FASTA index not found: ${REFERENCE_FASTA}.fai"
    exit 1
fi

# Check if input VCF exists
if [[ ! -f "$INPUT_VCF" ]]; then
    echo "ERROR: Input VCF not found: $INPUT_VCF"
    exit 1
fi

# Load required modules
module load bcftools
module load htslib

# Determine chromosome name based on build
if [[ "$BUILD" == "hg19" || "$BUILD" == "GRCh37" ]]; then
    # hg19/GRCh37: chromosomes without prefix (e.g., "1", "2", "22")
    CHR_NAME="$CHR"
    echo "Using hg19/GRCh37 chromosome naming: $CHR_NAME"
elif [[ "$BUILD" == "hg38" || "$BUILD" == "GRCh38" ]]; then
    # hg38/GRCh38: chromosomes with 'chr' prefix (e.g., "chr1", "chr2", "chr22")
    CHR_NAME="chr${CHR}"
    echo "Using hg38/GRCh38 chromosome naming: $CHR_NAME"
else
    echo "ERROR: Unknown genome build: $BUILD"
    exit 1
fi

# Count input variants
INPUT_VARIANTS=$(bcftools view -H "$INPUT_VCF" | wc -l)
echo "Input variants for chr${CHR}: $INPUT_VARIANTS"

# Create temporary file
TEMP_VCF="${OUT_DIR}/${DATASET}_chr${CHR}_temp.vcf.gz"

# Step 1: Remove indels and keep only SNPs
echo "Removing indels and keeping only SNPs..."
bcftools view \
    -v snps \
    -O z \
    -o "$TEMP_VCF" \
    "$INPUT_VCF"

# Index temporary file
tabix -p vcf "$TEMP_VCF"

# Step 2: Update chromosome names if needed and normalize against reference
echo "Normalizing and splitting multiallelic variants..."
if [[ "$BUILD" == "hg38" || "$BUILD" == "GRCh38" ]]; then
    # For hg38, ensure chromosome names have 'chr' prefix
    bcftools annotate \
        --rename-chrs <(echo -e "${CHR}\tchr${CHR}") \
        "$TEMP_VCF" \
    | bcftools norm \
        -f "$REFERENCE_FASTA" \
        -m -both \
        -c s \
        -O z \
        -o "$OUTPUT_VCF"
else
    # For hg19, chromosome names should be without prefix
    bcftools norm \
        -f "$REFERENCE_FASTA" \
        -m -both \
        -c s \
        -O z \
        -o "$OUTPUT_VCF" \
        "$TEMP_VCF"
fi

# Index output file
tabix -p vcf "$OUTPUT_VCF"

# Count output variants
OUTPUT_VARIANTS=$(bcftools view -H "$OUTPUT_VCF" | wc -l)
echo "Harmonized variants for chr${CHR}: $OUTPUT_VARIANTS"

# Clean up temporary files
rm -f "$TEMP_VCF"
rm -f "${TEMP_VCF}.tbi"

echo "=== Chromosome $CHR harmonization completed ==="
echo "Input variants: $INPUT_VARIANTS"
echo "Output variants: $OUTPUT_VARIANTS"
echo "Completed at: $(date)"