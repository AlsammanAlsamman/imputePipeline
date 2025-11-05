#!/bin/bash
#
# Remove ambiguous A/T and G/C SNPs script
# Usage: remove_ambiguous_filter.sh <reference_fasta> <input_bed> <input_bim> <input_fam> <dataset> <out_bed> <out_bim> <out_fam> <out_stats> <done_file>
#

set -e

# Parse arguments
REFERENCE_FASTA="$1"
INPUT_BED="$2"
INPUT_BIM="$3"
INPUT_FAM="$4"
DATASET="$5"
OUT_BED="$6"
OUT_BIM="$7"
OUT_FAM="$8"
OUT_STATS="$9"
DONE_FILE="${10}"

# Create output directory
OUT_DIR=$(dirname "$OUT_BED")
mkdir -p "$OUT_DIR"
mkdir -p "results/logs"

echo "=== Remove Ambiguous SNPs ==="
echo "Dataset: $DATASET"
echo "Input prefix: ${INPUT_BED%.bed}"
echo "Output directory: $OUT_DIR"
echo "Reference: $REFERENCE_FASTA"

# Get input prefix (remove .bed extension)
INPUT_PREFIX="${INPUT_BED%.bed}"
OUT_PREFIX="${OUT_BED%.bed}"

# Load PLINK module
module load plink2/1.90b3w

# Create list of ambiguous A/T and G/C SNPs to exclude
echo "Identifying ambiguous A/T and G/C SNPs..."
echo "Checking .bim file format..."
head -5 "$INPUT_BIM"

AMBIGUOUS_SNPS="${OUT_DIR}/${DATASET}_ambiguous_snps.txt"

# Check for A/T and G/C SNPs only (columns 5 and 6 in .bim file)
awk 'NF>=6 {
    allele1=$5; allele2=$6;
    if ((allele1=="A" && allele2=="T") || 
        (allele1=="T" && allele2=="A") || 
        (allele1=="G" && allele2=="C") || 
        (allele1=="C" && allele2=="G")) {
        print $2
    }
}' "$INPUT_BIM" > "$AMBIGUOUS_SNPS"

echo "Sample of ambiguous SNPs found:"
head -5 "$AMBIGUOUS_SNPS"

# Count ambiguous SNPs
AMBIGUOUS_COUNT=$(wc -l < "$AMBIGUOUS_SNPS")
echo "Found $AMBIGUOUS_COUNT ambiguous A/T and G/C SNPs to remove"

# Remove ambiguous SNPs using PLINK
echo "Removing ambiguous SNPs with PLINK..."
plink --bfile "$INPUT_PREFIX" \
      --exclude "$AMBIGUOUS_SNPS" \
      --make-bed \
      --out "$OUT_PREFIX" \
      --memory 24000 \
      --threads 2 \
      --allow-no-sex

# Generate remove ambiguous statistics
echo "Generating remove ambiguous statistics..."
{
    echo "=== Remove Ambiguous SNPs Statistics ==="
    echo "Dataset: $DATASET"
    echo "Date: $(date)"
    echo "Step: Remove ambiguous SNPs (A/T and G/C)"
    echo ""
    
    # Count SNPs before and after
    if [[ -f "$INPUT_BIM" ]] && [[ -f "$OUT_BIM" ]]; then
        INPUT_SNPS=$(wc -l < "$INPUT_BIM")
        OUTPUT_SNPS=$(wc -l < "$OUT_BIM")
        echo "SNPs before ambiguous removal: $INPUT_SNPS"
        echo "Ambiguous SNPs removed: $AMBIGUOUS_COUNT"
        echo "SNPs after ambiguous removal: $OUTPUT_SNPS"
        echo "SNPs retained: $((INPUT_SNPS - AMBIGUOUS_COUNT))"
    fi
    
    # Count samples (should be unchanged)
    if [[ -f "$OUT_FAM" ]]; then
        SAMPLES=$(wc -l < "$OUT_FAM")
        echo "Samples retained: $SAMPLES (unchanged)"
    fi
    
    echo ""
    echo "Output files:"
    echo "  - ${OUT_PREFIX}.bed"
    echo "  - ${OUT_PREFIX}.bim" 
    echo "  - ${OUT_PREFIX}.fam"
    echo "  - ${OUT_PREFIX}.log"
} > "$OUT_STATS"

# Clean up temporary files
rm -f "$AMBIGUOUS_SNPS"

# Create done file
echo "Remove ambiguous SNPs completed successfully for dataset: $DATASET" > "$DONE_FILE"
echo "Ambiguous SNPs removed: $AMBIGUOUS_COUNT" >> "$DONE_FILE"
echo "Completed at: $(date)" >> "$DONE_FILE"

echo "=== Remove ambiguous SNPs completed ==="
echo "Statistics saved to: $OUT_STATS"
echo "Done marker created: $DONE_FILE"