#!/bin/bash
#
# Input QC filtering script for PLINK files
# Usage: input_qc_filter.sh <geno> <maf> <hwe> <input_plink> <dataset> <out_bed> <out_bim> <out_fam> <out_stats> <done_file>
#

set -e

# Parse arguments
GENO="$1"
MAF="$2"
HWE="$3"
INPUT_PLINK="$4"
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

echo "=== Input QC Filtering ==="
echo "Dataset: $DATASET"
echo "Input: $INPUT_PLINK"
echo "Output directory: $OUT_DIR"
echo "Parameters: geno=$GENO, maf=$MAF, hwe=$HWE (SNP filtering only)"

# Check input files exist
if [[ ! -f "${INPUT_PLINK}.bed" ]] || [[ ! -f "${INPUT_PLINK}.bim" ]] || [[ ! -f "${INPUT_PLINK}.fam" ]]; then
    echo "ERROR: Input PLINK files not found: ${INPUT_PLINK}.{bed,bim,fam}"
    exit 1
fi

# Load PLINK module
module load plink2/1.90b3w

# Get output prefix (remove .bed extension)
OUT_PREFIX="${OUT_BED%.bed}"

# Run PLINK QC filtering
echo "Running PLINK QC filtering..."
plink --bfile "$INPUT_PLINK" \
      --geno "$GENO" \
      --maf "$MAF" \
      --hwe "$HWE" \
      --make-bed \
      --out "$OUT_PREFIX" \
      --memory 32000 \
      --threads 2 \
      --allow-no-sex

# Generate simple statistics
echo "Generating QC statistics..."
{
    echo "=== Input QC Statistics ==="
    echo "Dataset: $DATASET"
    echo "Date: $(date)"
    echo "Parameters used:"
    echo "  - SNP missing rate (--geno): $GENO"
    echo "  - Minor allele frequency (--maf): $MAF" 
    echo "  - Hardy-Weinberg equilibrium (--hwe): $HWE"
    echo "  - Sample filtering: NONE (all samples retained)"
    echo ""
    
    # Count samples and SNPs
    if [[ -f "${OUT_PREFIX}.fam" ]] && [[ -f "${OUT_PREFIX}.bim" ]]; then
        SAMPLES=$(wc -l < "${OUT_PREFIX}.fam")
        SNPS=$(wc -l < "${OUT_PREFIX}.bim")
        echo "Output summary:"
        echo "  - Samples after QC: $SAMPLES"
        echo "  - SNPs after QC: $SNPS"
    fi
    
    echo ""
    echo "Output files:"
    echo "  - ${OUT_PREFIX}.bed"
    echo "  - ${OUT_PREFIX}.bim" 
    echo "  - ${OUT_PREFIX}.fam"
    echo "  - ${OUT_PREFIX}.log"
} > "$OUT_STATS"

# Create done file
echo "Input QC filtering completed successfully for dataset: $DATASET" > "$DONE_FILE"
echo "Completed at: $(date)" >> "$DONE_FILE"

echo "=== Input QC filtering completed ==="
echo "Statistics saved to: $OUT_STATS"
echo "Done marker created: $DONE_FILE"