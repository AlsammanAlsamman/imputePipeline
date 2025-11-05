#!/bin/bash
# Helper script for preparing VCF files for PLINK conversion
# Removes problematic variants that cause PLINK failures

# Input parameters
INPUT_VCF="$1"      # Input filtered VCF file
OUTPUT_VCF="$2"     # Output PLINK-ready VCF file
SUMMARY_FILE="$3"   # Brief summary output file
CHROMOSOME="$4"     # Chromosome number
THREADS="$5"        # Number of threads to use

# Validate input parameters
if [[ $# -ne 5 ]]; then
    echo "Usage: $0 <input_vcf> <output_vcf> <summary_file> <chromosome> <threads>"
    echo "Example: $0 chr1.filtered.vcf.gz chr1.plink_ready.vcf.gz chr1.summary.txt 1 2"
    exit 1
fi

if [[ ! -f "$INPUT_VCF" ]]; then
    echo "ERROR: Input VCF file not found: $INPUT_VCF"
    exit 1
fi

# Start timer
START_TIME=$(date +%s)

echo "=== VCF Preparation for PLINK ==="
echo "Chromosome: $CHROMOSOME"
echo "Input: $INPUT_VCF"
echo "Output: $OUTPUT_VCF"
echo "Summary: $SUMMARY_FILE"
echo "================================="
echo ""

# Step 1: Get input statistics
echo "Step 1: Analyzing input VCF..."
INPUT_VARIANTS=$(bcftools view -H "$INPUT_VCF" | wc -l)
echo "Input variants: $INPUT_VARIANTS"

# Step 2: Prepare VCF for PLINK by removing problematic variants
echo "Step 2: Removing problematic variants for PLINK..."
echo "- Keeping only biallelic SNPs"
echo "- Removing multi-allelic variants"
echo "- Filtering out indels and complex variants"

# Use bcftools to create PLINK-ready VCF
bcftools view \
    -m2 -M2 \
    -v snps \
    --threads "$THREADS" \
    -Oz \
    -o "$OUTPUT_VCF" \
    "$INPUT_VCF"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to prepare VCF for PLINK"
    exit 1
fi

# Step 3: Index the output VCF
echo "Step 3: Indexing PLINK-ready VCF..."
bcftools index "$OUTPUT_VCF"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to index PLINK-ready VCF"
    exit 1
fi

# Step 4: Get output statistics
echo "Step 4: Generating statistics..."
OUTPUT_VARIANTS=$(bcftools view -H "$OUTPUT_VCF" | wc -l)
REMOVED_VARIANTS=$((INPUT_VARIANTS - OUTPUT_VARIANTS))

# Step 5: Create brief summary
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

echo "Creating summary: $SUMMARY_FILE"

cat > "$SUMMARY_FILE" << EOF
VCF Preparation for PLINK - Chromosome $CHROMOSOME
=================================================
Date: $(date)
Processing time: ${DURATION} seconds

Input variants: $INPUT_VARIANTS
Output variants: $OUTPUT_VARIANTS
Variants removed: $REMOVED_VARIANTS
Retention rate: $(echo "scale=1; $OUTPUT_VARIANTS * 100 / $INPUT_VARIANTS" | bc)%

Operations performed:
- Removed multi-allelic variants (kept biallelic only)
- Removed indels and complex variants (kept SNPs only)
- Applied bcftools filters: -m2 -M2 -v snps

Files:
- Input: $INPUT_VCF
- Output: $OUTPUT_VCF
- Indexed: ${OUTPUT_VCF}.csi
EOF

echo ""
echo "=== Brief Summary ==="
echo "Input variants: $INPUT_VARIANTS"
echo "Output variants: $OUTPUT_VARIANTS"
echo "Removed: $REMOVED_VARIANTS variants"
echo "Retention: $(echo "scale=1; $OUTPUT_VARIANTS * 100 / $INPUT_VARIANTS" | bc)%"
echo "Processing time: ${DURATION} seconds"
echo "===================="
echo ""

echo "SUCCESS: VCF prepared for PLINK conversion - chromosome $CHROMOSOME"
echo "Output file: $OUTPUT_VCF"
echo "Summary file: $SUMMARY_FILE"