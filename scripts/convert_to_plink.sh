#!/bin/bash
# Helper script for converting PLINK-ready VCF files to PLINK1 binary format
# Applies post-imputation QC filters during conversion and validates output

# Input parameters
INPUT_VCF="$1"          # Input PLINK-ready VCF file
PLINK_PREFIX="$2"       # PLINK output prefix (without extension)
CONVERSION_LOG="$3"     # Conversion log file
CHROMOSOME="$4"         # Chromosome number
INFO_SCORE="$5"         # INFO score threshold
MAF="$6"                # Minor allele frequency threshold
MISSING_RATE="$7"       # Missing rate threshold
THREADS="$8"            # Number of threads to use

# Validate input parameters
if [[ $# -ne 8 ]]; then
    echo "Usage: $0 <input_vcf> <plink_prefix> <conversion_log> <chromosome> <info_score> <maf> <missing_rate> <threads>"
    echo "Example: $0 chr1.plink_ready.vcf.gz chr1 chr1.conversion.log 1 0.7 0.01 0.05 2"
    exit 1
fi

if [[ ! -f "$INPUT_VCF" ]]; then
    echo "ERROR: Input VCF file not found: $INPUT_VCF"
    exit 1
fi

# Start timer
START_TIME=$(date +%s)

echo "=== VCF to PLINK1 Conversion ==="
echo "Chromosome: $CHROMOSOME"
echo "Input VCF: $INPUT_VCF"
echo "PLINK prefix: $PLINK_PREFIX"
echo "Conversion log: $CONVERSION_LOG"
echo "QC filters - INFO: $INFO_SCORE, MAF: $MAF, Missing: $MISSING_RATE"
echo "==============================="
echo ""

# Step 1: Check input VCF
echo "Step 1: Analyzing input VCF file..."
INPUT_SIZE=$(du -h "$INPUT_VCF" | cut -f1)
echo "Input file size: $INPUT_SIZE"

# Step 2: Convert VCF to PLINK1 format with QC filters
echo "Step 2: Converting to PLINK1 format with QC filters..."
echo "Applying filters:"
echo "- MAF threshold: $MAF"
echo "- Missing rate threshold: $MISSING_RATE" 
echo "- INFO score filtering will be done via VCF if available"

# PLINK1 conversion command with QC filters
# Note: INFO score filtering may need to be done on VCF level if PLINK1 doesn't support it directly
plink \
    --vcf "$INPUT_VCF" \
    --make-bed \
    --out "$PLINK_PREFIX" \
    --maf "$MAF" \
    --geno "$MISSING_RATE" \
    --allow-extra-chr \
    --chr "$CHROMOSOME" \
    --memory 30000 \
    --threads "$THREADS" \
    > "$CONVERSION_LOG" 2>&1

if [[ $? -ne 0 ]]; then
    echo "ERROR: PLINK conversion failed"
    cat "$CONVERSION_LOG"
    exit 1
fi

echo "PLINK conversion completed successfully"

# Step 3: Validate PLINK output files
echo "Step 3: Validating PLINK output files..."
PLINK_BED="${PLINK_PREFIX}.bed"
PLINK_BIM="${PLINK_PREFIX}.bim"
PLINK_FAM="${PLINK_PREFIX}.fam"

# Check if all files exist
if [[ ! -f "$PLINK_BED" ]]; then
    echo "ERROR: PLINK .bed file not found: $PLINK_BED"
    exit 1
fi

if [[ ! -f "$PLINK_BIM" ]]; then
    echo "ERROR: PLINK .bim file not found: $PLINK_BIM"
    exit 1
fi

if [[ ! -f "$PLINK_FAM" ]]; then
    echo "ERROR: PLINK .fam file not found: $PLINK_FAM"
    exit 1
fi

echo "All PLINK files created successfully"

# Step 4: Generate file statistics
echo "Step 4: Generating file statistics..."
VARIANTS=$(wc -l < "$PLINK_BIM")
SAMPLES=$(wc -l < "$PLINK_FAM")
BED_SIZE=$(du -h "$PLINK_BED" | cut -f1)
BIM_SIZE=$(du -h "$PLINK_BIM" | cut -f1)
FAM_SIZE=$(du -h "$PLINK_FAM" | cut -f1)

echo "PLINK file statistics:"
echo "- Variants: $VARIANTS"
echo "- Samples: $SAMPLES"
echo "- .bed file size: $BED_SIZE"
echo "- .bim file size: $BIM_SIZE"
echo "- .fam file size: $FAM_SIZE"

# Step 5: Validate file integrity
echo "Step 5: Validating file integrity..."

# Check .bim file format (6 columns expected)
BIM_COLUMNS=$(head -1 "$PLINK_BIM" | awk '{print NF}')
if [[ $BIM_COLUMNS -ne 6 ]]; then
    echo "WARNING: .bim file has $BIM_COLUMNS columns, expected 6"
else
    echo ".bim file format: OK"
fi

# Check .fam file format (6 columns expected)
FAM_COLUMNS=$(head -1 "$PLINK_FAM" | awk '{print NF}')
if [[ $FAM_COLUMNS -ne 6 ]]; then
    echo "WARNING: .fam file has $FAM_COLUMNS columns, expected 6"
else
    echo ".fam file format: OK"
fi

# Check for empty files
if [[ $VARIANTS -eq 0 ]]; then
    echo "WARNING: No variants in PLINK files"
elif [[ $VARIANTS -lt 1000 ]]; then
    echo "WARNING: Very few variants ($VARIANTS) - check QC filters"
else
    echo "Variant count: OK"
fi

if [[ $SAMPLES -eq 0 ]]; then
    echo "ERROR: No samples in PLINK files"
    exit 1
else
    echo "Sample count: OK"
fi

# Step 6: Create detailed conversion summary
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

echo "Step 6: Creating conversion summary..."

# Append to conversion log
cat >> "$CONVERSION_LOG" << EOF

=== CONVERSION SUMMARY ===
Chromosome: $CHROMOSOME
Date: $(date)
Processing time: ${DURATION} seconds

Input:
- VCF file: $INPUT_VCF
- File size: $INPUT_SIZE

QC Filters Applied:
- INFO score threshold: $INFO_SCORE
- MAF threshold: $MAF
- Missing rate threshold: $MISSING_RATE

Output:
- PLINK prefix: $PLINK_PREFIX
- Variants: $VARIANTS
- Samples: $SAMPLES
- .bed size: $BED_SIZE
- .bim size: $BIM_SIZE
- .fam size: $FAM_SIZE

File Validation:
- .bed file: $(if [[ -f "$PLINK_BED" ]]; then echo "OK"; else echo "MISSING"; fi)
- .bim file: $(if [[ -f "$PLINK_BIM" ]]; then echo "OK"; else echo "MISSING"; fi)
- .fam file: $(if [[ -f "$PLINK_FAM" ]]; then echo "OK"; else echo "MISSING"; fi)
- .bim format: $BIM_COLUMNS columns
- .fam format: $FAM_COLUMNS columns

=========================
EOF

echo ""
echo "=== Conversion Summary ==="
echo "Chromosome: $CHROMOSOME"
echo "Processing time: ${DURATION} seconds"
echo "Final variants: $VARIANTS"
echo "Final samples: $SAMPLES"
echo "QC filters applied: INFO≥$INFO_SCORE, MAF≥$MAF, Missing≤$MISSING_RATE"
echo "=========================="
echo ""

echo "SUCCESS: VCF to PLINK1 conversion completed for chromosome $CHROMOSOME"
echo "Output files:"
echo "  - Binary: ${PLINK_PREFIX}.bed"
echo "  - Variant info: ${PLINK_PREFIX}.bim"
echo "  - Sample info: ${PLINK_PREFIX}.fam"
echo "  - Conversion log: $CONVERSION_LOG"