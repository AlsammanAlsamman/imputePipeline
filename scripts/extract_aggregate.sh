#!/bin/bash

# Aggregate chromosome extractions and create final done file
# Usage: extract_aggregate.sh <dataset> <chromosomes> <done_file>

set -e
set -u

# Parse command line arguments
DATASET="$1"
CHROMOSOMES="$2"
DONE_FILE="$3"

# Setup directories
EXTRACT_DIR=$(dirname "$DONE_FILE")

echo "=== Aggregating extraction results ==="
echo "Dataset: $DATASET"
echo "Chromosomes: 1-$CHROMOSOMES"
echo "Directory: $EXTRACT_DIR"

# Check that ALL chromosomes were successfully extracted
MISSING_CHROMOSOMES=()
TOTAL_EXTRACTED_FILES=0

for CHR in $(seq 1 $CHROMOSOMES); do
    CHR_DIR="${EXTRACT_DIR}/chr${CHR}"
    CHR_DONE="${EXTRACT_DIR}/chr${CHR}_extract_done.txt"
    
    if [[ ! -d "$CHR_DIR" || ! -f "$CHR_DONE" ]]; then
        MISSING_CHROMOSOMES+=("$CHR")
        echo "ERROR: Missing extraction for chromosome $CHR"
    else
        FILE_COUNT=$(find "$CHR_DIR" -type f | wc -l)
        TOTAL_EXTRACTED_FILES=$((TOTAL_EXTRACTED_FILES + FILE_COUNT))
        echo "Chromosome $CHR: $FILE_COUNT files extracted"
    fi
done

# Fail if any chromosomes are missing
if [[ ${#MISSING_CHROMOSOMES[@]} -gt 0 ]]; then
    echo "FATAL ERROR: Missing extractions for chromosomes: ${MISSING_CHROMOSOMES[*]}"
    echo "All chromosomes must be successfully extracted"
    exit 1
fi

echo "SUCCESS: All $CHROMOSOMES chromosomes extracted successfully"
echo "Total extracted files: $TOTAL_EXTRACTED_FILES"

# Create final done file
{
    echo "Imputed files extraction completed successfully"
    echo "Dataset: $DATASET"
    echo "Chromosomes extracted: 1-$CHROMOSOMES"
    echo "Total extracted files: $TOTAL_EXTRACTED_FILES"
    echo "Completed at: $(date)"
} > "$DONE_FILE"

echo "=== Extraction aggregation completed ==="
echo "Final done marker created: $DONE_FILE"