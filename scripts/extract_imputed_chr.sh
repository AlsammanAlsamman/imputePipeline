#!/bin/bash

# Extract password-protected imputed file for a single chromosome
# Usage: extract_imputed_chr.sh <dataset> <chr> <input_dir> <output_dir> <pass_file> <done_file>

set -e
set -u

# Parse command line arguments
DATASET="$1"
CHR="$2"
INPUT_DIR="$3"
OUTPUT_DIR="$4"
PASS_FILE="$5"
DONE_FILE="$6"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "=== Starting extraction for chromosome $CHR ==="
echo "Dataset: $DATASET"
echo "Chromosome: $CHR"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Password file: $PASS_FILE"
echo "Started at: $(date)"

# Check if password file exists
if [[ ! -f "$PASS_FILE" ]]; then
    echo "ERROR: Password file not found: $PASS_FILE"
    exit 1
fi

# Read password from file (assuming password is on first line)
PASSWORD=$(head -n1 "$PASS_FILE" | tr -d '\r\n')
if [[ -z "$PASSWORD" ]]; then
    echo "ERROR: Password file is empty or invalid: $PASS_FILE"
    exit 1
fi

echo "Password file loaded successfully"

# Check if unzip is available
if ! command -v unzip &> /dev/null; then
    echo "ERROR: unzip command not found. Please install unzip utility"
    exit 1
fi

ZIP_FILE=""
# Find the ZIP file for this chromosome using a robust method that handles special characters
echo "Searching for ZIP files for chromosome $CHR..."

# Use find to locate files that match the pattern (handles special characters and CR/LF issues)
FOUND_FILES=$(find "$INPUT_DIR" -name "*chr_${CHR}.zip*" -o -name "*chr${CHR}.zip*" 2>/dev/null)

if [[ -n "$FOUND_FILES" ]]; then
    # Take the first match
    ZIP_FILE=$(echo "$FOUND_FILES" | head -n1)
    echo "Found ZIP file: $ZIP_FILE"
else
    # Try using ls and grep to find files with potential special characters
    DIRECT_MATCH=$(ls -1 "$INPUT_DIR" 2>/dev/null | grep -E "^chr_?${CHR}\.zip" | head -n1)
    if [[ -n "$DIRECT_MATCH" ]]; then
        ZIP_FILE="${INPUT_DIR}/${DIRECT_MATCH}"
        echo "Found ZIP file: $ZIP_FILE"
    fi
fi

if [[ -z "$ZIP_FILE" ]]; then
    echo "ERROR: No ZIP file found for chromosome $CHR"
    echo "Checked files: ${POSSIBLE_FILES[*]}"
    exit 1
fi

echo "Found ZIP file: $ZIP_FILE"

# Extract with password
echo "Extracting with password..."
if unzip -P "$PASSWORD" -j "$ZIP_FILE" -d "$OUTPUT_DIR" > /dev/null 2>&1; then
    # Check if extraction was successful by looking for expected files
    EXTRACTED_FILES=$(find "$OUTPUT_DIR" -type f | wc -l)
    if [[ $EXTRACTED_FILES -gt 0 ]]; then
        echo "SUCCESS: Extracted $EXTRACTED_FILES files for chromosome $CHR"
        
        # List extracted files
        echo "Extracted files:"
        find "$OUTPUT_DIR" -type f -exec basename {} \; | sed 's/^/  /'
        
        # Get file sizes
        TOTAL_SIZE=$(find "$OUTPUT_DIR" -type f -exec stat -f%z {} + 2>/dev/null | awk '{s+=$1} END {print s}' || echo "0")
        echo "Total extracted size: $(numfmt --to=iec $TOTAL_SIZE 2>/dev/null || echo "${TOTAL_SIZE} bytes")"
        
    else
        echo "ERROR: No files extracted for chromosome $CHR"
        exit 1
    fi
else
    echo "ERROR: Failed to extract chromosome $CHR (wrong password or corrupted file?)"
    exit 1
fi

# Create done file with details
{
    echo "Chromosome $CHR extraction completed successfully"
    echo "Dataset: $DATASET"
    echo "ZIP file: $ZIP_FILE"
    echo "Extracted files: $EXTRACTED_FILES"
    echo "Total size: $(numfmt --to=iec $TOTAL_SIZE 2>/dev/null || echo "${TOTAL_SIZE} bytes")"
    echo "Completed at: $(date)"
    echo ""
    echo "Extracted files list:"
    find "$OUTPUT_DIR" -type f -exec basename {} \; | sed 's/^/  /'
} > "$DONE_FILE"

echo "=== Chromosome $CHR extraction completed ==="
echo "Extracted files: $EXTRACTED_FILES"
echo "Done marker created: $DONE_FILE"