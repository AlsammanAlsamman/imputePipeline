#!/bin/bash

# Download imputed chromosome files in batches of 5
# Usage: download_imputed_batch.sh <dataset> <links_file> <results_dir> <stats_file> <done_file>

set -e
set -u

# Parse command line arguments
DATASET="$1"
LINKS_FILE="$2"
RESULTS_DIR="$3"
STATS_FILE="$4"
DONE_FILE="$5"

# Setup directories
OUT_DIR="${RESULTS_DIR}/05_imputed_data"
mkdir -p "$OUT_DIR"

echo "=== Starting batch download of imputed files ==="
echo "Dataset: $DATASET"
echo "Links file: $LINKS_FILE"
echo "Output directory: $OUT_DIR"
echo "Started at: $(date)"

# Check if links file exists
if [[ ! -f "$LINKS_FILE" ]]; then
    echo "ERROR: Links file not found: $LINKS_FILE"
    exit 1
fi

# Count total URLs
TOTAL_URLS=$(wc -l < "$LINKS_FILE")
echo "Total files to download: $TOTAL_URLS"

# Download files in batches of 5
BATCH_SIZE=5
BATCH_NUM=1
DOWNLOADED_COUNT=0
FAILED_DOWNLOADS=()
SUCCESSFUL_DOWNLOADS=()

echo "Starting downloads in batches of $BATCH_SIZE..."

# Process URLs in batches
while IFS= read -r url || [[ -n "$url" ]]; do
    # Skip empty lines
    [[ -z "$url" ]] && continue
    
    # Extract filename from URL (e.g., chr_1.zip)
    FILENAME=$(basename "${url}")
    OUTPUT_FILE="${OUT_DIR}/${FILENAME}"
    
    echo "Adding to batch $BATCH_NUM: $FILENAME"
    
    # Start download in background
    (
        echo "Downloading $FILENAME..."
        if wget -q --timeout=300 --tries=3 -O "$OUTPUT_FILE" "$url"; then
            # Validate downloaded file
            if [[ -f "$OUTPUT_FILE" && -s "$OUTPUT_FILE" ]]; then
                FILE_SIZE=$(stat -f%z "$OUTPUT_FILE" 2>/dev/null || stat -c%s "$OUTPUT_FILE" 2>/dev/null || echo "0")
                if [[ "$FILE_SIZE" -gt 1000 ]]; then
                    echo "SUCCESS: $FILENAME downloaded ($FILE_SIZE bytes)"
                    echo "$FILENAME" >> "${OUT_DIR}/successful_downloads.tmp"
                else
                    echo "ERROR: $FILENAME appears corrupted (size: $FILE_SIZE bytes)"
                    rm -f "$OUTPUT_FILE"
                    echo "$FILENAME" >> "${OUT_DIR}/failed_downloads.tmp"
                fi
            else
                echo "ERROR: $FILENAME download failed or file is empty"
                echo "$FILENAME" >> "${OUT_DIR}/failed_downloads.tmp"
            fi
        else
            echo "ERROR: wget failed for $FILENAME"
            echo "$FILENAME" >> "${OUT_DIR}/failed_downloads.tmp"
        fi
    ) &
    
    DOWNLOADED_COUNT=$((DOWNLOADED_COUNT + 1))
    
    # If we've reached batch size or end of file, wait for batch to complete
    if [[ $((DOWNLOADED_COUNT % BATCH_SIZE)) -eq 0 ]] || [[ $DOWNLOADED_COUNT -eq $TOTAL_URLS ]]; then
        echo "Waiting for batch $BATCH_NUM to complete..."
        wait
        echo "Batch $BATCH_NUM completed"
        BATCH_NUM=$((BATCH_NUM + 1))
        
        # Small delay between batches
        sleep 2
    fi
    
done < "$LINKS_FILE"

echo "All download batches completed"

# Collect results
SUCCESSFUL_COUNT=0
FAILED_COUNT=0

if [[ -f "${OUT_DIR}/successful_downloads.tmp" ]]; then
    SUCCESSFUL_COUNT=$(wc -l < "${OUT_DIR}/successful_downloads.tmp")
    while IFS= read -r filename; do
        SUCCESSFUL_DOWNLOADS+=("$filename")
    done < "${OUT_DIR}/successful_downloads.tmp"
    rm -f "${OUT_DIR}/successful_downloads.tmp"
fi

if [[ -f "${OUT_DIR}/failed_downloads.tmp" ]]; then
    FAILED_COUNT=$(wc -l < "${OUT_DIR}/failed_downloads.tmp")
    while IFS= read -r filename; do
        FAILED_DOWNLOADS+=("$filename")
    done < "${OUT_DIR}/failed_downloads.tmp"
    rm -f "${OUT_DIR}/failed_downloads.tmp"
fi

# Generate download statistics
{
    echo "=== Imputed Files Download Statistics ==="
    echo "Dataset: $DATASET"
    echo "Date: $(date)"
    echo "Links file: $LINKS_FILE"
    echo ""
    
    echo "Download Summary:"
    echo "  - Total files to download: $TOTAL_URLS"
    echo "  - Successfully downloaded: $SUCCESSFUL_COUNT"
    echo "  - Failed downloads: $FAILED_COUNT"
    echo "  - Success rate: $(( (SUCCESSFUL_COUNT * 100) / TOTAL_URLS ))%"
    
    echo ""
    echo "Batch processing:"
    echo "  - Batch size: $BATCH_SIZE files"
    echo "  - Total batches: $((((TOTAL_URLS - 1) / BATCH_SIZE) + 1))"
    
    if [[ ${#SUCCESSFUL_DOWNLOADS[@]} -gt 0 ]]; then
        echo ""
        echo "Successfully downloaded files:"
        for file in "${SUCCESSFUL_DOWNLOADS[@]}"; do
            FILE_PATH="${OUT_DIR}/${file}"
            if [[ -f "$FILE_PATH" ]]; then
                FILE_SIZE=$(stat -f%z "$FILE_PATH" 2>/dev/null || stat -c%s "$FILE_PATH" 2>/dev/null || echo "0")
                echo "  - $file (${FILE_SIZE} bytes)"
            fi
        done
    fi
    
    if [[ ${#FAILED_DOWNLOADS[@]} -gt 0 ]]; then
        echo ""
        echo "Failed downloads:"
        for file in "${FAILED_DOWNLOADS[@]}"; do
            echo "  - $file"
        done
    fi
    
    echo ""
    echo "Output directory: $OUT_DIR"
    echo "Files in output directory:"
    ls -la "$OUT_DIR"/*.zip 2>/dev/null || echo "  No .zip files found"
    
} > "$STATS_FILE"

# Create done file
{
    echo "Imputed files download completed"
    echo "Dataset: $DATASET"
    echo "Total files: $TOTAL_URLS"
    echo "Successfully downloaded: $SUCCESSFUL_COUNT"
    echo "Failed downloads: $FAILED_COUNT"
    echo "Completed at: $(date)"
} > "$DONE_FILE"

# Check if all downloads succeeded
if [[ $FAILED_COUNT -gt 0 ]]; then
    echo "WARNING: $FAILED_COUNT files failed to download"
    echo "Check $STATS_FILE for details"
    echo "You may need to manually retry failed downloads"
else
    echo "SUCCESS: All $SUCCESSFUL_COUNT files downloaded successfully"
fi

echo "=== Download process completed ==="
echo "Successfully downloaded: $SUCCESSFUL_COUNT/$TOTAL_URLS files"
echo "Statistics saved to: $STATS_FILE"
echo "Done marker created: $DONE_FILE"