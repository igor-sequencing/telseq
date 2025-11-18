#!/bin/bash

# Script to extract reads from telomeric regions using BED file
# Usage: ./extract_telomeric_reads.sh <input.bam> [output.bam] [threads]

set -e

# Check arguments
if [ $# -lt 1 ]; then
    echo "Usage: $0 <input.bam> [output.bam] [threads]"
    echo ""
    echo "Arguments:"
    echo "  input.bam   - Input BAM/CRAM file"
    echo "  output.bam  - Output BAM file (optional, default: stdout)"
    echo "  threads     - Number of threads for samtools (optional, default: 4)"
    echo ""
    echo "Description:"
    echo "  Extracts reads from telomeric regions (first and last 10kb of each chromosome)"
    echo "  Automatically detects if BAM uses 'chr' prefix or not"
    echo ""
    echo "Examples:"
    echo "  $0 input.bam telomeric_reads.bam 8"
    echo "  $0 input.bam | samtools view -c"
    exit 1
fi

INPUT_BAM="$1"
OUTPUT_BAM="${2:-}"
THREADS="${3:-4}"

# Get script directory to find BED files
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BED_WITH_CHR="${SCRIPT_DIR}/telomeric_regions_grch38.bed"
BED_NO_CHR="${SCRIPT_DIR}/telomeric_regions_grch38_nochr.bed"

# Check if input file exists
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input file not found: $INPUT_BAM" >&2
    exit 1
fi

# Check if BED files exist
if [ ! -f "$BED_WITH_CHR" ] || [ ! -f "$BED_NO_CHR" ]; then
    echo "Error: BED files not found in $SCRIPT_DIR" >&2
    echo "Expected files:" >&2
    echo "  $BED_WITH_CHR" >&2
    echo "  $BED_NO_CHR" >&2
    exit 1
fi

# Detect if BAM uses chr prefix
echo "Detecting chromosome naming convention..." >&2
FIRST_CHR=$(samtools view -H "$INPUT_BAM" | grep "^@SQ" | head -1 | awk '{print $2}' | cut -d: -f2)

if [[ "$FIRST_CHR" == chr* ]]; then
    echo "Detected 'chr' prefix in chromosome names" >&2
    BED_FILE="$BED_WITH_CHR"
else
    echo "Detected no 'chr' prefix in chromosome names" >&2
    BED_FILE="$BED_NO_CHR"
fi

echo "Using BED file: $BED_FILE" >&2
echo "Extracting reads from telomeric regions..." >&2

# Extract reads using samtools view with -L (BED file)
if [ -z "$OUTPUT_BAM" ]; then
    # Output to stdout
    samtools view -@ "$THREADS" -h -L "$BED_FILE" "$INPUT_BAM"
else
    # Output to file
    samtools view -@ "$THREADS" -h -b -L "$BED_FILE" "$INPUT_BAM" -o "$OUTPUT_BAM"
    echo "Indexing output BAM file..." >&2
    samtools index -@ "$THREADS" "$OUTPUT_BAM"

    # Show statistics
    TOTAL_READS=$(samtools view -@ "$THREADS" -c "$OUTPUT_BAM")
    echo "" >&2
    echo "Done!" >&2
    echo "Output file: $OUTPUT_BAM" >&2
    echo "Total telomeric reads: $TOTAL_READS" >&2

    # Show per-region statistics
    echo "" >&2
    echo "Reads per telomeric region:" >&2
    samtools idxstats "$OUTPUT_BAM" | awk 'BEGIN{sum=0} $3>0 {print $1": "$3" reads"; sum+=$3} END{print "Total mapped: "sum}' >&2
fi
