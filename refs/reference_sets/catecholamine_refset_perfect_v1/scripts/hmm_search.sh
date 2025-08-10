#!/bin/bash
# HMM-based search for catecholamine decarboxylases

if [ $# -ne 2 ]; then
    echo "Usage: $0 <query_file> <output_dir>"
    echo "Example: $0 my_proteins.faa results/"
    exit 1
fi

QUERY=$1
OUTPUT_DIR=$2
SCRIPT_DIR=$(dirname "$0")
HMM_DIR="$SCRIPT_DIR/../hmm"
THRESHOLD_FILE="$SCRIPT_DIR/../metadata/thresholds.tsv"

mkdir -p "$OUTPUT_DIR"

echo "Searching with HMM models..."

# TDC search
if [ -f "$HMM_DIR/TDC.hmm" ]; then
    echo "Searching for TDC..."
    hmmsearch --domtblout "$OUTPUT_DIR/TDC_hits.domtbl" "$HMM_DIR/TDC.hmm" "$QUERY"
    
    # Apply thresholds
    if [ -f "$THRESHOLD_FILE" ]; then
        TDC_THRESHOLD=$(grep "^TDC" "$THRESHOLD_FILE" | grep "recommended" | cut -f3)
        if [ ! -z "$TDC_THRESHOLD" ]; then
            awk -v thresh="$TDC_THRESHOLD" '$8 >= thresh' "$OUTPUT_DIR/TDC_hits.domtbl" > "$OUTPUT_DIR/TDC_filtered.domtbl"
            echo "TDC hits (filtered): $(wc -l < "$OUTPUT_DIR/TDC_filtered.domtbl")"
        fi
    fi
fi

# SadA/AADC search
if [ -f "$HMM_DIR/SadA_AADC.hmm" ]; then
    echo "Searching for SadA/AADC..."
    hmmsearch --domtblout "$OUTPUT_DIR/SadA_AADC_hits.domtbl" "$HMM_DIR/SadA_AADC.hmm" "$QUERY"
    
    # Apply thresholds
    if [ -f "$THRESHOLD_FILE" ]; then
        SADA_THRESHOLD=$(grep "^SadA_AADC" "$THRESHOLD_FILE" | grep "recommended" | cut -f3)
        if [ ! -z "$SADA_THRESHOLD" ]; then
            awk -v thresh="$SADA_THRESHOLD" '$8 >= thresh' "$OUTPUT_DIR/SadA_AADC_hits.domtbl" > "$OUTPUT_DIR/SadA_AADC_filtered.domtbl"
            echo "SadA/AADC hits (filtered): $(wc -l < "$OUTPUT_DIR/SadA_AADC_filtered.domtbl")"
        fi
    fi
fi

echo "Results saved to $OUTPUT_DIR"
echo "Use the filtered files for high-confidence hits"
