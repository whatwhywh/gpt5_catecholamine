#!/bin/bash
# BLAST search script for catecholamine decarboxylase identification

if [ $# -ne 2 ]; then
    echo "Usage: $0 <query_file> <output_dir>"
    echo "Example: $0 my_proteins.faa results/"
    exit 1
fi

QUERY=$1
OUTPUT_DIR=$2
SCRIPT_DIR=$(dirname "$0")
SEQ_DIR="$SCRIPT_DIR/../sequences"

mkdir -p "$OUTPUT_DIR"

echo "Creating BLAST databases..."
makeblastdb -in "$SEQ_DIR/TDC_sequences.faa" -dbtype prot -out "$OUTPUT_DIR/TDC_db" -title "TDC_database"
makeblastdb -in "$SEQ_DIR/SadA_AADC_sequences.faa" -dbtype prot -out "$OUTPUT_DIR/SadA_AADC_db" -title "SadA_AADC_database"
makeblastdb -in "$SEQ_DIR/other_decarboxylases.faa" -dbtype prot -out "$OUTPUT_DIR/other_db" -title "Other_decarboxylases"

echo "Searching against TDC sequences..."
blastp -query "$QUERY" -db "$OUTPUT_DIR/TDC_db" -out "$OUTPUT_DIR/TDC_hits.txt" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -evalue 1e-5 -max_target_seqs 10

echo "Searching against SadA/AADC sequences..."
blastp -query "$QUERY" -db "$OUTPUT_DIR/SadA_AADC_db" -out "$OUTPUT_DIR/SadA_AADC_hits.txt" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -evalue 1e-5 -max_target_seqs 10

echo "Searching against other decarboxylases..."
blastp -query "$QUERY" -db "$OUTPUT_DIR/other_db" -out "$OUTPUT_DIR/other_hits.txt" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -evalue 1e-5 -max_target_seqs 5

echo "Filtering results (≥70% identity, ≥100 bit score)..."
awk '$3 >= 70 && $12 >= 100' "$OUTPUT_DIR/TDC_hits.txt" > "$OUTPUT_DIR/TDC_filtered.txt"
awk '$3 >= 70 && $12 >= 100' "$OUTPUT_DIR/SadA_AADC_hits.txt" > "$OUTPUT_DIR/SadA_AADC_filtered.txt"
awk '$3 >= 70 && $12 >= 100' "$OUTPUT_DIR/other_hits.txt" > "$OUTPUT_DIR/other_filtered.txt"

echo "Results summary:"
echo "TDC hits: $(wc -l < "$OUTPUT_DIR/TDC_filtered.txt")"
echo "SadA/AADC hits: $(wc -l < "$OUTPUT_DIR/SadA_AADC_filtered.txt")"
echo "Other decarboxylase hits: $(wc -l < "$OUTPUT_DIR/other_filtered.txt")"
