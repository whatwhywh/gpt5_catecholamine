# Usage Guide - Perfect Edition

Comprehensive instructions for using the Catecholamine Decarboxylase Reference Set.

## 🎯 Identification Workflow

### Step 1: Choose Your Method

**HMM-based search (Recommended)**
- More sensitive and specific
- Uses validated statistical thresholds
- Better for distant homologs

**BLAST-based search**
- Faster for large datasets
- Good for close homologs
- Easier to interpret

### Step 2: Run the Search

#### HMM Search
```bash
cd scripts/
./hmm_search.sh your_proteins.faa results/
```

This will:
1. Search against TDC and SadA/AADC models
2. Apply validated thresholds automatically
3. Generate filtered results for high-confidence hits

#### BLAST Search
```bash
./blast_search.sh your_proteins.faa results/
```

This will:
1. Create BLAST databases for each enzyme type
2. Search against positive and negative sets
3. Filter results by identity and bit score

### Step 3: Analyze Results

```bash
python3 analyze_results.py results/
```

## 📊 Interpreting Results

### HMM Results

**Output files:**
- `TDC_hits.domtbl`: All TDC hits
- `TDC_filtered.domtbl`: High-confidence TDC hits
- `SadA_AADC_hits.domtbl`: All SadA/AADC hits
- `SadA_AADC_filtered.domtbl`: High-confidence SadA/AADC hits

**Key columns in domtbl files:**
1. Target name (your protein)
2. Target accession
3. Query name (HMM model)
4. Query accession
5. Full sequence E-value
6. Full sequence score
7. Domain E-value
8. Domain score ⭐ (main criterion)
9. Domain bias

**Interpretation:**
- **Filtered files**: Already meet quality thresholds
- **Domain score**: Higher = better match
- **E-value**: Lower = more significant

### BLAST Results

**Output files:**
- `TDC_hits.txt`: All TDC BLAST hits
- `TDC_filtered.txt`: High-quality TDC hits (≥70% identity, ≥100 bit score)
- `SadA_AADC_hits.txt`: All SadA/AADC BLAST hits
- `SadA_AADC_filtered.txt`: High-quality SadA/AADC hits
- `other_hits.txt`: Hits to negative controls

**Columns in BLAST output:**
1. Query ID (your protein)
2. Subject ID (reference protein)
3. % identity ⭐
4. Alignment length
5. Mismatches
6. Gap opens
7-10. Alignment coordinates
11. E-value ⭐
12. Bit score ⭐

## 🎯 Decision Guidelines

### High Confidence (Likely catecholamine decarboxylase)
- **HMM**: Appears in filtered results
- **BLAST**: ≥80% identity to TDC or SadA/AADC, E-value ≤1e-10
- **Length**: 400-700 amino acids
- **Negative check**: No strong hits to other decarboxylases

### Moderate Confidence (Possible candidate)
- **HMM**: Above threshold but lower score
- **BLAST**: 70-80% identity, E-value 1e-10 to 1e-5
- **Additional validation recommended**

### Low Confidence (Unlikely)
- **HMM**: Below threshold
- **BLAST**: <70% identity or E-value >1e-5
- **Strong hits to negative controls**

## 🔬 Advanced Analysis

### Manual Threshold Adjustment

If you need different sensitivity/specificity:

```bash
# More sensitive (lower threshold)
awk '$8 >= 50' TDC_hits.domtbl > TDC_sensitive.domtbl

# More specific (higher threshold)
awk '$8 >= 200' TDC_hits.domtbl > TDC_specific.domtbl
```

### Phylogenetic Analysis

For evolutionary analysis:

```bash
# Extract hit sequences
seqkit grep -f <(awk '{print $1}' TDC_filtered.domtbl) your_proteins.faa > tdc_candidates.faa

# Combine with reference sequences
cat tdc_candidates.faa sequences/TDC_sequences.faa > combined_tdc.faa

# Multiple sequence alignment
mafft combined_tdc.faa > tdc_aligned.faa

# Phylogenetic tree
fasttree tdc_aligned.faa > tdc_tree.newick
```

### Domain Analysis

Check for PLP-dependent domain:

```bash
# Search for PLP domains using Pfam
hmmsearch /path/to/Pfam/PF00282.hmm your_proteins.faa
hmmsearch /path/to/Pfam/PF01276.hmm your_proteins.faa
```

## 🚨 Troubleshooting

### No Hits Found
1. Check input file format (FASTA)
2. Verify sequence length (should be 300-800 aa)
3. Try more sensitive thresholds
4. Check for truncated sequences

### Too Many Hits
1. Use more stringent thresholds
2. Check for contamination in input
3. Verify hits against negative controls

### Conflicting Results
1. HMM and BLAST disagree: Trust HMM for distant homologs
2. Multiple enzyme types: Check domain architecture
3. Weak signals: Require additional validation

## 📋 Validation Checklist

For each candidate:

- [ ] Passes HMM threshold OR ≥70% BLAST identity
- [ ] Appropriate length (300-800 aa)
- [ ] Bacterial origin (if known)
- [ ] No stronger hits to negative controls
- [ ] Contains PLP-dependent domain (if possible to check)
- [ ] Phylogenetically clusters with known enzymes
