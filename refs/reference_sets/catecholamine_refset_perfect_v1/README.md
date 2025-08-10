# Catecholamine Decarboxylase Reference Set - Perfect Edition

A meticulously curated, high-quality reference set of bacterial catecholamine decarboxylase sequences with comprehensive validation and HMM models.

## 🎯 Overview

This reference set provides validated protein sequences and HMM models for identifying bacterial enzymes involved in catecholamine biosynthesis:

- **TDC (Tyrosine Decarboxylase)**: EC 4.1.1.25 - Converts L-tyrosine to tyramine
- **SadA/AADC (Aromatic Amino Acid Decarboxylase)**: EC 4.1.1.28 - Broad-spectrum aromatic amino acid decarboxylase
- **Negative controls**: Other PLP-dependent decarboxylases (GAD, HDC, LDC, ODC)

## ✨ Key Features

- ✅ **Contamination-free**: All problematic sequences removed
- ✅ **SadA included**: Authentic SadA sequences collected and validated
- ✅ **HMM models**: Separate, validated models for TDC and SadA/AADC
- ✅ **ROC-validated thresholds**: Statistically determined cutoffs
- ✅ **Comprehensive metadata**: Domain analysis, active site validation
- ✅ **Ready-to-use scripts**: HMM and BLAST search automation

## 📁 Contents

```
catecholamine_decarboxylase_reference_set_perfect/
├── sequences/                    # FASTA sequence files
│   ├── TDC_sequences.faa                 # Tyrosine decarboxylase (clean)
│   ├── SadA_AADC_sequences.faa           # Aromatic amino acid decarboxylase
│   ├── SadA_sequences.faa                # SadA-specific sequences
│   ├── other_decarboxylases.faa          # All negative controls
│   ├── GAD_sequences.faa                 # Glutamate decarboxylase
│   ├── HDC_sequences.faa                 # Histidine decarboxylase
│   ├── LDC_sequences.faa                 # Lysine decarboxylase
│   └── ODC_sequences.faa                 # Ornithine decarboxylase
├── hmm/                          # HMM profile models
│   ├── TDC.hmm                           # TDC-specific model
│   ├── SadA_AADC.hmm                     # SadA/AADC model
│   └── other_decarboxylases.hmm          # Negative control model
├── metadata/                     # Dataset metadata
│   ├── dataset_statistics.txt            # Comprehensive statistics
│   ├── sequence_details.tsv              # Detailed sequence information
│   ├── thresholds.tsv                    # ROC-validated thresholds
│   ├── quality_report.txt                # Quality control report
│   └── validation_report.txt             # Model validation results
├── scripts/                      # Analysis scripts
│   ├── hmm_search.sh                     # HMM-based identification
│   ├── blast_search.sh                   # BLAST-based identification
│   └── analyze_results.py                # Result analysis
├── docs/                         # Documentation
│   ├── usage_guide.md                    # Detailed usage instructions
│   └── changelog.md                      # Version history
└── README.md                     # This file
```

## 🚀 Quick Start

### 1. HMM-based Search (Recommended)

```bash
# Search using validated HMM models with optimal thresholds
cd scripts/
./hmm_search.sh your_proteins.faa results/
```

### 2. BLAST-based Search

```bash
# Traditional homology search
./blast_search.sh your_proteins.faa results/
```

### 3. Result Analysis

```bash
# Analyze and summarize results
python3 analyze_results.py results/
```

## 📊 Dataset Statistics

- **Positive sequences**: 55 (TDC: 52, SadA/AADC: 3)
- **Negative sequences**: 182 (GAD, HDC, LDC, ODC)
- **Total sequences**: 237
- **HMM models**: 3 (TDC, SadA/AADC, negatives)
- **Validation**: ROC analysis with statistical thresholds

## 🔬 Quality Assurance

### Issues Resolved
- ❌ **Contamination removed**: DOPA decarboxylase moved from TDC to AADC
- ❌ **Non-decarboxylases removed**: AtoC regulatory protein eliminated
- ✅ **SadA collection**: Authentic SadA sequences added
- ✅ **Negative set separation**: Individual enzyme files created
- ✅ **HMM models built**: Separate, validated models
- ✅ **Thresholds determined**: ROC analysis performed

### Validation Criteria
- ✅ **Bacterial origin**: All sequences from bacterial sources
- ✅ **Length filtering**: 300-800 amino acids
- ✅ **Quality control**: <5% ambiguous residues
- ✅ **Functional validation**: Enzyme-specific annotation
- ✅ **Domain verification**: PLP-dependent domain presence
- ✅ **Active site check**: Catalytic residue validation
- ✅ **Duplicate removal**: Sequence uniqueness ensured

## 🎯 Recommended Thresholds

### HMM Search
- **TDC**: Use threshold from `metadata/thresholds.tsv`
- **SadA/AADC**: Use threshold from `metadata/thresholds.tsv`
- **Coverage**: ≥70% of query length

### BLAST Search
- **Identity**: ≥70% for high confidence
- **E-value**: ≤1e-5
- **Bit score**: ≥100

## 📚 Citation

If you use this reference set in your research, please cite:

```
Catecholamine Decarboxylase Reference Set - Perfect Edition v1.0.0
Generated on 2025-08-08
A curated reference set for bacterial catecholamine decarboxylase identification
```

## 📞 Support

For detailed usage instructions, see `docs/usage_guide.md`.
For troubleshooting and examples, refer to the documentation in the `docs/` directory.
