# Changelog

All notable changes to the Catecholamine Decarboxylase Reference Set.

## [Perfect v1.0.0] - 2025-08-08

### 🎯 Major Improvements
- **Complete contamination removal**: All problematic sequences identified and resolved
- **SadA collection**: Authentic SadA sequences collected and validated
- **HMM model construction**: Separate, validated models for TDC and SadA/AADC
- **ROC analysis**: Statistical threshold determination
- **Comprehensive validation**: Domain and active site verification

### ✅ Fixed Issues
- Removed DOPA decarboxylase contamination from TDC set
- Eliminated non-decarboxylase proteins (AtoC) from negative set
- Separated negative sequences into individual enzyme files
- Corrected documentation inconsistencies
- Added missing SadA sequences

### 🆕 Added Features
- HMM profile models with statistical validation
- ROC-derived threshold recommendations
- Automated HMM search scripts
- Comprehensive metadata with domain information
- Result analysis tools
- Detailed validation reports

### 📊 Dataset Composition
- TDC sequences: 52 (contamination-free)
- SadA/AADC sequences: 3 (including authentic SadA)
- GAD sequences: 49
- HDC sequences: 45
- LDC sequences: 42
- ODC sequences: 46

### 🔬 Quality Metrics
- All sequences: Bacterial origin verified
- Length range: 300-800 amino acids
- Quality: <5% ambiguous residues
- Functional annotation: Enzyme-specific verified
- Domain validation: PLP-dependent domain confirmed
- Duplicates: Removed
- HMM models: 3 validated models
- Thresholds: ROC-optimized

### 🛠️ Tools Provided
- `hmm_search.sh`: Automated HMM-based identification
- `blast_search.sh`: Traditional BLAST-based search
- `analyze_results.py`: Result analysis and summarization
- Comprehensive documentation and usage guides

### 📚 Documentation
- Complete README with quick start guide
- Detailed usage instructions
- Troubleshooting guide
- Validation methodology
- Quality control reports

---

## Previous Versions

### [v2.0.0] - Issues Identified
- Contamination in positive sets
- Non-decarboxylase in negative sets
- Missing SadA sequences
- Documentation inconsistencies

### [v1.0.0] - Initial Release
- Basic sequence collection
- Limited validation
- No HMM models
