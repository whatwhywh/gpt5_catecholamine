#!/usr/bin/env python3
"""
Result analysis script for catecholamine decarboxylase identification
"""

import sys
import pandas as pd
from pathlib import Path

def analyze_hmm_results(results_dir):
    """Analyze HMM search results"""
    results_dir = Path(results_dir)
    
    print("=== HMM Search Results Analysis ===")
    
    # TDC results
    tdc_file = results_dir / 'TDC_filtered.domtbl'
    if tdc_file.exists():
        with open(tdc_file) as f:
            tdc_hits = [line for line in f if not line.startswith('#')]
        print(f"TDC hits: {len(tdc_hits)}")
        
        if tdc_hits:
            print("Top TDC hits:")
            for i, hit in enumerate(tdc_hits[:5]):
                fields = hit.strip().split()
                if len(fields) >= 8:
                    print(f"  {i+1}. {fields[0]} (score: {fields[7]})")
    
    # SadA/AADC results
    sada_file = results_dir / 'SadA_AADC_filtered.domtbl'
    if sada_file.exists():
        with open(sada_file) as f:
            sada_hits = [line for line in f if not line.startswith('#')]
        print(f"SadA/AADC hits: {len(sada_hits)}")
        
        if sada_hits:
            print("Top SadA/AADC hits:")
            for i, hit in enumerate(sada_hits[:5]):
                fields = hit.strip().split()
                if len(fields) >= 8:
                    print(f"  {i+1}. {fields[0]} (score: {fields[7]})")

def analyze_blast_results(results_dir):
    """Analyze BLAST search results"""
    results_dir = Path(results_dir)
    
    print("\n=== BLAST Search Results Analysis ===")
    
    blast_files = {
        'TDC': results_dir / 'TDC_filtered.txt',
        'SadA/AADC': results_dir / 'SadA_AADC_filtered.txt',
        'Other': results_dir / 'other_filtered.txt'
    }
    
    for category, file_path in blast_files.items():
        if file_path.exists():
            with open(file_path) as f:
                hits = [line for line in f if line.strip()]
            
            print(f"{category} hits: {len(hits)}")
            
            if hits:
                print(f"Top {category} hits:")
                for i, hit in enumerate(hits[:3]):
                    fields = hit.strip().split('\t')
                    if len(fields) >= 12:
                        print(f"  {i+1}. {fields[0]} -> {fields[1]} ({fields[2]}% identity, {fields[11]} bit score)")

def main():
    if len(sys.argv) != 2:
        print('Usage: python3 analyze_results.py <results_directory>')
        sys.exit(1)
    
    results_dir = sys.argv[1]
    
    analyze_hmm_results(results_dir)
    analyze_blast_results(results_dir)

if __name__ == '__main__':
    main()
