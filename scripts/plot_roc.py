#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits_dir", required=True)
    ap.add_argument("--thresholds", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # Placeholder ROC: if we had labeled positives/negatives we'd compute TPR/FPR.
    # For now, emit an informative PDF page.
    plt.figure(figsize=(6, 4))
    plt.text(0.5, 0.7, "ROC not computed in skeleton", ha="center")
    plt.text(0.5, 0.5, f"hits: {args.hits_dir}", ha="center")
    plt.text(0.5, 0.3, f"thresholds: {args.thresholds}", ha="center")
    plt.axis("off")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.out)
    plt.close()


if __name__ == "__main__":
    main()

