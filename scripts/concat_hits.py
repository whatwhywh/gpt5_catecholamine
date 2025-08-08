#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits_dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = []
    for f in sorted(Path(args.hits_dir).glob("*.hits.tsv")):
        try:
            df = pd.read_csv(f, sep="\t")
        except Exception:
            continue
        rows.append(df)
    if rows:
        all_df = pd.concat(rows, ignore_index=True)
    else:
        all_df = pd.DataFrame(columns=["sample","gene_id","method","model","bits","evalue","id","qcov","mcov","domain","active_site","placement","score_stub"])
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    all_df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()

