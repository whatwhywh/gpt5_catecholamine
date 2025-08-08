#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Tuple

import pandas as pd
import yaml


def read_tsv(path: Path) -> pd.DataFrame:
    cols = [
        "query",
        "target",
        "pident",
        "evalue",
        "bits",
        "alnlen",
        "qlen",
        "tlen",
        "qstart",
        "qend",
        "tstart",
        "tend",
        "qcov",
        "tcov",
    ]
    try:
        df = pd.read_csv(path, sep="\t", header=None, names=cols, dtype={"query": str, "target": str})
    except Exception as e:
        raise SystemExit(f"Failed to read {path}: {e}")
    return df


def best_hits(df: pd.DataFrame, key_query: str, key_target: str) -> pd.DataFrame:
    # Choose best by bits desc, then evalue asc
    return (
        df.sort_values(["bits", "evalue"], ascending=[False, True])
        .drop_duplicates(subset=[key_query])[ [key_query, key_target, "pident", "evalue", "bits", "qcov", "tcov", "qlen", "tlen", "alnlen"] ]
        .reset_index(drop=True)
    )


def compute_rbh(q2r: pd.DataFrame, r2q: pd.DataFrame) -> pd.DataFrame:
    b1 = best_hits(q2r, "query", "target")
    b2 = best_hits(r2q, "query", "target").rename(columns={"query": "target", "target": "query", "qlen": "tlen", "tlen": "qlen"})
    rbh = b1.merge(b2[["query", "target"]], on=["query", "target"], how="inner")
    return rbh


def apply_filters(rbh: pd.DataFrame, min_id: float, min_cov: float, max_e: float) -> pd.DataFrame:
    rbf = rbh[(rbh["pident"] >= min_id * 100.0) & (rbh["qcov"] >= min_cov * 100.0) & (rbh["tcov"] >= min_cov * 100.0) & (rbh["evalue"] <= max_e)].copy()
    rbf["method"] = "RBH"
    return rbf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q2r", required=True)
    ap.add_argument("--r2q", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    cfg = yaml.safe_load(Path(args.config).read_text())
    min_id = float(cfg.get("rbh", {}).get("min_id", 0.4))
    min_cov = float(cfg.get("rbh", {}).get("min_cov", 0.8))
    max_e = float(cfg.get("rbh", {}).get("max_evalue", 1e-10))

    q2r = read_tsv(Path(args.q2r))
    r2q = read_tsv(Path(args.r2q))
    rbh = compute_rbh(q2r, r2q)
    rbf = apply_filters(rbh, min_id=min_id, min_cov=min_cov, max_e=max_e)
    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    rbf.to_csv(outp, sep="\t", index=False)


if __name__ == "__main__":
    main()

