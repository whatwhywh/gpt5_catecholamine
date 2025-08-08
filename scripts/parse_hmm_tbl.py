#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import List, Dict, Tuple

import pandas as pd
import yaml


DOMTBL_COLS = [
    "target_name", "target_accession", "tlen", "query_name", "query_accession", "qlen",
    "full_evalue", "full_score", "full_bias",
    "dom_num", "dom_of", "c_evalue", "i_evalue", "dom_score", "dom_bias",
    "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc", "target_description"
]


def read_domtbl(path: Path) -> pd.DataFrame:
    rows = []
    with path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 22:
                # skip malformed
                continue
            rows.append(parts[:22])
    if not rows:
        return pd.DataFrame(columns=DOMTBL_COLS)
    df = pd.DataFrame(rows, columns=DOMTBL_COLS)
    # cast numeric
    num_cols = ["tlen", "qlen", "full_evalue", "full_score", "full_bias", "dom_num", "dom_of", "c_evalue", "i_evalue", "dom_score", "dom_bias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc"]
    for c in num_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def load_thresholds(thr_path: Path) -> Dict[str, Dict[str, float]]:
    thr = {}
    if thr_path.exists():
        tdf = pd.read_csv(thr_path, sep="\t")
        for _, r in tdf.iterrows():
            thr[str(r["model"])]= {
                "bitscore": float(r["bitscore_cut"]),
                "evalue": float(r["evalue_hint"]) if not pd.isna(r["evalue_hint"]) else 1e-5,
                "qcov": float(r.get("min_query_cov", 0.0)),
                "mcov": float(r.get("min_model_cov", 0.0)),
            }
    return thr


def parse_hmm_hits(domtbl_files: List[Path], thr: Dict[str, Dict[str, float]]) -> pd.DataFrame:
    hits = []
    for f in domtbl_files:
        if not f.exists():
            continue
        df = read_domtbl(f)
        if df.empty:
            continue
        # Keep best domain per target_name for this model
        df["model"] = df["query_name"].astype(str)
        df["query_id"] = df["target_name"].astype(str)
        df["bits"] = df["full_score"].astype(float)
        # coverage calcs
        df["model_cov"] = (df["hmm_to"] - df["hmm_from"] + 1.0) / df["qlen"].clip(lower=1)
        df["query_cov"] = (df["ali_to"] - df["ali_from"] + 1.0) / df["tlen"].clip(lower=1)
        df_best = df.sort_values(["query_id", "bits"], ascending=[True, False]).drop_duplicates("query_id")
        for _, r in df_best.iterrows():
            model = str(r["model"]).split()[0]
            cutoff = thr.get(model, {"bitscore": -1e9, "evalue": 1e6, "qcov": 0.0, "mcov": 0.0})
            passed = (r["bits"] >= cutoff["bitscore"]) and (r["full_evalue"] <= cutoff["evalue"]) and (r["query_cov"] >= cutoff.get("qcov", 0.0)) and (r["model_cov"] >= cutoff.get("mcov", 0.0))
            hits.append({
                "gene_id": r["query_id"],
                "method": "HMM",
                "model": model,
                "bits": r["bits"],
                "evalue": r["full_evalue"],
                "query_cov": r["query_cov"],
                "model_cov": r["model_cov"],
                "passed": passed,
            })
    if not hits:
        return pd.DataFrame(columns=["gene_id","method","model","bits","evalue","query_cov","model_cov","passed"])
    return pd.DataFrame(hits)


def load_rbh(path: Path) -> pd.DataFrame:
    if not path or not path.exists():
        return pd.DataFrame(columns=["query","target","pident","evalue","bits","qcov","tcov","method"])  # empty
    df = pd.read_csv(path, sep="\t")
    if "method" not in df.columns:
        df["method"] = "RBH"
    df = df.rename(columns={"query": "gene_id"})
    return df


def write_hits(out_path: Path, sample: str, faa: Path, hmm_df: pd.DataFrame, rbh_df: pd.DataFrame):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    if not hmm_df.empty:
        for _, r in hmm_df.iterrows():
            rows.append({
                "sample": sample,
                "gene_id": r["gene_id"],
                "method": r["method"],
                "model": r["model"],
                "bits": r["bits"],
                "evalue": r["evalue"],
                "id": None,
                "qcov": r["query_cov"],
                "mcov": r["model_cov"],
                "domain": None,
                "active_site": None,
                "placement": None,
                "score_stub": 0,
            })
    if not rbh_df.empty:
        for _, r in rbh_df.iterrows():
            rows.append({
                "sample": sample,
                "gene_id": r["gene_id"],
                "method": r.get("method", "RBH"),
                "model": "RBH",  # model not specified here
                "bits": r.get("bits", None),
                "evalue": r.get("evalue", None),
                "id": r.get("pident", None),
                "qcov": r.get("qcov", None),
                "mcov": r.get("tcov", None),
                "domain": None,
                "active_site": None,
                "placement": None,
                "score_stub": 0,
            })
    df = pd.DataFrame(rows)
    df.to_csv(out_path, sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--faa", required=True)
    ap.add_argument("--hmm", required=False, help=", separated domtbl files", default="")
    ap.add_argument("--rbh", required=False)
    ap.add_argument("--thresholds", required=False, default="config/thresholds.seed.tsv")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sample = Path(args.faa).stem
    thr = load_thresholds(Path(args.thresholds))
    domtbl_files = [Path(p) for p in args.hmm.split(",") if p]
    hmm_df = parse_hmm_hits(domtbl_files, thr)
    rbh_df = load_rbh(Path(args.rbh)) if args.rbh else pd.DataFrame()

    write_hits(Path(args.out), sample, Path(args.faa), hmm_df, rbh_df)


if __name__ == "__main__":
    main()

