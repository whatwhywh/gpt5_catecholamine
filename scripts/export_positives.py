#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Set, Dict

import pandas as pd
from Bio import SeqIO


def load_hits(hits_path: Path) -> Set[str]:
    if not hits_path.exists():
        return set()
    df = pd.read_csv(hits_path, sep="\t")
    accepted: Set[str] = set()
    for _, r in df.iterrows():
        method = str(r.get("method", ""))
        gene_id = str(r.get("gene_id", ""))
        if not gene_id:
            continue
        if method == "RBH":
            accepted.add(gene_id)
        elif method == "HMM" and bool(r.get("passed", False)):
            accepted.add(gene_id)
    return accepted


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--proteins_dir", required=True)
    ap.add_argument("--hits_dir", required=True)
    ap.add_argument("--out_faa", required=True)
    args = ap.parse_args()

    out_path = Path(args.out_faa)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    records = []
    for faa in sorted(Path(args.proteins_dir).glob("*.faa")):
        sample = faa.stem
        hits = Path(args.hits_dir) / f"{sample}.hits.tsv"
        accepted = load_hits(hits)
        if not accepted:
            continue
        id_to_rec: Dict[str, SeqIO.SeqRecord] = {rec.id: rec for rec in SeqIO.parse(str(faa), "fasta")}
        for gid in accepted:
            rec = id_to_rec.get(gid)
            if rec is None:
                continue
            # mark sample in id
            rec = rec[:]
            rec.id = f"{sample}|{rec.id}"
            rec.description = ""
            records.append(rec)

    with out_path.open("w") as fh:
        SeqIO.write(records, fh, "fasta")


if __name__ == "__main__":
    main()

