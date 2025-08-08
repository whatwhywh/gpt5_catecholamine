#!/usr/bin/env python3
"""
Translate putative CDS nucleotide FASTA to protein FASTA.

Heuristic:
- For each sequence, choose the best ORF across 6 frames (longest AA without internal stop "*")
- If multiple ORFs tie, pick the first.

This is a practical fallback when input FASTA contains CDS with unknown frame.
If contigs/genomes are given, upstream workflow should use Prodigal.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, Tuple

from Bio import SeqIO
from Bio.Seq import Seq


def iter_orfs_6frames(nt_seq: Seq) -> Iterable[Tuple[int, bool, str]]:
    """Yield tuples (frame_offset, is_reverse, aa_string) for 3+3 frames.
    frame_offset: 0,1,2
    is_reverse: False for + strand, True for - strand
    """
    seq = nt_seq.upper()
    frames = [seq[i:] for i in range(3)]
    rc = seq.reverse_complement()
    rc_frames = [rc[i:] for i in range(3)]
    for i, f in enumerate(frames):
        yield i, False, str(f.translate(to_stop=False))
    for i, f in enumerate(rc_frames):
        yield i, True, str(f.translate(to_stop=False))


def best_orf_aa(nt_seq: Seq, min_len: int = 30) -> Tuple[str, int, bool]:
    """Pick the longest AA segment without internal stop from 6 frames.
    Returns (aa, frame_offset, is_reverse). If all very short, returns longest anyway.
    """
    best = ("", 0, False)
    best_len = -1
    for frame, is_rev, aa in iter_orfs_6frames(nt_seq):
        # Split by stop codon and choose the longest subsegment
        parts = aa.split("*")
        longest = max(parts, key=len) if parts else ""
        if len(longest) > best_len:
            best = (longest, frame, is_rev)
            best_len = len(longest)
    if best_len < min_len:
        # fallback to naive translation frame 0
        aa0 = str(nt_seq.translate(to_stop=False))
        parts = aa0.split("*")
        best_seg = max(parts, key=len) if parts else ""
        return best_seg, 0, False
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--min_len", type=int, default=30)
    args = ap.parse_args()

    in_path = Path(args.inp)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    records = []
    for rec in SeqIO.parse(str(in_path), "fasta"):
        aa, frame, is_rev = best_orf_aa(rec.seq, args.min_len)
        new_id = f"{rec.id}|bestorf|frame={frame}{'-' if is_rev else '+'}|len={len(aa)}"
        records.append(SeqIO.SeqRecord(Seq(aa), id=new_id, description=""))

    with out_path.open("w") as fh:
        SeqIO.write(records, fh, "fasta")


if __name__ == "__main__":
    main()

