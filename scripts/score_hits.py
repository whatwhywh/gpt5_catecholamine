#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict

import pandas as pd
import yaml
import seaborn as sns
import matplotlib.pyplot as plt


def load_cfg(path: Path) -> Dict:
    return yaml.safe_load(path.read_text())


def classify(score: int) -> str:
    if score >= 8:
        return "High"
    if score >= 5:
        return "Medium"
    return "Low"


def length_rule_ok(length: int, model: str, cfg: Dict) -> bool:
    lr = cfg.get("length_rules", {}).get(model, None)
    if not lr:
        return False
    return lr.get("min", 0) <= length <= lr.get("max", 10**9)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits_dir", required=True)
    ap.add_argument("--out_scores", required=True)
    ap.add_argument("--out_heatmap", required=True)
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    cfg = load_cfg(Path(args.config))
    w = cfg.get("scoring", {})

    rows = []
    for p in Path(args.hits_dir).glob("*.hits.tsv"):
        df = pd.read_csv(p, sep="\t")
        sample = df["sample"].iloc[0] if not df.empty else p.stem
        # aggregate per gene
        gene_groups = df.groupby("gene_id") if not df.empty else []
        for gene_id, g in gene_groups:
            # evidence flags
            has_hmm = any((g["method"] == "HMM") & (g["qcov"].fillna(0) >= cfg.get("min_cov", 0.75)))
            has_rbh = any(g["method"] == "RBH")
            # length proxy from gene_id tag (translate_orf adds len=)
            aa_len = None
            if "|len=" in str(gene_id):
                try:
                    aa_len = int(str(gene_id).split("|len=")[-1])
                except Exception:
                    aa_len = None
            length_ok = False
            if has_hmm:
                # try to infer model by highest bits among HMM hits for this gene
                hmm_hits = g[g["method"] == "HMM"].sort_values("bits", ascending=False)
                top_model = hmm_hits["model"].iloc[0] if len(hmm_hits) else "TDC"
                if aa_len is not None:
                    length_ok = length_rule_ok(aa_len, str(top_model), cfg)
            # score
            score = 0
            score += w.get("hmm", 3) if has_hmm else 0
            score += w.get("rbh", 2) if has_rbh else 0
            score += w.get("domain_length", 1) if length_ok else 0
            # placement/active_site not yet computed in skeleton
            verdict = classify(score)
            rows.append({
                "sample": sample,
                "gene_id": gene_id,
                "hmm": int(has_hmm),
                "rbh": int(has_rbh),
                "length_ok": int(length_ok),
                "score": score,
                "verdict": verdict,
            })

    out_dir = Path(args.out_scores).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    df_all = pd.DataFrame(rows)
    if df_all.empty:
        # create empty with header
        df_all = pd.DataFrame(columns=["sample","gene_id","hmm","rbh","length_ok","score","verdict"])
    df_all.to_csv(args.out_scores, sep="\t", index=False)

    # Per-strain best score
    if not df_all.empty:
        best = df_all.groupby("sample")["score"].max().reset_index(name="best_score")
        # verdict per strain by best score
        best["verdict"] = best["best_score"].apply(classify)
        # heatmap
        heat = best.set_index("sample")["best_score"].to_frame()
        plt.figure(figsize=(max(6, 0.3 * len(heat)), 6))
        sns.heatmap(heat, annot=True, cmap="viridis", cbar=True)
        plt.tight_layout()
        Path(args.out_heatmap).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(args.out_heatmap, dpi=200)
        plt.close()
    else:
        # emit an empty placeholder heatmap
        plt.figure(figsize=(4, 3))
        plt.text(0.5, 0.5, "No hits to plot", ha="center", va="center")
        plt.axis("off")
        plt.savefig(args.out_heatmap, dpi=200)
        plt.close()


if __name__ == "__main__":
    main()

