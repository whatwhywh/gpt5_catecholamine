## AAAD (TDC/SadA) Screening Pipeline

목표: 장내 미생물 유래 뉴클레오타이드 FASTA에서 TDC/SadA 등 AAAD 보유 여부를 pHMM+RBH(+추후 phylo)로 판정하고, 균주별 점수/판정을 요약.

### 빠른 시작

```bash
mamba env create -f envs/base.yml
conda activate aaad

# refs/positives/*.faa 와 refs/hmm/*.hmm 제공 필요

snakemake -j 16 --use-conda
```

입력: `Fasta_files/*.fa|fna|fasta`

주요 산출:
- `results/hits/{sample}.hits.tsv`
- `results/summary/per_strain_scores.tsv`
- `results/vis/strain_pathway_heatmap.png`

참고: ROC/계통 배치 모듈은 스켈레톤에는 placeholder 입니다.

