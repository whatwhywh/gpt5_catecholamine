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

---

# 1) 참조 세트(curated reference) 구축 — 긍·부정 대조군을 함께

**왜:** TDC/SadA는 모두 PLP-의존성 decarboxylase라 GAD/HDC/LDC와 **교차오류**가 납니다. “긍정(정답)”과 “부정(혼동군)”을 함께 써야 경계가 선명해집니다.

* **양성(positives):**

  * *Tyrosine decarboxylase* (TyrDC; Enterococcus, Lactobacillus 등)
  * *SadA-like aromatic amino acid decarboxylase* (Staphylococcus 등; L-DOPA/5-HTP/Trp/Tyr/ Phe 넓은 기질폭)
  * (선택) DDC-like AAAD(세균 유래 보고 사례)
* **음성(negatives):**

  * *Glutamate decarboxylase* (GadA/GadB), *Histidine decarboxylase* (Hdc), *Lysine decarboxylase* (Ldc) 등
* **절차:** Swiss-Prot(리뷰드) 위주로 종 다양성 확보 → 중복 제거(90% 클러스터) → **MAFFT L-INS-i** 정렬 → **trimAl**로 양질 영역만 보존.

결과: `refs/tdc_pos.faa`, `refs/sada_pos.faa`, `refs/dec_neg.faa`, 그리고 병합 정렬.

---

# 2) 계열별 HMM 만들기 — “패밀리-특이” 확정

**왜:** BLAST류는 빠르지만 **특이성**이 떨어집니다. **HMMER** pHMM이 근거로 가장 탄탄합니다.

* `hmmbuild`로 **TDC용 pHMM**, **SadA용 pHMM**, **NEG decarboxylase용 pHMM** 생성
* 교차검증: `hmmsearch`를 positives/negatives에 모두 때려 **ROC 곡선**을 그리고, **bitscore cutoff**를 “TPR≥0.95 & FPR≤0.01” 지점으로 고정(모델별 별도 threshold).
* (보강) **KofamScan**/**eggNOG-mapper**로 orthology 태그 병행 → “패밀리+직교 orthology”의 이중증거.

---

# 3) 입력이 뉴클레오타이드일 때의 전처리

* \*\*만약 ‘유전자 단위 CDS fasta’\*\*라면: **transeq/EMBOSS**로 6-frame 검사보단 *주어진 프레임* 번역 유지(+ 내부 stop 검사).
* \*\*만약 ‘유전체/contig’\*\*라면: **Prodigal**(`-p meta`)로 CDS 예측 → 단백질 생성.
* 결과: `proteins/*.faa`

---

# 4) 1차 스크리닝 — HMMER 우선, MMseqs2/DIAMOND 보조

**왜:** HMMER로 *민감도* 확보 → MMseqs2로 *속도/대역폭* 확보.

1. **hmmsearch**

   * TDC.pHMM, SadA.pHMM 각각 실행 → **E-value**, **bit score**, **alignment coverage**(Ref/Query 기준 모두) 수집
   * **컷오프:** pHMM별 ROC에서 정한 bit score(모델-특이), Query/Ref coverage ≥ **0.7–0.8** 권장

2. **MMseqs2 easy-search** (참조 단백질 패널)

   * 상호검증용 **RBH(Reciprocal Best Hit)** 수행
   * **ID ≥35–40%**, **aln cov ≥70–80%**, **E ≤1e-10**를 가이드(보수적으로 설정)
   * GPU가 있으면 **MMseqs2-GPU**로 가속 (5070Ti 사용)

> 해석 원칙: **“HMM pass ∧ RBH pass”= 강한 양성**, “HMM pass ∧ RBH 약함”= 중간, “BLAST류만 pass”= 낮음.

---

# 5) 계통 배치(phylogenetic placement) — 최종 특이성 잠금

**왜:** decarboxylase끼리 *서열 보존이 높아서* 오판단이 남습니다. **참조 계통수에 쿼리를 배치**해 clade-level로 확정합니다.

* **IQ-TREE**/**FastTree**로 *참조* decarboxylase들(TDC/SadA/NEG)을 포함한 **백본 트리** 생성
* **EPA-ng** 또는 **SEPP**로 후보 서열 **placement** → **TDC clade** 또는 **SadA clade**에 고신뢰 배치되는지 확인(LWR, pplacer-like 지표).
* 트리 상의 **분지 길이/지지도**와 함께 보고 → “왜 TDC/SadA로 본다”의 **설명가능성** 확보.

---

# 6) 활성부위/길이·도메인 규칙 체크 — 기계적 근거 보강

**왜:** 기능 손실 변이, 비기능 동족체를 걸러야 합니다.

* **PLP-결합 Lys**(Schiff base) 등 **보존 잔기 패턴** 정렬 위치 확인(참조 MSA 좌표 기반)
* **길이 대역**: TDC 대략 \~600 aa, SadA 대략 \~450–500 aa (예외 허용폭 ±10–15%)
* **도메인 서명**: InterProScan/Pfam으로 **Group II PLP-DC** 서명 확인 + NEG 도메인과의 **동시 히트 방지**
* 결과를 \*\*Active-site integrity score(0–1)\*\*로 수치화 → 종합 점수에 가중치로 반영

---

# 7) (가능시) 유전자 이웃/오페론 문맥

**왜:** *tyrDC–tyrP* 같은 **공존 시그널**은 기능적 개연성을 키웁니다.

* **gff/contig**가 있다면: **clinker/clustermap.js**로 주변 **아미노산 퍼미에이스/조절인자** 공존 확인
* 문맥 점수(+0/ +1) 부여 → *증거 가중치*

---

# 8) 증거 가중치 스코어링 (최종 판정)

각 후보 단백질에 다음을 합산 (예시 가중치, 조정 가능):

| 증거     | 기준                        | 점수              |
| ------ | ------------------------- | --------------- |
| HMM 점수 | 모델 컷오프 초과                 | +3              |
| RBH    | ID≥40% & cov≥80%          | +2 (완화: +1)     |
| 계통 배치  | TDC/SadA clade에 LWR 상위 5% | +3 (상위 20%: +2) |
| 활성부위   | PLP-Lys 및 핵심 잔기 보존        | +2              |
| 길이/도메인 | 길이대·Pfam 서명 부합            | +1              |
| 이웃성    | 관련 운반체/조절 유전자 공존          | +1              |

* **≥8점**: High confidence
* **5–7점**: Medium
* **≤4점**: Low (추가 검증 권장)

---

# 9) 산출물(심사위원 친화형)

* **per-gene hit 테이블**: 파일/유전자ID, TDC/SadA 구분, HMM bit, E, cov, RBH ID/cov/E, clade, 활성부위 보존(yes/no), 총점
* **strain × pathway 히트맵**: 균주별 TDC/SadA 보유 여부·점수
* **annotated phylogeny PDF**: 후보에 색/아이콘 표시
* **QC 리포트**: 컷오프 ROC, 음성군 혼입률, 경계 사례 목록
* (선택) **FASTA of positives**: 후속 클로닝/합성용

---

# 10) 실행 가이드 (하드웨어 반영)

* 병렬화: HMMER/MMseqs2/InterProScan/MAFFT/IQ-TREE 모두 **–threads 16**
* **MMseqs2-GPU**: 5070Ti로 검색 가속 (DB가 크면 효과 큼)
* 워크플로: **Snakemake/Nextflow**로 파이프라인화 → 재현성·체크포인트·로그 확보
* 컨테이너: `Dockerfile`/`environment.yml` 제공(심사/이식성↑)

---

## 요약 — 왜 이 방식이 “가장 타당”한가

1. **pHMM 기반(계열특이) + RBH(직교 검증)** → *민감도와 특이성*을 동시에 잡음
2. **계통 배치**로 “TDC vs GAD 등” *최종 구분*을 계통학적으로 설명 가능
3. **활성부위/도메인/길이/오페론**까지 체크 → *기능적 개연성* 확보
4. **가중치 스코어링** → 판정 기준이 **명시적**이고 재현 가능
5. **컨테이너·워크플로** → 제3자 재현성(심사위원 설득력)

---

# 카테콜아민 생합성 유전자 검출 파이프라인 설계도 (v0.1)

> **목표**: 장내 미생물 유래 뉴클레오타이드 FASTA (\~2,400개, `./Fasta_files/`)에서 **TDC/SadA 등 방향족 아미노산 탈탄산효소(AAAD)** 보유 여부를 다중 증거 기반(pHMM+RBH+계통배치+도메인/활성부위)으로 판정하여, **균주 후보** 및 **근거 지표**를 산출.

---

## 1. 입·출력 I/O 규격

**입력**

* 디렉토리: `./Fasta_files/` (확장자: `.fa`, `.fasta`, `.fna`)
* 내용: 뉴클레오타이드 서열 (유전자 단일 CDS 또는 contig/유전체 혼재 가능)

**출력 (핵심 파일)**

* `results/hits/per_gene_hits.tsv`: 후보 유전자별 hit (HMM/BLAST류/도메인/배치 정보)
* `results/summary/per_strain_scores.tsv`: 균주별 종합 점수 및 판정 (High/Medium/Low)
* `results/phylogeny/aaad_backbone_tree.nwk`, `results/phylogeny/placements.tsv` (계통수 및 배치)
* `results/qc/roc_curves.pdf`, `results/qc/thresholds.tsv` (모델별 컷오프 검증)
* `results/vis/strain_pathway_heatmap.png` (균주×경로 히트맵)
* `results/fasta/positives.faa` (양성 후보 단백질 서열)

---

## 2. 파이프라인 개요 (Mermaid)

```mermaid
flowchart TD
    A[입력 FASTA (nucleotide)] --> B{입력 유형 판별}
    B -->|CDS 유전자| C[프레임 유지 번역]
    B -->|유전체/contig| D[Prodigal -p meta (CDS 예측)]
    C --> E[단백질 서열 세트]
    D --> E[단백질 서열 세트]

    subgraph REF[참조 구축]
      R1[Positive 패널 (TDC/SadA)] --> R3[MSA (MAFFT L-INS-i)] --> R4[trimAl] --> R5[hmmbuild]
      R2[Negative 패널 (GAD/HDC/LDC 등)] --> R3
    end

    E --> H1[hmmsearch (TDC/SadA pHMM)]
    E --> H2[MMseqs2 GPU easy-search (RBH)]
    H1 --> J[Hit 통합]
    H2 --> J[Hit 통합]

    J --> K[InterPro/Pfam 도메인 & 길이 규칙]
    J --> L[백본 계통수(IQ-TREE/FastTree)]
    L --> M[EPA-ng/SEPP 배치]

    K --> N[가중치 스코어링]
    M --> N
    J --> N

    N --> O[per_gene_hits.tsv]
    N --> P[per_strain_scores.tsv]
    N --> Q[히트맵/트리/리포트]
```

---

## 3. 툴체인 & 자원 설정 (16C/128GB/Rtx 5070Ti + 3070Ti)

* **Gene prediction/translation**: Prodigal(유전체), EMBOSS transeq/맞춤 번역(CDS)
* **MSA/HMM**: MAFFT(L-INS-i), trimAl, HMMER (hmmbuild/hmmsearch)
* **동질성 검색**: MMseqs2 (GPU 가속, 5070Ti 우선), DIAMOND(옵션)
* **도메인**: InterProScan(풀) 또는 Pfam-A hmmsearch(경량)
* **계통**: IQ-TREE2(Fast) 또는 FastTree + **EPA-ng/SEPP** 배치
* **오페론/이웃성(선택)**: clinker, clustermap.js (contig 메타데이터 있을 때)
* **워크플로**: Snakemake (권장) 또는 Nextflow DSL2
* **컨테이너/환경**: mamba/conda `environment.yml` + 선택적 Docker/Singularity

**스레드/메모리 권장**

* `--threads 16`, HMMER/MAFFT/IQ-TREE 병렬화
* MMseqs2 GPU: 5070Ti 사용, 3070Ti는 보조 큐로 분할 가능
* InterProScan은 메모리 부하가 큼 → Pfam-only 경량 모드 제공

---

## 4. 참조 세트 구성 (양성·음성 동시)

* **양성(positives)**: TDC (Enterococcus/Lactobacillus 다종), SadA-like (Staphylococcus/Clostridia 등)
* **음성(negatives)**: GAD(H+/PLP), HDC, LDC 등 decarboxylase 패밀리
* **중복제거**: `seqkit rmdup -s` 또는 MMseqs2 `linclust --min-seq-id 0.9`
* **정렬**: `mafft --maxiterate 1000 --localpair` (또는 L-INS-i)
* **트리밍**: `trimal -automated1`
* **pHMM 생성**: `hmmbuild TDC.hmm tdc.aln`, `hmmbuild SadA.hmm sada.aln`

> 선택: 기존 구축된 `refs/*.faa`, `refs/*.hmm`를 제공하거나, Snakemake 모듈로 자동 구축.

---

## 5. 컷오프(ROC) 보정

* positives/negatives 모두에 대해 `hmmsearch` 수행 → bit score 분포 획득
* **ROC 곡선** 그려 `TPR≥0.95 & FPR≤0.01` 지점의 **모델별 bit cutoff** 도출
* `results/qc/thresholds.tsv` 저장(모델명, bitscore\_cut, evalue\_hint, ref/qry coverage 제한)

---

## 6. 입력 전처리

* **CDS 유전자 FASTA**: 기록된 프레임 유지 번역(내부 stop 검증) → 단백질
* **유전체/contig**: `prodigal -p meta -i sample.fna -a sample.faa -d sample.cds.fna`
* 산출: `proteins/{sample}.faa`

---

## 7. 1차 스크리닝 (HMMER 우선, MMseqs2 보조)

**HMMER**

* `hmmsearch --cpu 16 --tblout sample.TDC.tbl TDC.hmm proteins/sample.faa`
* 커버리지 계산(Ref/Query 기준 모두) 및 필터: `bit ≥ cutoff`, `cov ≥ 0.75–0.8`

**MMseqs2 (GPU)**

* 참조 단백질 패널 DB 생성 → `mmseqs easy-search proteins/sample.faa refs/aaad_refdb results/mmseqs.tsv tmp --gpu 1`
* Reciprocal Best Hit(RBH): 후보 ↔ 참조 상호 최고 hit 확인

**태깅(선택)**

* KofamScan/eggNOG-mapper로 ortholog 태그 추가 (네트 제한 시 스킵 가능)

---

## 8. 계통학적 배치

* **백본 트리**: positives(다양종) + 대표 negatives 혼합으로 IQ-TREE2 또는 FastTree 구성
* **배치**: EPA-ng 또는 SEPP로 candidate 서열을 백본에 placement → **LWR 등 지표** 산출
* 산출: `aaad_backbone_tree.nwk`, `placements.tsv`, 배치 시각화 PDF/PNG

---

## 9. 기능 보강 규칙 (활성부위/도메인/길이/이웃성)

* **활성부위**: PLP-결합 Lys 등 보존 잔기 좌표 체크 (MSA 기준에 매핑)
* **도메인**: Pfam/InterPro 서명 존재 확인 (Group II decarboxylases 등)
* **길이**: TDC \~600aa, SadA \~450–500aa (허용 ±10–15%)
* **이웃성(옵션)**: tyrP 등 운반체/조절자 공존 시 +1 가점

---

## 10. 가중치 스코어링 (판정 기준)

| 증거     | 기준                            | 점수              |
| ------ | ----------------------------- | --------------- |
| HMM    | 모델 cutoff 초과 (bit) & cov≥0.75 | +3              |
| RBH    | ID≥40% & cov≥80% & E≤1e-10    | +2 (완화: +1)     |
| 배치     | TDC/SadA clade에 LWR 상위 5%     | +3 (상위 20%: +2) |
| 활성부위   | 핵심 잔기 보존                      | +2              |
| 도메인/길이 | 규칙 부합                         | +1              |
| 이웃성    | 관련 유전자 공존                     | +1              |

**총점 판정**: ≥8 High / 5–7 Medium / ≤4 Low

---

## 11. 산출물 명세

* `per_gene_hits.tsv`: sample, gene\_id, method(HMM/MMseqs/RBH), model, bit/e/id/cov, domain, active\_site, placement, score
* `per_strain_scores.tsv`: strain/species, n\_hits(TDC/SadA), best\_score, verdict(H/M/L), notes
* `positives.faa`: High/Medium 후보 단백질 FASTA
* `qc/`: thresholds.tsv, roc\_curves.pdf, boundary\_cases.tsv (경계 사례 목록)
* `vis/`: heatmap.png, tree\_annot.pdf

---

## 12. 프로젝트 디렉토리 스켈레톤

```
project/
├─ config/
│  ├─ config.yaml               # 경로/스레드/컷오프/모드
│  └─ thresholds.seed.tsv       # 초기 컷오프(없으면 자동 학습)
├─ refs/
│  ├─ positives/*.faa           # TDC/SadA 대표 패널
│  ├─ negatives/*.faa           # GAD/HDC/LDC 등
│  └─ hmm/*.hmm                 # TDC.hmm, SadA.hmm …
├─ workflows/
│  ├─ Snakefile
│  └─ rules/*.smk               # 모듈별 규칙
├─ envs/
│  ├─ base.yml                  # mamba 환경(필수 도구)
│  └─ interpro.yml              # 선택: InterProScan
├─ scripts/
│  ├─ translate_orf.py
│  ├─ parse_hmm_tbl.py
│  ├─ rbh_mmseqs.py
│  ├─ score_hits.py
│  └─ plot_roc.py
├─ results/                     # 자동 생성
└─ logs/                        # 실행 로그
```

---

## 13. Snakemake 파이프라인 (핵심 규칙 스케치)

```python
# workflows/Snakefile (요약 스케치)
configfile: "config/config.yaml"
SAMPLES = [f.replace('.fna','').replace('.fa','').replace('.fasta','') for f in glob_wildcards("Fasta_files/{sample}.{{fa,fna,fasta}}" ).sample]
THREADS = config.get("threads", 16)

rule all:
    input:
        expand("results/hits/{sample}.hits.tsv", sample=SAMPLES),
        "results/summary/per_strain_scores.tsv",
        "results/qc/roc_curves.pdf",
        "results/phylogeny/placements.tsv"

rule predict_or_translate:
    input: "Fasta_files/{sample}.fasta"
    output: "proteins/{sample}.faa"
    threads: THREADS
    conda: "../envs/base.yml"
    shell:
        """
        if grep -qi ">.*cds\|gene" {input}; then \
            python scripts/translate_orf.py --in {input} --out {output}; \
        else \
            prodigal -p meta -i {input} -a {output} -q; \
        fi
        """

rule hmmsearch:
    input: faa="proteins/{sample}.faa", hmm="refs/hmm/TDC.hmm"
    output: tbl="results/hmm/{sample}.TDC.tbl"
    threads: THREADS
    shell:
        "hmmsearch --cpu {threads} --tblout {output.tbl} {input.hmm} {input.faa} > logs/{wildcards.sample}.tdc.hmm.log 2>&1"

rule mmseqs_rbh:
    input: faa="proteins/{sample}.faa", db="refs/mmseqs/aaad_refdb"
    output: tsv="results/mmseqs/{sample}.tsv"
    threads: THREADS
    shell:
        "mmseqs easy-search {input.faa} {input.db} {output.tsv} tmp --gpu 1 --threads {threads} > logs/{wildcards.sample}.mmseqs.log 2>&1"

rule integrate_hits:
    input:
        hmm_tdc="results/hmm/{sample}.TDC.tbl",
        mmseqs="results/mmseqs/{sample}.tsv"
    output: "results/hits/{sample}.hits.tsv"
    conda: "../envs/base.yml"
    shell:
        "python scripts/parse_hmm_tbl.py --hmm {input.hmm_tdc} --mmseqs {input.mmseqs} --out {output} --thresholds config/thresholds.seed.tsv"

rule score_and_summarize:
    input: expand("results/hits/{sample}.hits.tsv", sample=SAMPLES)
    output:
        scores="results/summary/per_strain_scores.tsv",
        heatmap="results/vis/strain_pathway_heatmap.png"
    conda: "../envs/base.yml"
    shell:
        "python scripts/score_hits.py --hits results/hits --out_scores {output.scores} --out_heatmap {output.heatmap}"
```

> 실제 규칙은 TDC/SadA/NEG 모델 병렬 처리, ROC 학습, 계통 배치(EPA-ng) 등을 포함하도록 확장합니다.

---

## 14. `config/config.yaml` 예시

```yaml
threads: 16
use_gpu: true
min_cov: 0.75
rbh:
  min_id: 0.40
  min_cov: 0.80
  max_evalue: 1e-10
scoring:
  hmm: 3
  rbh: 2
  placement: 3
  active_site: 2
  domain_length: 1
  synteny: 1
cutoffs:
  TDC_bitscore: auto   # ROC로 산출, 없으면 seed 사용
  SadA_bitscore: auto
paths:
  input_dir: ./Fasta_files
  refs_dir: ./refs
  results_dir: ./results
```

---

## 15. 실행 예시

```bash
# 1) 환경 준비
mamba env create -f envs/base.yml
mamba activate aaad

# 2) (선택) 참조 구축
snakemake -j 16 refs/hmm/TDC.hmm refs/hmm/SadA.hmm

# 3) 본 분석 실행
snakemake -j 16 --use-conda

# 4) 결과 확인
open results/summary/per_strain_scores.tsv
open results/vis/strain_pathway_heatmap.png
```

---

## 16. 성능 및 스케일링 노트

* 2,400 샘플 × 평균 단백질 2–3k 서열 가정 시, HMMER는 16스레드에서 수 시간대; MMseqs2 GPU로 RBH는 수십 분\~수 시간.
* InterProScan(풀)은 대기열 큼 → **Pfam-only** 경량 모드 권장. 필요 시 하이콘피던스 후보에만 후처리.
* 계통 배치(EPA-ng): 백본 사이즈 1–2k 시 합리적. 더 크면 FastTree + subsampling.

---

## 17. QC/감사 추적성

* 모든 규칙에 `logs/*.log`와 `software_versions.txt` 기록
* `results/qc/roc_curves.pdf`와 `thresholds.tsv`로 컷오프 근거 제시
* `boundary_cases.tsv`에 경계 사례(낮은 여유·상충 증거) 명시

---

## 18. 확장 포인트

* **대사 네이버후드 점수** 추가 (tyrP/antiporter 공존)
* **메타지놈 모드**: 중복 제거(clstr) 후 대표 서열만 검색
* **Nextflow DSL2 포팅**: 모듈화·클라우드 배치 대응

---

### 결론

이 설계는 \*\*pHMM(민감도) + RBH(특이성) + 계통 배치(설명가능성)\*\*를 중심으로, **활성부위/도메인/길이/이웃성**을 결합한 다중 증거 판정 체계를 제공합니다. 재현성(워크플로/환경 고정)과 심사 친화적 산출물(QC·ROC·트리·히트맵)을 기본 내장하여, 제3자가 동일 하드웨어에서 곧바로 검증·확장할 수 있습니다.

---

