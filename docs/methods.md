# CPTAC ccRCC Drug-Target Discovery — Methods

**Dataset:** CPTAC clear cell renal cell carcinoma (ccRCC)
discovery cohort — Clark DJ *et al.*, *Cell* 2019, DOI
[10.1016/j.cell.2019.10.007](https://doi.org/10.1016/j.cell.2019.10.007)
**Disease ontology:** MONDO_0005005 · **Organism:** *Homo sapiens* (taxon 9606)

## In plain terms

1. **Find what the tumour makes too much of.** Using the CPTAC ccRCC data, we find genes whose
   **protein** is much higher in tumour than in normal kidney, and keep only those where the
   **RNA** agrees. (Protein is the lead evidence; RNA is corroboration.)
2. **Keep only the reachable ones.** A drug like an antibody, ADC, CAR-T or radioligand can only
   hit a protein on the cell surface or one that is secreted — so we discard everything that lives
   inside the cell.
3. **Throw out immune-cell decoys.** In bulk tumour data, some "tumour-high" genes are really just
   markers of infiltrating immune cells. We detect and remove those from the data itself.
4. This leaves a shortlist of **12 candidates** (this is *n* — see the funnel below).
5. **Score each candidate 0–5 on four things** — how strongly it is overexpressed, how reachable
   it is by a drug, how safe the target looks, and how selective it is — using public annotation
   databases (UniProt, DepMap, gnomAD, HPA). Each score's exact calculation and source is in
   Step 5.
6. **Combine the four scores into one 0–100 composite** using fixed weights, and rank the 12.
7. **Only at the very end** do we look up drug/clinical-trial status — information deliberately
   kept out of every step above — to ask: did this blind method rank the genes that are *already*
   drug targets near the top? (It did: both land in the top 3.)

## Overview

This is a **single, drug-blind, expression-driven discovery method** for surface/secreted
drug-target candidates in ccRCC. It takes as input only the CPTAC ccRCC proteogenomic matrices
and public gene-level annotation (subcellular localisation, paralog structure, genetic-dependency
and constraint metrics). It never consumes drug or clinical-trial information during selection or
scoring. Drug/trial status is brought in **only after the ranking is fixed**, as held-out
annotation, to test what the blind method recovered.

### Candidate funnel (how *n* shrinks to 12)

| Stage | Filter | Genes remaining |
|---|---|---:|
| 0 | Proteins quantified in both tumour and normal proteomes | 11,710 |
| 1 | Protein-overexpressed (log2FC ≥ 1 **and** BH-FDR < 0.05) | 364 |
| 2 | **+** RNA-concordant (RNA log2FC > 0) → candidate pool | **310** |
| 3 | **+** surface/secreted localisation (drug-reachable) | **101** |
| 4 | **+** immune-decoy filter, then top 12 by protein fold-change | **12** |

So the scored candidate list is **n = 12**. Steps 1–4 below define each filter; Step 5 scores the
12 survivors.

## Data

Processed CPTAC ccRCC matrices:

- `HS_CPTAC_CCRCC_proteome_Tumor.cct` — protein abundance, 110 tumours
- `HS_CPTAC_CCRCC_proteome_Normal.cct` — protein abundance, 84 normal adjacent tissues
- `HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Tumor.cct` — RNA-seq log2 FPKM, tumours
- `HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Normal.cct` — RNA-seq log2 FPKM, normals
- `HS_CPTAC_CCRCC_CLI.tsi` — clinical annotation

## Step 1 — Differential overexpression (primary, protein)

For every gene present in both tumour and normal proteomes we compute tumour-vs-normal log2
fold change and a two-sided Welch *t*-test, correcting across all genes with the
Benjamini–Hochberg procedure. **Primary selection is protein-level**: a gene must show
`protein log2FC ≥ 1` and `BH-FDR < 0.05` to enter the candidate pool. Protein is the primary
filter because it is the modality a surface-targeting therapeutic actually engages.

## Step 2 — RNA concordance (directional corroboration)

The protein-selected pool is required to be **directionally concordant** at the RNA level:
`RNA log2FC > 0` (tumour-up). This is a corroboration filter, not a second independent
significance test — it removes genes whose protein signal has no transcriptional support,
without demanding dual significance (which would over-penalise post-transcriptionally regulated
surface proteins).

## Step 3 — Surface/secreted restriction

Candidates are restricted to gene products plausibly accessible to an antibody, ADC, CAR-T or
radioligand: cell-surface (single-/multi-pass membrane, GPI-anchored), broad plasma-membrane
transporters, or secreted/matrix proteins. Localisation is taken from UniProt subcellular-location
annotation.

## Step 4 — Immune filter (data-driven)

A recurring failure mode in bulk-tumour overexpression screens is that "tumour-up" genes are in
fact markers of infiltrating leukocytes, not of tumour cells. We identify and remove these
directly from the expression data:

1. Build a 16-marker leukocyte signature (pan-immune and lineage markers: *PTPRC, CD3D, CD3E,
   CD2, CD8A, CD4, MS4A1, CD79A, CD79B, CD68, CD14, LYZ, ITGAM, CSF1R, CD163, NKG7* — of 20
   attempted; markers absent from the matrix are dropped).
2. z-score each marker across the 110 tumours and average to a per-tumour immune score.
3. For each candidate compute the Spearman correlation (ρ) of its tumour abundance with the
   immune score.
4. **Exclude any candidate with ρ ≥ 0.40.**

The threshold sits in a clean gap in the data (`immune_filter.csv`, `immune_filter.png`): known
leukocyte markers fall at ρ 0.48–0.83, every retained surface candidate falls below 0.30, and
tumour-intrinsic control genes (*NNMT, VIM, VEGFA, ENO2, ANGPTL4*) fall at ρ ≤ 0. *CD70* is
tumour-intrinsic by this measure (ρ = 0.118, not significant) and is retained.

## Step 5 — Scoring rubric (four dimensions, drug-blind)

Each of the 12 surviving candidates is scored on four 0–5 dimensions. **No dimension uses drug
approval, clinical-trial status, or Open Targets overall/known-drug/literature evidence** — those
carry publication and drug-development bias and would make the method partly a lookup of what is
already drugged. The exact calculation and data source for each dimension:

**1. Overexpression / association** — *source: CPTAC only (this study's own data).*
Rank the 12 candidates by protein log2FC and, separately, by RNA log2FC; average the two rank
positions and rescale to 0–5:
`assoc = 5 × (rank(protein log2FC) + rank(RNA log2FC)) / (2 × 12)`.
Reading the formula piece by piece: each gene gets a **rank** from 1 (lowest overexpression) to
12 (highest) on each of the two measures — we use ranks, not the raw log2FC values, so that one
gene with an extreme fold-change cannot dominate the dimension. The two ranks are **summed**, which
rewards a candidate that is high on *both* protein and RNA (concordance) rather than on one alone.
The denominator `2 × 12 = 24` is the largest the sum can be (rank 12 on both measures), so the
ratio is a 0–1 position; the leading `× 5` puts it on the same 0–5 scale as the other three
dimensions. A candidate top-ranked on both protein and RNA approaches 5; the weakest approaches 0.
Association is defined **entirely** from CPTAC overexpression — no external evidence. This is the
one dimension that uses percentile-style (rank-based) scaling, which is appropriate here because
the raw inputs are directly comparable log2 fold-changes on a common scale.

**2. Tractability** — *source: UniProt subcellular localisation.*
Intrinsic reachability by a surface-directed drug. Cell-surface / GPI-anchored / broad
plasma-membrane transporter (antibody-, ADC-, CAR-accessible) → **3.0**; secreted-only,
lipid-droplet or matrix-embedded product (harder to hit cleanly) → **1.5**. No clinical-status
or "how advanced is the drug" bucket enters this score.

**3. Safety** — *source: DepMap, gnomAD, Human Protein Atlas.* Start at 3.0 and adjust for how
tolerable inhibiting the gene looks:
- DepMap mean gene-effect (essentiality across cell lines): strongly essential (≤ −1.0) **−1.5**;
  moderately essential (≤ −0.5) **−0.7**; non-essential (> −0.5) **+0.7**.
- gnomAD LOEUF (loss-of-function constraint): highly constrained (< 0.35) **−0.8**; tolerant
  (> 1.0) **+0.5**.
- HPA vital-organ enrichment: enriched in a vital organ **−1.0**, otherwise **+0.5**.
- Each annotated safety liability: **−0.5** (capped at 2). Result clamped to 0–5.

**4. Selectivity** — *source: paralog analysis + Human Protein Atlas RNA tissue-specificity.*
Selectivity asks how likely a drug aimed at this protein is to hit something it was not meant to.
Two independent off-target risks are combined. Start at 3.5 and adjust:
- **Close paralogs** (other human proteins with similar sequence — structural cousins a binder
  could cross-react with). Count of paralogs at ≥ 30 % identity: none **+0.8** (a monoclonal or
  CAR binder is unlikely to have a near-relative to stray onto); two to four **−0.8**; five or more
  **−1.5** (a large, tightly related family is hard to bind cleanly).
- **Maximum paralog identity ≥ 70 %** **−1.0** — one very close homolog is difficult to
  distinguish structurally even if the rest of the family is distant.
- **RNA tissue pattern** (HPA): tissue-enriched **+0.8** — expression concentrated in a few
  tissues means less chance of hitting the same antigen on healthy tissue elsewhere;
  low tissue-specificity (broadly expressed) **−0.8** — the antigen is present across many normal
  tissues, raising on-target/off-tumour toxicity risk. Clamped to 0–5.

Worked example: **CA9** = 3.5 − 1.5 (9 paralogs in the carbonic-anhydrase family) + 0.8
(tissue-enriched) = **2.8**. The method is deliberately harsh here — CA9 is a validated drug
target, but its large paralog family is a genuine, well-known selectivity challenge, and the
blind score reflects that without knowing CA9 is drugged.

**Were the safety and selectivity inputs normalised?** No — and this is deliberate. The four
inputs to each (essentiality, genetic constraint, vital-organ expression, liability count for
safety; paralog count/identity and tissue-specificity for selectivity) are **not** z-scored,
min-maxed, or otherwise rescaled against the other 11 candidates. Each raw value is compared to a
**fixed, absolute threshold** (e.g. LOEUF < 0.35, paralog identity ≥ 70 %) and converted to a
hand-set point adjustment; the running total is then clamped to 0–5. Two reasons: (i) with only
n = 12 candidates, continuous batch-relative normalisation would be unstable — a single outlier
would rescale everyone — whereas a threshold keeps each gene's score meaning the same regardless
of who else is in the set; (ii) the thresholds encode expert judgement about what counts as
"constrained" or "a close homolog," so they are a chosen rubric, not parameters fitted to the data.
The trade-off is honesty about that: the point values (−1.5, −0.8, +0.8 …) are expert-set, not
learned, and a different rubric would shift absolute scores — which is exactly why Step 9 resamples
the dimension weights to show the *ranking* is robust to that choice. Only the association
dimension uses rank-based (percentile-style) scaling, because there the raw inputs are directly
comparable log2 fold-changes.

**Combining the four into a composite.** Weights are renormalised from a 30/25/20/15 allocation
(association weighted most, selectivity least):

```
association 0.3333   tractability 0.2778   safety 0.2222   selectivity 0.1667
```

`composite (0–100) = 20 × (weighted mean of the four 0–5 scores).`

The 20× factor rescales the 0–5 weighted mean onto 0–100. Full per-gene scores and the raw
annotation values behind them are in `results/scorecards/scorecard_clean.csv` and
`results/evidence/evidence_detail.json`.

A "novelty" dimension (how un-crowded a target is) is **not** scored: it is a statement about
the existing drug landscape, the very information held out of the method. It is addressed
qualitatively in the Discussion.

## Step 6 — Development tiers

- **T1 fast-follow** — a therapeutic already engages the antigen directly (held-out label).
- **T2 discovery** — composite ≥ 70, no direct drug.
- **T3 watch** — composite 55–70.
- **T4 deprioritise** — composite < 55.

## Step 7 — Mutation drivers as context (unscored)

The five canonical ccRCC loss-of-function drivers (*VHL, PBRM1, SETD2, BAP1, KDM5C*) are
reported as an **unscored context panel** with mutation frequency, consequence, and the adjacent
clinical approach where one exists (e.g. HIF-2α inhibition downstream of *VHL*; synthetic-lethal
strategies for *SETD2*/*BAP1*). They are tumour-suppressor lesions, not overexpressed surface
proteins, so they lie structurally outside the surface/secreted target space. Their overlap with
the scored set is 0/5 — a factual divergence, not a competing ranking.

## Step 8 — Blind-recovery test

To ask whether the drug-blind ranking captured actionable biology, we bring in the held-out
annotation. Two ccRCC antigens are the subject of a **direct** therapeutic (a drug that binds the
antigen itself): CA9 (girentuximab radioconjugate, NCT05663710) and CD70 (allogeneic anti-CD70
CAR-T ALLO-316, NCT04696731). We test:

- **Rank recovery** — a hypergeometric test for both direct-drug antigens landing within the
  observed top-*k* of the ranked list.
- **Score separation** — a one-sided Mann–Whitney *U* test comparing the composite scores of the
  direct-drug antigens against the rest.

## Step 9 — Weight-sensitivity (robustness)

Rubric weights are a modelling choice, so we resample them: 2000 Dirichlet draws over the four
dimensions (seed 42), recomputing the composite and ranks each time, and report each gene's mean
rank, rank SD and P(rank = 1).

## Software

Python 3.11 (numpy, pandas, scipy, statsmodels, matplotlib). Full pipeline: `code/pipeline.py`;
environment snapshot: `environment_snapshot.txt`.

## What is deliberately *not* used in selection or scoring

- Drug approval / development status
- Clinical-trial existence
- Open Targets overall association, known-drug or literature evidence
- Any hand-curated gene blocklist

Drug/trial data are held out for the recovery test; the immune exclusion is computed from the
expression data.
