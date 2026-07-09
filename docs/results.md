# CPTAC ccRCC Drug-Target Discovery — Results

## In one paragraph

Starting from 11,710 proteins measured in CPTAC ccRCC, we kept those over-produced in tumour at
both the protein and RNA level (→ 310), narrowed to the ones a drug can physically reach on the
cell surface or in secretions (→ 101), removed immune-cell decoys, and took the **top 12** as the
candidate list. Each of the 12 was scored 0–5 on four qualities — overexpression, drug-reachability
(tractability), safety and selectivity — from public databases (UniProt, DepMap, gnomAD, HPA);
those four combine, by fixed weights, into a single **0–100 composite** that ranks them. The whole
pipeline was run **blind to drug and trial information**. Only at the end did we cross-reference the
ranking against clinical-trial data — and the two ccRCC antigens that *already* have a direct
therapeutic (CA9, CD70) both landed in the top 3.

## Headline

A single expression-driven method, run **blind to all drug and clinical-trial information**,
ranks 12 surface/secreted candidates in ccRCC. When drug status is revealed afterwards, the two
ccRCC antigens that already have a direct therapeutic — **CA9** and **CD70** — both land in the
**top 3 of 12** (CD70 #1, CA9 #3). Under a hypergeometric null this is unlikely by chance
(**p = 0.0455**); the direct-drug antigens also score higher than the rest as a group
(one-sided Mann–Whitney **p = 0.0303**). The method recovers clinically-validated biology
without ever being told what is drugged.

## Ranked candidates (drug-blind composite)

| Rank | Gene | Protein log2FC | RNA log2FC | Assoc. | Tract. | Safety | Select. | Composite | Tier | Direct drug |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| 1 | **CD70** | 2.36 | 4.27 | 4.17 | 3.0 | 4.7 | 4.3 | **79.7** | T1 fast-follow | anti-CD70 CAR-T (ALLO-316) |
| 2 | HILPDA | 2.74 | 4.83 | 4.79 | 1.5 | 4.7 | 4.3 | 75.5 | T2 discovery | — |
| 3 | **CA9** | 2.23 | 6.02 | 3.96 | 3.0 | 4.7 | 2.8 | **73.3** | T1 fast-follow | girentuximab radioconjugate |
| 4 | COL23A1 | 2.46 | 3.95 | 3.96 | 1.5 | 4.7 | 4.3 | 70.0 | T2 discovery | — |
| 5 | SLC16A3 | 2.24 | 3.30 | 3.12 | 3.0 | 4.7 | 2.0 | 65.0 | T3 watch | — |
| 6 | POSTN | 2.35 | 1.03 | 2.08 | 3.0 | 4.2 | 3.5 | 60.9 | T3 watch | — |
| 7 | NPTX2 | 1.73 | 4.03 | 2.08 | 1.5 | 4.7 | 4.3 | 57.4 | T3 watch | — |
| 8 | SCARB1 | 1.80 | 2.89 | 2.08 | 3.0 | 2.7 | 3.5 | 54.2 | T4 deprioritise | — |
| 9 | SLC2A1 | 2.00 | 2.30 | 2.08 | 3.0 | 3.4 | 2.0 | 52.3 | T4 deprioritise | — |
| 10 | HAPLN1 | 2.10 | 1.44 | 1.67 | 1.5 | 4.2 | 3.5 | 49.8 | T4 deprioritise | — |
| 11 | HAVCR1 | 1.80 | 1.86 | 1.46 | 1.5 | 3.2 | 5.0 | 49.0 | T4 deprioritise | — |
| 12 | SLC2A3 | 1.73 | 1.69 | 1.04 | 3.0 | 4.2 | 1.0 | 45.6 | T4 deprioritise | — |

Full table with all evidence columns: `results/scorecards/scorecard_clean.csv`.

## Blind recovery of drugged antigens

- **CA9** is a fully mechanical recovery: it survived every selection filter with no manual
  intervention and carries the strongest RNA overexpression in the set (log2FC = 6.02, RNA rank
  1 of 12). It is #3 on the composite.
- **CD70** is retained on its own data — it is tumour-intrinsic by the immune filter
  (ρ = 0.118, not significant), not a leukocyte marker — and tops the composite at 79.7.
- Both direct-drug antigens fall in the top 3 → hypergeometric **p = 0.0455**. (Both-in-top-2 is
  not met: HILPDA is #2.)
- Composite scores of the direct-drug antigens vs the rest: one-sided Mann–Whitney **U = 19,
  p = 0.0303**.

Statistics: `results/scorecards/recovery_stats.json`.

## Immune filter

The data-driven immune filter cleanly separates leukocyte markers (ρ 0.48–0.83) from tumour-cell
surface candidates (all < 0.30), with tumour-intrinsic controls at ρ ≤ 0. This removes the main
false-positive mode of bulk-tumour overexpression screens without a hand-curated blocklist. See
`immune_filter.png` and `results/enrichment` for the marker-level table.

## Robustness to rubric weights

Across 2000 random weightings of the four dimensions, **CD70 is a robust #1** (P(rank 1) = 0.69,
mean rank 1.3). CA9 has mean rank 3.7 (SD 1.8). The ordering does not depend on the specific
weight vector chosen (`weight_sensitivity.png`, `results/enrichment/weight_sensitivity.csv`).

## Mutation drivers (context, unscored)

| Driver | Freq | Consequence | Adjacent clinical approach |
|---|---:|---|---|
| VHL | 85 % | HIF-2α/VEGF axis stabilised | HIF-2α inhibitor (belzutifan) |
| PBRM1 | 40 % | SWI/SNF chromatin remodelling | — (IO-response association) |
| SETD2 | 12 % | H3K36me3 loss | WEE1 inhibitor (synthetic lethal) |
| BAP1 | 10 % | nuclear deubiquitinase loss | PARP inhibitor (synthetic lethal) |
| KDM5C | 7 % | H3K4 demethylase loss | — |

**0 of 5 drivers appear in the expression-derived surface/secreted target set.** ccRCC's
genetic drivers are loss-of-function tumour suppressors — not directly druggable as surface
antigens — so the expression method and the mutation landscape point at disjoint target spaces.
Both are actionable, by different modalities.

## Discussion

**Novelty / crowding.** The method does not score how novel a target is, because novelty is a
property of the existing drug landscape — exactly what was held out. Read post-hoc: CA9 and CD70
are the crowded, validated end (they anchor the recovery test); HILPDA, COL23A1 and SLC16A3 are
the less-explored high-composite candidates worth a closer look.

**HILPDA caveat.** HILPDA ranks #2 on expression, but its dominant biology is intracellular
lipid-droplet regulation; UniProt annotates it as lipid-droplet / secreted / membrane. Its
surface-drug tractability is genuinely low (intrinsic tractability 1.5), which the rubric already
reflects — its #2 rank is driven by overexpression magnitude, not by being an easy surface
target. Treat it as a biology lead, not a ready-made antigen.

**Grade/stage association (supporting).** Among the candidates, only HAVCR1 and SLC16A3 track
tumour grade, and both only marginally (BH-FDR ≈ 0.047, ρ ≈ 0.27, ~7 % of variance). CA9 and CD70
do not track grade. This is a weak signal and is not overstated.

**Limitations.** Trial existence is not efficacy — the anchor trials are small, early-phase
(n ≈ 18–120). The recovery p-values are nominal and modest (0.03–0.05) on a 12-gene list. The
method finds tumour-vs-normal overexpressed surface proteins; it is blind to targets that are not
overexpressed (e.g. the driver suppressors) and to intracellular or fusion-driven targets.
