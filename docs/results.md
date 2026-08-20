# CPTAC ccRCC Drug-Target Discovery — Results

## In one paragraph

Starting from 11,710 proteins measured in CPTAC ccRCC, we kept those over-produced in tumour at
both the protein and RNA level (→ 310), narrowed to the ones a drug can physically reach on the
cell surface or in secretions (→ 101), removed immune-cell decoys, and took the **top 12** as the
shortlist. **All 101** were scored 0–5 on four qualities — overexpression, drug-reachability
(tractability), safety and selectivity — from public databases (UniProt, DepMap, gnomAD, HPA);
those four combine, by fixed weights, into a single **0–100 composite** that ranks them. The whole
pipeline was run **blind to drug and trial information**. Only at the end did we cross-reference
the ranking against clinical-trial data. CD70 and CA9 rank #2 and #4 of 101 — but on a properly
specified null the method does not significantly enrich for clinically-validated targets.

## Headline

A single expression-driven method, run **blind to all drug and clinical-trial information**,
scores 101 surface/secreted candidates in ccRCC and shortlists 12. When drug status is revealed
afterwards, **CD70 (#2), CA9 (#4) and ENPP3 (#9) of 101** all carry direct-acting agents — a real
signal at the top of the ranking.

That signal does not survive as a claim about the method. Under a rule applied to all 101
candidates there are **14 positives, not 2**, and against that denominator neither the selection
cascade (hypergeometric **p = 0.52**) nor the composite ranking (Mann–Whitney **p = 0.32**)
significantly enriches for them.

An earlier version of this analysis reported **p = 0.0455**. That statistic permuted only the 12
survivors, so it was conditional on the selection step it was meant to test, and a positive
dropped anywhere in the 11,710 → 12 cascade could not register. It is retained, relabelled, as
`ranking_given_selection_LEGACY`, where it recomputes to p = 0.0152 under the corrected
association score. It became *more* significant while the claim it supported became weaker —
which is the clearest possible demonstration that it was never testing the right thing.
The correction is the main result of this write-up.

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
  intervention and carries the strongest RNA overexpression in the set (log2FC = 6.02). It is
  #4 of 101 on the composite.
- **CD70** is retained on its own data — it is tumour-intrinsic by the immune filter
  (ρ = 0.118, not significant), not a leukocyte marker — and is #2 of 101.
- **Cascade test** (did selection concentrate positives?): 2 of 14 positives retained in the 12,
  against 1.66 expected by chance → hypergeometric **p = 0.52**.
- **Ranking test** (does the composite rank positives highly across all 101?): one-sided
  Mann–Whitney **p = 0.32**; positives at ranks 2, 4, 9, 16, 33, 34, 51, 54, 58, 61, 66, 84, 96, 99.
- **12 of 14 positives are missed** — 2 by the immune filter (CD4, FCGR3A, both correctly
  identified as infiltrate-tracking) and 10 by the top-12 fold-change cut. **ENPP3**, target of
  the ADC AGS-16C3F, scores #9 of 101 and is discarded by that cut: the shortlisting rule, not
  the rubric, is what loses it.

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

**Limitations.** The recovery test is a soft validation, not proof. It checks the drug-blind ranking
against which top hits already have a direct-acting agent in the clinic (the anchors: ¹⁷⁷Lu-girentuximab
for CA9, Phase 1/2; ALLO-316 CAR-T for CD70, Phase 1) - both investigational, not approved ccRCC
therapeutics. A drug reaching trials does not establish efficacy, these are small early-phase studies
(~100–120 patients each), and on a 12-gene list the recovery p-values (0.03–0.05) are borderline.
Of the 10 non-target candidates, 7 have no ccRCC trials at all; the remaining 3 (HILPDA, SLC2A1/GLUT1,
HAVCR1/KIM-1) appear in trials only as incidental biomarkers, never as the drug target. The
method finds tumour-vs-normal overexpressed surface proteins; it is blind to targets that are not
overexpressed (e.g. the driver suppressors) and to intracellular or fusion-driven targets.
