# Drug-target scoring rubric — CPTAC ccRCC (Clark et al. 2019)

The candidate set is a single pool of **tumour-overexpressed cell-surface / secreted
proteins** derived de novo from the CPTAC ccRCC proteome + transcriptome. Each candidate
is scored on **four 0–5 dimensions**, all computed from CPTAC expression and public
gene-level annotation. **No dimension uses drug approval, clinical-trial status, or Open
Targets overall / known-drug / literature evidence** — those carry publication and
drug-development bias and would turn the score into a partial lookup of what is already
drugged. Drug/trial status is held out entirely and revealed only afterwards, to test
what the blind ranking recovered.

The five recurrently mutated ccRCC drivers (*VHL, PBRM1, BAP1, KDM5C, SETD2*) are **not
scored**. They are loss-of-function tumour suppressors, not overexpressed surface
proteins, so they sit outside this target space; they are reported as an unscored context
panel (see `drivers_context.csv`).

---

## Dimension 1 — Overexpression (association)
*How strongly is the protein over-abundant in tumour vs normal, at both protein and RNA level?*

Association is defined **entirely from CPTAC overexpression** — not from any external
disease-association database. Score = mean of the protein- and RNA-log2FC percentile ranks
within the candidate pool, scaled to 0–5.

| Score | Anchor |
|---|---|
| 5 | Top of the pool on both protein and RNA overexpression |
| 4 | High on both modalities |
| 3 | Moderate, or high on one modality and mid on the other |
| 2 | Modest overexpression |
| 1 | Weak overexpression (bottom of the surviving pool) |
| 0 | No net overexpression |

**Data sources:** CPTAC ccRCC proteome (tumour vs NAT log2FC) and RNA-seq (log2 FPKM
tumour vs NAT), 110 tumours / 84 normals.

---

## Dimension 2 — Tractability (intrinsic surface access)
*How accessible is the protein to a surface-directed modality — antibody, ADC, CAR, radioligand?*

Scored on **intrinsic localisation evidence only** — the clinical-status buckets
(`Approved Drug`, `Advanced Clinical`, `Phase 1 Clinical`) are explicitly removed before
scoring, as are ligand/PDB counts and drug precedent, all of which leak drug knowledge.

| Score | Anchor |
|---|---|
| 3.0 | Localisation corroborated by **two orthogonal high-confidence sources** — UniProt curated localisation *and* GO cellular component — i.e. a clean, well-evidenced cell-surface epitope |
| 1.5 | Supported by only one high-confidence source, or by medium-confidence calls / signal-peptide-and-topology prediction alone — engageable, but the surface evidence is thinner |
| 0 | No surface or secreted access |

**Data sources:** Open Targets antibody-tractability buckets, which aggregate UniProt
subcellular localisation, GO cellular component, and SignalP/TMHMM topology prediction.

**Note on an earlier version.** This dimension was previously implemented as a hardcoded
five-gene list (`HILPDA, COL23A1, NPTX2, HAPLN1, HAVCR1` → 1.5, everything else → 3.0) inside a
function documented as localisation-derived; the localisation variable it computed was never
used. It was also load-bearing — flattening it moved the then-headline p-value from 0.046 to
0.091. The rule above replaces it with a derived criterion and reproduces the same twelve
assignments. That agreement was **verified after** the rule was chosen, and the rule uses a
pre-existing external evidence ladder rather than thresholds tuned here, but it is recorded
plainly so the reader can discount it as they see fit.

---

## Dimension 3 — Safety (tolerance to inhibition)
*How risky is systemic engagement for normal tissue?*

| Score | Anchor (higher = safer) |
|---|---|
| 5 | Non-essential, LoF-tolerant, narrow normal-tissue expression, low vital-organ expression |
| 4 | Mostly tolerant; limited essentiality or expression concern |
| 3 | Mixed signals (e.g. constrained but not pan-essential) |
| 2 | Some essentiality or broad vital-tissue expression |
| 1 | Pan-essential or high constraint + broad vital expression |
| 0 | Severe predicted on-target toxicity |

**Data sources:** DepMap CRISPR mean gene-effect (essentiality); gnomAD LOEUF constraint;
Human Protein Atlas vital-organ (heart / liver / CNS / kidney) expression; count of
annotated safety liabilities.

---

## Dimension 4 — Selectivity (off-target / cross-reactivity risk)
*Off-target risk from close paralogs and breadth of normal expression.*

| Score | Anchor (higher = more selective / lower risk) |
|---|---|
| 5 | No close paralogs; tumour-restricted; low heart/liver/CNS expression |
| 4 | Few distant paralogs; manageable selectivity |
| 3 | Moderate family size or moderate vital-organ expression |
| 2 | Several close paralogs (high identity) or notable vital-organ expression |
| 1 | Large family, high-identity paralogs, broad vital expression |
| 0 | Extreme cross-reactivity risk |

**Data sources:** Ensembl/UniProt paralog lists (count + % identity); RNA tissue-specificity
(HPA); GTEx/HPA vital-organ expression.

---

## Composite score
`composite (0–100) = 20 × weighted mean of the four 0–5 dimensions.`

Weights are given as a 30/25/20/15 ratio and normalised to sum to 1:

| Dimension | Ratio | Normalised weight |
|---|---|---|
| Overexpression (association) | 30 | 0.3333 |
| Tractability (surface access) | 25 | 0.2778 |
| Safety | 20 | 0.2222 |
| Selectivity | 15 | 0.1667 |

Overexpression is weighted highest — the method is expression-driven by design. The weights
are the main subjective lever; `weight_sensitivity.csv` reports how the ranking moves over
2000 random (Dirichlet) reweightings.

A **novelty** dimension (how un-crowded a target's drug landscape is) is deliberately **not
scored**: crowdedness is a statement about the existing drug landscape, the very information
held out of the method. It is discussed qualitatively in the write-up instead.

## Development tiers
- **T1 fast-follow** — a therapeutic already engages the antigen directly (held-out label).
- **T2 discovery** — composite ≥ 70, no direct drug.
- **T3 watch** — composite 55–70.
- **T4 deprioritise** — composite < 55.
