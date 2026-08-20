"""Regenerate the embedded DATA blob and the recovery panel inside index.html.

index.html is a self-contained report with its data inlined as one `const DATA = {...}`
line. Rather than hand-editing it (which is how it drifted out of sync with results/ in
the first place), this script rebuilds that line — and the recovery-panel markup that
reads it — from the current contents of results/.

    python code/update_report.py
"""
import json
import os
import re

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
HTML = os.path.join(ROOT, "index.html")

sc = pd.read_csv(f"{RES}/scorecards/scorecard_clean.csv").sort_values("rank")
sc101 = pd.read_csv(f"{RES}/scorecards/scorecard_101.csv").sort_values("rank")
imm = pd.read_csv(f"{RES}/enrichment/immune_filter.csv")
mc = pd.read_csv(f"{RES}/enrichment/weight_sensitivity.csv")
drv = pd.read_csv(f"{RES}/scorecards/drivers_context.csv")
rec = json.load(open(f"{RES}/scorecards/recovery_stats.json"))
ev = json.load(open(f"{RES}/evidence/evidence_detail.json"))

AGENT_LABEL = {"CD70": "anti-CD70 CAR-T (ALLO-316)",
               "CA9": "girentuximab radioconjugate"}


def _num(v, nd=3):
    try:
        return round(float(v), nd)
    except (TypeError, ValueError):
        return None


genes = []
for r in sc.itertuples():
    e = ev.get(r.gene, {})
    genes.append({
        "gene": r.gene, "prot_log2fc": round(r.prot_log2fc, 2),
        "rna_log2fc": round(r.rna_log2fc, 2), "immune_rho": round(r.immune_rho, 3),
        "dims": {"assoc_expr": r.assoc_expr, "tractability": r.tractability,
                 "safety": r.safety, "selectivity": r.selectivity},
        "composite": r.composite_0_100, "tier": r.dev_tier,
        "direct_drug": AGENT_LABEL.get(r.gene, "") if r.validation_heldout != "none" else "",
        "locations": "", "rna_tissue_spec": e.get("rna_tissue_spec") or "",
        "depmap": _num(e.get("depmap_mean_geneeffect"), 4),
        "loeuf": _num(e.get("loeuf"), 4)})

# The immune tab shows the same subset the figure does: strongest exclusions,
# the retained shortlist, and the tumour-intrinsic controls.
keep = set(sc.gene) | set(imm[imm.group == "tumour_intrinsic_control"].gene) \
    | set(imm[imm.group == "immune_filtered"].nlargest(10, "immune_rho").gene)
immune = [{"gene": r.gene, "rho": float(r.immune_rho),
           "group": ("leukocyte_marker" if r.group == "immune_filtered"
                     else ("tumour_intrinsic_control"
                           if r.group == "tumour_intrinsic_control" else "candidate"))}
          for r in imm[imm.gene.isin(keep)].sort_values(
              "immune_rho", ascending=False).itertuples()]

cd = rec["ccrcc_direct"]
top101 = {r.gene: int(r.rank) for r in sc101.itertuples()}
DATA = {
    "genes": genes,
    "drivers": [{"symbol": r.symbol, "mutation_rate_pct": r.mutation_rate_pct,
                 "lesion_type": r.lesion_type,
                 "in_surface_target_set": r.in_surface_target_set} for r in drv.itertuples()],
    "immune": immune,
    "mc": {r.gene: {"mean_rank": r.mc_mean_rank, "std": r.mc_rank_std, "p1": r.p_rank1}
           for r in mc.itertuples()},
    "recovery": {
        "method": rec["method"],
        "n_surface_candidates": rec["n_surface_candidates"],
        "n_positives": cd["n_positives"],
        "positives": cd["positives"],
        "ranks_of_101": {g: top101.get(g) for g in ("CD70", "CA9", "ENPP3")},
        "cascade_retained": cd["cascade_test"]["positives_retained"],
        "cascade_expected": cd["cascade_test"]["expected_by_chance"],
        "cascade_p": cd["cascade_test"]["hypergeom_p"],
        "ranking_mw_p": cd["ranking_test"]["mannwhitney_p"],
        "ranking_top12_p": cd["ranking_test"]["top_k_p"]["12"],
        "n_missed": cd["misses"]["n"],
        "missed_immune": cd["misses"]["dropped_by_immune_filter"],
        "missed_fc": cd["misses"]["dropped_by_fc_cutoff"],
        "legacy_p": rec.get("ranking_given_selection_LEGACY", {}).get("hypergeom_p"),
    },
    "weights": {"assoc_expr": 0.3333, "tractability": 0.2778,
                "safety": 0.2222, "selectivity": 0.1667},
}

html = open(HTML).read()

# 1. swap the data blob
html = re.sub(r"const DATA = \{.*?\n", "const DATA = " + json.dumps(DATA) + "\n",
              html, count=1, flags=re.S)

# 2. rewrite the recovery panel to report the unconditional tests
new_render = '''function renderRecovery(){
  const r=DATA.recovery;const kv=document.getElementById("rec-kv");
  const R=r.ranks_of_101;
  kv.innerHTML=`
    <div>Held-out positives</div><div>${r.n_positives} of ${r.n_surface_candidates} surface candidates</div>
    <div>CD70 / CA9 / ENPP3 rank</div><div>#${R.CD70} / #${R.CA9} / #${R.ENPP3} of ${r.n_surface_candidates}</div>
    <div>Cascade test (selection)</div><div>${r.cascade_retained} retained vs ${r.cascade_expected} expected &mdash; <b>p = ${r.cascade_p}</b></div>
    <div>Ranking test (all ${r.n_surface_candidates})</div><div>Mann&ndash;Whitney <b>p = ${r.ranking_mw_p}</b>; top-12 <b>p = ${r.ranking_top12_p}</b></div>
    <div>Positives missed</div><div>${r.n_missed} of ${r.n_positives} &mdash; ${r.missed_immune.length} by the immune filter, ${r.missed_fc.length} by the fold-change cut</div>
    <div>Legacy statistic</div><div>p = ${r.legacy_p} &mdash; <i>conditional on selection; not the headline</i></div>`;
}'''
html = re.sub(r"function renderRecovery\(\)\{.*?\n\}", new_render, html, count=1, flags=re.S)

# 3. the two summary pills under the recovery heading
html = html.replace(
    'document.querySelectorAll(".pill")[0].textContent="hypergeometric p = "+DATA.recovery.both_in_top3_hypergeom_p;',
    'document.querySelectorAll(".pill")[0].textContent="cascade p = "+DATA.recovery.cascade_p;')
html = html.replace(
    'document.querySelectorAll(".pill")[1].textContent="Mann\\u2013Whitney p = "+DATA.recovery.mannwhitney_p;',
    'document.querySelectorAll(".pill")[1].textContent="ranking Mann\\u2013Whitney p = "+DATA.recovery.ranking_mw_p;')

# 4. the narrative paragraph that asserted the retracted result
old_claim = re.search(r"ranks 12 surface/secreted candidates.*?without ever being told what is drugged\.",
                      html, flags=re.S)
if old_claim:
    html = html.replace(old_claim.group(0), (
        "scores all 101 surface/secreted candidates and shortlists 12. When drug status is "
        "revealed afterwards, <b>CD70 (#2), CA9 (#4) and ENPP3 (#9) of 101</b> all carry "
        "direct-acting agents &mdash; a real signal at the top of the ranking. It does not, "
        "however, survive as a claim about the method: under a rule applied to all 101 "
        "candidates there are <b>14 positives, not 2</b>, and against that denominator neither "
        "selection (hypergeometric <b>p = 0.52</b>) nor the composite ranking "
        "(Mann&ndash;Whitney <b>p = 0.32</b>) significantly enriches for them. An earlier "
        "version reported p = 0.0455 from a null that permuted only the 12 survivors &mdash; "
        "conditional on the very selection step it was meant to test."))

open(HTML, "w").write(html)
print(f"index.html updated ({len(html)//1024} KB); "
      f"{len(genes)} genes, {len(immune)} immune rows, {cd['n_positives']} positives")
