"""
Regenerate figures/triage_scorecard.png from the committed result files.

Self-contained plotting stage (the original run built the figure inline; this
script makes it reproducible). Reads only results/*, writes figures/triage_scorecard.png.

Four panels:
  A  Rubric heatmap        — 12 candidates x 4 drug-blind dimensions (0-5)
  B  Composite & dev tier  — ranked composite (0-100), coloured by tier, drugged antigens marked
  C  Immune filter         — Spearman rho vs leukocyte signature; leukocyte markers vs targets
  D  Blind-recovery test   — the two direct-drug antigens (CA9, CD70) vs the rest on the
                             composite axis; visualises what the drug-blind ranking recovered.

Run:  python code/make_figure.py
"""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib import cm

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
def rp(*p): return os.path.join(ROOT, *p)

# ---- data -------------------------------------------------------------------
sc = pd.read_csv(rp("results", "scorecards", "scorecard_clean.csv")).sort_values("rank")
sc101 = pd.read_csv(rp("results", "scorecards", "scorecard_101.csv")).sort_values("rank")
imm = pd.read_csv(rp("results", "enrichment", "immune_filter.csv"))
rec = json.load(open(rp("results", "scorecards", "recovery_stats.json")))

DRUG_LABEL = {"CD70": "anti-CD70 CAR-T", "CA9": "girentuximab (radioconj.)"}
DIRECT = sc.loc[sc["validation_heldout"] == "ccrcc-direct-agent", "gene"].tolist()

TIER_COLOR = {"T1 fast-follow": "#2f9e57", "T2 discovery": "#e8a838",
              "T3 watch": "#9aa0a6", "T4 deprioritize": "#d1495b"}
GROUP_COLOR = {"immune_filtered": "#d1495b", "surface_candidate": "#3b7dd8",
               "tumour_intrinsic_control": "#9aa0a6"}
GROUP_LABEL = {"immune_filtered": "excluded by immune filter",
               "surface_candidate": "retained surface candidate",
               "tumour_intrinsic_control": "tumour-intrinsic control"}

fig = plt.figure(figsize=(13.6, 9.4))
gs = fig.add_gridspec(2, 2, width_ratios=[1.12, 1.0], height_ratios=[1.0, 1.0],
                      hspace=0.30, wspace=0.24, left=0.07, right=0.965, top=0.90, bottom=0.075)
fig.suptitle("CPTAC ccRCC drug-target discovery — drug-blind expression-driven method",
             fontsize=16, fontweight="bold", y=0.965)

# ---- Panel A: rubric heatmap ------------------------------------------------
axA = fig.add_subplot(gs[0, 0])
dims = [("assoc_expr", "Overexpression\n(protein+RNA)"), ("tractability", "Tractability\n(surface access)"),
        ("safety", "Safety"), ("selectivity", "Selectivity")]
M = sc[[d[0] for d in dims]].to_numpy(dtype=float)
genes = sc["gene"].tolist()
cmap = cm.get_cmap("viridis")
im = axA.imshow(M, cmap=cmap, vmin=0, vmax=5, aspect="auto")
axA.set_xticks(range(len(dims))); axA.set_xticklabels([d[1] for d in dims], fontsize=9.5)
axA.set_yticks(range(len(genes)))
axA.set_yticklabels([f"$\\it{{{g}}}$" for g in genes], fontsize=9.5)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        r, g, b, _ = cmap(M[i, j] / 5.0)
        lum = 0.299 * r + 0.587 * g + 0.114 * b
        axA.text(j, i, f"{M[i, j]:.1f}", ha="center", va="center", fontsize=9.5,
                 color="black" if lum > 0.55 else "white")
axA.set_title("A   Clean rubric — expression-derived surface/secreted targets",
              loc="left", fontsize=11.5, fontweight="bold", pad=8)
cb = fig.colorbar(im, ax=axA, fraction=0.045, pad=0.02)
cb.set_label("dimension score (0–5)", fontsize=9)

# ---- Panel B: composite & development tier ----------------------------------
axB = fig.add_subplot(gs[0, 1])
y = np.arange(len(sc))[::-1]
axB.barh(y, sc["composite_0_100"], color=[TIER_COLOR[t] for t in sc["dev_tier"]],
         edgecolor="white", height=0.72)
for yi, (_, row) in zip(y, sc.iterrows()):
    axB.text(row["composite_0_100"] + 1, yi, f"{row['composite_0_100']:.0f}",
             va="center", fontsize=9)
    if row["gene"] in DRUG_LABEL:
        axB.text(row["composite_0_100"] + 6.5, yi, f"● {DRUG_LABEL[row['gene']]}",
                 va="center", fontsize=8.5)
axB.set_yticks(y); axB.set_yticklabels([f"$\\it{{{g}}}$" for g in sc["gene"]], fontsize=9.5)
axB.set_xlim(0, 122)   # headroom for the held-out drug labels and the legend
axB.set_xlabel("composite (0–100), drug-blind", fontsize=9.5)
axB.set_title("B   Composite & development tier", loc="left", fontsize=11.5, fontweight="bold", pad=8)
tier_handles = [Patch(facecolor=c, label=t) for t, c in TIER_COLOR.items()]
tier_handles.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="black",
                           markersize=8, label="antigen with a direct drug"))
axB.legend(handles=tier_handles, fontsize=8, loc="lower right", frameon=False)

# ---- Panel C: immune filter -------------------------------------------------
axC = fig.add_subplot(gs[1, 0])
# All 101 candidates are diagnosed, but plotting 101 bars is unreadable: show the
# 10 strongest exclusions, the 12 retained candidates, and the controls.
_keep = set(sc["gene"]) | set(imm.loc[imm.group == "tumour_intrinsic_control", "gene"]) \
    | set(imm[imm.group == "immune_filtered"].nlargest(10, "immune_rho")["gene"])
imm_s = imm[imm.gene.isin(_keep)].sort_values("immune_rho", ascending=True).reset_index(drop=True)
yc = np.arange(len(imm_s))
axC.barh(yc, imm_s["immune_rho"], color=[GROUP_COLOR[g] for g in imm_s["group"]],
         edgecolor="white", height=0.74)
axC.axvline(0.40, ls="--", color="#333", lw=1)
axC.text(0.40, len(imm_s) * 0.52, "exclude ≥ 0.40", rotation=90,
         va="center", ha="right", fontsize=8.5, color="#333")
axC.set_yticks(yc); axC.set_yticklabels([f"$\\it{{{g}}}$" for g in imm_s["gene"]], fontsize=7.8)
axC.set_xlabel("Spearman ρ vs leukocyte signature (110 tumours)", fontsize=9.5)
axC.set_title("C   Immune filter separates infiltrate-tracking genes from tumour-cell targets",
              loc="left", fontsize=11.5, fontweight="bold", pad=8)
axC.legend(handles=[Patch(facecolor=GROUP_COLOR[k], label=GROUP_LABEL[k])
                    for k in ["surface_candidate", "immune_filtered", "tumour_intrinsic_control"]],
           fontsize=8, loc="lower right", frameon=False)

# ---- Panel D: blind-recovery test (Option B) --------------------------------
axD = fig.add_subplot(gs[1, 1])
# The recovery test is evaluated over ALL 101 surface candidates, not the 12
# survivors: restricting to the survivors conditions on the very selection step
# being tested, and cannot record a positive that selection threw away.
r = rec["ccrcc_direct"]
pos = set(r["positives"])
d101 = sc101.copy()
d101["is_pos"] = d101["gene"].isin(pos)
rows = {f"No ccRCC-direct\nagent (n={len(d101)-len(pos)})": 0.0,
        f"ccRCC-direct agent\n(n={len(pos)})": 1.0}
rng = np.random.default_rng(0)
for _, row in d101.iterrows():
    yb = (1.0 if row["is_pos"] else 0.0) + rng.uniform(-0.09, 0.09)  # jitter for overplotting
    if row["is_pos"]:
        axD.scatter(row["composite_0_100"], yb, s=95, color="#2f9e57", zorder=5,
                    edgecolor="black", linewidth=0.7)
    else:
        axD.scatter(row["composite_0_100"], yb, s=34, color="#9aa0a6", alpha=0.65, zorder=3)
# Label the three antibody/ADC-relevant hits downward, into the empty band between
# the two rows, so they clear the statistics box above.
# These three sit within ~10 composite points of each other, so the labels are
# staggered vertically rather than placed at a common offset.
for g, dy in (("CD70", -18), ("CA9", -34), ("ENPP3", -50)):
    if g in set(d101.gene):
        rr = d101[d101.gene == g].iloc[0]
        axD.annotate(f"{g} (#{int(rr['rank'])})", (rr["composite_0_100"], 0.9),
                     xytext=(0, dy), textcoords="offset points", ha="center",
                     fontsize=8.5, fontweight="bold",
                     arrowprops=dict(arrowstyle="-", lw=0.7, color="#555"))
axD.set_yticks([0.0, 1.0]); axD.set_yticklabels(list(rows.keys()), fontsize=9.5)
axD.set_ylim(-0.45, 1.95)
axD.set_xlim(d101.composite_0_100.min() - 4, d101.composite_0_100.max() + 5)
axD.set_xlabel("composite (0–100), drug-blind", fontsize=9.5)
axD.set_title("D   Held-out recovery, evaluated over all 101 candidates", loc="left",
              fontsize=11.5, fontweight="bold", pad=8)
txt = (f"{len(pos)} of {len(d101)} surface candidates have a ccRCC-direct agent.\n"
       f"Selection retained {r['cascade_test']['positives_retained']} of them "
       f"(chance: {r['cascade_test']['expected_by_chance']}) — p = "
       f"{r['cascade_test']['hypergeom_p']}\n"
       f"Ranking, all 101:  Mann–Whitney p = {r['ranking_test']['mannwhitney_p']}\n"
       f"CD70 and CA9 rank top-4, but {r['misses']['n']} positives are missed —\n"
       "ENPP3 (anti-ENPP3 ADC) scores #9 yet fails the fold-change cut.")
axD.text(0.025, 0.995, txt, transform=axD.transAxes, fontsize=8.2, va="top", ha="left",
         bbox=dict(boxstyle="round,pad=0.5", facecolor="#fdf1f1", edgecolor="#e8c9c9"))

out = rp("figures", "triage_scorecard.png")
fig.savefig(out, dpi=200, bbox_inches="tight")
print("wrote", out)
