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
imm = pd.read_csv(rp("results", "enrichment", "immune_filter.csv"))
rec = json.load(open(rp("results", "scorecards", "recovery_stats.json")))

DRUG_LABEL = {"CD70": "anti-CD70 CAR-T", "CA9": "girentuximab (radioconj.)"}
DIRECT = sc.loc[sc["validation_heldout"] == "direct-drug", "gene"].tolist()

TIER_COLOR = {"T1 fast-follow": "#2f9e57", "T2 discovery": "#e8a838",
              "T3 watch": "#9aa0a6", "T4 deprioritize": "#d1495b"}
GROUP_COLOR = {"leukocyte_marker": "#d1495b", "surface_candidate": "#3b7dd8",
               "tumour_intrinsic_control": "#9aa0a6"}
GROUP_LABEL = {"leukocyte_marker": "leukocyte marker", "surface_candidate": "surface candidate",
               "tumour_intrinsic_control": "tumour-intrinsic control"}

fig = plt.figure(figsize=(13.6, 9.4))
gs = fig.add_gridspec(2, 2, width_ratios=[1.12, 1.0], height_ratios=[1.0, 1.0],
                      hspace=0.30, wspace=0.24, left=0.07, right=0.965, top=0.90, bottom=0.075)
fig.suptitle("CPTAC ccRCC drug-target discovery — single blind expression-driven method",
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
axB.set_xlim(0, 100); axB.set_xlabel("composite (0–100), drug-blind", fontsize=9.5)
axB.set_title("B   Composite & development tier", loc="left", fontsize=11.5, fontweight="bold", pad=8)
tier_handles = [Patch(facecolor=c, label=t) for t, c in TIER_COLOR.items()]
tier_handles.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="black",
                           markersize=8, label="antigen with a direct drug"))
axB.legend(handles=tier_handles, fontsize=8, loc="lower right", frameon=False)

# ---- Panel C: immune filter -------------------------------------------------
axC = fig.add_subplot(gs[1, 0])
imm_s = imm.sort_values("immune_rho", ascending=True).reset_index(drop=True)
yc = np.arange(len(imm_s))
axC.barh(yc, imm_s["immune_rho"], color=[GROUP_COLOR[g] for g in imm_s["group"]],
         edgecolor="white", height=0.74)
axC.axvline(0.40, ls="--", color="#333", lw=1)
axC.text(0.40, len(imm_s) * 0.52, "exclude ≥ 0.40", rotation=90,
         va="center", ha="right", fontsize=8.5, color="#333")
axC.set_yticks(yc); axC.set_yticklabels([f"$\\it{{{g}}}$" for g in imm_s["gene"]], fontsize=7.8)
axC.set_xlabel("Spearman ρ vs leukocyte signature (110 tumours)", fontsize=9.5)
axC.set_title("C   Immune filter separates leukocyte markers from tumour-cell targets",
              loc="left", fontsize=11.5, fontweight="bold", pad=8)
axC.legend(handles=[Patch(facecolor=GROUP_COLOR[k], label=GROUP_LABEL[k])
                    for k in ["surface_candidate", "leukocyte_marker", "tumour_intrinsic_control"]],
           fontsize=8, loc="lower right", frameon=False)

# ---- Panel D: blind-recovery test (Option B) --------------------------------
axD = fig.add_subplot(gs[1, 1])
sc2 = sc.copy()
sc2["is_direct"] = sc2["gene"].isin(DIRECT)
rng = np.random.default_rng(0)
rows = {"Other candidates\n(n=10)": 0.0, "Direct-drug antigens\n(n=2)": 1.0}
for _, row in sc2.iterrows():
    grp = "Direct-drug antigens\n(n=2)" if row["is_direct"] else "Other candidates\n(n=10)"
    yb = rows[grp] + (rng.uniform(-0.14, 0.14) if not row["is_direct"] else 0.0)
    if row["is_direct"]:
        axD.scatter(row["composite_0_100"], yb, s=160, color="#2f9e57", zorder=5,
                    edgecolor="black", linewidth=0.8)
        # splay the two labels outward so they don't collide (composites 73 vs 80)
        dx, ha = ((-9, "right") if row["gene"] == "CA9" else (9, "left"))
        axD.annotate(f"{row['gene']}  (rank {int(row['rank'])})",
                     (row["composite_0_100"], yb), xytext=(dx, 16), textcoords="offset points",
                     ha=ha, fontsize=9, fontweight="bold")
    else:
        axD.scatter(row["composite_0_100"], yb, s=70, color="#9aa0a6", alpha=0.85, zorder=3)
axD.set_yticks([0.0, 1.0]); axD.set_yticklabels(list(rows.keys()), fontsize=9.5)
axD.set_ylim(-0.6, 1.55); axD.set_xlim(40, 92)
axD.set_xlabel("composite (0–100), drug-blind", fontsize=9.5)
axD.set_title("D   Blind recovery of the drugged antigens", loc="left",
              fontsize=11.5, fontweight="bold", pad=8)
hp = rec["both_in_top3_hypergeom_p"]; mw = rec["mannwhitney_p"]
txt = ("Run blind to drug knowledge, both direct-drug antigens land\n"
       "in the top 3 of 12 (CD70 #1, CA9 #3).\n"
       f"Rank recovery:  hypergeometric p = {hp}\n"
       f"Score separation:  Mann–Whitney p = {mw}\n"
       "Driver overlap:  0 / 5 mutation drivers in the target set")
axD.text(0.025, 0.60, txt, transform=axD.transAxes, fontsize=8.5, va="top", ha="left",
         bbox=dict(boxstyle="round,pad=0.5", facecolor="#eef4fb", edgecolor="#cfe0f2"))

out = rp("figures", "triage_scorecard.png")
fig.savefig(out, dpi=200, bbox_inches="tight")
print("wrote", out)
