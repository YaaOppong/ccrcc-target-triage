#!/usr/bin/env python3
"""
Build the Open Graph / social preview card for the Pages site.

OUTPUT: figures/og_card.png  — 1200x630 (1.91:1), the size LinkedIn, Slack,
        X and Facebook all read from <meta property="og:image">.

WHY A PURPOSE-BUILT CARD: none of the analysis figures survive the ~350px a
LinkedIn "Featured" card is rendered at — triage_scorecard.png is a dense
multi-panel figure that becomes unreadable at that size. A preview card is a
*hero figure*, not a chart: one headline, three numbers, one claim.

Palette is inherited from index.html (accent #2166ac) so the card and the page
it opens read as one piece.

Regenerate with:  python code/make_og_card.py
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

# --- palette: same tokens as index.html -------------------------------------
INK    = "#1a1d21"
MUTED  = "#6b7280"
ACCENT = "#2166ac"
ARROW  = "#c3cbd4"
BG     = "#ffffff"
TINT   = "#eef4fb"   # the .headline gradient top stop from index.html
TINTED = "#cfe0f2"   # .headline border

# Font stack ordered by whether a REAL bold face exists. Helvetica / Helvetica
# Neue ship only weight 400 here, and matplotlib does not synthesise a bold —
# it silently renders regular, which is what makes a card look washed out.
# Every family below exposes a genuine 700.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Source Sans Pro", "Arial", "Segoe UI",
                                   "Roboto", "DejaVu Sans"]

W, H, DPI = 1200, 630, 100
fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI, facecolor=BG)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off()
ax.set_xlim(0, 1); ax.set_ylim(0, 1)

def text(x, y, s, size, color=INK, weight="normal", ha="left", va="baseline", **kw):
    return ax.text(x, y, s, size=size, color=color, weight=weight,
                   ha=ha, va=va, transform=ax.transAxes, **kw)

# --- left accent keyline ----------------------------------------------------
ax.add_patch(Rectangle((0, 0), 0.013, 1, color=ACCENT, transform=ax.transAxes))

L, R = 0.070, 0.960          # text margins

# --- masthead: context left, source right, one baseline ---------------------
text(L, 0.882, "CPTAC PROTEOGENOMICS  ·  110 TUMOURS  ·  DRUG-BLIND",
     13, MUTED, "bold")
text(R, 0.882, "YaaOppong.github.io/ccrcc-target-triage", 13, MUTED, ha="right")

# --- headline ---------------------------------------------------------------
text(L, 0.735, "Surface target triage", 49, INK, "bold")
text(L, 0.612, "in clear cell RCC", 49, ACCENT, "bold")

# --- funnel: three numbers, the last one is the hero ------------------------
# Reading left to right IS the method; the arrows carry the narrowing. Spread
# across the full measure so the row balances the headline above it.
stages = [(0.070, "11,710", "proteins\nmeasured",                 36, MUTED),
          (0.375, "101",    "surface / secreted\ndrug-reachable", 36, MUTED),
          (0.715, "12",     "scored\ncandidates",                 60, ACCENT)]

for x, num, cap, size, col in stages:
    text(x, 0.360, num, size, col, "bold")
    text(x, 0.298, cap, 14.5, MUTED, va="top", linespacing=1.35)

# Drawn, not typed: the arrow glyph is absent from several of these families and
# matplotlib drops missing glyphs silently rather than falling back.
for x in (0.276, 0.578):
    ax.annotate("", xy=(x + 0.056, 0.381), xytext=(x, 0.381),
                xycoords=ax.transAxes, textcoords=ax.transAxes,
                arrowprops=dict(arrowstyle="-|>,head_width=0.22,head_length=0.45",
                                color=ARROW, linewidth=2.6, shrinkA=0, shrinkB=0))

# --- the payoff -------------------------------------------------------------
ax.add_patch(FancyBboxPatch((L - 0.024, 0.055), (R - L) + 0.024, 0.135,
                            boxstyle="round,pad=0,rounding_size=0.012",
                            facecolor=TINT, edgecolor=TINTED,
                            linewidth=1.2, transform=ax.transAxes))
text(L, 0.103, "Both antigens with agents in trials ranked #1 and #3 — blind",
     20, INK, "bold")

fig.savefig("figures/og_card.png", dpi=DPI, facecolor=BG)
print(f"wrote figures/og_card.png ({W}x{H})")
