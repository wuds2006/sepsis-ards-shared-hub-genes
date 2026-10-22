#!/usr/bin/env python3
# Graphical Abstract for AJRCCM submission
# Nature-family palette: #6A5C95 #D99AAE #E7A889 #A9CBB8 #CFE3E6 #8FA2D6 #5B608C
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# palette
PURPLE = "#6A5C95"
PINK = "#D99AAE"
PEACH = "#E7A889"
GREEN = "#A9CBB8"
LIGHTBLUE = "#CFE3E6"
PERIWINKLE = "#8FA2D6"
DARKBLUE = "#5B608C"
GREY = "#C9CCD3"
INK = "#2A2A3A"

fig, ax = plt.subplots(figsize=(12.5, 9.0), dpi=300)
ax.set_xlim(0, 100)
ax.set_ylim(0, 72)
ax.axis("off")

def box(x, y, w, h, text, fc, ec="none", fs=11, tc=INK, weight="normal", style="round,pad=0.6,rounding_size=1.2", italic=False):
    p = FancyBboxPatch((x, y), w, h, boxstyle=style, fc=fc, ec=ec, lw=0)
    ax.add_patch(p)
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs,
            color=tc, weight=weight, style="italic" if italic else "normal")
    return (x + w/2, y + h/2)

def arrow(p1, p2, color=INK, lw=1.8):
    a = FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=18, color=color, lw=lw)
    ax.add_patch(a)

# Title
ax.text(50, 69, "Shared Hub Genes and Immune Signatures of\nSepsis-Associated Acute Respiratory Distress Syndrome",
        ha="center", va="center", fontsize=15.5, weight="bold", color=PURPLE)

# --- Row 1: datasets ---
d1 = box(20, 58, 26, 6.5, "GSE76293\nARDS (12 vs 12)\nblood PMNs", LIGHTBLUE, fs=9.5, weight="bold")
d2 = box(56, 58, 26, 6.5, "GSE69528\nSepsis (83 vs 28)\nwhole blood", LIGHTBLUE, fs=9.5, weight="bold")

# --- Row 2: DEG + WGCNA ---
m1 = box(20, 47, 26, 6, "DEG + WGCNA\n557 / 1,458 DEGs", PERIWINKLE, fs=9.5)
m2 = box(56, 47, 26, 6, "DEG + WGCNA\n13 / 11 modules", PERIWINKLE, fs=9.5)
arrow(d1, m1); arrow(d2, m2)

# --- Row 3: intersection ---
inter = box(50, 37.5, 44, 5.5, "Intersection: 179 shared candidate genes", PEACH, fs=11, weight="bold")
arrow(m1, inter); arrow(m2, inter)

# --- Row 4: ML ---
ml = box(50, 28.5, 44, 5.5, "11 machine-learning classifiers + SHAP", PINK, fs=11, weight="bold")
arrow(inter, ml)

# --- Row 5: 6 hub genes ---
ax.text(50, 22.0, "Six shared hub genes", ha="center", va="center", fontsize=10.5, weight="bold", color=DARKBLUE)
genes = ["CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1"]
gw, gh, gap = 14.0, 5.2, 1.2
total = len(genes)*gw + (len(genes)-1)*gap
x0 = 50 - total/2
for i, g in enumerate(genes):
    bx = x0 + i*(gw+gap)
    box(bx, 15.2, gw, gh, g, PURPLE, fs=11, tc="white", weight="bold", italic=True)
arrow(ml, (50, 22.0))

# --- Bottom: three result cards ---
cards = [
    (8,  3, 26, 9, "M2 macrophage\npolarization", "shared immune\nsuppression", GREEN),
    (37, 3, 26, 9, "Complement\nactivation", "sole shared top-20\nimmune pathway", LIGHTBLUE),
    (66, 3, 26, 9, "Glucocorticoid\nrepositioning", "VNN1 as druggable\ntarget", PEACH),
]
for x, y, w, h, t, s, c in cards:
    box(x, y, w, h, "", "white", ec=c, style="round,pad=0.6,rounding_size=1.2")
    ax.text(x+w/2, y+h-3.0, t, ha="center", va="center", fontsize=10.5, weight="bold", color=DARKBLUE)
    ax.text(x+w/2, y+1.7, s, ha="center", va="center", fontsize=8.3, color=INK)
# arrows from hub genes to cards
arrow((50, 15.2), (21, 12.4))
arrow((50, 15.2), (50, 12.4))
arrow((50, 15.2), (79, 12.4))

# footer
ax.text(50, 0.6, "Integrated multi-omics · machine learning · single-cell · drug repositioning & docking",
        ha="center", va="center", fontsize=8.5, color=GREY)

out = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission/04_graphical_abstract/graphical_abstract"
fig.savefig(out + ".png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(out + ".pdf", bbox_inches="tight", facecolor="white")
fig.savefig(out + ".tiff", dpi=300, bbox_inches="tight", facecolor="white")
print("Saved graphical abstract (png/pdf/tiff)")
