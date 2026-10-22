#!/usr/bin/env python3
"""Figure S15 assembly: ssGSEA immune infiltration, four cohorts.
A-D: per-cohort boxplots (17-18 cell types, Wilcoxon stars), full width.
E-H: per-cohort hub-gene x cell-type Spearman heatmaps, 2 per row.
178 mm double-column page, Helvetica bold letters."""
import pymupdf

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission"
PDIR = f"{BASE}/03_figures/Figures_ggplot2_NPG/supp_panels"
OUT = f"{BASE}/03_figures/Figures_ggplot2_NPG/composites/FigureS15_composite.pdf"

MM = 72 / 25.4
MARGIN_X = 8.0
TOP = 8.0
ROW_GAP = 4.0
COL_GAP = 4.0
FW = 178 - 2 * MARGIN_X          # 162 mm usable
CW = (FW - COL_GAP) / 2          # half-width column

ORDER = ["GSE76293", "GSE69528", "GSE243066", "GSE100159"]
BOX = [f"{PDIR}/FigS15_box_{c}.pdf" for c in ORDER]
COR = [f"{PDIR}/FigS15_cor_{c}.pdf" for c in ORDER]


def panel_h(path, w_mm):
    src = pymupdf.open(path)
    r = src[0].rect
    src.close()
    return w_mm * r.height / r.width


def mm(v):
    return v * MM


# ---- pass 1: total height ----
h_box = [panel_h(p, FW) for p in BOX]
h_cor = [panel_h(p, CW) for p in COR]
total_h = (TOP + sum(h_box) + 4 * ROW_GAP
           + max(h_cor[0], h_cor[1]) + ROW_GAP
           + max(h_cor[2], h_cor[3]) + 6)

# ---- pass 2: place ----
doc = pymupdf.open()
page = doc.new_page(width=178 * MM, height=total_h * MM)
letters = "ABCDEFGH"
li = 0
y = TOP


def place(path, x_mm, y_mm, w_mm, letter):
    src = pymupdf.open(path)
    r = src[0].rect
    h_mm = w_mm * r.height / r.width
    dest = pymupdf.Rect(mm(x_mm), mm(y_mm), mm(x_mm + w_mm), mm(y_mm + h_mm))
    page.show_pdf_page(dest, src, 0)
    page.insert_text(pymupdf.Point(mm(x_mm), mm(y_mm + 3.2)), letter,
                     fontname="hebo", fontsize=9, color=(0, 0, 0))
    src.close()


for i, p in enumerate(BOX):
    place(p, MARGIN_X, y, FW, letters[li]); li += 1
    y += h_box[i] + ROW_GAP

for row in range(2):
    for col in range(2):
        idx = row * 2 + col
        place(COR[idx], MARGIN_X + col * (CW + COL_GAP), y, CW, letters[li]); li += 1
    y += max(h_cor[row * 2], h_cor[row * 2 + 1]) + ROW_GAP

doc.save(OUT, deflate=True, garbage=3)
doc.close()
print("saved", OUT, f"(178x{total_h:.0f}mm)")
