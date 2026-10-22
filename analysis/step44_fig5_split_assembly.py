#!/usr/bin/env python3
"""Split revised Figure 5 into two main figures:
- Figure 5: A (hub6 ssGSEA) + B-E (CIBERSORT stacked bars, 4 datasets)
- Figure 6 (new): A (Spearman correlation heatmaps) + B (differential cell violins)
Original Figures 6/7/8 are renumbered to 7/8/9 separately (file renames)."""
import pymupdf

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission"
P = f"{BASE}/03_figures/Figures_ggplot2_NPG/main_panels"
OUT = f"{BASE}/03_figures/Figures_ggplot2_NPG/composites"

MM = 72 / 25.4
PAGE_W = 178 * MM
MARGIN_X = 8.0
ROW_GAP = 6.0
COL_GAP = 4.0
TOP = 8.0


def mm(v):
    return v * MM


def pdf_h(path, w_mm):
    d = pymupdf.open(path)
    r = d[0].rect
    d.close()
    return w_mm * r.height / r.width


class Sheet:
    def __init__(self, height_mm):
        self.doc = pymupdf.open()
        self.page = self.doc.new_page(width=PAGE_W, height=mm(height_mm))
        self.cache = {}

    def _pdf(self, path):
        if path not in self.cache:
            self.cache[path] = pymupdf.open(path)
        return self.cache[path]

    def letter(self, x_mm, y_mm, text):
        self.page.insert_text((mm(x_mm), mm(y_mm) + 3.4), text,
                              fontsize=9, fontname="hebo", color=(0, 0, 0))

    def pdf_panel(self, path, x_mm, y_mm, w_mm, letter=None):
        src = self._pdf(path)
        r = src[0].rect
        h_mm = w_mm * r.height / r.width
        dest = pymupdf.Rect(mm(x_mm), mm(y_mm), mm(x_mm + w_mm), mm(y_mm + h_mm))
        self.page.show_pdf_page(dest, src, 0)
        if letter:
            self.letter(x_mm - 4.5, y_mm, letter)
        return h_mm

    def save(self, name):
        h_pt = self.page.rect.height
        self.doc.save(f"{OUT}/{name}", deflate=True, garbage=3)
        for d in self.cache.values():
            d.close()
        self.doc.close()
        print("saved", name, f"({PAGE_W / MM:.0f}x{h_pt / MM:.0f}mm)")


CW = (178 - 2 * MARGIN_X - COL_GAP) / 2
X1 = MARGIN_X
X2 = MARGIN_X + CW + COL_GAP
FW = 178 - 2 * MARGIN_X  # full-width usable

## ---------------- Figure 5 (split: A-E only) ----------------
fA = f"{P}/Fig5A_hub6_GSVA.pdf"
fB = f"{P}/Fig5B_CIBERSORT_stacked_GSE76293.pdf"
fC = f"{P}/Fig5C_CIBERSORT_stacked_GSE69528.pdf"
fD = f"{P}/Fig5D_CIBERSORT_stacked_GSE243066.pdf"
fE = f"{P}/Fig5E_CIBERSORT_stacked_GSE100159.pdf"

hA = pdf_h(fA, FW)
hB = pdf_h(fB, CW)
hC = pdf_h(fC, CW)
hD = pdf_h(fD, CW)
hE = pdf_h(fE, CW)

H = TOP + hA + ROW_GAP + max(hB, hC) + ROW_GAP + max(hD, hE) + 6
s = Sheet(H)
y = TOP
s.pdf_panel(fA, MARGIN_X, y, FW, "A")
y += hA + ROW_GAP
s.pdf_panel(fB, X1, y, CW, "B")
s.pdf_panel(fC, X2, y, CW, "C")
y += max(hB, hC) + ROW_GAP
s.pdf_panel(fD, X1, y, CW, "D")
s.pdf_panel(fE, X2, y, CW, "E")
s.save("Figure5_composite.pdf")

## ---------------- Figure 6 (new: former Fig5 F+G, re-lettered A/B) ----------------
gA = f"{P}/Fig5F_spearman_corrplot.pdf"
gB = f"{P}/Fig5G_celltype_violins.pdf"

hGa = pdf_h(gA, FW)
hGb = pdf_h(gB, FW)

H = TOP + hGa + ROW_GAP + hGb + 6
s = Sheet(H)
y = TOP
s.pdf_panel(gA, MARGIN_X, y, FW, "A")
y += hGa + ROW_GAP
s.pdf_panel(gB, MARGIN_X, y, FW, "B")
s.save("Figure6_composite.pdf")
print("step44 done.")
