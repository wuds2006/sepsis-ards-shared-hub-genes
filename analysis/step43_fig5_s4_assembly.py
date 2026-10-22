#!/usr/bin/env python3
"""Assemble revised Figure 5 (CIBERSORT immune infiltration) and Figure S4
(hub-gene x CIBERSORT correlation heatmaps) composites, vector-preserving."""
import pymupdf

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission"
P = f"{BASE}/03_figures/Figures_ggplot2_NPG/main_panels"
S = f"{BASE}/03_figures/Figures_ggplot2_NPG/supp_panels"
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

## ---------------- Figure 5 (v2) ----------------
fA = f"{P}/Fig5A_hub6_ssGSEA.pdf"
fB = f"{P}/Fig5B_CIBERSORT_stacked_GSE76293.pdf"
fC = f"{P}/Fig5C_CIBERSORT_stacked_GSE69528.pdf"
fD = f"{P}/Fig5D_CIBERSORT_stacked_GSE243066.pdf"
fE = f"{P}/Fig5E_CIBERSORT_stacked_GSE100159.pdf"
fF = f"{P}/Fig5F_spearman_corrplot.pdf"
fG = f"{P}/Fig5G_celltype_violins.pdf"

hA = pdf_h(fA, FW)
hB = pdf_h(fB, CW)
hC = pdf_h(fC, CW)
hD = pdf_h(fD, CW)
hE = pdf_h(fE, CW)
hF = pdf_h(fF, FW)
hG = pdf_h(fG, FW)

H = (TOP + hA + ROW_GAP + max(hB, hC) + ROW_GAP + max(hD, hE) + ROW_GAP
     + hF + ROW_GAP + hG + 6)
s = Sheet(H)
y = TOP
s.pdf_panel(fA, MARGIN_X, y, FW, "A")
y += hA + ROW_GAP
s.pdf_panel(fB, X1, y, CW, "B")
s.pdf_panel(fC, X2, y, CW, "C")
y += max(hB, hC) + ROW_GAP
s.pdf_panel(fD, X1, y, CW, "D")
s.pdf_panel(fE, X2, y, CW, "E")
y += max(hD, hE) + ROW_GAP
s.pdf_panel(fF, MARGIN_X, y, FW, "F")
y += hF + ROW_GAP
s.pdf_panel(fG, MARGIN_X, y, FW, "G")
s.save("Figure5_composite.pdf")

## ---------------- Figure S4 (v2) ----------------
gA = f"{S}/FigS4A_hub_ciber_cor_GSE76293.pdf"
gB = f"{S}/FigS4B_hub_ciber_cor_GSE69528.pdf"
gC = f"{S}/FigS4C_hub_ciber_cor_GSE243066.pdf"
gD = f"{S}/FigS4D_hub_ciber_cor_GSE100159.pdf"
hS = max(pdf_h(gA, CW), pdf_h(gB, CW), pdf_h(gC, CW), pdf_h(gD, CW))
H = TOP + hS + ROW_GAP + hS + 6
s = Sheet(H)
y = TOP
s.pdf_panel(gA, X1, y, CW, "A")
s.pdf_panel(gB, X2, y, CW, "B")
y += hS + ROW_GAP
s.pdf_panel(gC, X1, y, CW, "C")
s.pdf_panel(gD, X2, y, CW, "D")
s.save("FigureS4_composite.pdf")
print("step43 done.")
