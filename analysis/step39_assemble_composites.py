#!/usr/bin/env python3
"""Assemble Figure 1-8 composites from ggplot2 NPG panels (vector-preserving).

PDF panels are placed with show_pdf_page (vectors and text stay editable in
Illustrator). PNG panels (workflow schematic, PyMOL poses) are embedded as
high-resolution rasters, as in the original submission.
Layout is cursor-based: panel heights are measured before placement, so rows
never overlap and page height is exact.
"""
import pymupdf

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission"
P = f"{BASE}/03_figures/Figures_ggplot2_NPG/main_panels"
EXP = f"{BASE}/07_figure_exports"
SC = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/08_singlecell"
DK = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/09_docking_drug"
OUT = f"{BASE}/03_figures/Figures_ggplot2_NPG/composites"

MM = 72 / 25.4
PAGE_W = 178 * MM
MARGIN_X = 8.0
ROW_GAP = 6.0
COL_GAP = 4.0


def mm(v):
    return v * MM


def pdf_h(path, w_mm):
    d = pymupdf.open(path)
    r = d[0].rect
    d.close()
    return w_mm * r.height / r.width


def png_h(path, w_mm):
    pix = pymupdf.Pixmap(path)
    h = w_mm * pix.height / pix.width
    pix = None
    return h


class Sheet:
    def __init__(self, height_mm):
        self.doc = pymupdf.open()
        self.page = self.doc.new_page(width=PAGE_W, height=mm(height_mm))
        self.page_width_pt = PAGE_W
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

    def png_panel(self, path, x_mm, y_mm, w_mm, letter=None):
        pix = pymupdf.Pixmap(path)
        h_mm = w_mm * pix.height / pix.width
        dest = pymupdf.Rect(mm(x_mm), mm(y_mm), mm(x_mm + w_mm), mm(y_mm + h_mm))
        self.page.insert_image(dest, filename=path)
        if letter:
            self.letter(x_mm - 4.5, y_mm, letter)
        return h_mm

    def save(self, name):
        h_pt = self.page.rect.height
        self.doc.save(f"{OUT}/{name}", deflate=True, garbage=3)
        for d in self.cache.values():
            d.close()
        self.doc.close()
        print("saved", name, f"({self.page_width_pt:.0f}x{h_pt:.0f}pt)")


CW = (178 - 2 * MARGIN_X - COL_GAP) / 2
X1 = MARGIN_X
X2 = MARGIN_X + CW + COL_GAP
THIRD = (178 - 2 * MARGIN_X - 2 * COL_GAP) / 3
XA = MARGIN_X
XB = MARGIN_X + THIRD + COL_GAP
XC = MARGIN_X + 2 * (THIRD + COL_GAP)
TOP = 8.0

## ---------------- Figure 1 ----------------
wf = f"{EXP}/Figure 1/parts/A_workflow_schematic.png"
hA = png_h(wf, CW)
hB = pdf_h(f"{P}/Fig1B1_volcano_GSE76293.pdf", CW - 5)
hC = pdf_h(f"{P}/Fig1C1_module_trait_GSE76293.pdf", CW - 5)
hD = pdf_h(f"{P}/Fig1D_venn_4sets.pdf", 62)
H = TOP + hA + ROW_GAP + hB + ROW_GAP + hC + ROW_GAP + hD + 6
s = Sheet(H)
y = TOP
s.png_panel(wf, X1, y, CW, "A")
y += hA + ROW_GAP
s.pdf_panel(f"{P}/Fig1B1_volcano_GSE76293.pdf", X1, y, CW - 5, "B")
s.pdf_panel(f"{P}/Fig1B2_volcano_GSE69528.pdf", X2, y, CW - 5)
y += hB + ROW_GAP
s.pdf_panel(f"{P}/Fig1C1_module_trait_GSE76293.pdf", X1, y, CW - 5, "C")
s.pdf_panel(f"{P}/Fig1C2_module_trait_GSE69528.pdf", X2, y, CW - 5)
y += hC + ROW_GAP
s.pdf_panel(f"{P}/Fig1D_venn_4sets.pdf", X1, y, 62, "D")
s.save("Figure1_composite.pdf")

## ---------------- Figure 2 ----------------
w3 = THIRD
hb = pdf_h(f"{P}/Fig2A_boxplot_GSE76293.pdf", w3)
hr = pdf_h(f"{P}/Fig2E1_ROC_GSE76293.pdf", w3)
wF = (178 - 2 * MARGIN_X - COL_GAP) / 2 - 6
hF = pdf_h(f"{P}/Fig2F1_ML_AUC_barplot.pdf", wF)
H = TOP + hb + ROW_GAP + max(hb, hr) + ROW_GAP + hr + ROW_GAP + hF + 6
s = Sheet(H)
y = TOP
s.pdf_panel(f"{P}/Fig2A_boxplot_GSE76293.pdf", XA, y, w3, "A")
s.pdf_panel(f"{P}/Fig2B_boxplot_GSE243066.pdf", XB, y, w3, "B")
s.pdf_panel(f"{P}/Fig2C_boxplot_GSE69528.pdf", XC, y, w3, "C")
y += hb + ROW_GAP
s.pdf_panel(f"{P}/Fig2D_boxplot_GSE100159.pdf", XA, y, w3, "D")
s.pdf_panel(f"{P}/Fig2E1_ROC_GSE76293.pdf", XB, y, w3, "E")
s.pdf_panel(f"{P}/Fig2E2_ROC_GSE243066.pdf", XC, y, w3)
y += max(hb, hr) + ROW_GAP
s.pdf_panel(f"{P}/Fig2E3_ROC_GSE69528.pdf", XB, y, w3)
s.pdf_panel(f"{P}/Fig2E4_ROC_GSE100159.pdf", XC, y, w3)
y += hr + ROW_GAP
s.pdf_panel(f"{P}/Fig2F1_ML_AUC_barplot.pdf", X1, y, wF, "F")
s.pdf_panel(f"{P}/Fig2F2_SHAP_barplot.pdf", X2, y, wF)
s.save("Figure2_composite.pdf")

## ---------------- Figure 3 ----------------
hA3 = pdf_h(f"{P}/Fig3A_PPI_network.pdf", 84)
hB3 = pdf_h(f"{P}/Fig3B_GO_BP_dotplot.pdf", 74)
hC3 = pdf_h(f"{P}/Fig3C_GO_sankey_bubble_179.pdf", 100)
hD3 = pdf_h(f"{P}/Fig3D_KEGG_chord_179.pdf", 60)
hE3 = pdf_h(f"{P}/Fig3E_GO_sankey_bubble_6hub.pdf", 96)
hF3 = pdf_h(f"{P}/Fig3F_KEGG_chord_6hub.pdf", 60)
H = TOP + max(hA3, hB3) + ROW_GAP + max(hC3, hD3) + ROW_GAP + max(hE3, hF3) + 6
s = Sheet(H)
y = TOP
s.pdf_panel(f"{P}/Fig3A_PPI_network.pdf", X1, y, 84, "A")
s.pdf_panel(f"{P}/Fig3B_GO_BP_dotplot.pdf", X1 + 88, y, 74, "B")
y += max(hA3, hB3) + ROW_GAP
s.pdf_panel(f"{P}/Fig3C_GO_sankey_bubble_179.pdf", X1, y, 100, "C")
s.pdf_panel(f"{P}/Fig3D_KEGG_chord_179.pdf", X1 + 104, y, 60, "D")
y += max(hC3, hD3) + ROW_GAP
s.pdf_panel(f"{P}/Fig3E_GO_sankey_bubble_6hub.pdf", X1, y, 96, "E")
s.pdf_panel(f"{P}/Fig3F_KEGG_chord_6hub.pdf", X1 + 100, y, 60, "F")
s.save("Figure3_composite.pdf")

## ---------------- Figure 4 ----------------
hN = pdf_h(f"{P}/Fig4A_nomogram_ARDS.pdf", CW)
hS = pdf_h(f"{P}/Fig4C_ROC_cal_DCA_ARDS.pdf", 162)
H = TOP + hN + ROW_GAP + hS + ROW_GAP + hS + 6
s = Sheet(H)
y = TOP
s.pdf_panel(f"{P}/Fig4A_nomogram_ARDS.pdf", X1, y, CW, "A")
s.pdf_panel(f"{P}/Fig4B_nomogram_Sepsis.pdf", X2, y, CW, "B")
y += hN + ROW_GAP
s.pdf_panel(f"{P}/Fig4C_ROC_cal_DCA_ARDS.pdf", X1, y, 162, "C")
y += hS + ROW_GAP
s.pdf_panel(f"{P}/Fig4D_ROC_cal_DCA_sepsis.pdf", X1, y, 162, "D")
s.save("Figure4_composite.pdf")

## ---------------- Figure 5 ----------------
hBx = pdf_h(f"{P}/Fig5A_ssGSEA_boxplot_GSE76293.pdf", 96)
hHm = pdf_h(f"{P}/Fig5C_corr_heatmap_GSE76293.pdf", 62)
H = TOP + hBx + ROW_GAP + hBx + 6
s = Sheet(H)
s.pdf_panel(f"{P}/Fig5A_ssGSEA_boxplot_GSE76293.pdf", X1, TOP, 96, "A")
s.pdf_panel(f"{P}/Fig5B_ssGSEA_boxplot_GSE69528.pdf", X1, TOP + hBx + ROW_GAP, 96, "B")
s.pdf_panel(f"{P}/Fig5C_corr_heatmap_GSE76293.pdf", X1 + 100, TOP, 62, "C")
s.pdf_panel(f"{P}/Fig5D_corr_heatmap_GSE69528.pdf", X1 + 100, TOP + hHm + ROW_GAP, 62, "D")
s.save("Figure5_composite.pdf")

## ---------------- Figure 6 ----------------
hG = pdf_h(f"{P}/Fig6A_GSVA_top20_GSE76293.pdf", CW)
H = TOP + hG + 6
s = Sheet(H)
s.pdf_panel(f"{P}/Fig6A_GSVA_top20_GSE76293.pdf", X1, TOP, CW, "A")
s.pdf_panel(f"{P}/Fig6B_GSVA_top20_GSE69528.pdf", X2, TOP, CW, "B")
s.save("Figure6_composite.pdf")

## ---------------- Figure 7 ----------------
hU = pdf_h(f"{P}/Fig7A_UMAP_GSE151263.pdf", 76)
hCC = pdf_h(f"{SC}/cellchat/diff_heatmap.pdf", 84)
hT = pdf_h(f"{P}/Fig7C_trajectory_AUCell_GSE151263.pdf", 162)
hO = pdf_h(f"{P}/Fig7D_CellOracle_GSE151263.pdf", 84)
H = (TOP + max(hU, hCC) + ROW_GAP + hT + ROW_GAP + hO + ROW_GAP
     + max(hU, hCC) + ROW_GAP + hT + ROW_GAP + hO + 6)
s = Sheet(H)
y = TOP
s.pdf_panel(f"{P}/Fig7A_UMAP_GSE151263.pdf", X1, y, 76, "A")
s.pdf_panel(f"{SC}/cellchat/diff_heatmap.pdf", X1 + 80, y, 84, "B")
y += max(hU, hCC) + ROW_GAP
s.pdf_panel(f"{P}/Fig7C_trajectory_AUCell_GSE151263.pdf", X1, y, 162, "C")
y += hT + ROW_GAP
s.pdf_panel(f"{P}/Fig7D_CellOracle_GSE151263.pdf", X1, y, 84, "D")
y += hO + ROW_GAP
s.pdf_panel(f"{P}/Fig7E_UMAP_GSE167363.pdf", X1, y, 76, "E")
s.pdf_panel(f"{SC}/GSE167363/cellchat/diff_heatmap.pdf", X1 + 80, y, 84, "F")
y += max(hU, hCC) + ROW_GAP
s.pdf_panel(f"{P}/Fig7G_trajectory_AUCell_GSE167363.pdf", X1, y, 162, "G")
y += hT + ROW_GAP
s.pdf_panel(f"{P}/Fig7H_CellOracle_GSE167363.pdf", X1, y, 84, "H")
s.save("Figure7_composite.pdf")

## ---------------- Figure 8 ----------------
hDk = pdf_h(f"{P}/Fig8A_docking_site_heatmap.pdf", 100)
p1 = f"{DK}/poses_pymol/VNN1_luteolin_pymol.png"
p2 = f"{DK}/poses_pymol/GRB10_EGCG_pymol.png"
hPose = png_h(p1, 31)
hGc = pdf_h(f"{P}/Fig8C_glucocorticoid_repositioning.pdf", 92)
H = TOP + max(hDk, hPose) + ROW_GAP + hGc + 6
s = Sheet(H)
y = TOP
s.pdf_panel(f"{P}/Fig8A_docking_site_heatmap.pdf", X1, y, 100, "A")
s.png_panel(p1, X1 + 104, y, 31, "B")
s.png_panel(p2, X1 + 137, y, 31)
y += max(hDk, hPose) + ROW_GAP
s.pdf_panel(f"{P}/Fig8C_glucocorticoid_repositioning.pdf", X1, y, 92, "C")
s.save("Figure8_composite.pdf")

print("All composites done.")
