#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 87 —— RESULTS 图组 composite 拼装（10 张主图，对应 RESULTS 1-10）
输出: AJRCCM_submission/03_figures_RESULTS_NPG/composites/
      FigureN_composite.{pdf,png(300dpi),tiff(1000dpi LZW)}
布局单位：全程 pt（PyMuPDF 原生），杜绝 mm/pt 混用。
"""
import os, io
import pymupdf
from PIL import Image

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission"
P  = f"{BASE}/03_figures/Figures_ggplot2_NPG/main_panels"
NP = f"{BASE}/03_figures_RESULTS_NPG/panels"
OUT = f"{BASE}/03_figures_RESULTS_NPG/composites"
SC = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/08_singlecell"
DOCK = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/09_docking_drug"

MM = 72 / 25.4
PAGE_W = 178 * MM          # 504.57 pt
MARGIN = 8 * MM
CW = PAGE_W - 2 * MARGIN   # 459 pt
GAP_X = 4 * MM
GAP_Y = 5 * MM
TOP = 13 * MM
os.makedirs(OUT, exist_ok=True)

def pdf_size(path):
    d = pymupdf.open(path); r = d[0].rect; d.close()
    return r.width, r.height

def aspect_h(path, w):
    pw, ph = pdf_size(path)
    return w * ph / pw

class Sheet:
    def __init__(self, height_pt):
        self.doc = pymupdf.open()
        self.page = self.doc.new_page(width=PAGE_W, height=height_pt)
    def letter(self, x, y, ch):
        self.page.insert_text((x - 10, y + 3.4), ch, fontsize=9,
                              fontname="hebo", color=(0, 0, 0))
    def pdf_panel(self, path, x, y, w, ch):
        src = pymupdf.open(path)
        pw, ph = pdf_size(path)
        h = w * ph / pw
        self.page.show_pdf_page(pymupdf.Rect(x, y, x + w, y + h), src, 0)
        src.close()
        if ch: self.letter(x, y, ch)
        return h
    def img_panel(self, path, x, y, w, ch):
        im = Image.open(path)
        h = w * im.size[1] / im.size[0]
        self.page.insert_image(pymupdf.Rect(x, y, x + w, y + h), filename=path)
        if ch: self.letter(x, y, ch)
        return h
    def save(self, name):
        out = os.path.join(OUT, f"{name}_composite.pdf")
        self.doc.save(out, deflate=True)
        self.doc.close()
        return out

def export_rasters(pdf_path, name):
    d = pymupdf.open(pdf_path)
    pg = d[0]
    pix = pg.get_pixmap(matrix=pymupdf.Matrix(300/72, 300/72), alpha=False)
    Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save(
        os.path.join(OUT, f"{name}_composite.png"), dpi=(300, 300))
    pix = pg.get_pixmap(matrix=pymupdf.Matrix(1000/72, 1000/72), alpha=False)
    Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save(
        os.path.join(OUT, f"{name}_composite.tiff"),
        format="TIFF", compression="tiff_lzw", dpi=(1000, 1000))
    d.close()

half = (CW - GAP_X) / 2

# ---------- Figure 1: shared candidate genes ----------
hB = aspect_h(f"{P}/Fig1B1_volcano_GSE76293.pdf", half)
hC = aspect_h(f"{P}/Fig1C1_module_trait_GSE76293.pdf", half)
D_W = 106 * MM
hD = aspect_h(f"{P}/Fig1D_venn_4sets.pdf", D_W)
H = TOP + hB + GAP_Y + hD + GAP_Y + hC + 6 * MM
s = Sheet(H)
y = TOP
s.pdf_panel(f"{P}/Fig1B1_volcano_GSE76293.pdf", MARGIN, y, half, "A")
s.pdf_panel(f"{P}/Fig1B2_volcano_GSE69528.pdf", MARGIN + half + GAP_X, y, half, "B")
y += hB + GAP_Y
xD = (PAGE_W - D_W) / 2
s.pdf_panel(f"{P}/Fig1D_venn_4sets.pdf", xD, y, D_W, "E")
y += hD + GAP_Y
s.pdf_panel(f"{P}/Fig1C1_module_trait_GSE76293.pdf", MARGIN, y, half, "C")
s.pdf_panel(f"{P}/Fig1C2_module_trait_GSE69528.pdf", MARGIN + half + GAP_X, y, half, "D")
out = s.save("Figure1"); export_rasters(out, "Figure1")
print("Figure1", round(H/72, 2), "in")

# ---------- Figure 2: GO/KEGG of 179 ----------
hGO = aspect_h(f"{P}/Fig3B_GO_BP_dotplot.pdf", 118 * MM)
hSK = aspect_h(f"{P}/Fig3C_GO_sankey_bubble_179.pdf", CW)
hCH = aspect_h(f"{P}/Fig3D_KEGG_chord_179.pdf", 120 * MM)
H = TOP + hGO + GAP_Y + hSK + GAP_Y + hCH + 6 * MM
s = Sheet(H)
y = TOP
s.pdf_panel(f"{P}/Fig3B_GO_BP_dotplot.pdf", (PAGE_W - 118*MM)/2, y, 118 * MM, "A")
y += hGO + GAP_Y
s.pdf_panel(f"{P}/Fig3C_GO_sankey_bubble_179.pdf", MARGIN, y, CW, "B")
y += hSK + GAP_Y
s.pdf_panel(f"{P}/Fig3D_KEGG_chord_179.pdf", (PAGE_W - 120*MM)/2, y, 120 * MM, "C")
out = s.save("Figure2"); export_rasters(out, "Figure2")
print("Figure2", round(H/72, 2), "in")

# ---------- Figure 3: ML prioritization ----------
hML = aspect_h(f"{P}/Fig2F1_ML_AUC_barplot.pdf", half)
hSH = aspect_h(f"{P}/Fig2F2_SHAP_barplot.pdf", half)
hX  = aspect_h(f"{NP}/R3C_crossdisease_AUC.pdf", CW)
hPI = aspect_h(f"{P}/Fig3A_PPI_network.pdf", 100 * MM)
H = TOP + max(hML, hSH) + GAP_Y + hX + GAP_Y + hPI + 6 * MM
s = Sheet(H)
y = TOP
s.pdf_panel(f"{P}/Fig2F1_ML_AUC_barplot.pdf", MARGIN, y, half, "A")
s.pdf_panel(f"{P}/Fig2F2_SHAP_barplot.pdf", MARGIN + half + GAP_X, y, half, "B")
y += max(hML, hSH) + GAP_Y
s.pdf_panel(f"{NP}/R3C_crossdisease_AUC.pdf", MARGIN, y, CW, "C")
y += hX + GAP_Y
s.pdf_panel(f"{P}/Fig3A_PPI_network.pdf", (PAGE_W - 100*MM)/2, y, 100 * MM, "D")
out = s.save("Figure3"); export_rasters(out, "Figure3")
print("Figure3", round(H/72, 2), "in")

# ---------- Figure 4: expression + single-gene GSEA ----------
hBX = aspect_h(f"{P}/Fig2A_boxplot_GSE76293.pdf", half)
hGS = aspect_h(f"{NP}/R4E_gsea_NES_heatmap.pdf", CW)
H = TOP + hBX + GAP_Y + hBX + GAP_Y + hGS + 6 * MM
s = Sheet(H)
y = TOP
s.pdf_panel(f"{P}/Fig2A_boxplot_GSE76293.pdf", MARGIN, y, half, "A")
s.pdf_panel(f"{P}/Fig2B_boxplot_GSE243066.pdf", MARGIN + half + GAP_X, y, half, "B")
y += hBX + GAP_Y
s.pdf_panel(f"{P}/Fig2C_boxplot_GSE69528.pdf", MARGIN, y, half, "C")
s.pdf_panel(f"{P}/Fig2D_boxplot_GSE100159.pdf", MARGIN + half + GAP_X, y, half, "D")
y += hBX + GAP_Y
s.pdf_panel(f"{NP}/R4E_gsea_NES_heatmap.pdf", MARGIN, y, CW, "E")
out = s.save("Figure4"); export_rasters(out, "Figure4")
print("Figure4", round(H/72, 2), "in")

# ---------- Figure 5: disease prediction models ----------
hN = aspect_h(f"{P}/Fig4A_nomogram_ARDS.pdf", half)
hRC = aspect_h(f"{P}/Fig4C_ROC_cal_DCA_ARDS.pdf", 162 * MM)
H = TOP + hN + GAP_Y + hRC + GAP_Y + hRC + 6 * MM
s = Sheet(H)
y = TOP
s.pdf_panel(f"{P}/Fig4A_nomogram_ARDS.pdf", MARGIN, y, half, "A")
s.pdf_panel(f"{P}/Fig4B_nomogram_Sepsis.pdf", MARGIN + half + GAP_X, y, half, "B")
y += hN + GAP_Y
s.pdf_panel(f"{P}/Fig4C_ROC_cal_DCA_ARDS.pdf", MARGIN, y, 162 * MM, "C")
y += hRC + GAP_Y
s.pdf_panel(f"{P}/Fig4D_ROC_cal_DCA_sepsis.pdf", MARGIN, y, 162 * MM, "D")
out = s.save("Figure5"); export_rasters(out, "Figure5")
print("Figure5", round(H/72, 2), "in")

# ---------- Figure 6: external confirmation ----------
h6A = aspect_h(f"{NP}/R6A_validation_AUC.pdf", half)
h6B = aspect_h(f"{NP}/R6B_meta_forest.pdf", half)
hRO = aspect_h(f"{P}/Fig2E2_ROC_GSE243066.pdf", half)
H = TOP + max(h6A, h6B) + GAP_Y + hRO + 6 * MM
s = Sheet(H)
y = TOP
s.pdf_panel(f"{NP}/R6A_validation_AUC.pdf", MARGIN, y, half, "A")
s.pdf_panel(f"{NP}/R6B_meta_forest.pdf", MARGIN + half + GAP_X, y, half, "B")
y += max(h6A, h6B) + GAP_Y
s.pdf_panel(f"{P}/Fig2E2_ROC_GSE243066.pdf", MARGIN, y, half, "C")
s.pdf_panel(f"{P}/Fig2E4_ROC_GSE100159.pdf", MARGIN + half + GAP_X, y, half, "D")
out = s.save("Figure6"); export_rasters(out, "Figure6")
print("Figure6", round(H/72, 2), "in")

# ---------- Figure 7: MR null under limited power ----------
h7A = aspect_h(f"{NP}/R7A_mr_forest_forward.pdf", half)
h7B = aspect_h(f"{NP}/R7B_mr_forest_reverse.pdf", half)
h7C = aspect_h(f"{NP}/R7C_mr_power.pdf", half)
h7D = aspect_h(f"{NP}/R7D_smr_forest.pdf", half)
H = TOP + max(h7A, h7B) + GAP_Y + max(h7C, h7D) + 6 * MM
s = Sheet(H)
y = TOP
s.pdf_panel(f"{NP}/R7A_mr_forest_forward.pdf", MARGIN, y, half, "A")
s.pdf_panel(f"{NP}/R7B_mr_forest_reverse.pdf", MARGIN + half + GAP_X, y, half, "B")
y += max(h7A, h7B) + GAP_Y
s.pdf_panel(f"{NP}/R7C_mr_power.pdf", MARGIN, y, half, "C")
s.pdf_panel(f"{NP}/R7D_smr_forest.pdf", MARGIN + half + GAP_X, y, half, "D")
out = s.save("Figure7"); export_rasters(out, "Figure7")
print("Figure7", round(H/72, 2), "in")

# ---------- Figure 8: immune microenvironment ----------
W3 = (CW - 2 * GAP_X) / 3
hGV = aspect_h(f"{P}/Fig5A_hub6_GSVA.pdf", half)
hCB = aspect_h(f"{P}/Fig5B_CIBERSORT_stacked_GSE76293.pdf", half)
hCO = aspect_h(f"{P}/Fig5F_spearman_corrplot.pdf", W3)
hT20 = aspect_h(f"{P}/Fig6A_GSVA_top20_GSE76293.pdf", (CW - GAP_X) / 2)
H = TOP + hGV + GAP_Y + hCB + GAP_Y + hCB + GAP_Y + max(hCO, hT20) + 6 * MM
s = Sheet(H)
y = TOP
s.pdf_panel(f"{P}/Fig5A_hub6_GSVA.pdf", MARGIN, y, half, "A")
# top-20 GSVA 放 A 右侧（complement 共享通路证据）
s.pdf_panel(f"{P}/Fig6A_GSVA_top20_GSE76293.pdf", MARGIN + half + GAP_X, y, half, "G")
y += hGV + GAP_Y
s.pdf_panel(f"{P}/Fig5B_CIBERSORT_stacked_GSE76293.pdf", MARGIN, y, half, "B")
s.pdf_panel(f"{P}/Fig5C_CIBERSORT_stacked_GSE69528.pdf", MARGIN + half + GAP_X, y, half, "C")
y += hCB + GAP_Y
s.pdf_panel(f"{P}/Fig5D_CIBERSORT_stacked_GSE243066.pdf", MARGIN, y, half, "D")
s.pdf_panel(f"{P}/Fig5E_CIBERSORT_stacked_GSE100159.pdf", MARGIN + half + GAP_X, y, half, "E")
y += hCB + GAP_Y
s.pdf_panel(f"{P}/Fig5F_spearman_corrplot.pdf", MARGIN, y, half, "F")
s.pdf_panel(f"{P}/Fig6B_GSVA_top20_GSE69528.pdf", MARGIN + half + GAP_X, y, half, "H")
out = s.save("Figure8"); export_rasters(out, "Figure8")
print("Figure8", round(H/72, 2), "in")

# ---------- Figure 9: single-cell ----------
hU = aspect_h(f"{P}/Fig7A_UMAP_GSE151263.pdf", 76 * MM)
hCC = aspect_h(f"{SC}/cellchat/diff_heatmap.pdf", 84 * MM)
hT = aspect_h(f"{P}/Fig7C_trajectory_AUCell_GSE151263.pdf", 162 * MM)
hO = aspect_h(f"{P}/Fig7D_CellOracle_GSE151263.pdf", 84 * MM)
H = (TOP + max(hU, hCC) + GAP_Y + hT + GAP_Y + hO + GAP_Y
     + max(hU, hCC) + GAP_Y + hT + GAP_Y + hO + 6 * MM)
s = Sheet(H)
y = TOP
s.pdf_panel(f"{P}/Fig7A_UMAP_GSE151263.pdf", MARGIN, y, 76 * MM, "A")
s.pdf_panel(f"{SC}/cellchat/diff_heatmap.pdf", MARGIN + 80 * MM, y, 84 * MM, "B")
y += max(hU, hCC) + GAP_Y
s.pdf_panel(f"{P}/Fig7C_trajectory_AUCell_GSE151263.pdf", MARGIN, y, 162 * MM, "C")
y += hT + GAP_Y
s.pdf_panel(f"{P}/Fig7D_CellOracle_GSE151263.pdf", MARGIN, y, 84 * MM, "D")
y += hO + GAP_Y
s.pdf_panel(f"{P}/Fig7E_UMAP_GSE167363.pdf", MARGIN, y, 76 * MM, "E")
s.pdf_panel(f"{SC}/GSE167363/cellchat/diff_heatmap.pdf", MARGIN + 80 * MM, y, 84 * MM, "F")
y += max(hU, hCC) + GAP_Y
s.pdf_panel(f"{P}/Fig7G_trajectory_AUCell_GSE167363.pdf", MARGIN, y, 162 * MM, "G")
y += hT + GAP_Y
s.pdf_panel(f"{P}/Fig7H_CellOracle_GSE167363.pdf", MARGIN, y, 84 * MM, "H")
out = s.save("Figure9"); export_rasters(out, "Figure9")
print("Figure9", round(H/72, 2), "in")

# ---------- Figure 10: protein validation + therapeutic ----------
hHP = aspect_h(f"{NP}/R10A_hpa_infection.pdf", half)
hCM = aspect_h(f"{P}/Fig8C_glucocorticoid_repositioning.pdf", half)
hDK = aspect_h(f"{P}/Fig8A_docking_site_heatmap.pdf", half)
pose_png = f"{DOCK}/poses_pymol/VNN1_luteolin_pymol_flat.png"
hPS = None
im = Image.open(pose_png); hPS = half * im.size[1] / im.size[0]
hPL = aspect_h(f"{DOCK}/plip_vnn1_luteolin_2D.pdf", half)
H = TOP + max(hHP, hCM) + GAP_Y + hDK + GAP_Y + max(hPS, hPL) + 6 * MM
s = Sheet(H)
y = TOP
s.pdf_panel(f"{NP}/R10A_hpa_infection.pdf", MARGIN, y, half, "A")
s.pdf_panel(f"{P}/Fig8C_glucocorticoid_repositioning.pdf", MARGIN + half + GAP_X, y, half, "B")
y += max(hHP, hCM) + GAP_Y
s.pdf_panel(f"{P}/Fig8A_docking_site_heatmap.pdf", MARGIN, y, half, "C")
s.img_panel(pose_png, MARGIN + half + GAP_X, y, half, "D")
y += max(hDK, hPS) + GAP_Y
s.pdf_panel(f"{DOCK}/plip_vnn1_luteolin_2D.pdf", MARGIN, y, half, "E")
out = s.save("Figure10"); export_rasters(out, "Figure10")
print("Figure10", round(H/72, 2), "in")

print("all composites done ->", OUT)
