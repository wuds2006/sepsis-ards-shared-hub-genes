#!/usr/bin/env python3
"""Assemble supplementary Figure S1-S14 composites.

New ggplot2 NPG panels are combined with original pipeline renders where the
underlying per-level data were not saved (CellChat S5, running-ES GSEA S8,
monocle3 trajectory graphs S11, WGCNA dendrograms, MR SNP-level plots).
"""
import pymupdf

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
SP = f"{BASE}/AJRCCM_submission/03_figures/Figures_ggplot2_NPG/supp_panels"
OUT = f"{BASE}/AJRCCM_submission/03_figures/Figures_ggplot2_NPG/composites"
ORIG = f"{BASE}/AJRCCM_submission/05_supplementary/Supplementary_Figures_Tables/Figures"
WG = f"{BASE}/02_WGCNA/figS1_parts"
MR = f"{BASE}/04_MR_SMR"
SC = f"{BASE}/08_singlecell"

MM = 72 / 25.4
PAGE_W = 178 * MM
MX = 8.0
TOP = 8.0
RG = 5.0
CG = 4.0


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

    def pan(self, path, x_mm, y_mm, w_mm, letter=None):
        src = self._pdf(path)
        r = src[0].rect
        h_mm = w_mm * r.height / r.width
        dest = pymupdf.Rect(mm(x_mm), mm(y_mm), mm(x_mm + w_mm), mm(y_mm + h_mm))
        self.page.show_pdf_page(dest, src, 0)
        if letter:
            self.letter(x_mm - 4.5, y_mm, letter)
        return h_mm

    def save(self, name):
        self.doc.save(f"{OUT}/{name}", deflate=True, garbage=3)
        for d in self.cache.values():
            d.close()
        self.doc.close()
        print("saved", name)


CW = (178 - 2 * MX - CG) / 2
X1 = MX
X2 = MX + CW + CG

## S1: soft-threshold (new) + dendrograms (original)
h1 = pdf_h(f"{SP}/FigS1A_soft_threshold_GSE76293.pdf", 140)
h2 = pdf_h(f"{WG}/c_dendro_GSE76293.pdf", CW)
s = Sheet(TOP + h1 + RG + h2 + RG + h2 + 6)
s.pan(f"{SP}/FigS1A_soft_threshold_GSE76293.pdf", X1 + 10, TOP, 140, "A")
y = TOP + h1 + RG
s.pan(f"{WG}/c_dendro_GSE76293.pdf", X1, y, CW, "B")
s.pan(f"{WG}/f_dendro_GSE69528.pdf", X2, y, CW, "C")
s.save("FigureS1_composite.pdf")

## S2: OOF ROC + confusion (new)
hA = pdf_h(f"{SP}/FigS2A_OOF_ROC_GSE76293.pdf", CW)
hC = pdf_h(f"{SP}/FigS2C_confusion_GSE76293.pdf", 70)
s = Sheet(TOP + hA + RG + hC + 6)
s.pan(f"{SP}/FigS2A_OOF_ROC_GSE76293.pdf", X1, TOP, CW, "A")
s.pan(f"{SP}/FigS2B_OOF_ROC_GSE69528.pdf", X2, TOP, CW, "B")
y = TOP + hA + RG
s.pan(f"{SP}/FigS2C_confusion_GSE76293.pdf", X1 + 10, y, 70, "C")
s.pan(f"{SP}/FigS2D_confusion_GSE69528.pdf", X2 + 10, y, 70, "D")
s.save("FigureS2_composite.pdf")

## S3: MR method forest (new) + SNP-level scatter/forest (original renders)
hF = pdf_h(f"{SP}/FigS3_MR_method_forest.pdf", 120)
hS = pdf_h(f"{MR}/MR_sepsis4980_scatter.pdf", CW)
s = Sheet(TOP + hF + RG + hS + 6)
s.pan(f"{SP}/FigS3_MR_method_forest.pdf", X1 + 20, TOP, 120, "A")
y = TOP + hF + RG
s.pan(f"{MR}/MR_sepsis4980_scatter.pdf", X1, y, CW, "B")
s.pan(f"{MR}/MR_sepsis4980_forest.pdf", X2, y, CW, "C")
s.save("FigureS3_composite.pdf")

## S4: ssGSEA heatmap GSE69528 (new)
s = Sheet(TOP + pdf_h(f"{SP}/FigS4_ssGSEA_heatmap_GSE69528.pdf", 150) + 6)
s.pan(f"{SP}/FigS4_ssGSEA_heatmap_GSE69528.pdf", X1 + 5, TOP, 150, None)
s.save("FigureS4_composite.pdf")

## S5: CellChat pathway comparison (kept original render)
s = Sheet(TOP + pdf_h(f"{ORIG}/FigureS5_CellChat_pathway_comparison.pdf", 162) + 6)
s.pan(f"{ORIG}/FigureS5_CellChat_pathway_comparison.pdf", X1, TOP, 162, None)
s.save("FigureS5_composite.pdf")

## S6: KO top10 (new) + stratified GO_BP (original)
hK = pdf_h(f"{SP}/FigS6A_KO_UPP1_top10.pdf", CW)
hG = pdf_h(f"{SC}/knockout_stratified/GSE151263_KO_stratified_GO_BP_UPP1.pdf", CW)
s = Sheet(TOP + hK + RG + hG + 6)
s.pan(f"{SP}/FigS6A_KO_UPP1_top10.pdf", X1, TOP, CW, "A")
s.pan(f"{SP}/FigS6B_KO_EXOSC4_top10.pdf", X2, TOP, CW, "B")
y = TOP + hK + RG
s.pan(f"{SC}/knockout_stratified/GSE151263_KO_stratified_GO_BP_UPP1.pdf", X1, y, CW, "C")
s.pan(f"{SC}/knockout_stratified/GSE151263_KO_stratified_GO_BP_EXOSC4.pdf", X2, y, CW, "D")
s.save("FigureS6_composite.pdf")

## S7: HPA infection logFC (new)
s = Sheet(TOP + pdf_h(f"{SP}/FigS7_HPA_infection_logFC.pdf", 150) + 6)
s.pan(f"{SP}/FigS7_HPA_infection_logFC.pdf", X1 + 5, TOP, 150, None)
s.save("FigureS7_composite.pdf")

## S8: single-gene GSEA (kept original render)
s = Sheet(TOP + pdf_h(f"{ORIG}/FigureS8_single_gene_gsea.pdf", 162) + 6)
s.pan(f"{ORIG}/FigureS8_single_gene_gsea.pdf", X1, TOP, 162, None)
s.save("FigureS8_composite.pdf")

## S9: single-cell atlas GSE167363 (new hub heatmap + original marker dotplot)
hH = pdf_h(f"{SP}/FigS9_hub_expr_heatmap_GSE167363.pdf", 130)
hD = pdf_h(f"{SC}/GSE167363/GSE167363_marker_dotplot.pdf", 120)
s = Sheet(TOP + hH + RG + hD + 6)
s.pan(f"{SP}/FigS9_hub_expr_heatmap_GSE167363.pdf", X1 + 15, TOP, 130, "A")
y = TOP + hH + RG
s.pan(f"{SC}/GSE167363/GSE167363_marker_dotplot.pdf", X1 + 20, y, 120, "B")
s.save("FigureS9_composite.pdf")

## S10: AUCell heatmap + CellOracle TF bars (new)
hA10 = pdf_h(f"{SP}/FigS10A_AUCell_heatmap_GSE167363.pdf", 140)
hB10 = pdf_h(f"{SP}/FigS10B_CellOracle_TF_GSE167363.pdf", 162)
s = Sheet(TOP + hA10 + RG + hB10 + 6)
s.pan(f"{SP}/FigS10A_AUCell_heatmap_GSE167363.pdf", X1 + 10, TOP, 140, "A")
y = TOP + hA10 + RG
s.pan(f"{SP}/FigS10B_CellOracle_TF_GSE167363.pdf", X1, y, 162, "B")
s.save("FigureS10_composite.pdf")

## S11: monocle3 trajectories (kept original render)
s = Sheet(TOP + pdf_h(f"{ORIG}/FigureS11_monocle3_trajectories.pdf", 162) + 6)
s.pan(f"{ORIG}/FigureS11_monocle3_trajectories.pdf", X1, TOP, 162, None)
s.save("FigureS11_composite.pdf")

## S12: metabolic delta (new) + hub correlation (original)
hM = pdf_h(f"{SP}/FigS12_metabolic_delta.pdf", 120)
hC12 = pdf_h(f"{SC}/metabolism/GSE151263_metab_hub_correlation.pdf", CW)
s = Sheet(TOP + hM + RG + hC12 + 6)
s.pan(f"{SP}/FigS12_metabolic_delta.pdf", X1 + 20, TOP, 120, "A")
y = TOP + hM + RG
s.pan(f"{SC}/metabolism/GSE151263_metab_hub_correlation.pdf", X1, y, CW, "B")
s.pan(f"{SC}/metabolism/GSE167363_metab_hub_correlation.pdf", X2, y, CW, "C")
s.save("FigureS12_composite.pdf")

## S13: stratified KO (new A,B) + stratified top10 (original)
hA13 = pdf_h(f"{SP}/FigS13A_KO_concordance.pdf", 110)
hB13 = pdf_h(f"{SP}/FigS13B_CEBPB_UPP1_gradient.pdf", 90)
hT13 = pdf_h(f"{SC}/knockout_stratified/GSE151263_KO_stratified_top10_UPP1.pdf", CW)
s = Sheet(TOP + max(hA13, hB13) + RG + hT13 + 6)
s.pan(f"{SP}/FigS13A_KO_concordance.pdf", X1, TOP, 110, "A")
s.pan(f"{SP}/FigS13B_CEBPB_UPP1_gradient.pdf", X1 + 116, TOP, 56, "B")
y = TOP + max(hA13, hB13) + RG
s.pan(f"{SC}/knockout_stratified/GSE151263_KO_stratified_top10_UPP1.pdf", X1, y, CW, "C")
s.pan(f"{SC}/knockout_stratified/GSE151263_KO_stratified_top10_EXOSC4.pdf", X2, y, CW, "D")
s.save("FigureS13_composite.pdf")

## S14: donor sensitivity (new)
s = Sheet(TOP + pdf_h(f"{SP}/FigS14_donor_sensitivity.pdf", 162) + 6)
s.pan(f"{SP}/FigS14_donor_sensitivity.pdf", X1, TOP, 162, None)
s.save("FigureS14_composite.pdf")

print("All supplementary composites done.")
