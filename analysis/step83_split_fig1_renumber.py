#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 83 —— 将 Figure1（GA + B/C/D）拆分为两张图并全链重编号：
  新 Figure 1 = Graphical abstract 单图（全宽）
  新 Figure 2 = 旧 Figure 1 的 B（火山图）/ D（Venn，居中）/ C（WGCNA）
  旧 Figure 2-9 → 新 Figure 3-10（文件名、图注、正文引用同步后移）
产出：
  Figures_final/Figure{1,2}_composite.{pdf,png,tiff}
  Figures_final/Figure{3..10}_composite.* （自旧 2-9 重命名）
  Figure_legends.md / manuscript.md 引用更新
"""
import os, re, shutil, glob
import pymupdf
from PIL import Image

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission"
FIN = f"{BASE}/03_figures/Figures_final"
P = f"{BASE}/03_figures/Figures_ggplot2_NPG/main_panels"
GA = f"{BASE}/04_graphical_abstract/graphical_abstract.pdf"
MM = 72 / 25.4
PAGE_W = 178 * MM
MARGIN = 23.0
GAP = 12.0
CW = PAGE_W - 2 * MARGIN

def raster(pdf_path, tag):
    d = pymupdf.open(pdf_path)
    pg = d[0]
    pix = pg.get_pixmap(matrix=pymupdf.Matrix(300 / 72, 300 / 72), alpha=False)
    Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save(
        f"{FIN}/{tag}_composite.png", dpi=(300, 300))
    dpi = 1000
    pix = pg.get_pixmap(matrix=pymupdf.Matrix(dpi / 72, dpi / 72), alpha=False)
    img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    img.save(f"{FIN}/{tag}_composite.tiff", format="TIFF",
             compression="tiff_lzw", dpi=(dpi, dpi))
    d.close()
    for ext in ["pdf", "png", "tiff"]:
        print(f"  {tag}.{ext}", round(os.path.getsize(f"{FIN}/{tag}_composite.{ext}") / 1e6, 2), "MB")

def ph(path, w):
    d = pymupdf.open(path); r = d[0].rect; d.close()
    return w * r.height / r.width

# ---------- 1) 旧 Figure 2-9 → 3-10 重命名（降序避免覆盖） ----------
for n in range(9, 1, -1):
    for ext in ["pdf", "png", "tiff"]:
        src = f"{FIN}/Figure{n}_composite.{ext}"
        if os.path.exists(src):
            shutil.move(src, f"{FIN}/Figure{n+1}_composite.{ext}")
print("renamed Figure2-9 -> Figure3-10")

# ---------- 2) 新 Figure 1：GA 单图 ----------
shutil.copy2(f"{FIN}/Figure1_composite.pdf", "/tmp/Figure1_composite_pre_split.pdf.bak")
ga = pymupdf.open(GA)
gr = ga[0].rect
GA_W = CW
GA_H = GA_W * gr.height / gr.width
doc = pymupdf.open()
page = doc.new_page(width=PAGE_W, height=20 + GA_H + 18)
page.show_pdf_page(pymupdf.Rect(MARGIN, 20, MARGIN + GA_W, 20 + GA_H), ga, 0)
doc.save(f"{FIN}/Figure1_composite.pdf", deflate=True)
doc.close(); ga.close()
print(f"Figure1 (GA only): {PAGE_W/72:.2f} x {(20+GA_H+18)/72:.2f} in")
raster(f"{FIN}/Figure1_composite.pdf", "Figure1")

# ---------- 3) 新 Figure 2：B（上）/ D（中，居中）/ C（下） ----------
panels = {
    "B1": f"{P}/Fig1B1_volcano_GSE76293.pdf",
    "B2": f"{P}/Fig1B2_volcano_GSE69528.pdf",
    "C1": f"{P}/Fig1C1_module_trait_GSE76293.pdf",
    "C2": f"{P}/Fig1C2_module_trait_GSE69528.pdf",
    "D":  f"{P}/Fig1D_venn_4sets.pdf",
}
half = (CW - GAP) / 2
hB = ph(panels["B1"], half)
hC = ph(panels["C1"], half)
D_W = 300.0
hD = ph(panels["D"], D_W)
TOP = 36.0
PAGE_H = TOP + hB + GAP + hD + GAP + hC + 18
doc = pymupdf.open()
page = doc.new_page(width=PAGE_W, height=PAGE_H)
def letter(x, y, ch):
    page.insert_text((x - 12, y + 3.4), ch, fontsize=9,
                     fontname="hebo", color=(0, 0, 0))
# Row 1: B1 B2
y = TOP
page.show_pdf_page(pymupdf.Rect(MARGIN, y, MARGIN + half, y + hB),
                   pymupdf.open(panels["B1"]), 0)
page.show_pdf_page(pymupdf.Rect(MARGIN + half + GAP, y, PAGE_W - MARGIN, y + hB),
                   pymupdf.open(panels["B2"]), 0)
letter(MARGIN, y, "B")
# Row 2: D centered
y += hB + GAP
xD = (PAGE_W - D_W) / 2
page.show_pdf_page(pymupdf.Rect(xD, y, xD + D_W, y + hD),
                   pymupdf.open(panels["D"]), 0)
letter(xD, y, "D")
# Row 3: C1 C2
y += hD + GAP
page.show_pdf_page(pymupdf.Rect(MARGIN, y, MARGIN + half, y + hC),
                   pymupdf.open(panels["C1"]), 0)
page.show_pdf_page(pymupdf.Rect(MARGIN + half + GAP, y, PAGE_W - MARGIN, y + hC),
                   pymupdf.open(panels["C2"]), 0)
letter(MARGIN, y, "C")
doc.save(f"{FIN}/Figure2_composite.pdf", deflate=True)
doc.close()
print(f"Figure2 (B/D/C): {PAGE_W/72:.2f} x {PAGE_H/72:.2f} in | "
      f"B {half:.0f}x{hB:.0f} D {D_W:.0f}x{hD:.0f} C {half:.0f}x{hC:.0f}")
raster(f"{FIN}/Figure2_composite.pdf", "Figure2")

# ---------- 4) manuscript.md 引用重编号 ----------
MP = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission/01_manuscript/manuscript.md"
md = open(MP, encoding="utf-8").read()
assert "and 6B" in md and md.count("and 6B") == 1
counts = {}
for n in range(9, 1, -1):
    for word in ["Figure", "Figures"]:
        pat = rf"{word} {n}(?![0-9])"
        md, k = re.subn(pat, f"{word} {n+1}", md)
        counts[f"{word} {n}"] = counts.get(f"{word} {n}", 0) + k
md = md.replace("(Figure 1; Supplementary Table S5)", "(Figure 2; Supplementary Table S5)")
md = md.replace("and 6B", "and 7B")
open(MP, "w", encoding="utf-8").write(md)
print("manuscript.md renumbered:", {k: v for k, v in counts.items() if v})

# ---------- 5) Figure_legends.md 重写 ----------
LG = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/AJRCCM_submission/03_figures/Figure_legends.md"
lg = open(LG, encoding="utf-8").read()
# 去掉旧 Figure 1 段
m = re.search(r"Figure 1\..*?(?=Figure 2\.)", lg, re.S)
assert m
lg = lg[:m.start()] + lg[m.end():]
# 降序重编号标题
for n in range(9, 1, -1):
    lg = re.sub(rf"^Figure {n}\.", f"Figure {n+1}.", lg, flags=re.M)
new_legends = """Figure 1. Graphical abstract.
Graphical abstract summarizing the study design and main findings: identification of six shared hub genes from the two discovery cohorts by DEG\u2013WGCNA intersection and machine-learning prioritization with external validation, mechanistic exploration by immune-infiltration profiling, single-gene GSEA, and single-cell virtual-knockout analyses, and in silico drug prediction by CMap-based repositioning and molecular docking.

Figure 2. Identification of shared candidate genes.
(B) Volcano plots of differentially expressed genes (DEGs) in GSE76293 (ARDS vs. control; 557 DEGs: 291 upregulated, 266 downregulated) and GSE69528 (sepsis vs. control; 1,458 DEGs: 658 upregulated, 800 downregulated). Red/pink points indicate upregulated genes and blue points indicate downregulated genes (|log2 fold change| > 1 and adjusted P < 0.05). (D) Four-set Venn diagram showing the intersection of DEGs and WGCNA module genes, yielding 179 shared candidate genes. (C) Weighted gene co-expression network analysis (WGCNA) module\u2013trait relationships in the two discovery cohorts, identifying 13 (ARDS) and 11 (sepsis) disease-associated modules.

"""
open(LG, "w", encoding="utf-8").write("FIGURE LEGENDS\n\n" + new_legends + lg.replace("FIGURE LEGENDS\n\n", "", 1).lstrip("\n"))
print("Figure_legends.md rewritten")
print("done")
