#!/usr/bin/env python
"""Assemble composite figures 1-5 and 7 (Figure 6/8 composites already exist).
Layout follows the panel letters of 03_figures/Figure_legends.md.
Style: same as scripts/assemble_figure6.py (PIL, Nature palette, 300dpi triple export).
"""
from PIL import Image, ImageDraw, ImageFont
import os

ROOT = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
FX = os.path.join(ROOT, "AJRCCM_submission/07_figure_exports")

GAP = 50
MARGIN = 60
LETTER_SZ = 80
font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", LETTER_SZ)

def scaled(img, target_w):
    w, h = img.size
    return img.resize((target_w, int(round(h * target_w / w))), Image.LANCZOS)

def load(src, target_w):
    im = Image.open(src)
    if im.mode in ("RGBA", "LA", "P"):
        im = im.convert("RGBA")
        bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
        im = Image.alpha_composite(bg, im).convert("RGB")
    else:
        im = im.convert("RGB")
    return scaled(im, target_w)

def assemble(rows_spec, out_stem):
    """rows_spec: list of rows; each row = list of (PIL-image-or-list-of-images, letter).
    Cells within a letter that have multiple images are pasted side by side first."""
    # determine full width from widest row after composing
    composed = []
    for row in rows_spec:
        cells = []
        for imgs, letter, cw in row:
            if isinstance(imgs, list):
                parts = [load(p, (cw - 30 * (len(imgs) - 1)) // len(imgs)) for p in imgs]
                h = max(p.size[1] for p in parts)
                c = Image.new("RGB", (cw, h), "white")
                x = 0
                for p in parts:
                    c.paste(p, (x, (h - p.size[1]) // 2)); x += p.size[0] + 30
                cells.append((c, letter))
            else:
                cells.append((load(imgs, cw), letter))
        composed.append(cells)
    W = max(sum(c.size[0] for c, _ in row) + GAP * (len(row) - 1) for row in composed)
    Hh = MARGIN * 2 + sum(max(c.size[1] for c, _ in row) for row in composed) + GAP * (len(composed) - 1)
    canvas = Image.new("RGB", (W + MARGIN * 2, Hh), "white")
    d = ImageDraw.Draw(canvas)
    y = MARGIN
    for row in composed:
        h = max(c.size[1] for c, _ in row)
        row_w = sum(c.size[0] for c, _ in row) + GAP * (len(row) - 1)
        x = MARGIN + (W - row_w) // 2  # center narrower rows
        for c, letter in row:
            canvas.paste(c, (x, y + (h - c.size[1]) // 2))
            d.rectangle([x, y, x + 105, y + 105], fill="white")
            d.text((x + 12, y + 6), letter, fill="black", font=font)
            x += c.size[0] + GAP
        y += h + GAP
    for ext, kw in [("png", dict(dpi=(300, 300))),
                    ("tiff", dict(compression="tiff_lzw", dpi=(300, 300))),
                    ("pdf", dict(resolution=300))]:
        out = f"{out_stem}.{ext}"
        canvas.save(out, **kw)
        print("saved", out, canvas.size)

FULL = 3050          # 2-cell full width: 2*1500 + 50
CELL2 = (FULL - GAP) // 2            # 1500
CELL4 = (FULL - 3 * GAP) // 4        # 725

F1 = os.path.join(FX, "Figure 1")
assemble([
    [(os.path.join(F1, "parts/A_workflow_schematic.png"), "A", FULL)],
    [([os.path.join(F1, "A_volcano_ARDS.png"), os.path.join(F1, "B_volcano_Sepsis.png")], "B", FULL)],
    [([os.path.join(F1, "C_moduleTrait_ARDS.png"), os.path.join(F1, "D_moduleTrait_Sepsis.png")], "C", FULL)],
    [(os.path.join(F1, "E_venn4sets.png"), "D", 1400)],
], os.path.join(F1, "Figure1_composite"))

F2E = os.path.join(FX, "Figure2_expression_ROC")
F2 = os.path.join(FX, "Figure2")
assemble([
    [(os.path.join(F2E, "Figure2_boxplot_GSE76293.png"), "A", CELL4),
     (os.path.join(F2E, "Figure2_boxplot_GSE243066.png"), "B", CELL4),
     (os.path.join(F2E, "Figure2_boxplot_GSE69528.png"), "C", CELL4),
     (os.path.join(F2E, "Figure2_boxplot_GSE100159.png"), "D", CELL4)],
    [(os.path.join(F2E, "Figure2_ROC_AUC_CI.png"), "E", FULL)],
    [(os.path.join(F2, "parts", "F_ML_SHAP_dual.png"), "F", FULL)],
], os.path.join(F2, "Figure2_composite"))

F3 = os.path.join(FX, "Figure3")
# PPI network: crop the embedded in-image title (redundant with legend; avoids letter-box overlap)
_ppi = Image.open(os.path.join(F3, "A_PPI_network.png")).convert("RGB")
os.makedirs(os.path.join(F3, "parts"), exist_ok=True)
_ppi.crop((0, 260, _ppi.size[0], _ppi.size[1])).save(
    os.path.join(F3, "parts/A_PPI_network_notitle.png"), dpi=(1000, 1000))
F3S = os.path.join(FX, "Figure 2")
assemble([
    [(os.path.join(F3, "parts/A_PPI_network_notitle.png"), "A", CELL2),
     (os.path.join(F3, "B_GO_BP_dotplot.png"), "B", CELL2)],
    [(os.path.join(F3S, "179CGs_GO_sankey_bubble.png"), "C", CELL2),
     (os.path.join(F3S, "179CGs_KEGG_chord.png"), "D", CELL2)],
    [(os.path.join(F3S, "6hub_GO_sankey_bubble.png"), "E", CELL2),
     (os.path.join(F3S, "6hub_KEGG_chord.png"), "F", CELL2)],
], os.path.join(F3, "Figure3_composite"))

F4 = os.path.join(FX, "Figure4")
assemble([
    [(os.path.join(F4, "A_immune_box_ARDS.png"), "A", CELL2),
     (os.path.join(F4, "B_immune_box_Sepsis.png"), "B", CELL2)],
    [(os.path.join(F4, "C_hubHeat_ARDS.png"), "C", CELL2),
     (os.path.join(F4, "D_hubHeat_Sepsis.png"), "D", CELL2)],
], os.path.join(F4, "Figure4_composite"))

F5 = os.path.join(FX, "Figure5")
assemble([
    [(os.path.join(F5, "A_GSVA_top20_ARDS.png"), "A", CELL2),
     (os.path.join(F5, "B_GSVA_top20_Sepsis.png"), "B", CELL2)],
], os.path.join(F5, "Figure5_composite"))

F7 = os.path.join(FX, "Figure7")
assemble([
    [(os.path.join(F7, "A_docking_heatmap.png"), "A", FULL)],
    [([os.path.join(F7, "B_VNN1_luteolin.png"), os.path.join(F7, "D_GRB10_EGCG.png")], "B", FULL)],
    [(os.path.join(F7, "parts/C_cmap_glucocorticoids.png"), "C", 1800)],
], os.path.join(F7, "Figure7_composite"))
