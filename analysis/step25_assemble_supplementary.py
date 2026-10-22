#!/usr/bin/env python
"""Assemble supplementary composite figures S4, S5, S9, S10 (and S1 from parts).
Style follows scripts/assemble_figures.py: PIL, white bg, letter labels, 300dpi triple export.
PDF panels are rendered to 300dpi PNG via Ghostscript before compositing.
"""
from PIL import Image, ImageDraw, ImageFont
import os, subprocess, tempfile

ROOT = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
SC = os.path.join(ROOT, "08_singlecell")
G167 = os.path.join(SC, "GSE167363")
SUB = os.path.join(ROOT, "AJRCCM_submission/05_supplementary/Supplementary_Figures_Tables/Figures")

GAP = 50
MARGIN = 60
LETTER_SZ = 80
font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", LETTER_SZ)

def render_pdf(pdf_path):
    """Render first page of a PDF to a 300dpi RGB PIL image via gs."""
    fd, tmp = tempfile.mkstemp(suffix=".png")
    os.close(fd)
    subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=png16m", "-r300",
                    "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
                    "-sOutputFile=" + tmp, pdf_path], check=True)
    im = Image.open(tmp).convert("RGB")
    os.remove(tmp)
    return im

def load(src, target_w):
    if src.lower().endswith(".pdf"):
        im = render_pdf(src)
    else:
        im = Image.open(src)
        if im.mode in ("RGBA", "LA", "P"):
            im = im.convert("RGBA")
            bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
            im = Image.alpha_composite(bg, im).convert("RGB")
        else:
            im = im.convert("RGB")
    w, h = im.size
    return im.resize((target_w, int(round(h * target_w / w))), Image.LANCZOS)

def assemble(rows_spec, out_stem):
    composed = []
    for row in rows_spec:
        cells = []
        for src, letter, cw in row:
            cells.append((load(src, cw), letter))
        composed.append(cells)
    W = max(sum(c.size[0] for c, _ in row) + GAP * (len(row) - 1) for row in composed)
    Hh = MARGIN * 2 + sum(max(c.size[1] for c, _ in row) for row in composed) + GAP * (len(composed) - 1)
    canvas = Image.new("RGB", (W + MARGIN * 2, Hh), "white")
    d = ImageDraw.Draw(canvas)
    y = MARGIN
    for row in composed:
        h = max(c.size[1] for c, _ in row)
        row_w = sum(c.size[0] for c, _ in row) + GAP * (len(row) - 1)
        x = MARGIN + (W - row_w) // 2
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

FULL = 3050
CELL2 = (FULL - GAP) // 2            # 1500
CELL3 = (FULL - 2 * GAP) // 3        # ~983
CELL4 = (FULL - 3 * GAP) // 4        # 725

# ---------- Figure S4: GSE69528 immune landscape ----------
IM = os.path.join(ROOT, "07_hub_immune")
assemble([
    [(os.path.join(IM, "GSE69528_immune_boxplot.pdf"), "A", FULL)],
    [(os.path.join(IM, "GSE69528_hub_immune_heatmap.pdf"), "B", FULL)],
], os.path.join(SUB, "FigureS4_immune_infiltration_GSE69528"))

# ---------- Figure S5: CellChat pathway comparison ----------
CC263 = os.path.join(SC, "cellchat")            # GSE151263 (PBMC: ARDS vs Sepsis)
CC363 = os.path.join(G167, "cellchat")          # GSE167363 (CD14+: Sepsis vs HC)
assemble([
    [(os.path.join(CC263, "rankNet_pathways.png"), "A", CELL2),
     (os.path.join(CC363, "rankNet_pathways.pdf"), "B", CELL2)],
    [(os.path.join(CC263, "diff_LR_bubble.png"), "C", CELL2),
     (os.path.join(CC363, "diff_LR_bubble.pdf"), "D", CELL2)],
], os.path.join(SUB, "FigureS5_CellChat_pathway_comparison"))

# ---------- Figure S9: GSE167363 single-cell atlas ----------
assemble([
    [(os.path.join(G167, "GSE167363_UMAP.pdf"), "A", CELL3),
     (os.path.join(G167, "GSE167363_celltype_prop.pdf"), "B", CELL3),
     (os.path.join(G167, "GSE167363_marker_dotplot.pdf"), "C", CELL3)],
    [(os.path.join(G167, "GSE167363_hub_dotplot.pdf"), "D", CELL3),
     (os.path.join(G167, "GSE167363_hub_featureplot.pdf"), "E", CELL3),
     (os.path.join(G167, "GSE167363_hub_violin_group.pdf"), "F", CELL3)],
], os.path.join(SUB, "FigureS9_single_cell_atlas_GSE167363"))

# ---------- Figure S10: pseudotime + AUCell + CellOracle ----------
SL = os.path.join(G167, "slingshot")
AU = os.path.join(G167, "aucell")
CO = os.path.join(G167, "celloracle")
assemble([
    [(os.path.join(SL, "GSE167363_slingshot_overview.pdf"), "A", CELL3),
     (os.path.join(SL, "GSE167363_slingshot_trajectory_PCA.pdf"), "B", CELL3),
     (os.path.join(SL, "GSE167363_slingshot_gene_trend.pdf"), "C", CELL3)],
    [(os.path.join(AU, "GSE167363_AUCell_violin.pdf"), "D", CELL3),
     (os.path.join(AU, "GSE167363_AUCell_heatmap.pdf"), "E", CELL3),
     (os.path.join(AU, "GSE167363_AUCell_umap.pdf"), "F", CELL3)],
    [(os.path.join(CO, "CellOracle_KO_effect_UPP1.pdf"), "G", CELL4),
     (os.path.join(CO, "CellOracle_KO_effect_CDK5RAP2.pdf"), "H", CELL4),
     (os.path.join(CO, "CellOracle_KO_effect_EXOSC4.pdf"), "I", CELL4),
     (os.path.join(CO, "CellOracle_KO_effect_VNN1.pdf"), "J", CELL4)],
], os.path.join(SUB, "FigureS10_pseudotime_AUCell_CellOracle"))

# ---------- Figure S1: WGCNA soft-threshold + dendrograms ----------
WG = os.path.join(ROOT, "02_WGCNA/figS1_parts")
if os.path.exists(os.path.join(WG, "a_sft_GSE76293.png")):
    assemble([
        [(os.path.join(WG, "a_sft_GSE76293.png"), "A", CELL3),
         (os.path.join(WG, "b_conn_GSE76293.png"), "B", CELL3),
         (os.path.join(WG, "c_dendro_GSE76293.pdf"), "C", CELL3)],
        [(os.path.join(WG, "d_sft_GSE69528.png"), "D", CELL3),
         (os.path.join(WG, "e_conn_GSE69528.png"), "E", CELL3),
         (os.path.join(WG, "f_dendro_GSE69528.pdf"), "F", CELL3)],
    ], os.path.join(SUB, "FigureS1_WGCNA_soft_threshold_dendrograms"))

print("DONE")
