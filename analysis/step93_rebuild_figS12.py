# Step 93 —— 重建 FigureS12 (CellChat) 为合规双栏矢量 composite
# 旧版整体拷贝自光栅化的 FigureS10_CellChat_pathway_comparison (10.57in, 非矢量)
# 新布局: A/B rankNet 信息流 (half+half) → C/D diff_LR_bubble (full width)
import pymupdf
from PIL import Image

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
CC1  = f"{BASE}/08_singlecell/cellchat"               # GSE151263 ARDS vs Sepsis
CC2  = f"{BASE}/08_singlecell/GSE167363/cellchat"     # GSE167363 HC vs Sepsis
OUT  = f"{BASE}/AJRCCM_submission/05_supplementary_NPG/composites"

PAGE_W = 7.01 * 72          # 504.7 pt
MARGIN = 14
GAP_X  = 16
GAP_Y  = 16
LETTER = 14

doc = pymupdf.open()

def pdf_size(path):
    d = pymupdf.open(path); r = d[0].rect; d.close()
    return r.width, r.height

def aspect_h(path, w):
    pw, ph = pdf_size(path)
    return w * ph / pw

panels = [
    # (path, letter, width_mode)  width_mode: 'half' | number(pt)
    (f"{CC1}/rankNet_pathways.pdf",  "A", "half"),
    (f"{CC2}/rankNet_pathways.pdf",  "B", "half"),
    (f"{CC1}/diff_LR_bubble.pdf",    "C", 440),
    (f"{CC2}/diff_LR_bubble.pdf",    "D", 440),
]

half_w = (PAGE_W - 2 * MARGIN - GAP_X) / 2
# row heights
hA = aspect_h(panels[0][0], half_w)
hB = aspect_h(panels[1][0], half_w)
row1 = max(hA, hB)
hC = aspect_h(panels[2][0], panels[2][2])
hD = aspect_h(panels[3][0], panels[3][2])
PAGE_H = MARGIN + row1 + GAP_Y + hC + GAP_Y + hD + 6

page = doc.new_page(width=PAGE_W, height=PAGE_H)

def place(path, x, y, w, letter):
    r = pymupdf.Rect(x, y, x + w, y + aspect_h(path, w))
    page.show_pdf_page(r, pymupdf.open(path), 0)
    page.insert_text(pymupdf.Point(x, y - 3), letter,
                     fontname="hebo", fontsize=9, color=(0, 0, 0))

y = MARGIN + LETTER
place(panels[0][0], MARGIN, y, half_w, "A")
place(panels[1][0], MARGIN + half_w + GAP_X, y, half_w, "B")
y += row1 + GAP_Y
place(panels[2][0], (PAGE_W - panels[2][2]) / 2, y, panels[2][2], "C")
y += hC + GAP_Y
place(panels[3][0], (PAGE_W - panels[3][2]) / 2, y, panels[3][2], "D")

pdf_path = f"{OUT}/FigureS12_composite.pdf"
doc.save(pdf_path, deflate=True)
doc.close()
print("saved", pdf_path, f"{PAGE_W/72:.2f} x {PAGE_H/72:.2f} in")

# rasters
d = pymupdf.open(pdf_path)
for dpi, ext in [(300, "png"), (1000, "tiff")]:
    pix = d[0].get_pixmap(matrix=pymupdf.Matrix(dpi / 72, dpi / 72), alpha=False)
    im = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    if ext == "png":
        im.save(f"{OUT}/FigureS12_composite.png", dpi=(dpi, dpi))
    else:
        im.save(f"{OUT}/FigureS12_composite.tiff", format="TIFF",
                compression="tiff_lzw", dpi=(dpi, dpi))
d.close()
print("rasters exported")
