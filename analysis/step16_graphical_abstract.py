# ============================================================
# 重做 Graphical Abstract (按参考图的内容与模块)
# 三个圆角矩形 + 中心枢纽 + 详细子模块
# ============================================================
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Ellipse, Rectangle, FancyArrow
import numpy as np

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT = f"{BASE}/AJRCCM_submission/04_graphical_abstract"

# 配色
COL_GREEN = "#B5D7B4"; COL_ORANGE = "#F0C99A"; COL_BLUE = "#B6C8DD"
COL_CENTER = "#F5E0D3"; COL_TEXT = "#2A2A3A"; COL_EDGE = "#7A7A8A"
COL_ARROW = "#5B608C"; COL_BAR = "#6A5C95"; COL_PINK = "#D99AAE"

fig, ax = plt.subplots(figsize=(20, 14))
ax.set_xlim(0, 20); ax.set_ylim(0, 14); ax.set_aspect("equal"); ax.axis("off")

# 标题
ax.text(10, 13.55, "Shared Hub Genes and Immune Signatures of",
        ha="center", va="center", fontsize=20, fontweight="bold", color=COL_TEXT)
ax.text(10, 13.0, "Sepsis and Acute Respiratory Distress Syndrome",
        ha="center", va="center", fontsize=20, fontweight="bold", color=COL_TEXT)

# ============================================================
# 中心椭圆 "Six shared hub genes" (在绿/橙框之间)
# ============================================================
ax.add_patch(Ellipse((10.0, 8.0), 4.5, 1.4,
    facecolor=COL_CENTER, edgecolor="#6A5C95", linewidth=2))
ax.text(10.0, 8.50, "Six shared", ha="center", fontsize=14, fontweight="bold", color=COL_TEXT)
ax.text(10.0, 8.18, "hub genes", ha="center", fontsize=14, fontweight="bold", color=COL_TEXT)
ax.text(10.0, 7.80, "*CDK5RAP2*   *CD177*   *GRB10*", ha="center", fontsize=10.5, color=COL_TEXT)
ax.text(10.0, 7.50, "*UPP1*   *EXOSC4*   *VNN1*", ha="center", fontsize=10.5, color=COL_TEXT)

# ============================================================
# 左上绿框 "Identification of potential core hub genes" (y=5.8-12.5)
# ============================================================
ax.add_patch(FancyBboxPatch((0.2, 5.8), 9.0, 6.7,
    boxstyle="round,pad=0.08,rounding_size=0.30",
    facecolor=COL_GREEN, edgecolor=COL_EDGE, linewidth=1.5))
ax.text(4.7, 12.20, "Identification of potential core hub genes",
        ha="center", va="center", fontsize=15, fontweight="bold", color=COL_TEXT)

# 数据源 (左上：图标 + 标签)
ax.text(0.8, 11.55, "GSE76293", ha="left", fontsize=11, fontweight="bold", color=COL_TEXT)
ax.text(0.8, 11.25, "ARDS (n=12 vs 12)", ha="left", fontsize=9, color=COL_TEXT)
ax.text(0.8, 10.95, "blood PMNs", ha="left", fontsize=8.5, color=COL_EDGE, style="italic")
# 数据库图标 (肺 + 细胞)
ax.add_patch(plt.Circle((0.7, 11.10), 0.30, facecolor="#FF8C8C", edgecolor="white", linewidth=1.2))
ax.text(0.7, 11.10, "L", ha="center", va="center", fontsize=7, color="white", fontweight="bold")

ax.text(3.3, 11.55, "GSE69528", ha="left", fontsize=11, fontweight="bold", color=COL_TEXT)
ax.text(3.3, 11.25, "Sepsis (n=83 vs 28)", ha="left", fontsize=9, color=COL_TEXT)
ax.text(3.3, 10.95, "whole blood", ha="left", fontsize=8.5, color=COL_EDGE, style="italic")
ax.add_patch(plt.Circle((3.2, 11.10), 0.30, facecolor="#8FA2D6", edgecolor="white", linewidth=1.2))
ax.text(3.2, 11.10, "B", ha="center", va="center", fontsize=7, color="white", fontweight="bold")

# DEG 统计
ax.text(1.0, 10.45, "557 DEGs", ha="left", fontsize=10, color=COL_TEXT)
ax.text(1.0, 10.15, "(291 up / 266 down)", ha="left", fontsize=8, color=COL_EDGE, style="italic")
ax.text(1.0, 9.80, "13 WGCNA modules", ha="left", fontsize=9, color=COL_TEXT)
ax.text(3.5, 10.45, "1,458 DEGs", ha="left", fontsize=10, color=COL_TEXT)
ax.text(3.5, 10.15, "(658 up / 800 down)", ha="left", fontsize=8, color=COL_EDGE, style="italic")
ax.text(3.5, 9.80, "11 WGCNA modules", ha="left", fontsize=9, color=COL_TEXT)

# DEGs & WGCNA → 交集
ax.text(2.4, 9.30, "DEGs & WGCNA", ha="center", fontsize=11, fontweight="bold", color=COL_TEXT,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=COL_EDGE, linewidth=0.8))
ax.annotate("", xy=(2.4, 8.80), xytext=(2.4, 9.10),
            arrowprops=dict(arrowstyle="->", color=COL_ARROW, lw=1.2))

# 交集
ax.text(2.4, 8.55, "179 shared candidate genes", ha="center", fontsize=10.5,
        fontweight="bold", color="#6A5C95",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#6A5C95", linewidth=1.2))
ax.annotate("", xy=(2.4, 8.05), xytext=(2.4, 8.35),
            arrowprops=dict(arrowstyle="->", color=COL_ARROW, lw=1.2))

# ML 柱状图 (左下)
ax.text(2.4, 7.75, "Model selection", ha="center", fontsize=10, fontweight="bold", color=COL_TEXT)
# 简化双柱 (Accuracy / AUC) - 6 个分类器
xpos = np.arange(0.5, 4.4, 0.6)
acc = [0.85, 0.90, 0.92, 0.95, 0.97, 1.00]
auc = [0.80, 0.85, 0.88, 0.92, 0.95, 1.00]
for i, (x, a, u) in enumerate(zip(xpos, acc, auc)):
    ax.add_patch(Rectangle((x, 6.20), 0.20, (a-0.6)*4, facecolor="#8FA2D6", edgecolor="white", linewidth=0.5))
    ax.add_patch(Rectangle((x+0.22, 6.20), 0.20, (u-0.6)*4, facecolor=COL_BAR, edgecolor="white", linewidth=0.5))
# 图例
ax.add_patch(Rectangle((3.7, 6.95), 0.15, 0.15, facecolor="#8FA2D6", edgecolor="white"))
ax.text(3.9, 7.03, "Accuracy", va="center", fontsize=7, color=COL_TEXT)
ax.add_patch(Rectangle((3.7, 7.15), 0.15, 0.15, facecolor=COL_BAR, edgecolor="white"))
ax.text(3.9, 7.23, "AUC", va="center", fontsize=7, color=COL_TEXT)
ax.text(2.4, 6.05, "(11 ML classifiers: LASSO / RF / XGBoost / etc., AUC 0.75-1.00)",
        ha="center", fontsize=7.5, color=COL_EDGE, style="italic")

# GO & KEGG 气泡图 (右下)
ax.text(6.7, 7.75, "GO & KEGG enrichment", ha="center", fontsize=10, fontweight="bold", color=COL_TEXT)
# 气泡 (Hit Ratio vs -log10 Pvalue)
bubbles = [(0.05, 8, 0.30, "#6A5C95"), (0.10, 6, 0.20, "#D99AAE"), (0.08, 5, 0.18, "#E7A889"),
           (0.15, 4, 0.15, "#A9CBB8"), (0.20, 3.5, 0.12, "#CFE3E6"), (0.25, 3, 0.10, "#8FA2D6")]
for hr, lp, sz, c in bubbles:
    ax.add_patch(plt.Circle((5.7 + hr*4, 6.20 + lp*0.10), sz, facecolor=c, edgecolor="white", linewidth=0.8, alpha=0.85))
# 坐标轴
ax.plot([5.5, 8.9], [6.20, 6.20], "k-", linewidth=0.5)
ax.plot([5.5, 5.5], [6.20, 7.20], "k-", linewidth=0.5)
ax.text(5.5, 6.05, "Hit Ratio", fontsize=7, color=COL_EDGE)
ax.text(8.7, 6.05, "1.0", fontsize=7, color=COL_EDGE, ha="right")
ax.text(5.35, 7.30, "−log10 P", fontsize=7, color=COL_EDGE, rotation=90, va="center")
ax.text(5.5, 7.18, "Count", fontsize=7, color=COL_EDGE, style="italic")
# 中心椭圆箭头
ax.annotate("", xy=(9.0, 7.80), xytext=(5.0, 7.80),
            arrowprops=dict(arrowstyle="->", color=COL_ARROW, lw=1.5,
                            connectionstyle="arc3,rad=-0.2"))

# ============================================================
# 右上橙框 "Exploration of molecular mechanisms" (y=5.8-12.5)
# ============================================================
ax.add_patch(FancyBboxPatch((11.0, 5.8), 8.8, 6.7,
    boxstyle="round,pad=0.08,rounding_size=0.30",
    facecolor=COL_ORANGE, edgecolor=COL_EDGE, linewidth=1.5))
ax.text(15.4, 12.20, "Exploration of molecular mechanisms",
        ha="center", va="center", fontsize=15, fontweight="bold", color=COL_TEXT)

# 左列: Immune Infiltration & Single-cell
ax.text(13.0, 11.80, "Immune infiltration & Single-cell",
        ha="center", fontsize=11, fontweight="bold", color=COL_TEXT)
# ssGSEA 折线 (Immunity)
ax.text(11.6, 11.40, "ssGSEA", fontsize=9, fontweight="bold", color=COL_TEXT)
xs = np.linspace(11.4, 12.8, 30)
ys = 11.20 + 0.15*np.sin(xs*5) + 0.05*np.random.RandomState(42).randn(30)
ax.plot(xs, ys, color="#6A5C95", linewidth=1.2)
ax.plot(xs, 11.20+0.05*np.cos(xs*4), color="#D99AAE", linewidth=1.2)
ax.text(12.0, 10.90, "M2 macrophages", fontsize=7, color="#D99AAE")
ax.text(11.5, 10.75, "(sepsis whole blood)", fontsize=7, color=COL_EDGE, style="italic")
ax.text(11.5, 10.60, "P=5.8e-14 (Sepsis)", fontsize=7, color=COL_EDGE, style="italic")

# 单细胞方法
ax.text(11.6, 10.30, "Single-cell:", fontsize=9, fontweight="bold", color=COL_TEXT)
methods = [("CellChat (MIF)", 10.00), ("Slingshot (mono→M2)", 9.80),
           ("scTenifoldKnk", 9.60), ("CellOracle (CEBPB KO)", 9.40)]
for txt, y in methods:
    ax.text(11.6, y, txt, fontsize=8, color=COL_TEXT)

# 中列: Pathway enrichment (GSEA-style 曲线)
ax.text(15.4, 11.80, "Pathway enrichment",
        ha="center", fontsize=11, fontweight="bold", color=COL_TEXT)
# GSEA 曲线
xs2 = np.linspace(14.0, 16.8, 100)
ax.plot(xs2, 10.80 + 0.30*np.exp(-(xs2-15.0)**2/0.3), color="#D99AAE", linewidth=2)
ax.text(15.0, 11.10, "Single-gene GSEA", fontsize=8.5, color=COL_TEXT, style="italic")
ax.text(15.4, 10.45, "(complement activation)", fontsize=8.5, color=COL_EDGE, style="italic")
# TF-mRNA 网络符号
ax.add_patch(plt.Circle((14.5, 10.00), 0.18, facecolor="#6A5C95", edgecolor="white"))
ax.text(14.5, 10.00, "TF", ha="center", va="center", fontsize=7, color="white", fontweight="bold")
ax.text(15.0, 10.00, "→", ha="center", va="center", fontsize=14, color=COL_ARROW)
ax.add_patch(plt.Circle((15.6, 10.00), 0.18, facecolor="#D99AAE", edgecolor="white"))
ax.text(15.6, 10.00, "mR", ha="center", va="center", fontsize=6, color="white", fontweight="bold")
ax.add_patch(plt.Circle((16.7, 10.00), 0.18, facecolor="#E7A889", edgecolor="white"))
ax.text(16.7, 10.00, "mr", ha="center", va="center", fontsize=6, color="white", fontweight="bold")
ax.plot([15.78, 16.52], [10.00, 10.00], color=COL_ARROW, linewidth=1.5)
ax.text(15.6, 10.30, "TF-mRNA network", fontsize=8, color=COL_TEXT, style="italic")

# 信号通路卡通
ax.text(15.4, 9.65, "Signaling pathways:", fontsize=8.5, fontweight="bold", color=COL_TEXT)
ax.text(14.2, 9.40, "p53 axis", fontsize=8, color=COL_TEXT,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="#CFE3E6", edgecolor="#8FA2D6"))
ax.text(15.4, 9.40, "PI3K/Akt", fontsize=8, color=COL_TEXT,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="#CFE3E6", edgecolor="#8FA2D6"))
ax.text(16.6, 9.40, "MAPK", fontsize=8, color=COL_TEXT,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="#CFE3E6", edgecolor="#8FA2D6"))

# 转录翻译 + 免疫细胞
ax.text(15.4, 9.10, "Transcription & Translation", fontsize=8.5, fontweight="bold", color=COL_TEXT, ha="center")
ax.text(14.0, 8.80, "Tryptophan", fontsize=8, color=COL_TEXT)
ax.text(14.0, 8.60, "metabolism", fontsize=8, color=COL_TEXT)
ax.text(15.0, 8.70, "→", fontsize=14, color=COL_ARROW, ha="center", va="center")
ax.text(15.5, 8.70, "M2 polarization", fontsize=8, color=COL_PINK, fontweight="bold", ha="left", va="center")
# 免疫细胞卡通
ax.text(17.0, 8.90, "T cell", fontsize=8, color=COL_TEXT, fontweight="bold")
ax.text(17.0, 8.70, "(exhausted)", fontsize=7.5, color=COL_EDGE, style="italic")
ax.text(17.0, 8.45, "Macrophage", fontsize=8, color=COL_TEXT, fontweight="bold")
ax.text(17.0, 8.25, "(M2 polarized)", fontsize=7.5, color=COL_EDGE, style="italic")

# 共享免疫微环境
ax.text(15.4, 7.80, "shared immune microenvironment",
        ha="center", fontsize=12, fontweight="bold", color="#6A5C95",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#6A5C95", linewidth=1.2))
ax.annotate("", xy=(11.0, 7.80), xytext=(15.0, 8.00),
            arrowprops=dict(arrowstyle="->", color=COL_ARROW, lw=1.5,
                            connectionstyle="arc3,rad=0.15"))

# ============================================================
# 下蓝框 "Drug prediction" (y=0.3-5.0)
# ============================================================
ax.add_patch(FancyBboxPatch((3.0, 0.3), 14.0, 4.7,
    boxstyle="round,pad=0.08,rounding_size=0.30",
    facecolor=COL_BLUE, edgecolor=COL_EDGE, linewidth=1.5))
ax.text(10.0, 4.70, "Drug prediction",
        ha="center", va="center", fontsize=15, fontweight="bold", color=COL_TEXT)

# 数据源: CMap 立方体 + AutoDock Vina
ax.add_patch(FancyBboxPatch((3.8, 3.50), 1.3, 0.9,
    boxstyle="round,pad=0.05", facecolor="#CFE3E6", edgecolor="#8FA2D6", linewidth=0.8))
ax.text(4.45, 3.95, "CMap", ha="center", va="center", fontsize=11, fontweight="bold", color="#2A2A3A")
ax.text(4.45, 3.65, "L1000", ha="center", va="center", fontsize=7, color=COL_EDGE, style="italic")
ax.add_patch(FancyBboxPatch((6.5, 3.50), 1.3, 0.9,
    boxstyle="round,pad=0.05", facecolor="#CFE3E6", edgecolor="#8FA2D6", linewidth=0.8))
ax.text(7.15, 3.95, "AutoDock", ha="center", va="center", fontsize=11, fontweight="bold", color="#2A2A3A")
ax.text(7.15, 3.65, "Vina", ha="center", va="center", fontsize=7, color=COL_EDGE, style="italic")
ax.text(3.5, 4.55, "databases", ha="left", fontsize=9, color=COL_EDGE, style="italic")
ax.annotate("", xy=(4.45, 3.50), xytext=(4.45, 3.20), arrowprops=dict(arrowstyle="->", color=COL_ARROW, lw=1.0))
ax.annotate("", xy=(7.15, 3.50), xytext=(7.15, 3.20), arrowprops=dict(arrowstyle="->", color=COL_ARROW, lw=1.0))

# CMap → Glucocorticoids
ax.text(4.45, 2.95, "Glucocorticoids", ha="center", fontsize=11, fontweight="bold", color="#6A5C95")
ax.text(4.45, 2.55, "(prednisolone, dexamethasone,\nhydrocortisone)",
        ha="center", fontsize=8, color=COL_TEXT)
ax.text(4.45, 1.95, "standard-of-care therapy", ha="center", fontsize=8.5, color=COL_EDGE, style="italic")

# Docking: 卡通蛋白 + 分子
ax.add_patch(plt.Circle((7.15, 2.50), 0.55, facecolor="#A9CBB8", edgecolor="white", linewidth=1.5))
ax.text(7.15, 2.50, "VNN1", ha="center", va="center", fontsize=9, fontweight="bold", color="#2A2A3A")
# 分子结构 (简化)
ax.text(8.20, 2.50, "+", ha="center", va="center", fontsize=18, color=COL_ARROW)
ax.add_patch(plt.Circle((9.10, 2.50), 0.55, facecolor="#D99AAE", edgecolor="white", linewidth=1.5))
ax.text(9.10, 2.50, "lig", ha="center", va="center", fontsize=9, fontweight="bold", color="#2A2A3A")
ax.annotate("", xy=(9.65, 2.50), xytext=(9.50, 2.50),
            arrowprops=dict(arrowstyle="->", color=COL_ARROW, lw=1.5))
ax.text(10.7, 2.50, "−10.01 kcal/mol", ha="left", va="center", fontsize=10, fontweight="bold", color="#6A5C95")
ax.text(7.15, 1.95, "VNN1 druggable target", ha="center", fontsize=9.5, color=COL_TEXT)
ax.text(7.15, 1.70, "flavonoid natural products", ha="center", fontsize=8.5, color=COL_EDGE, style="italic")
ax.text(7.15, 1.45, "(luteolin / quercetin / EGCG / etc.)", ha="center", fontsize=7.5, color=COL_EDGE, style="italic")

# 关键结果右侧
ax.text(12.0, 4.30, "Key results:", ha="left", fontsize=10, fontweight="bold", color=COL_TEXT)
ax.text(12.0, 4.00, "• VNN1–luteolin: −10.01 kcal/mol", ha="left", fontsize=9, color=COL_TEXT)
ax.text(12.0, 3.75, "• GRB10–EGCG: −9.61 kcal/mol", ha="left", fontsize=9, color=COL_TEXT)
ax.text(12.0, 3.50, "• 6-gene signature: AUC 0.984 (ARDS)", ha="left", fontsize=9, color=COL_TEXT)
ax.text(12.0, 3.25, "• M2 (sepsis) + complement axis", ha="left", fontsize=9, color=COL_TEXT)
ax.text(12.0, 3.00, "• Glucocorticoid repositioning", ha="left", fontsize=9, color=COL_TEXT)
ax.text(12.0, 2.75, "• VNN1 / 6KOO as drug nodes", ha="left", fontsize=9, color=COL_TEXT)
# 卡通蛋白
ax.add_patch(plt.Circle((14.5, 1.80), 0.30, facecolor="#A9CBB8", edgecolor="white", linewidth=1.2))
ax.add_patch(plt.Circle((15.3, 1.80), 0.30, facecolor="#A9CBB8", edgecolor="white", linewidth=1.2))
ax.add_patch(plt.Circle((16.1, 1.80), 0.30, facecolor="#A9CBB8", edgecolor="white", linewidth=1.2))
ax.add_patch(plt.Circle((15.0, 1.20), 0.30, facecolor="#D99AAE", edgecolor="white", linewidth=1.2))
ax.add_patch(plt.Circle((15.8, 1.20), 0.30, facecolor="#D99AAE", edgecolor="white", linewidth=1.2))
ax.text(15.8, 0.75, "docking poses (PDB / AlphaFold)", ha="center", fontsize=7, color=COL_EDGE, style="italic")

# Therapeutic relevance
ax.text(10.0, 0.65, "↑ Therapeutic relevance", ha="center", fontsize=11, fontweight="bold", color="#5B608C")

# ============================================================
# 箭头连接 (弯曲/直线)
# ============================================================
# 绿框 → 中心椭圆 (弧向左)
ax.annotate("", xy=(8.0, 7.80), xytext=(9.2, 7.80),
            arrowprops=dict(arrowstyle="->", color=COL_ARROW, lw=1.5,
                            connectionstyle="arc3,rad=-0.2"))
# 橙框 → 中心椭圆 (弧向右)
ax.annotate("", xy=(12.0, 7.80), xytext=(10.8, 7.80),
            arrowprops=dict(arrowstyle="->", color=COL_ARROW, lw=1.5,
                            connectionstyle="arc3,rad=0.2"))
# 中心 → 蓝框 (直线)
ax.annotate("", xy=(10.0, 5.00), xytext=(10.0, 7.00),
            arrowprops=dict(arrowstyle="->", color=COL_ARROW, lw=1.8))

# 底部署名
ax.text(10.0, 0.10, "Transcriptomic · Machine learning · Single-cell · Drug repositioning & docking",
        ha="center", fontsize=10, style="italic", color=COL_EDGE)

plt.tight_layout()
plt.savefig(f"{OUT}/graphical_abstract.pdf", bbox_inches="tight", dpi=300, facecolor="white")
plt.savefig(f"{OUT}/graphical_abstract.png", bbox_inches="tight", dpi=300, facecolor="white")
plt.savefig(f"{OUT}/graphical_abstract.tiff", bbox_inches="tight", dpi=500, facecolor="white",
            pil_kwargs={"compression": "tiff_lzw"})
plt.close()
print("Graphical abstract 已保存（PDF/PNG 300dpi + TIFF 500dpi）")