#!/usr/bin/env python
"""Step3n-1: TCR/BCR 克隆性分析（转录本水平）— GSE151263 (PBMC) & GSE167363 (CD14+ Mono)。

两个数据集均为 5' GEX 文库、无 VDJ 富集文库（GEO supplementary 仅有计数矩阵），
因此无法做 CDR3 水平真正克隆型。本脚本采用既定的替代口径：
  受体链转录本组合（TRAC/TRBC1/TRBC2/IGKC/IGLC2/... 的表达组合）作为
  "伪克隆型"（pseudo-clonotype），量化克隆性扩增代理指标：
  1. 受体转录本检出率（细胞类型 × 链）
  2. 伪克隆型丰富度 / Shannon 熵 / Simpson / Gini / top-10 克隆占比（样本×分组）
  3. 分组间克隆性扩增比较（ARDS vs Sepsis；Sepsis vs HC）
  4. top 伪克隆型构成与扩张

输出（08_singlecell/clonality/）:
  GSE151263_receptor_dotplot.{pdf,png,tiff}      受体链转录本 × 细胞类型
  GSE151263_clonality_diversity.{pdf,png,tiff}    多样性指数组间比较（4 指数面板）
  GSE151263_clonality_rank_abundance.{pdf,png,tiff}  rank-abundance 扩张曲线
  GSE151263_top_clonotypes.{pdf,png,tiff}        top 伪克隆型分组堆积
  GSE167363 同名文件
  clonality_indices_{dataset}.csv / top_clonotypes_{dataset}.csv / receptor_detection_{dataset}.csv
  TCR_BCR_clonality_methods.md                   方法学与局限性声明
"""
import os
import csv
from collections import Counter, OrderedDict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import scipy.sparse as sp
from matplotlib import font_manager
from scipy.stats import mannwhitneyu

ROOT = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
SC = os.path.join(ROOT, "08_singlecell")
OUT = os.path.join(SC, "clonality")
os.makedirs(OUT, exist_ok=True)

PURPLE = "#6A5C95"; PINK = "#D99AAE"; ORANGE = "#E7A889"
GREEN = "#A9CBB8"; MINT = "#CFE3E6"; BLUE = "#8FA2D6"; NAVY = "#5B608C"
PAL7 = [PURPLE, PINK, ORANGE, GREEN, MINT, BLUE, NAVY]

for fp in ("/System/Library/Fonts/Helvetica.ttc",):
    if os.path.exists(fp):
        font_manager.fontManager.addfont(fp)
plt.rcParams["font.family"] = ["Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42

TCR_GENES = ["TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGV9"]
BCR_GENES = ["IGHM", "IGHG1", "IGHG2", "IGHG3", "IGKC", "IGLC2", "IGLC3", "IGLL1", "IGLL5"]


def save3(fig, name, line_art=False):
    out = os.path.join(OUT, name)
    fig.savefig(out + ".pdf", facecolor="white")
    fig.savefig(out + ".png", dpi=300, facecolor="white")
    fig.savefig(out + ".tiff", dpi=500 if not line_art else 1000, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("[saved]", name)


def gini(x):
    x = np.sort(np.asarray(x, dtype=float))
    if x.sum() <= 0:
        return np.nan
    n = len(x)
    return (2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum())


def diversity_stats(counts_by_combo):
    """counts_by_combo: Counter of pseudo-clonotype -> cell number"""
    v = np.array(list(counts_by_combo.values()), dtype=float)
    p = v / v.sum()
    shannon = -(p[p > 0] * np.log(p[p > 0])).sum()
    simpson = 1 - (p ** 2).sum()
    return dict(n_clonotypes=len(v), shannon=shannon, simpson=simpson,
                gini=gini(v), top10_frac=(np.sort(v)[::-1][:10].sum() / v.sum()))


def load_dataset(which):
    if which == "GSE151263":
        cd = os.path.join(SC, "cellchat_input")
    else:
        cd = os.path.join(SC, "GSE167363/input")
    cnt = sio.mmread(os.path.join(cd, "counts.mtx")).tocsr()
    genes = [l.strip().split("\t")[0] for l in open(os.path.join(cd, "genes.tsv"))]
    rows = list(csv.DictReader(open(os.path.join(cd, "meta.csv"))))
    cell_key = "cell" if rows and "cell" in rows[0] else ("" if rows else "cell")
    meta = {r[cell_key]: r for r in rows}
    cells = list(meta.keys())
    assert cnt.shape[0] == len(genes) and cnt.shape[1] == len(cells), (cnt.shape, len(genes), len(cells))
    return cnt, genes, cells, meta


def analyze(which, group_case, group_ctrl):
    print("\n=====", which, "=====")
    cnt, genes, cells, meta = load_dataset(which)
    gidx = {g: i for i, g in enumerate(genes)}
    tcr_idx = [gidx[g] for g in TCR_GENES if g in gidx]
    bcr_idx = [gidx[g] for g in BCR_GENES if g in gidx]
    tcr_names = [g for g in TCR_GENES if g in gidx]
    bcr_names = [g for g in BCR_GENES if g in gidx]
    ridx = tcr_idx + bcr_idx
    rnames = tcr_names + bcr_names

    dens = cnt[ridx, :].toarray().astype(float)  # chains x cells
    detected = dens > 0

    groups = np.array([meta[c]["group"] for c in cells])
    ctypes = np.array([meta[c]["celltype"] for c in cells])
    samples = np.array([meta[c]["sample"] for c in cells])
    uniq_ct = sorted(set(ctypes))

    # ---------- 1. 受体链检出率（细胞类型 × 链） ----------
    det_pct = np.zeros((len(uniq_ct), len(rnames)))
    det_mean = np.zeros((len(uniq_ct), len(rnames)))
    for i, ct in enumerate(uniq_ct):
        m = ctypes == ct
        det_pct[i] = detected[:, m].mean(axis=1) * 100
        det_mean[i] = dens[:, m].mean(axis=1)
    with open(os.path.join(OUT, f"receptor_detection_{which}.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["celltype"] + rnames)
        for i, ct in enumerate(uniq_ct):
            w.writerow([ct] + [f"{det_pct[i, j]:.2f}|{det_mean[i, j]:.3f}" for j in range(len(rnames))])

    fig, ax = plt.subplots(figsize=(9.2, 0.55 * len(uniq_ct) + 2.2))
    X, Y = np.meshgrid(np.arange(len(rnames)), np.arange(len(uniq_ct)))
    sc_ = ax.scatter(X.ravel(), Y.ravel(), s=det_pct.ravel() * 5.5 + (det_pct.ravel() > 0) * 8,
                     c=det_mean.ravel(), cmap="BuPu", vmin=0,
                     vmax=np.percentile(det_mean, 99) or 1, edgecolor="#4D4D4D", linewidth=0.3)
    ax.set_xticks(range(len(rnames)), rnames, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(uniq_ct)), uniq_ct, fontsize=8)
    cb = fig.colorbar(sc_, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label("Mean transcript count", fontsize=8)
    ax.set_title(f"{which}: TCR/BCR receptor-chain transcript detection\n(size = % cells detected; VDJ library not available)",
                 fontsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xlim(-0.7, len(rnames) - 0.3); ax.set_ylim(len(uniq_ct) - 0.3, -0.7)
    fig.tight_layout()
    save3(fig, f"{which}_receptor_dotplot")

    # ---------- 2. 伪克隆型（受体链组合） ----------
    tcr_only = detected[:len(tcr_idx)]
    bcr_only = detected[len(tcr_idx):]
    combos = []
    for j in range(len(cells)):
        # TCR 组合（T/B/NK 细胞用）；BCR 组合（B 细胞用）
        if ctypes[j] in ("B", "Plasma") or (bcr_only[:, j].any() and not tcr_only[:, j].any()):
            chains = tuple(bcr_names[k] for k in np.where(bcr_only[:, j])[0])
        else:
            chains = tuple(tcr_names[k] for k in np.where(tcr_only[:, j])[0])
        combos.append("+".join(chains) if chains else "none")

    combos = np.array(combos)
    per_sample = OrderedDict()
    for smp in sorted(set(samples)):
        m = samples == smp
        c = Counter(combos[m])
        c.pop("none", None)
        if not c:
            c = {"none": int(m.sum())}
        per_sample[smp] = c

    sample2group = {}
    for c in cells:
        sample2group[meta[c]["sample"]] = meta[c]["group"]

    rows_idx = []
    for smp, c in per_sample.items():
        st = diversity_stats(c) if "none" not in c or len(c) > 1 else dict(
            n_clonotypes=0, shannon=0, simpson=0, gini=0, top10_frac=0)
        grp = sample2group[smp]
        n_recep = int(sum(v for k, v in c.items() if k != "none"))
        rows_idx.append(dict(sample=smp, group=grp, n_receptor_cells=n_recep, **st))
    with open(os.path.join(OUT, f"clonality_indices_{which}.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_idx[0].keys())); w.writeheader()
        for r in rows_idx: w.writerow(r)

    # ---------- 3. 多样性指数组间比较（4 面板） ----------
    metrics = [("n_clonotypes", "Richness (# pseudo-clonotypes)"),
               ("shannon", "Shannon entropy"), ("gini", "Gini coefficient"),
               ("top10_frac", "Top-10 clonotype fraction")]
    fig, axes = plt.subplots(1, 4, figsize=(12.5, 3.6))
    case_col, ctrl_col = PINK, BLUE
    for ax, (key, lab) in zip(axes, metrics):
        vals = {g: [r[key] for r in rows_idx if r["group"] == g] for g in (group_case, group_ctrl)}
        for k, (grp, col) in enumerate(((group_case, case_col), (group_ctrl, ctrl_col))):
            x = np.random.normal(k, 0.06, len(vals[grp]))
            ax.scatter(x, vals[grp], s=42, color=col, edgecolor="k", linewidth=0.4, zorder=3,
                       label=f"{grp} (n={len(vals[grp])})")
        bp = ax.boxplot([vals[group_case], vals[group_ctrl]], positions=[0, 1], widths=0.5,
                        patch_artist=True, showfliers=False)
        for patch, col in zip(bp["boxes"], (case_col, ctrl_col)):
            patch.set_facecolor(col); patch.set_alpha(0.45); patch.set_edgecolor("#4D4D4D")
        try:
            u, pv = mannwhitneyu(vals[group_case], vals[group_ctrl])
            ax.set_title(f"{lab}\nMW p = {pv:.2g}", fontsize=9)
        except ValueError:
            ax.set_title(lab, fontsize=9)
        ax.set_xticks([0, 1], [group_case, group_ctrl], fontsize=9)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].legend(fontsize=7, loc="upper right", frameon=False)
    fig.suptitle(f"{which}: receptor-transcript clonality, {group_case} vs {group_ctrl}\n(pseudo-clonotypes from chain-combination; per-sample)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    save3(fig, f"{which}_clonality_diversity")

    # ---------- 4. rank-abundance ----------
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    for grp, col in ((group_case, case_col), (group_ctrl, ctrl_col)):
        agg = Counter()
        for smp in per_sample:
            if sample2group[smp] == grp:
                for k, v in per_sample[smp].items():
                    if k != "none":
                        agg[k] += v
        if not agg:
            continue
        v = np.array(sorted(agg.values())[::-1], dtype=float)
        p = v / v.sum()
        n = min(50, len(p))
        ax.plot(np.arange(1, n + 1), p[:n], "-o", ms=3.5, lw=1.4, color=col,
                label=f"{grp} ({len(v)} pseudo-clonotypes)")
    ax.set_yscale("log")
    ax.set_xlabel("Pseudo-clonotype rank"); ax.set_ylabel("Frequency (log)")
    ax.set_title(f"{which}: clonotype rank-abundance (top 50)", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    save3(fig, f"{which}_clonality_rank_abundance")

    # ---------- 5. top 伪克隆型构成 ----------
    agg_all = Counter()
    for smp in per_sample:
        for k, v in per_sample[smp].items():
            if k != "none":
                agg_all[k] += v
    top_combos = [k for k, _ in agg_all.most_common(8)]
    grp_mat = np.zeros((len(top_combos), 2))
    for j, grp in enumerate((group_case, group_ctrl)):
        m = groups == grp
        for i, combo in enumerate(top_combos):
            grp_mat[i, j] = (combos[m] == combo).mean() * 100
    with open(os.path.join(OUT, f"top_clonotypes_{which}.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["pseudo_clonotype"] + [group_case, group_ctrl, "n_cells_total"])
        for i, c_ in enumerate(top_combos):
            w.writerow([c_, f"{grp_mat[i,0]:.2f}", f"{grp_mat[i,1]:.2f}", agg_all[c_]])

    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    x = np.arange(len(top_combos)); wdt = 0.38
    ax.bar(x - wdt / 2, grp_mat[:, 0], wdt, color=case_col, edgecolor="k", linewidth=0.4, label=group_case)
    ax.bar(x + wdt / 2, grp_mat[:, 1], wdt, color=ctrl_col, edgecolor="k", linewidth=0.4, label=group_ctrl)
    ax.set_xticks(x, [c_ if len(c_) < 26 else c_[:24] + "…" for c_ in top_combos],
                  rotation=35, ha="right", fontsize=7.5)
    ax.set_ylabel("% of group cells")
    ax.set_title(f"{which}: top-8 pseudo-clonotypes (receptor-chain combinations)", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    save3(fig, f"{which}_top_clonotypes")


analyze("GSE151263", "ARDS", "Sepsis")
analyze("GSE167363", "Sepsis", "HC")

md = """# TCR/BCR 克隆性分析（转录本水平）方法学说明

## 背景与局限（如实声明）
- GSE151263 与 GSE167363 的 GEO supplementary 均为 **基因表达（GEX）计数矩阵**，
  不含 VDJ 富集文库（TCR/BCR CDR3 序列）。
- 因此无法定义 CDR3 水平的真实克隆型、也无法计算正规克隆性指标（如
  Scirpy/Immcantation 的 clone size / CDR3 谱系重叠）。
- 本分析采用既定的替代口径：**受体链转录本表达组合 = 伪克隆型（pseudo-clonotype）**，
  可比较组间相对克隆扩增与多样性，但同一"伪克隆型"内的不同真实克隆无法区分；
  结果应解释为"受体谱系层面的克隆性代理"，不能作为真正的克隆型鉴定。

## 分析口径
1. TCR 链：TRAC/TRBC1/TRBC2/TRDC/TRGC1/TRGV9；BCR 链：IGHM/IGHG1-3/IGKC/IGLC2/IGLC3/IGLL1/IGLL5。
2. B 系细胞（B/浆细胞或仅 BCR 链阳性）用 BCR 链组合，其余用 TCR 链组合；无任一链转录 → "none"。
3. 指标（每样本）：丰富度（不同组合数）、Shannon 熵、Simpson、Gini、top-10 占比；
   组间 Mann-Whitney U 检验。
4. rank-abundance 曲线按组合并后取 top 50。

## 预期生物学
- GSE151263（PBMC）：T/B 细胞占约 40%，TRAC/TRBC1/IGKC 检出预期丰富；
  脓毒症常见 T 细胞克隆性扩张倾向。
- GSE167363（CD14+ 单核）：非 T/B 细胞群，受体链转录本预期接近零——
  该数据集的"克隆性"结论以"未检出/背景水平"如实报告，不强行解读。
"""
with open(os.path.join(OUT, "TCR_BCR_clonality_methods.md"), "w") as fh:
    fh.write(md)
print("\nall done.")
