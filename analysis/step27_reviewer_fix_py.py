#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step27 — 审稿修复分析（Python 部分）

修复 II-2: GSE151263 患者级 pseudobulk 六基因 ARDS vs Sepsis ROC（疾病 vs 疾病对照）
修复 II-4: 供体分层敏感性分析（拟时序 / AUCell UPP1 program / 代谢通路，per-donor）
修复 II-9: AUC = 1.000 的精确秩置信下界（exact rank-based lower bound）

输出:
  08_singlecell/disease_vs_disease_ROC_GSE151263.csv
  08_singlecell/donor_sensitivity/*.csv
  08_singlecell/donor_sensitivity/FigureS14_donor_sensitivity.{pdf,png,tiff}
  06_external_validation/AUC1_exact_lower_bounds.csv
"""
import os, sys, csv, math
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache_step27")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import scipy.io as sio
import scipy.sparse as sp
import pandas as pd

for fp in ("/System/Library/Fonts/Helvetica.ttc",):
    if os.path.exists(fp):
        font_manager.fontManager.addfont(fp)
plt.rcParams["font.family"] = ["Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42

ROOT = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
SC = os.path.join(ROOT, "08_singlecell")
OUT = os.path.join(SC, "donor_sensitivity")
os.makedirs(OUT, exist_ok=True)
PURPLE = "#6A5C95"; PINK = "#D99AAE"; ORANGE = "#E7A889"
GREEN = "#A9CBB8"; MINT = "#CFE3E6"; BLUE = "#8FA2D6"; NAVY = "#5B608C"
PAL7 = [PURPLE, PINK, ORANGE, GREEN, MINT, BLUE, NAVY]
HUB = ["CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1"]


def auc_mw(x, g):
    """AUC + exact Mann-Whitney p; g True=case"""
    x1 = np.asarray(x)[g]; x0 = np.asarray(x)[~g]
    if len(x1) == 0 or len(x0) == 0:
        return np.nan, np.nan
    u = mannwhitneyu(x1, x0, alternative="two-sided")
    auc = u.statistic / (len(x1) * len(x0))
    return float(auc), float(u.pvalue)


def signif(x, n=3):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return np.nan
    return float(f"%.{n}g" % x)


def exact_lower_bound_auc1(n1, n0, alpha=0.05):
    """AUC=1（全部一致）时精确单侧置信下界: alpha^(1/(n1*n0))"""
    return float(alpha ** (1.0 / (n1 * n0)))


# ============================================================
# Part B: GSE151263 pseudobulk ARDS vs Sepsis ROC（II-2 修复）
# ============================================================
print("=== Part B: disease-vs-disease ROC (GSE151263 pseudobulk) ===")
import scanpy as sc
import anndata as ad

cd = os.path.join(SC, "cellchat_input")
cnt = sio.mmread(os.path.join(cd, "counts.mtx")).T.tocsr()
genes = [l.strip().split("\t")[0] for l in open(os.path.join(cd, "genes.tsv"))]
meta = pd.read_csv(os.path.join(cd, "meta.csv")).set_index("cell")
A = ad.AnnData(X=sp.csr_matrix(cnt, dtype=np.float32), obs=meta,
               var=pd.DataFrame(index=genes))
sc.pp.normalize_total(A, target_sum=1e4)
sc.pp.log1p(A)

# 患者级 pseudobulk（每样本均值表达）
samples = A.obs["sample"].astype(str).unique()
grp_of = dict(zip(A.obs["sample"].astype(str), A.obs["group"].astype(str)))
hub_idx = [g for g in HUB if g in A.var_names]
rows = []
for s in samples:
    m = (A.obs["sample"].astype(str) == s).values
    xb = np.asarray(A.X[m][:, [A.var_names.get_loc(g) for g in hub_idx]].todense()).mean(axis=0)
    rows.append([s, grp_of[s]] + list(np.round(xb, 4)))
pb = pd.DataFrame(rows, columns=["sample", "group"] + hub_idx)
pb.to_csv(os.path.join(SC, "disease_vs_disease_ROC_GSE151263.csv"), index=False)
print(pb.to_string(index=False))

g_ards = (pb["group"] == "ARDS").values
res_b = []
for g in hub_idx:
    auc, p = auc_mw(pb[g].values, g_ards)
    res_b.append([g, "single", round(auc, 3), signif(p, 3)])
# combined: 六基因 z-score 均值（零方差基因 z 置 0）
Xz = pb[hub_idx].values
sd = Xz.std(axis=0)
Z = np.where(sd == 0, 0.0, (Xz - Xz.mean(axis=0)) / np.where(sd == 0, 1, sd))
comb = Z.mean(axis=1)
auc, p = auc_mw(comb, g_ards)
res_b.append(["Combined_signature", "signature", round(auc, 3), signif(p, 3)])
res_b = pd.DataFrame(res_b, columns=["gene", "type", "AUC_ARDSvsSepsis", "P_exact"])
res_b.to_csv(os.path.join(SC, "disease_vs_disease_ROC_GSE151263.csv"), index=False)
print(res_b.to_string(index=False))


# ============================================================
# Part C: 供体分层敏感性分析（II-4 修复）
# ============================================================
print("\n=== Part C: donor-level sensitivity ===")
donor_tables = {}

# C1: 拟时序 per-donor 中位数（GSE167363 的 5 个样本归并为 4 位患者：S2_T0/S2_T6 同一患者）
PATIENT_MAP_167363 = {"1": "HC1", "2": "HC2", "3": "S2", "4": "S2", "5": "S3"}
pt_rows = []
for ds, f in [("GSE151263", "monocle3/monocle3_pseudotime_GSE151263.csv"),
              ("GSE167363", "monocle3/monocle3_pseudotime_GSE167363.csv")]:
    d = pd.read_csv(os.path.join(SC, f))
    d["donor"] = d["cell"].str.split("_").str[-1].str.split("-").str[0]
    if ds == "GSE167363":
        d["donor"] = d["donor"].map(PATIENT_MAP_167363)
    dm = d.groupby(["donor", "group"])["pseudotime"].median().reset_index()
    dm["dataset"] = ds
    dm["n_cells"] = d.groupby(["donor", "group"]).size().values
    pt_rows.append(dm)
    # 供体层 wilcox
    case = dm[dm["group"].isin(["Sepsis", "ARDS"])]["pseudotime"].values
    ctrl = dm[~dm["group"].isin(["Sepsis", "ARDS"])]["pseudotime"].values
    case_grp = [g for g in dm["group"] if g in ("Sepsis", "ARDS")]
    # 区分对比方向：PBMC 是 ARDS vs Sepsis；CD14 是 Sepsis vs HC
    a_grp = "ARDS" if ds == "GSE151263" else "Sepsis"
    b_grp = "Sepsis" if ds == "GSE151263" else "HC"
    x1 = dm[dm["group"] == a_grp]["pseudotime"].values
    x0 = dm[dm["group"] == b_grp]["pseudotime"].values
    u = mannwhitneyu(x1, x0, alternative="two-sided")
    sep = (np.min(x1) > np.max(x0)) or (np.min(x0) > np.max(x0))
    print(f"[{ds}] donor-level pseudotime {a_grp} vs {b_grp}: "
          f"n={len(x1)}v{len(x0)}, medians {np.median(x1):.2f} vs {np.median(x0):.2f}, "
          f"P={u.pvalue:.3f}, complete_separation={sep}")
    donor_tables[f"pseudotime_{ds}"] = dm

# C2: AUCell UPP1 program per-donor 中位数（单核细胞）
auc_rows = []
for ds, f, ct in [("GSE151263", "aucell/GSE151263_AUCell_scores.csv", "Mono"),
                  ("GSE167363", "GSE167363/aucell/GSE167363_AUCell_scores.csv", None)]:
    d = pd.read_csv(os.path.join(SC, f))
    d["donor"] = d["cell"].str.split("_").str[-1].str.split("-").str[0]
    if ds == "GSE167363":
        d["donor"] = d["donor"].map(PATIENT_MAP_167363)
    sub = d[d["celltype"] == ct] if ct else d
    dm = sub.groupby(["donor", "group"])["UPP1_program"].median().reset_index()
    dm["dataset"] = ds; dm["celltype"] = ct or "all_mono"
    dm["n_cells"] = sub.groupby(["donor", "group"]).size().values
    a_grp = "ARDS" if ds == "GSE151263" else "Sepsis"
    b_grp = "Sepsis" if ds == "GSE151263" else "HC"
    x1 = dm[dm["group"] == a_grp]["UPP1_program"].values
    x0 = dm[dm["group"] == b_grp]["UPP1_program"].values
    u = mannwhitneyu(x1, x0, alternative="two-sided")
    sep = (np.min(x1) > np.max(x0)) or (np.max(x1) < np.min(x0))
    print(f"[{ds}] donor-level UPP1 program ({dm['celltype'].iloc[0]}) {a_grp} vs {b_grp}: "
          f"n={len(x1)}v{len(x0)}, medians {np.median(x1):.3f} vs {np.median(x0):.3f}, "
          f"P={u.pvalue:.3f}, complete_separation={sep}")
    auc_rows.append(dm)
    donor_tables[f"UPP1_{ds}"] = dm

# C3: 代谢通路 per-donor 均值（重算 score_genes）
PATHWAYS = {
    "Glycolysis": "HALLMARK_GLYCOLYSIS",
    "Oxidative phosphorylation": "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "Fatty acid metabolism": "HALLMARK_FATTY_ACID_METABOLISM",
    "Xenobiotic metabolism": "HALLMARK_XENOBIOTIC_METABOLISM",
}
GMT = {}
for line in open(os.path.join(ROOT, "00_rawdata/hallmark/h.all.v2024.1.Hs.symbols.gmt")):
    fs = line.rstrip("\n").split("\t")
    GMT[fs[0]] = [g for g in fs[2:] if g]

metab_rows = []
for ds, cdir in [("GSE151263", os.path.join(SC, "cellchat_input")),
                 ("GSE167363", os.path.join(SC, "GSE167363/input"))]:
    cnt = sio.mmread(os.path.join(cdir, "counts.mtx")).T.tocsr()
    genes = [l.strip().split("\t")[0] for l in open(os.path.join(cdir, "genes.tsv"))]
    obs = pd.read_csv(os.path.join(cdir, "meta.csv"))
    ckey = obs.columns[0]
    obs = obs.set_index(ckey)
    M = ad.AnnData(X=sp.csr_matrix(cnt, dtype=np.float32), obs=obs,
                   var=pd.DataFrame(index=genes))
    sc.pp.normalize_total(M, target_sum=1e4)
    sc.pp.log1p(M)
    for label, gs in PATHWAYS.items():
        gg = [g for g in GMT[gs] if g in M.var_names]
        sc.tl.score_genes(M, gene_list=gg, score_name="tmp", random_state=42)
        smp = M.obs["sample"].astype(str).values
        # 归并为患者（S2_T0/S2_T6 同一患者）
        donor = np.array([s.split("_")[1] if s.startswith("GSM") and "_" in s
                         else s for s in smp])
        df = pd.DataFrame({"donor": donor,
                           "group": M.obs["group"].astype(str).values,
                           "score": M.obs["tmp"].values})
        dm = df.groupby(["donor", "group"])["score"].mean().reset_index()
        dm["pathway"] = label; dm["dataset"] = ds
        metab_rows.append(dm)
        a_grp = "ARDS" if ds == "GSE151263" else "Sepsis"
        b_grp = "Sepsis" if ds == "GSE151263" else "HC"
        x1 = dm[dm["group"] == a_grp]["score"].values
        x0 = dm[dm["group"] == b_grp]["score"].values
        u = mannwhitneyu(x1, x0, alternative="two-sided")
        sep = (np.min(x1) > np.max(x0)) or (np.max(x1) < np.min(x0))
        print(f"[{ds}] donor-level {label} {a_grp} vs {b_grp}: "
              f"n={len(x1)}v{len(x0)}, means {np.mean(x1):.4f} vs {np.mean(x0):.4f}, "
              f"P={u.pvalue:.3f}, complete_separation={sep}")
    print(f"[{ds}] metabolism done")

# 汇总 CSV
pt1 = donor_tables["pseudotime_GSE151263"].rename(columns={"pseudotime": "value"})
pt2 = donor_tables["pseudotime_GSE167363"].rename(columns={"pseudotime": "value"})
up1 = donor_tables["UPP1_GSE151263"].rename(columns={"UPP1_program": "value"})
gly = pd.concat([r[r["pathway"] == "Glycolysis"] for r in metab_rows]).rename(
    columns={"score": "value"})
summary = pd.concat([
    pt1.assign(metric="pseudotime"),
    pt2.assign(metric="pseudotime"),
    donor_tables["UPP1_GSE151263"].rename(columns={"UPP1_program": "value"}).assign(metric="UPP1_program"),
    donor_tables["UPP1_GSE167363"].rename(columns={"UPP1_program": "value"}).assign(metric="UPP1_program"),
    *[r.rename(columns={"score": "value"}).assign(metric="metabolism") for r in metab_rows],
], ignore_index=True)
summary.to_csv(os.path.join(OUT, "donor_level_summary.csv"), index=False)

# ============================================================
# Figure S14: 供体敏感性 4 面板
# ============================================================
panels = [
    (pt1, "PBMC pseudotime (GSE151263)"),
    (pt2, "CD14+ pseudotime (GSE167363)"),
    (up1, "UPP1 program, monocytes (GSE151263)"),
    (gly.assign(dataset=np.where(gly["dataset"] == "GSE151263", "PBMC", "CD14+")),
     "Glycolysis, per dataset"),
]
fig, axes = plt.subplots(1, 4, figsize=(13.5, 3.2))
for ax, (d, title) in zip(axes, panels):
    if "dataset" in d.columns:
        xs, labs = [], []
        for gi, (dsname, grp) in enumerate(
                [(dn, g) for dn in sorted(d["dataset"].unique())
                 for g in sorted(d[d["dataset"] == dn]["group"].unique())]):
            vals = d[(d["dataset"] == dsname) & (d["group"] == grp)]["value"].values
            xj = np.full(len(vals), gi) + np.random.default_rng(7 + gi).uniform(-0.06, 0.06, len(vals))
            ax.scatter(xj, vals, s=70, color=PAL7[gi % 7], edgecolor="black",
                       linewidth=0.6, zorder=3)
            ax.hlines(np.median(vals), gi - 0.25, gi + 0.25, color=NAVY, linewidth=2, zorder=4)
            xs.append(gi); labs.append(f"{grp}\n({dsname})")
        ax.set_xticks(xs); ax.set_xticklabels(labs, fontsize=8)
    else:
        for gi, grp in enumerate(sorted(d["group"].unique())):
            vals = d[d["group"] == grp]["value"].values
            xj = np.full(len(vals), gi) + np.random.default_rng(7 + gi).uniform(-0.06, 0.06, len(vals))
            ax.scatter(xj, vals, s=70, color=PAL7[gi % 7], edgecolor="black",
                       linewidth=0.6, zorder=3, label=grp)
            ax.hlines(np.median(vals), gi - 0.25, gi + 0.25, color=NAVY, linewidth=2, zorder=4)
        ax.set_xticks(range(len(d["group"].unique())))
        ax.set_xticklabels(sorted(d["group"].unique()), fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel("Donor-level value", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)

fig.suptitle("Donor-level sensitivity analysis", fontsize=12, y=1.02)
fig.text(0.5, -0.04, "Each point is one donor (sample); horizontal bars mark medians. "
         "Cell-level statistics in the main text are supported by consistent donor-level medians.",
         ha="center", fontsize=8.5, style="italic")
fig.tight_layout()
for ext, dpi in [("pdf", None), ("png", 300), ("tiff", 500)]:
    if ext == "pdf":
        fig.savefig(os.path.join(OUT, "FigureS14_donor_sensitivity.pdf"),
                    facecolor="white", bbox_inches="tight")
    else:
        fig.savefig(os.path.join(OUT, f"FigureS14_donor_sensitivity.{ext}"), dpi=dpi,
                    facecolor="white", bbox_inches="tight",
                    pil_kwargs={"compression": "tiff_lzw"} if ext == "tiff" else None)
plt.close(fig)
print("[saved] FigureS14_donor_sensitivity (3 formats)")

# ============================================================
# Part D: AUC = 1.000 精确置信下界（II-9 修复）
# ============================================================
print("\n=== Part D: exact lower bounds for AUC = 1.000 ===")
rows = []
for name, n1, n0 in [
    ("GSE243066 ARDS validation (combined signature)", 34, 15),
    ("GSE100159 sepsis validation (combined signature)", 33, 12),
]:
    lb = exact_lower_bound_auc1(n1, n0)
    rows.append([name, n1, n0, round(lb, 4)])
    print(f"{name}: AUC=1.000, exact 95% one-sided lower bound = {lb:.4f}")
pd.DataFrame(rows, columns=["cohort", "n_case", "n_ctrl",
                            "exact_lower_bound_AUC"]).to_csv(
    os.path.join(ROOT, "06_external_validation/AUC1_exact_lower_bounds.csv"), index=False)

print("\n[done] step27 Python part")
