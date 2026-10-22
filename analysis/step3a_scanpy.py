#!/usr/bin/env python3
"""
Step 3a —— GSE151263 单细胞 QC + 聚类 + 注释 + Hub基因表达 (scanpy)
"""
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import os, gzip

sc.settings.verbosity = 1
BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
PAPER = "/Users/mac/Desktop/2026年撰写的文章/2026论文撰写一Acute Respiratory Distress Syndrome and Sepsis"
OUT = f"{BASE}/08_singlecell"
os.makedirs(OUT, exist_ok=True)

raw_dir = f"{PAPER}/ARDS/单细胞/GSE151263_RAW"
files = sorted([f for f in os.listdir(raw_dir) if f.endswith("_processed_UMI.txt.gz")])
print("样本文件:", files)

hub_genes = ["CDK5RAP2","CD177","GRB10","UPP1","EXOSC4","VNN1"]

adatas = []
for f in files:
    sample = f.replace("_processed_UMI.txt.gz", "")
    group = "ARDS" if "ARDS" in sample else "Sepsis"
    print(f"读取 {sample} ...")
    # 读 基因x细胞 矩阵
    df = pd.read_csv(f"{raw_dir}/{f}", sep="\t", index_col=0)
    # 聚合重复基因
    df = df.groupby(level=0).sum()
    # 下采样细胞
    if df.shape[1] > 2000:
        np.random.seed(42)
        keep = np.random.choice(df.shape[1], 2000, replace=False)
        df = df.iloc[:, sorted(keep)]
    # 转 AnnData (细胞 x 基因)
    X = sp.csr_matrix(df.values.T)   # 稀疏
    ad = sc.AnnData(X, var=pd.DataFrame(index=df.index))
    ad.obs_names = df.columns.astype(str)
    ad.obs["sample"] = sample
    ad.obs["group"] = group
    adatas.append(ad)
    del df; 

# 合并
adata = adatas[0].concatenate(adatas[1:], batch_key="sample")
print(f"合并后: {adata.shape[0]} 细胞 x {adata.shape[1]} 基因")

# QC
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.n_genes_by_counts < 6000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]
print(f"QC 后: {adata.shape[0]} 细胞")

# 归一化 + HVG + PCA
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")

# 批次校正(harmony) - 尝试导入
try:
    import harmonypy as hm
    ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, "sample")
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] == adata.n_obs:
        adata.obsm["X_pca_harmony"] = Z
    else:
        adata.obsm["X_pca_harmony"] = Z.T
    use_rep = "X_pca_harmony"
    print("使用 harmony 批次校正, shape:", adata.obsm["X_pca_harmony"].shape)
except Exception as e:
    use_rep = "X_pca"
    print("无 harmony, 用 PCA:", e)

# 邻居 + UMAP + Leiden
sc.pp.neighbors(adata, use_rep=use_rep, n_pcs=20)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

# 细胞注释(标记基因打分)
markers = {
    "CD4_T": ["CD3D","CD3E","CD4","IL7R"],
    "CD8_T": ["CD3D","CD8A","CD8B"],
    "NK": ["NKG7","GNLY","NCAM1","KLRD1"],
    "B": ["MS4A1","CD79A","CD19"],
    "Mono": ["CD14","LYZ","FCN1","S100A8"],
    "DC": ["FCER1A","CST3","CLEC9A"],
    "Platelet": ["PPBP","PF4"],
}
for ct, genes in markers.items():
    sc.tl.score_genes(adata, [g for g in genes if g in adata.var_names], score_name=f"score_{ct}")
score_cols = [f"score_{ct}" for ct in markers]
adata.obs["celltype"] = adata.obs[score_cols].idxmax(axis=1).str.replace("score_","")
print("细胞类型分布:\n", adata.obs["celltype"].value_counts())

# 保存
adata.write(f"{OUT}/GSE151263_scanpy.h5ad")
adata.obs.to_csv(f"{OUT}/GSE151263_cell_annot.csv")

# Hub 基因表达(按细胞类型 + 按组)
hub_avail = [g for g in hub_genes if g in adata.var_names]
print("可验证 Hub 基因:", hub_avail)
if hub_avail:
    # 按细胞类型均值
    expr_ct = pd.DataFrame(index=adata.obs["celltype"].unique(), columns=hub_avail)
    for g in hub_avail:
        expr_ct[g] = adata.obs.groupby("celltype").apply(lambda x: adata[x.index, g].X.mean())
    expr_ct.to_csv(f"{OUT}/GSE151263_hub_celltype_mean.csv")

    # 按组(ARDS vs Sepsis)均值
    expr_grp = pd.DataFrame(index=["ARDS","Sepsis"], columns=hub_avail)
    for g in hub_avail:
        expr_grp[g] = adata.obs.groupby("group").apply(lambda x: adata[x.index, g].X.mean())
    expr_grp.to_csv(f"{OUT}/GSE151263_hub_group_mean.csv")
    print("\nHub 基因按组均值:\n", expr_grp.round(3))

print("\n===== Step3a scanpy 完成 =====")
