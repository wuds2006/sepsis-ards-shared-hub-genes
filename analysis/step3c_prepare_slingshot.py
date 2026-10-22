#!/usr/bin/env python3
"""
Step3c 数据准备: Mono 细胞亚聚类 + 巨噬细胞极化注释, 导出 slingshot 输入
"""
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import os

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT = f"{BASE}/08_singlecell/slingshot"
os.makedirs(OUT, exist_ok=True)

# 用 step3b 导出的 raw counts 提取 Mono 细胞
counts = sp.load_npz(f"{BASE}/08_singlecell/cellchat_input/counts.npz")  # 基因 x 细胞
genes = pd.read_csv(f"{BASE}/08_singlecell/cellchat_input/genes.csv")['gene'].values
meta = pd.read_csv(f"{BASE}/08_singlecell/cellchat_input/meta.csv", index_col=0)

mono_mask = (meta['celltype'] == 'Mono').values
mono_cells = meta.index[mono_mask]
mono_counts = counts[:, mono_mask].T.tocsr()  # 细胞 x 基因
print(f"Mono 细胞数: {mono_counts.shape[0]}")

mono = sc.AnnData(X=mono_counts, var=pd.DataFrame(index=genes))
mono.obs = meta.loc[mono_cells].copy()
mono.obs_names = mono_cells

# 完整流程(raw counts)
sc.pp.normalize_total(mono, target_sum=1e4)
sc.pp.log1p(mono)
sc.pp.highly_variable_genes(mono, n_top_genes=1500)
sc.pp.scale(mono, max_value=10)
sc.tl.pca(mono, svd_solver='arpack')
# 用 harmony 校正(样本批次)
try:
    import harmonypy as hm
    ho = hm.run_harmony(mono.obsm['X_pca'], mono.obs, 'sample')
    Z = np.asarray(ho.Z_corr)
    mono.obsm['X_pca_harmony'] = Z if Z.shape[0] == mono.n_obs else Z.T
    use_rep = 'X_pca_harmony'
except Exception as e:
    use_rep = 'X_pca'
    print("无 harmony:", e)

sc.pp.neighbors(mono, use_rep=use_rep, n_pcs=15)
sc.tl.umap(mono)
sc.tl.leiden(mono, resolution=0.6)

# 巨噬细胞极化标记
m1_genes = ["CD86","CD80","TNF","IL1B","CXCL10","NOS2"]
m2_genes = ["CD163","MRC1","MSR1","TGM2","IL10","CCL18"]
hub_genes = ["CDK5RAP2","CD177","GRB10","UPP1","EXOSC4","VNN1"]

for g in [g for g in m1_genes if g in mono.var_names]:
    sc.tl.score_genes(mono, [g], score_name=f"tmp_{g}")
for g in [g for g in m2_genes if g in mono.var_names]:
    sc.tl.score_genes(mono, [g], score_name=f"tmp_{g}")

# 计算 M1/M2 打分
m1_avail = [g for g in m1_genes if g in mono.var_names]
m2_avail = [g for g in m2_genes if g in mono.var_names]
if m1_avail:
    sc.tl.score_genes(mono, m1_avail, score_name="M1_score")
if m2_avail:
    sc.tl.score_genes(mono, m2_avail, score_name="M2_score")

# 导出
mono.write(f"{OUT}/mono_subcluster.h5ad")

# 导出 slingshot 输入: 降维坐标 + 聚类 + 注释
coords = pd.DataFrame(mono.obsm['X_umap'], columns=['UMAP1','UMAP2'], index=mono.obs_names)
pca = pd.DataFrame(mono.obsm['X_pca'][:, :20], index=mono.obs_names,
                   columns=[f'PC{i+1}' for i in range(20)])
meta = mono.obs[['group','leiden']].copy()
meta.columns = ['group','cluster']

# Hub 基因 + M1/M2 打分
expr = pd.DataFrame(index=mono.obs_names)
for g in [g for g in hub_genes if g in mono.var_names]:
    expr[g] = np.asarray(mono[:, g].X.toarray()).ravel() if hasattr(mono[:, g].X, 'toarray') else np.asarray(mono[:, g].X).ravel()
if 'M1_score' in mono.obs:
    expr['M1_score'] = mono.obs['M1_score'].values
if 'M2_score' in mono.obs:
    expr['M2_score'] = mono.obs['M2_score'].values

coords.to_csv(f"{OUT}/mono_umap.csv")
pca.to_csv(f"{OUT}/mono_pca.csv")
meta.to_csv(f"{OUT}/mono_meta.csv")
expr.to_csv(f"{OUT}/mono_gene_expr.csv")

print("\n=== Mono 亚群分布 ===")
print(mono.obs['leiden'].value_counts().sort_index())
print("\n=== 各亚群 M1/M2 打分 ===")
if 'M1_score' in mono.obs and 'M2_score' in mono.obs:
    print(mono.obs.groupby('leiden')[['M1_score','M2_score']].mean().round(3))
print("\n===== Step3c 数据准备完成 =====")
