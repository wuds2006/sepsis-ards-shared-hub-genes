#!/usr/bin/env python3
"""
Step 3b 数据准备 —— 从原始 UMI 提取 raw counts + 细胞注释, 供 CellChat(R) 使用
复现 step3a 的过滤逻辑(同 seed=42 下采样 + 同 QC), 保证与已有 celltype 注释一致
"""
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import os

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
PAPER = "/Users/mac/Desktop/2026年撰写的文章/2026论文撰写一Acute Respiratory Distress Syndrome and Sepsis"
OUT = f"{BASE}/08_singlecell/cellchat_input"
os.makedirs(OUT, exist_ok=True)

raw_dir = f"{PAPER}/ARDS/单细胞/GSE151263_RAW"
files = sorted([f for f in os.listdir(raw_dir) if f.endswith("_processed_UMI.txt.gz")])
print("样本文件:", files)

# 已有注释(基于 step3a)
annot = pd.read_csv(f"{BASE}/08_singlecell/GSE151263_cell_annot.csv", index_col=0)
print("已有注释细胞数:", len(annot), "| celltype:", sorted(annot['celltype'].unique()))

parts = []
for i, f in enumerate(files):
    sample = f.replace("_processed_UMI.txt.gz", "")
    group = "ARDS" if "ARDS" in sample else "Sepsis"
    print(f"读取 {sample} (batch={i}) ...")
    df = pd.read_csv(f"{raw_dir}/{f}", sep="\t", index_col=0)
    df = df.groupby(level=0).sum()
    if df.shape[1] > 2000:
        np.random.seed(42)
        keep = np.random.choice(df.shape[1], 2000, replace=False)
        df = df.iloc[:, sorted(keep)]
    # 复现 scanpy concatenate: barcode 加 -{batch_index} 后缀
    df.columns = [f"{bc}-{i}" for bc in df.columns.astype(str)]
    parts.append((sample, group, df))

# 合并所有样本到统一基因空间(inner join)
all_genes = None
for sample, group, df in parts:
    all_genes = df.index if all_genes is None else all_genes.intersection(df.index)
all_genes = sorted(all_genes)
print("统一基因数:", len(all_genes))

# 构造合并 counts 矩阵(基因 x 细胞)
cell_barcodes = []
gene_matrices = []
sample_meta = []
for sample, group, df in parts:
    sub = df.loc[all_genes]
    gene_matrices.append(sp.csr_matrix(sub.values))
    cell_barcodes.extend(list(sub.columns.astype(str)))
    sample_meta.extend([(sample, group)] * sub.shape[1])
    del df

counts = sp.hstack(gene_matrices).tocsr()  # 基因 x 细胞
print("合并 counts:", counts.shape, "(基因 x 细胞)")

# 构造 meta
meta = pd.DataFrame(sample_meta, columns=["sample", "group"], index=cell_barcodes)
meta.index.name = "cell"

# QC 过滤(与 step3a 一致): n_genes > 200, < 6000, pct_mt < 20
n_genes = np.asarray(counts.sum(axis=0)).ravel()  # 每细胞表达基因数(需用 counts>0)
n_genes_per_cell = np.asarray((counts > 0).sum(axis=0)).ravel()
# mt 比例
mt_mask = np.array([g.startswith("MT-") for g in all_genes])
mt_counts = np.asarray(counts[mt_mask, :].sum(axis=0)).ravel()
total_counts = np.asarray(counts.sum(axis=0)).ravel()
pct_mt = mt_counts / np.maximum(total_counts, 1) * 100

qc_keep = (n_genes_per_cell > 200) & (n_genes_per_cell < 6000) & (pct_mt < 20)
print("QC 后细胞数:", qc_keep.sum())

counts = counts[:, qc_keep]
meta = meta.iloc[np.where(qc_keep)[0]]

# 匹配已有注释(用 cell barcode)
annot_match = annot.reindex(meta.index)
matched = annot_match['celltype'].notna()
print(f"注释匹配率: {matched.sum()}/{len(meta)}")

# 只保留有注释的细胞
counts = counts[:, np.where(matched)[0]]
meta = meta.iloc[np.where(matched)[0]]
meta["celltype"] = annot.reindex(meta.index)['celltype'].values
meta["group"] = annot.reindex(meta.index)['group'].values

print("最终: counts", counts.shape, "| meta", meta.shape)
print("celltype 分布:\n", meta['celltype'].value_counts())
print("group 分布:\n", meta['group'].value_counts())

# 导出
sp.save_npz(f"{OUT}/counts.npz", counts)           # 基因 x 细胞, 稀疏
pd.Series(all_genes).to_csv(f"{OUT}/genes.csv", index=False, header=["gene"])
meta.to_csv(f"{OUT}/meta.csv")
print("已导出 counts.npz + genes.csv + meta.csv")

# 同时导出 10x 风格 mtx(供 Seurat Read10X 或 R Matrix::readMM)
from scipy.io import mmwrite
mmwrite(f"{OUT}/counts.mtx", counts)               # 基因 x 细胞
with open(f"{OUT}/barcodes.tsv", "w") as f:
    f.write("\n".join(meta.index) + "\n")
with open(f"{OUT}/genes.tsv", "w") as f:
    f.write("\n".join(all_genes) + "\n")
print("已导出 10x 风格 mtx")
print("\n===== Step3b 数据准备完成 =====")
