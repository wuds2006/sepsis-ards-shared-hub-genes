#!/usr/bin/env python3
"""
Step3d: CellOracle 虚拟敲除（本地完整环境运行脚本）
======================================================
用途: 对 Hub 基因做 CellOracle 基序-based GRN + in silico perturbation 虚拟敲除

⚠️ 为何无法在云端沙箱运行:
  1. pycistarget 不在 PyPI, 需从 GitHub 安装(沙箱 GitHub 访问受限)
  2. velocyto 需编译 Cython 扩展
  3. 需下载 cisTarget 基序数据库(hg38, 数百MB)
  4. 14万细胞的 GRN 构建需 >16GB 内存

✅ 请在你的本地 macOS 环境(有 R + 可访问 GitHub + 足够内存)按以下步骤运行。
"""

# ============ 0. 环境准备（一次性） ============
# 建议用 conda 建独立环境:
#   conda create -n celloracle python=3.9 -y
#   conda activate celloracle
#   pip install celloracle
#   # 若 velocyto 报错, 可先: pip install numpy cython && pip install velocyto --no-build-isolation
#   pip install pycistarget  # 需 R + RcisTarget: Rscript -e 'BiocManager::install(c("RcisTarget","AUCell"))'

import os, sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
H5AD = os.path.join(BASE, "08_singlecell", "GSE151263_scanpy.h5ad")
OUT  = os.path.join(BASE, "08_singlecell", "celloracle")
os.makedirs(OUT, exist_ok=True)

HUB_GENES = ["CDK5RAP2","CD177","GRB10","UPP1","EXOSC4","VNN1"]
# 注意: 仅 UPP1(5336细胞) 和 EXOSC4(668细胞) 在 PBMC 单核细胞有可检测表达;
#       CD177 是中性粒标志(PBMC不含中性粒, 表达=0), 其余为低表达胞内蛋白。
KNOCKOUT_TARGETS = ["UPP1", "EXOSC4"]  # 实际可敲除的 Hub 基因

def main():
    import celloracle as co
    from celloracle.applications import Oracle_development_module as dev

    # ============ 1. 加载数据 ============
    adata = sc.read(H5AD)
    print("载入细胞数:", adata.n_obs, "基因数:", adata.n_vars)

    # ============ 2. 下载 base GRN (hg38, 按需选一种) ============
    # 完整 base GRN (需 ~2GB 下载 + 内存):
    #   https://celloracle.org/ 的 base_GRN 下载页
    # 这里用官方示例(需先手动下载 base_GRN 文件放到 OUT/)
    base_grn_path = os.path.join(OUT, "base_GRN.txt")
    if not os.path.exists(base_grn_path):
        print("⚠️ 请先从 https://celloracle.org/ 下载 hg38 base GRN 到:", base_grn_path)
        print("   (或运行 co.data.load_base_grn() 在线获取)")
        # 也可在线加载(需网络):
        # base_grn = co.data.load_base_grn()  # 交互式下载
        return

    base_grn = pd.read_parquet(base_grn_path) if base_grn_path.endswith(".parquet") \
               else pd.read_csv(base_grn_path, sep="\t")

    # ============ 3. 初始化 Oracle 对象 ============
    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(
        adata=adata,
        cluster_column_name="celltype",
        embedding_name="X_umap")
    oracle.import_TF_data(TF_info_matrix=base_grn)

    # ============ 4. 构建 GRN (基序-based) ============
    oracle.perform_PCA()
    n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))) > 0.002)[0][0]
    n_comps = min(n_comps, 50)
    oracle.knn_imputation(n_pca_dims=n_comps, k=30, balanced=True, b_sight=3000,
                          b_maxl=3000, n_jobs=8)
    oracle.fit_GRN_for_simulation(alpha=10, verbose_level=1)

    links = oracle.get_links()
    links.to_csv(os.path.join(OUT, "celloracle_GRN_links.csv"), index=False)
    print("GRN 边数:", links.shape[0])

    # ============ 5. 虚拟敲除 ============
    results = {}
    for g in KNOCKOUT_TARGETS:
        if g not in oracle.gene_expression_matrix.columns:
            print(f"⚠️ {g} 不在表达矩阵, 跳过")
            continue
        print(f"===== 虚拟敲除 {g} =====")
        oracle.simulate_shift(
            perturb_condition={g: 0.0},  # 敲低到 0
            n_propagation=3)
        oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
        oracle.calculate_embedding_shift(sigma_corr=0.05)
        shift = oracle.get_embedding_shift()
        results[g] = shift
        shift.to_csv(os.path.join(OUT, f"CellOracle_KO_{g}_embedding_shift.csv"), index=False)
        # 可视化
        fig = oracle.plot_whole_embedding(scale=0.5, show=False)
        fig.savefig(os.path.join(OUT, f"CellOracle_KO_{g}_shift.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"{g} 敲除完成")

    print("\n===== CellOracle 双敲除完成 =====")
    for g, s in results.items():
        print(f"{g}: 嵌入位移向量维度 {s.shape}")

if __name__ == "__main__":
    main()
