#!/usr/bin/env python3
"""
coloc 数据对齐：eQTL(目标基因) × GWAS(sepsis/ARDS)
- 按 SNP 匹配, 对齐等位基因(翻转 beta), 计算 MAF
- 输出每个基因的 eQTL + GWAS 对齐数据
"""
import pandas as pd
import numpy as np
import os

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
EQTL_DIR = f"{BASE}/04_MR_SMR/eqtl_query"
OUT_DIR = f"{BASE}/04_MR_SMR/coloc_input"
os.makedirs(OUT_DIR, exist_ok=True)

GENES = {
    "CDK5RAP2": "ENSG00000136861",
    "CD177": "ENSG00000204936",
    "GRB10": "ENSG00000106070",
    "UPP1": "ENSG00000183696",
    "VNN1": "ENSG00000112299",
}
EQTL_N = 30596  # eQTLGen 血液样本量

def load_eqtl(gene, ensg):
    """读目标基因的 cis-eQTL (Probe == ensg)"""
    df = pd.read_csv(f"{EQTL_DIR}/eqtl_{gene}.txt", sep="\t")
    df = df[df['Probe'] == ensg].copy()
    # SNP Chr BP A1 A2 Freq Probe ... b SE p
    df['maf'] = df['Freq'].apply(lambda x: min(x, 1-x))
    df['e_allele'] = df['A1']
    df['e_beta'] = df['b']
    df['e_se'] = df['SE']
    df['e_maf'] = df['maf']
    return df[['SNP','Chr','BP','e_allele','A2','e_beta','e_se','e_maf','p']].rename(
        columns={'A2':'e_other','p':'e_p'})

def match_gwas(eqtl, gwas_ma_path):
    """从 GWAS .ma 分块匹配 SNP, 对齐等位基因"""
    snps = set(eqtl['SNP'])
    gwas_cols = ['SNP','A1','A2','freq','b','se','p']
    matched = []
    for c in pd.read_csv(gwas_ma_path, sep='\t', chunksize=3000000, usecols=gwas_cols):
        c = c[c['SNP'].isin(snps)]
        if len(c):
            matched.append(c)
    if not matched:
        return pd.DataFrame()
    gwas = pd.concat(matched, ignore_index=True)
    # 对齐等位基因
    merged = eqtl.merge(gwas, on='SNP', how='inner')
    # 判断是否需要翻转: e_allele(A1) vs GWAS A1
    flip = merged['e_allele'] != merged['A1']
    merged['g_beta'] = np.where(flip, -merged['b'], merged['b'])
    merged['g_maf'] = np.where(flip, 1-merged['freq'], merged['freq'])
    merged['g_maf'] = merged['g_maf'].apply(lambda x: min(x, 1-x))
    merged['g_se'] = merged['se']
    merged['g_p'] = merged['p']
    return merged[['SNP','e_beta','e_se','e_maf','e_p','g_beta','g_se','g_maf','g_p']]

# 处理
for gene, ensg in GENES.items():
    eqtl = load_eqtl(gene, ensg)
    print(f"=== {gene}: {len(eqtl)} 个 cis-eQTL ===")
    for gwas_name, gwas_path in [
        ("sepsis", f"{BASE}/04_MR_SMR/sepsis_4980.ma"),
        ("ards", f"{BASE}/04_MR_SMR/ards_J10.ma"),
    ]:
        aligned = match_gwas(eqtl, gwas_path)
        out = f"{OUT_DIR}/{gene}_{gwas_name}_coloc.txt"
        aligned.to_csv(out, sep="\t", index=False)
        print(f"  {gene}-{gwas_name}: 对齐 {len(aligned)} 个 SNP -> {out}")

print("\n===== coloc 数据对齐完成 =====")
