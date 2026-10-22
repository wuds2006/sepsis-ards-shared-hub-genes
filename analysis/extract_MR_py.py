#!/usr/bin/env python3
"""MR 数据提取: 按 chrom:pos 精确匹配暴露显著SNP 与结局GWAS"""
import pandas as pd
import gzip, sys

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT = f"{BASE}/04_MR_SMR"
EXPO = f"{BASE}/00_rawdata/gwas/finngen_R11_AB1_OTHER_BACTERIAL.gz"
OUTC = f"{BASE}/00_rawdata/gwas/finngen_R11_J10_ARDS.gz"

cols = ['#chrom','pos','ref','alt','rsids','pval','beta','sebeta','af_alt']

# 1) 暴露显著 SNP (p<5e-8)
print("读取暴露显著SNP...")
expo = pd.read_csv(EXPO, sep='\t', compression='gzip', usecols=cols)
expo.columns = ['chrom','pos','ref','alt','rsids','pval','beta','sebeta','af_alt']
expo_sig = expo[expo['pval'] < 5e-8].copy()
print(f"暴露显著SNP: {len(expo_sig)}")
expo_sig.to_csv(f"{OUT}/expo_sig.tsv", sep='\t', index=False)

# 2) 结局匹配(按 chrom:pos 分块流式)
print("匹配结局...")
keys = set(zip(expo_sig['chrom'], expo_sig['pos']))
matched = []
header = None
for chunk in pd.read_csv(OUTC, sep='\t', compression='gzip', usecols=cols, chunksize=2000000):
    chunk.columns = ['chrom','pos','ref','alt','rsids','pval','beta','sebeta','af_alt']
    m = chunk[chunk.apply(lambda r: (r['chrom'], r['pos']) in keys, axis=1)]
    if len(m) > 0:
        matched.append(m)
    if header is None:
        header = True
outc = pd.concat(matched, ignore_index=True)
print(f"结局匹配SNP: {len(outc)}")
outc.to_csv(f"{OUT}/outc_match.tsv", sep='\t', index=False)
print("DONE")
