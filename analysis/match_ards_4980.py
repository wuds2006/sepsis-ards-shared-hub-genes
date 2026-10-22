#!/usr/bin/env python3
"""用 sepsis_sig_4980.tsv 的 chrom:pos 匹配 ARDS 结局(向量化 merge)"""
import pandas as pd

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT = f"{BASE}/04_MR_SMR"
ARDS = f"{BASE}/00_rawdata/gwas/finngen_R11_J10_ARDS.gz"
ards_cols = ['#chrom','pos','ref','alt','rsids','pval','beta','sebeta','af_alt']

sig = pd.read_csv(f"{OUT}/sepsis_sig_4980.tsv", sep='\t')
keys = sig[['chrom','pos']].drop_duplicates()
print("sepsis 工具变量:", len(keys))

matched = []
for c in pd.read_csv(ARDS, sep='\t', compression='gzip', usecols=ards_cols, chunksize=3000000):
    c.columns = ['chrom','pos','ref','alt','rsids','pval','beta','sebeta','af_alt']
    m = c.merge(keys, on=['chrom','pos'], how='inner')
    if len(m): matched.append(m)
ards_match = pd.concat(matched, ignore_index=True)
print("ARDS 结局匹配:", len(ards_match))
ards_match.to_csv(f"{OUT}/ards_match_4980.tsv", sep='\t', index=False)
print("DONE")
