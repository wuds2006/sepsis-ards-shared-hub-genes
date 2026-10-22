#!/usr/bin/env python3
"""MR 数据提取: 用 Hamilton sepsis GWAS (ieu-b-4980 = GCST90281174) 作暴露, FinnGen ARDS 作结局
按 chrom:pos 精确匹配"""
import pandas as pd

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT = f"{BASE}/04_MR_SMR"
SEPSIS = f"{BASE}/00_rawdata/gwas/sepsis_ieu-b-4980_GCST90281174.tsv.gz"
ARDS   = f"{BASE}/00_rawdata/gwas/finngen_R11_J10_ARDS.gz"

sepsis_cols = ['chromosome','base_pair_location','effect_allele','other_allele',
               'beta','standard_error','effect_allele_frequency','p_value','rs_id']
ards_cols = ['#chrom','pos','ref','alt','rsids','pval','beta','sebeta','af_alt']

# 1) sepsis 工具变量 SNP (p<5e-6, 因脓毒症无全基因组显著位点)
print("读取 sepsis 工具变量SNP (p<5e-6)...")
chunks = []
for c in pd.read_csv(SEPSIS, sep='\t', compression='gzip', usecols=sepsis_cols, chunksize=2000000):
    c = c[c['p_value'] < 5e-6]
    if len(c): chunks.append(c)
sepsis_sig = pd.concat(chunks, ignore_index=True)
sepsis_sig.columns = ['chrom','pos','effect_allele','other_allele','beta','se','eaf','pval','rsid']
print(f"sepsis 显著SNP: {len(sepsis_sig)}")

# 2) ARDS 结局匹配(按 chrom:pos)
print("匹配 ARDS 结局...")
keys = set(zip(sepsis_sig['chrom'], sepsis_sig['pos']))
matched = []
for c in pd.read_csv(ARDS, sep='\t', compression='gzip', usecols=ards_cols, chunksize=2000000):
    c.columns = ['chrom','pos','ref','alt','rsids','pval','beta','sebeta','af_alt']
    m = c[c.apply(lambda r: (r['chrom'], r['pos']) in keys, axis=1)]
    if len(m): matched.append(m)
ards_match = pd.concat(matched, ignore_index=True)
print(f"ARDS 结局匹配SNP: {len(ards_match)}")

sepsis_sig.to_csv(f"{OUT}/sepsis_sig_4980.tsv", sep='\t', index=False)
ards_match.to_csv(f"{OUT}/ards_match_4980.tsv", sep='\t', index=False)
print("DONE")
