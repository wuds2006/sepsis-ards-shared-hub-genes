"""反向 MR 数据提取: FinnGen ARDS (J10_ARDS) 作暴露, Hamilton sepsis (ieu-b-4980) 作结局
工具变量: p<5e-6 (与正向 MR 同一放宽阈值), rsID 匹配
输出: 04_MR_SMR/ards_sig.tsv (暴露), 04_MR_SMR/sepsis_match_ards.tsv (结局)
"""
import pandas as pd

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
OUT = f"{BASE}/04_MR_SMR"
ARDS = f"{BASE}/00_rawdata/gwas/finngen_R11_J10_ARDS.gz"
SEPSIS = f"{BASE}/00_rawdata/gwas/sepsis_ieu-b-4980_GCST90281174.tsv.gz"

ards_cols = ['#chrom', 'pos', 'ref', 'alt', 'rsids', 'pval', 'beta', 'sebeta', 'af_alt']
sepsis_cols = ['chromosome', 'base_pair_location', 'effect_allele', 'other_allele',
               'beta', 'standard_error', 'effect_allele_frequency', 'p_value', 'rs_id']

# 1) ARDS 工具变量 SNP (p<5e-6)
print("读取 ARDS 工具变量 SNP (p<5e-6)...")
sig = []
for c in pd.read_csv(ARDS, sep='\t', compression='gzip', usecols=ards_cols, chunksize=2000000):
    sig.append(c[c['pval'] < 5e-6])
ards_sig = pd.concat(sig, ignore_index=True)
ards_sig = ards_sig[ards_sig['rsids'].notna() & (ards_sig['rsids'] != '')]
print(f"ARDS 工具变量 SNP: {len(ards_sig)}")
ards_sig.to_csv(f"{OUT}/ards_sig.tsv", sep='\t', index=False)

# 2) sepsis 结局匹配(按 rsid)
print("匹配 sepsis 结局...")
rsids = set(ards_sig['rsids'])
matched = []
for c in pd.read_csv(SEPSIS, sep='\t', compression='gzip', usecols=sepsis_cols, chunksize=2000000):
    matched.append(c[c['rs_id'].isin(rsids)])
sepsis_match = pd.concat(matched, ignore_index=True)
print(f"sepsis 结局匹配 SNP: {len(sepsis_match)}")
sepsis_match.to_csv(f"{OUT}/sepsis_match_ards.tsv", sep='\t', index=False)
print("完成")
