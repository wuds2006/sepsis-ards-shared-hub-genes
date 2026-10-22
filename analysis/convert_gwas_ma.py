#!/usr/bin/env python3
"""
GWAS → SMR MA 格式转换 (SNP A1 A2 freq b se p n)
- sepsis: Hamilton ieu-b-4980 / GCST90281174 (GRCh37)
- ARDS:   FinnGen R11 J10_ARDS (GRCh38)
向量化分块处理, 避免 iterrows 慢
"""
import pandas as pd
import numpy as np

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
GWAS = f"{BASE}/00_rawdata/gwas"
OUT = f"{BASE}/04_MR_SMR"

def convert_sepsis():
    src = f"{GWAS}/sepsis_ieu-b-4980_GCST90281174.tsv.gz"
    dst = f"{OUT}/sepsis_4980.ma"
    usecols = ['rs_id','effect_allele','other_allele','effect_allele_frequency',
               'beta','standard_error','p_value','n']
    header_written = False
    n_rows = 0
    for c in pd.read_csv(src, sep='\t', compression='gzip', usecols=usecols, chunksize=3000000):
        c = c.dropna(subset=['rs_id','beta','standard_error','effect_allele','other_allele'])
        c = c[c['rs_id'].astype(str).str.startswith('rs')]
        c = c.rename(columns={'rs_id':'SNP','effect_allele':'A1','other_allele':'A2',
                              'effect_allele_frequency':'freq','beta':'b',
                              'standard_error':'se','p_value':'p','n':'n'})
        c = c[['SNP','A1','A2','freq','b','se','p','n']]
        c['n'] = c['n'].fillna(462927).astype(int)
        c.to_csv(dst, sep='\t', index=False, header=not header_written, mode='a')
        header_written = True
        n_rows += len(c)
    print(f"sepsis_4980.ma 完成: {n_rows} 行")

def convert_ards():
    src = f"{GWAS}/finngen_R11_J10_ARDS.gz"
    dst = f"{OUT}/ards_J10.ma"
    usecols = ['rsids','alt','ref','af_alt','beta','sebeta','pval']
    N = 493732
    header_written = False
    n_rows = 0
    for c in pd.read_csv(src, sep='\t', compression='gzip', usecols=usecols, chunksize=3000000):
        c = c.dropna(subset=['rsids','beta','sebeta'])
        c['rsids'] = c['rsids'].astype(str).str.split(',').str[0]
        c = c[c['rsids'].str.startswith('rs')]
        c['n'] = N
        c = c.rename(columns={'rsids':'SNP','alt':'A1','ref':'A2','af_alt':'freq',
                              'beta':'b','sebeta':'se','pval':'p'})
        c = c[['SNP','A1','A2','freq','b','se','p','n']]
        c.to_csv(dst, sep='\t', index=False, header=not header_written, mode='a')
        header_written = True
        n_rows += len(c)
    print(f"ards_J10.ma 完成: {n_rows} 行")

convert_sepsis()
convert_ards()
print("===== GWAS 转换完成 =====")
