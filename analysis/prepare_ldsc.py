#!/usr/bin/env python3
"""
LDSC 数据准备:
1. sepsis GWAS (hg37) -> liftover hg38 + 转 LDSC 格式
2. ARDS GWAS (hg38) -> 转 LDSC 格式
输出: LDSC munge 就绪的输入 (SNP A1 A2 BETA P N CHR BP)
"""
import pandas as pd
import numpy as np
from pyliftover import LiftOver

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
GWAS = f"{BASE}/00_rawdata/gwas"
LD = f"{BASE}/00_rawdata/ld_ref"
OUT = f"{BASE}/04_MR_SMR/ldsc_input"
import os
os.makedirs(OUT, exist_ok=True)

def prepare_ards():
    """ARDS (FinnGen, hg38) -> LDSC 格式"""
    src = f"{GWAS}/finngen_R11_J10_ARDS.gz"
    dst = f"{OUT}/ards_ldsc.txt"
    cols = ['#chrom','pos','alt','ref','rsids','beta','sebeta','pval']
    N = 493732
    w = open(dst, 'w')
    w.write("SNP\tA1\tA2\tBETA\tP\tN\tCHR\tBP\n")
    n = 0
    for c in pd.read_csv(src, sep='\t', compression='gzip', usecols=cols, chunksize=3000000):
        c = c.dropna(subset=['rsids','beta','sebeta','pval'])
        c['rsids'] = c['rsids'].astype(str).str.split(',').str[0]
        c = c[c['rsids'].str.startswith('rs')]
        c['N'] = N
        c[['rsids','alt','ref','beta','pval','N','#chrom','pos']].to_csv(
            w, sep='\t', index=False, header=False)
        n += len(c)
    w.close()
    print(f"ARDS LDSC 输入: {n} 行 -> {dst}")

def prepare_sepsis():
    """sepsis (Hamilton, hg37) -> liftover hg38 + LDSC 格式"""
    src = f"{GWAS}/sepsis_ieu-b-4980_GCST90281174.tsv.gz"
    dst = f"{OUT}/sepsis_ldsc.txt"
    cols = ['chromosome','base_pair_location','effect_allele','other_allele',
            'beta','standard_error','p_value','rs_id']
    N = 486484
    lo = LiftOver(f"{LD}/liftover/hg19ToHg38.over.chain.gz")
    w = open(dst, 'w')
    w.write("SNP\tA1\tA2\tBETA\tP\tN\tCHR\tBP\n")
    n = 0; lifted = 0; failed = 0
    cache = {}  # (chr,pos) -> (chr,pos) 缓存
    for c in pd.read_csv(src, sep='\t', compression='gzip', usecols=cols, chunksize=1000000):
        c = c.dropna(subset=['rs_id','beta','standard_error','p_value'])
        c = c[c['rs_id'].astype(str).str.startswith('rs')]
        new_chr = []; new_bp = []
        for chrom, pos in zip(c['chromosome'], c['base_pair_location']):
            key = (chrom, pos)
            if key in cache:
                nc, nb = cache[key]
            else:
                try:
                    res = lo.convert_coordinate(chrom, pos)
                except Exception:
                    res = []
                if res and len(res) > 0:
                    nc = int(res[0][0].replace('chr',''))
                    nb = res[0][1]
                else:
                    nc = -1; nb = -1
                cache[key] = (nc, nb)
            if nc > 0:
                new_chr.append(nc); new_bp.append(nb); lifted += 1
            else:
                new_chr.append(-1); new_bp.append(-1); failed += 1
        c['CHR'] = new_chr
        c['BP'] = new_bp
        c = c[c['CHR'] > 0]
        c['N'] = N
        c[['rs_id','effect_allele','other_allele','beta','p_value','N','CHR','BP']].to_csv(
            w, sep='\t', index=False, header=False)
        n += len(c)
    w.close()
    print(f"sepsis LDSC 输入: {n} 行 (liftover 成功 {lifted}, 失败 {failed}) -> {dst}")

prepare_ards()
prepare_sepsis()
print("===== LDSC 数据准备完成 =====")
