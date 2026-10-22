#!/usr/bin/env python3
"""
从 GEO 原始 series_matrix 解析探针级表达矩阵 + 正确分组。
不采信之前 Windows 预处理文件。
输出:
  01_DEG/GSE76293_expr.txt  (探针级, rows=探针, cols=样本)
  01_DEG/GSE76293_groups.csv (样本, source_name, 分组)
  01_DEG/GSE69528_expr.txt
  01_DEG/GSE69528_groups.csv
"""
import gzip, os, csv

BASE = "/Users/mac/Desktop/2026年撰写的文章/2026论文撰写一Acute Respiratory Distress Syndrome and Sepsis"
OUT  = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/01_DEG"
os.makedirs(OUT, exist_ok=True)

def parse_series_matrix(path):
    """返回 (meta_dict, header, data_rows)。meta: field->[values]; header: 样本GSM列表; data: {probe: [values]}"""
    meta = {}
    header = None
    data = []
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
        # 先读元数据 "!" 行
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('!'):
                parts = line.split('\t')
                field = parts[0]
                vals = [p.strip('"') for p in parts[1:]]
                if field.startswith('!Sample_'):
                    meta[field] = vals
                continue
            if line.startswith('"ID_REF"') or line.startswith('ID_REF'):
                header = [p.strip('"') for p in line.split('\t')]
                if header and header[0] in ('ID_REF', '"ID_REF"'):
                    header = header[1:]  # 去掉 ID_REF 列
                break
        # 读数据行
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('!') or line.strip() == '':
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            probe = parts[0].strip('"')
            vals = [p.strip('"') for p in parts[1:]]
            data.append((probe, vals))
    return meta, header, data

def write_expr_and_groups(ds, meta, header, data):
    n_sample = len(header)
    gsm = meta['!Sample_geo_accession'] if '!Sample_geo_accession' in meta else header
    source = meta.get('!Sample_source_name_ch1', ['']*n_sample)
    title  = meta.get('!Sample_title', ['']*n_sample)

    # 分组规则
    group = []
    for s in source:
        if ds == 'GSE76293':
            if s == 'PMNs from ARDS patient':
                group.append('ARDS')
            elif s == 'PMNs from HVT':
                group.append('Control')
            else:
                group.append('EXCLUDE_cultured')  # 体外培养样本排除
        else:  # GSE69528
            if s == 'Uninfected healthy':
                group.append('Control')
            elif 'sepsis caused by' in s:
                group.append('Sepsis')
            else:
                group.append('EXCLUDE_t2dm')  # 未感染T2DM排除

    # 写表达矩阵
    expr_path = os.path.join(OUT, f'{ds}_expr.txt')
    with open(expr_path, 'w') as f:
        f.write('ID_REF\t' + '\t'.join(header) + '\n')
        for probe, vals in data:
            f.write(probe + '\t' + '\t'.join(vals) + '\n')

    # 写分组
    grp_path = os.path.join(OUT, f'{ds}_groups.csv')
    with open(grp_path, 'w') as f:
        f.write('sample,source_name,title,group\n')
        for i in range(n_sample):
            f.write(f'{header[i]},{source[i]},{title[i]},{group[i]}\n')

    # 统计
    from collections import Counter
    cnt = Counter(group)
    print(f'  {ds}: 总样本 {n_sample}, 分组统计 {dict(cnt)}')

for ds, path in [
    ('GSE76293', f'{BASE}/ARDS/GSE76293/GSE76293_series_matrix.txt.gz'),
    ('GSE69528', f'{BASE}/脓毒症/GSE69528/GSE69528_series_matrix.txt.gz'),
]:
    print(f'===== {ds} =====')
    meta, header, data = parse_series_matrix(path)
    print(f'  探针数: {len(data)}, 样本列数: {len(header)}')
    write_expr_and_groups(ds, meta, header, data)

print('\nDONE')
