#!/bin/bash
# MR 数据提取: 流式提取暴露显著SNP + 结局匹配SNP, 写小文件
set -u
RAW="/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/00_rawdata/gwas"
OUT="/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/04_MR_SMR"
mkdir -p "$OUT"
EXPO="$RAW/finngen_R11_AB1_OTHER_BACTERIAL.gz"
OUTC="$RAW/finngen_R11_J10_ARDS.gz"

echo "[1/2] 提取暴露显著SNP (p<5e-8)..."
gunzip -c "$EXPO" | awk -F'\t' 'NR==1 || $7+0 < 5e-8' > "$OUT/expo_sig.tsv"
echo "暴露显著SNP行数: $(($(wc -l < "$OUT/expo_sig.tsv") - 1))"

echo "[2/2] 提取结局匹配SNP..."
# 提取rsids列(第5列), 构建grep模式
RS=$(tail -n +2 "$OUT/expo_sig.tsv" | cut -f5 | grep -v '^$' | tr '\n' '|' | sed 's/|$//')
gunzip -c "$OUTC" | head -1 > "$OUT/outc_match.tsv"
gunzip -c "$OUTC" | grep -E "$RS" >> "$OUT/outc_match.tsv"
echo "结局匹配行数: $(($(wc -l < "$OUT/outc_match.tsv") - 1))"
echo "DONE"
