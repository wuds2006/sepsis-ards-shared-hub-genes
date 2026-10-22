#!/bin/bash
set -u
RAW="/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/00_rawdata"
mkdir -p "$RAW/eqtlgen" "$RAW/validation" "$RAW/gwas"

echo "[1/3] 下载 eQTLGen cis-eQTL（约2.5GB，支持续传）..."
curl -sL -C - -m 7200 \
  -o "$RAW/eqtlgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz" \
  "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
echo "eQTLGen 大小: $(ls -la "$RAW/eqtlgen/" | grep txt.gz | awk '{print $5}')"

echo "[2/3] 下载 GSE243066 RAW（ARDS 外验 RNA-seq）..."
curl -sL -C - -m 3600 \
  -o "$RAW/validation/GSE243066_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243066/suppl/GSE243066_RAW.tar"
echo "GSE243066 大小: $(ls -la "$RAW/validation/GSE243066_RAW.tar" 2>/dev/null | awk '{print $5}')"

echo "ALL_DOWNLOADS_DONE"
