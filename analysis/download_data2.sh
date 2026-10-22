#!/bin/bash
set -u
RAW="/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis/00_rawdata"
mkdir -p "$RAW/eqtlgen" "$RAW/validation" "$RAW/gwas"

# 并行下载(各自断点续传)
curl -sL -C - -m 7200 -o "$RAW/eqtlgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz" \
  "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz" &
P1=$!

curl -sL -C - -m 1800 -o "$RAW/gwas/finngen_R11_J10_ARDS.gz" \
  "https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_J10_ARDS.gz" &
P2=$!

curl -sL -C - -m 1800 -o "$RAW/validation/GSE243066_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243066/suppl/GSE243066_RAW.tar" &
P3=$!

wait $P1 $P2 $P3
echo "=== 下载完成状态 ==="
for f in "$RAW/eqtlgen/"*.gz "$RAW/gwas/"*.gz "$RAW/validation/GSE243066_RAW.tar"; do
  echo "$(ls -la "$f" 2>/dev/null | awk '{print $5}')  $f"
done
echo "ALL_DOWNLOADS_DONE"
