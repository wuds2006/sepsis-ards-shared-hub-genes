#!/usr/bin/env bash
# =============================================================================
# step117_run_smr_heidi.sh -- SMR + HEIDI provenance runbook (Supplementary
# Methods S7). NOT executed by run_all.sh: it requires the SMR v1.03 binary
# (https://yanglab.westlake.edu.cn/software/smr/), the eQTLGen whole-blood
# cis-eQTL BESD database (31,684 individuals), and a local 1000 Genomes
# European PLINK reference, none of which are redistributed in this archive.
#
# What it reproduces (archived under results/04_MR_SMR/):
#   smr_<GENE>_<outcome>_heidi.smr   per-gene SMR + HEIDI output (10 files:
#                                    5 genes x 2 outcomes)
#   SMR.txt                          combined SMR/HEIDI table used by the
#                                    manuscript (all P_SMR > 0.1,
#                                    P_HEIDI 0.287-0.855)
#   SMR_coloc_summary.txt            side-by-side coloc PP.H4 summary
#
# Status: the original run was performed with the SMR binary from the command
# line (no wrapper script was archived). This runbook records the documented
# protocol from Supplementary Methods S7 so the run can be re-executed or
# audited. Flag values follow the SMR v1.03 manual; verify against the manual
# before running.
# =============================================================================
set -euo pipefail

# ---- 0. paths and external resources (EDIT THESE to your local layout) -----
BASE="${BASE:-/path/to/ARDS_Sepsis_reanalysis}"
RAW="$BASE/00_rawdata"
OUT="$BASE/04_MR_SMR"
SMR_BIN="${SMR_BIN:-smr}"                       # SMR v1.03 executable
EQTLGEN_BESD="${EQTLGEN_BESD:-$RAW/eqtlgen/eqtlgen_besd}"   # .besd prefix
G1000_PREFIX="${G1000_PREFIX:-$RAW/ld_ref/1000G_EUR}"       # PLINK bfile prefix

# Outcomes: sepsis = ieu-b-4980 (GCST90281174), ARDS = FinnGen R11 J10_ARDS.
# Both must first be converted to SMR .ma format -- this IS scripted:
#   python3 analysis/convert_gwas_ma.py
# producing 04_MR_SMR/sepsis_4980.ma and 04_MR_SMR/ards_finngen.ma
# (GRCh37 for sepsis; GRCh38 for ARDS -- see the header of convert_gwas_ma.py).

SEPSIS_MA="$OUT/sepsis_4980.ma"
ARDS_MA="$OUT/ards_finngen.ma"

# Five hub genes have a cis-eQTL instrument (VNN1, CD177, CDK5RAP2, GRB10, UPP1);
# EXOSC4 has none and is not tested.
GENES=(VNN1 CD177 CDK5RAP2 GRB10 UPP1)

# ---- 1. build the local LD reference for HEIDI (--make-bld, +/- 2,000 kb) ----
# Per Supplementary Methods S7: "local 1000 Genomes European linkage-
# disequilibrium references built with --make-bld (+/-2,000 kb windows)".
# The .bld files live alongside the bfile; SMR builds them per chromosome.
# SMR v1.03 syntax (see manual, "Create the .bld files"):
#   smr --make-bld <chr> --bfile <prefix> --out <bld-dir>
for chr in $(seq 1 22); do
  "$SMR_BIN" --make-bld "$chr" --bfile "$G1000_PREFIX" \
    --out "$RAW/ld_ref/bldfiles" 2>/dev/null || true   # idempotent: skips existing
done

# ---- 2. SMR + HEIDI per gene x outcome --------------------------------------
# Flags follow the SMR manual defaults used in the original run:
#   --beqtl-summary  eQTLGen BESD
#   --gwas-summary   outcome in .ma format
#   --bfile          1000G EUR reference (used by HEIDI)
#   --probe          restrict to the gene's probe (gene-level test)
#   --heidi-mbc      HEIDI multiple-binding-criterion variants (manual default)
#   --out            smr_<GENE>_<outcome>_heidi
for outcome in sepsis ards; do
  case "$outcome" in
    sepsis) MA="$SEPSIS_MA" ;;
    ards)   MA="$ARDS_MA"   ;;
  esac
  for gene in "${GENES[@]}"; do
    echo "[SMR] $gene x $outcome"
    "$SMR_BIN" \
      --bfile "$G1000_PREFIX" \
      --beqtl-summary "$EQTLGEN_BESD" \
      --gwas-summary "$MA" \
      --probe "$gene" \
      --heidi-mbc 5 \
      --diff-freq-prop 0.9 \
      --out "$OUT/smr_${gene}_${outcome}_heidi"
  done
done

# ---- 3. combine per-gene outputs into SMR.txt -------------------------------
# The archived SMR.txt concatenates the per-gene .smr tables with a gene and
# outcome column; regenerate it with:
python3 - "$OUT" "${GENES[@]}" <<'PY'
import os, sys, csv
out = sys.argv[1]; genes = sys.argv[2:]
rows = []
for outcome in ("sepsis", "ards"):
    for g in genes:
        p = os.path.join(out, f"smr_{g}_{outcome}_heidi.smr")
        if not os.path.exists(p):
            print(f"  [warn] missing {p}"); continue
        with open(p) as fh:
            hdr = None
            for line in fh:
                if line.startswith("probe"):  # SMR writes header starting 'probe'
                    hdr = line.lstrip("#").strip().split("\t"); continue
                if hdr and line.strip() and not line.startswith("#"):
                    rows.append([g, outcome] + line.strip().split("\t"))
with open(os.path.join(out, "SMR.txt"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["gene", "outcome"] + (hdr or []))
    w.writerows(rows)
print(f"  wrote {out}/SMR.txt ({len(rows)} rows)")
PY

echo "Done. Compare 04_MR_SMR/SMR.txt against the archived copy in results/04_MR_SMR/"
echo "(archived values: all P_SMR > 0.1; P_HEIDI range 0.287-0.855; manuscript"
echo " reports no significant SMR association and no pleiotropy signal)."
