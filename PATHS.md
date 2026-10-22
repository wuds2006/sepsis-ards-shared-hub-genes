# PATHS.md — absolute-path inventory and porting guide

The archived scripts were run on the original workstation and reference absolute
paths. This file inventories every path root they use and explains how to point
the pipeline at a different machine. Facts below were extracted mechanically
from all 213 scripts in `analysis/` (counts are script counts, 2026-10-22).

## 1. Path roots referenced by the scripts

| Root (as baked into the scripts) | Scripts | What it must contain |
|---|---|---|
| `/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis` | 129 define it as `BASE` (R) or `BASE` (Python); 13 `setwd()` into it or its `04_MR_SMR/` | the project root: `00_rawdata/`, `01_DEG/` ... `13_survival_prognosis/`, `scripts/` |
| `.../ARDS_Sepsis_reanalysis/AJRCCM_submission` | 28 | submission tree: `03_figures*/`, `05_supplementary*/`, used by figure-assembly scripts (not needed to reproduce numbers) |
| `.../AJRCCM_submission/09_review_deAI_2026-09-25` | 18 | the four master manuscript files; used by manuscript-editing scripts only |
| `.../ARDS_Sepsis_reanalysis/09_docking_drug` (and `/redocking`) | docking / MM-GBSA scripts define it as `BASE` | docking outputs, structures under `cocrystal/`, `ligands/` |
| `/Users/mac/Desktop/2026年撰写的文章/2026论文撰写一Acute Respiratory Distress Syndrome and Sepsis` | 5 | an external reference copy of the paper project (superseded by the archive itself; see section 3) |
| `source("/…/ARDS_Sepsis_reanalysis/scripts/{ajrccm_theme.R, common.R}")` | 25 | the only cross-script dependencies — both files are archived in `analysis/` |

Only one username appears anywhere: `mac`. No other machines, credentials, or
secrets are embedded.

## 2. Porting to your machine

Option A — recreate the expected layout (zero script edits):

```bash
mkdir -p /path/to/parent/2026-08-23-21-32-45
git clone <repo-url> /path/to/parent/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis
# then symlink /Users/mac -> your clone location, OR use Option B
```

Option B — one-time rewrite (recommended). Replace the roots from section 1
with your local equivalents, in that order (longest first):

```bash
ROOT=/your/local/ARDS_Sepsis_reanalysis
find analysis -type f \( -name '*.R' -o -name '*.py' -o -name '*.sh' \) -exec sed -i.bak \
  -e "s|/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis|$ROOT|g" \
  -e "s|/Users/mac/Desktop/2026年撰写的文章/2026论文撰写一Acute Respiratory Distress Syndrome and Sepsis|$ROOT|g" \
  {} +
```

Keep `analysis/ajrccm_theme.R` and `analysis/common.R` at `$ROOT/scripts/`
(the 25 `source()` calls resolve there). Then set `BASE`/`PROJ` in section 0 of
`analysis/run_all.sh` and run stage S1 to populate `$ROOT/00_rawdata/`.

Raw inputs land under `00_rawdata/` with this layout (create empty dirs for
manual downloads — see README "Data sources"):

```
00_rawdata/
  gwas/          finngen_R11_J10_ARDS.gz, sepsis_ieu-b-4980_GCST90281174.tsv.gz,
                 finngen_R11_AB1_OTHER_BACTERIAL.gz
  eqtlgen/       2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-...txt.gz (+ BESD for SMR)
  validation/    GSE100159/, GSE243066/ (+ GSE243066_RAW.tar), GSE65682/
  hallmark/      h.all.v2024.1.Hs.symbols.gmt
  msigdb/        c5.all.v7.5.1.symbols.gmt, immune_GO_BP.txt
  ld_ref/        1000 Genomes EUR PLINK files (+ bldfiles/ built by step117)
  ldsc_ref/      LDSC European LD scores / weights
  GSE167363/     single-cell matrices (also see 08_singlecell/ inputs)
```

## 3. Cosmetic vs. blocking references

- **Blocking** (must resolve for the pipeline to run): the project root and its
  numbered output dirs, `00_rawdata/` inputs, `scripts/ajrccm_theme.R`,
  `scripts/common.R`, the docking `BASE` variants.
- **Non-blocking** (figure assembly and manuscript editing only; not executed
  by `run_all.sh`): `AJRCCM_submission/`, `09_review_deAI_2026-09-25/`, and the
  Desktop paper-project path used by 5 early scripts — the archive's own
  `manuscript/` copy supersedes them.
- The Desktop path above contains Chinese characters; if your shell or sed
  build mishandles them, delete the five referencing references with the
  Option B command anyway (they are single-line constants).

## 4. Verification after porting

```bash
grep -rn "/Users/mac" analysis | grep -v "\.bak" | wc -l   # expect 0
bash analysis/run_all.sh --stage S8_microarray_preprocessing   # smoke test
```

S8 (`step1a_DEG.R`) is a good smoke test: it needs only the two discovery
series matrices and platform annotations, runs in minutes, and regenerates
`01_DEG/GSE76293_DEG.csv` / `GSE69528_DEG.csv`, which can be diffed against
`results/01_DEG/` in this archive.
