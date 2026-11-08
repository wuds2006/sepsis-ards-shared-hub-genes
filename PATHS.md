# PATHS.md — path conventions

Every script resolves paths from the repository root; nothing is hard-coded to a
user directory. Set the working directory to the repository root and run:

```bash
cd /path/to/repo
Rscript R/figures/06_fig6.R
bash run_all.sh --stage 1
```

## 1. Variables

| Variable | Default | What it must point at |
|---|---|---|
| working directory | repository root | required by every script (`R/00_common.R`, `R/00_setup.R`, `run_all.sh`) |
| `FIG_ROOT` | `../Figures` | where rendered figures go; `main/` and `supplementary/` are created on demand |
| `DOCK_ROOT` | `09_docking_drug` | docking pose, ligand and receptor files for the Figure S25 contact-map scripts; not distributed with the repository |
| `RAW_ROOT` | `00_rawdata` | downloaded public inputs (GEO series matrices, GWAS summary statistics, annotation); see `FETCH_PUBLIC_DATA.md` |

## 2. Where things live

| Path | Tracked? | Content |
|---|---|---|
| `data/frozen/` | yes | 13 modules of derived result tables — the only input the figure scripts need |
| `data/derived/hub6_expression_long.csv` | yes | per-sample expression of the six hub genes, extracted once from the raw matrices |
| `00_rawdata/`, `01_DEG/` … `13_survival_prognosis/` | no (`.gitignore`) | downloaded inputs and re-run outputs; recreated by `run_all.sh` |
| `R/figures/_figS25_work/` | no | intermediate 2D/3D coordinates for Figure S25, rebuilt by scripts 26–29 |

## 3. Porting notes

- The analysis scripts write into module-named directories at the repository
  root (`01_DEG/`, `05_ML_SHAP/`, …). These are ignored by git, so re-running a
  stage never dirties the repository.
- Paths inside `R/47_prognosis_extended_models.R` and `R/48_prognosis_figure.R`
  reach the raw survival cohort through `00_rawdata/validation/GSE65682/`
  (variable `SM` / `ANN`) and copy their figures to `FIG_ROOT/supplementary`.
- Two older scripts keep their original step numbers in output filenames
  (`step106_data.rds`, `step109_*.csv`). Those names are referenced by the
  frozen tables and by the figure scripts, so they are kept verbatim.
