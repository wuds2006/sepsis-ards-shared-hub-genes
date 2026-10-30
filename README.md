# sepsis-ards-shared-hub-genes

Analysis code and derived results for:

> **Shared Hub Genes and Immune Signatures of Sepsis and ARDS: A Multi-Omics Integration**

## Quick start

```bash
# 1. clone, then point the pipeline at your local layout (run_all.sh section 0 / PATHS.md)
# 2. fetch the raw public data
bash analysis/download_data.sh && bash analysis/download_data2.sh
#    (plus the manual downloads listed under "Data sources" below)
# 3. run the full pipeline (S1-S14), or pick stages with --stage / --list
bash analysis/run_all.sh --list
bash analysis/run_all.sh --stage S1_data_acquisition,S8_microarray_preprocessing,S9_deg_wgcna
bash analysis/run_all.sh
```

`run_all.sh` is the master entry point: it reproduces every RESULTS number of
the manuscript, stage by stage, in the order of the Supplementary Methods
sections S1-S14 (see "Pipeline map" below). Derived outputs land in the
numbered project directories `01_DEG` ... `13_survival_prognosis`; this
repository's `results/` folder holds the archived copies of those outputs.

## Where to look first

| File | Purpose |
|---|---|
| `RESULTS_TO_CODE.md` | **start here** — every RESULTS number of the manuscript mapped to the scripts and `results/` files that produce it, keyed to Supplementary Methods sections S1-S14 |
| `DATA_INVENTORY.md` / `DATA_INVENTORY.csv` | **data availability** — every public input (GEO accession, GWAS ID, eQTL, protein/structure/pharmacology resource) with platform, version, sample counts, exclusions, access date, and role |
| `FETCH_PUBLIC_DATA.md` | retrieval runbook: the exact commands to rebuild the raw input tree from GEO, IEU OpenGWAS, FinnGen, and eQTLGen |
| `run_all.sh` | master pipeline (stage registry; one stage per Supplementary Methods section) |
| `PATHS.md` | every absolute path the scripts reference, and how to point them at your machine |
| `METHODS_TO_CODE.md` | sentence-level map from every Methods step to scripts and outputs |
| `analysis/` | all scripts (R / Python / bash), described in the pipeline map below |
| `results/` | archived derived outputs; every manuscript number traces here |
| `manuscript/` | submitted manuscript text, supplementary methods, figure legends, and all 27 supplementary tables as submitted (synced 2026-10-30) |
| `environment/` | interpreter and tool versions (`tool_versions.txt`) plus R and Python package inventories |
| `MANIFEST.txt` | byte size and SHA-256 checksum for every file |

## Pipeline map (stage = Supplementary Methods section)

`run_all.sh` executes the stages below in this order. Within a stage, scripts
run in the listed order. Figure-assembly and manuscript-editing scripts are
archived for provenance but deliberately not run (see `METHODS_TO_CODE.md`).

| Stage | Supplementary Methods section | Key scripts (relative to `analysis/`) |
|---|---|---|
| S1 | Data acquisition and sample curation | `download_data.sh`, `download_data2.sh`, `parse_geo_matrix.py` |
| S8 | Microarray preprocessing and probe-to-gene mapping | `step1a_DEG.R` |
| S9 | Differential expression and WGCNA | `step1b_WGCNA.R`, `step1c_intersection.R`, `step55_wgcna_preservation.R`, `step27_reviewer_fix_R.R` |
| S2 | Machine-learning and SHAP | `step2a_ML.R`, `step2a2_ML_sepsis.R`, `step50_ml_stability_overfit.R`, `step2e_nomogram.R` (disease-prediction models), `step27_reviewer_fix_py.py` |
| S6 | External validation | `step2b_GSE100159_validate.R`, `step2b_validation_ARDS.R`, `step2b2_*`, `step2b3_*`, `step2b4_*`, `step2b_AUC_CI.R`, `step2d_expression_ROC.R`, `step49_cd177_ablation.R`, `step75_muc1_sensitivity.R`, `step52_validation_diagnostics.R`, `step59_cohort_pca.R` |
| S3 | Immune-cell deconvolution | `step41_cibersort_ssgsea.R`, `step45_hub6_gsva.R`, `step2c_hub_immune.R`, `step47_ssgsea_immune_4cohorts.R` |
| S10 | GSVA of immune pathways | `step2c2_GSVA_immune.R`, `step2c3_GSVA_immune_top20.R` |
| S4 | Single-cell processing | `step3a_scanpy.py`, `step3a_singlecell.R`, `step3a_annotate.R`, `step3e_GSE167363_annotate.R` (+ fix scripts), `step3b_*` (CellChat), `step3c_*` / `step3h_*` (slingshot), `step3o_monocle3.R`, `step3g_aucell.R`, `step3n_metabolism.py`, `step3n_seurat_integration.R`, `step54_donor_level.R`, `step53_celloracle_per_donor.py`, `step110_cellchat_netP_recalc.R`, `step111_cellchat_GSE167363_netP_recalc.R` (quantitative CellChat) |
| S5 | Virtual knockout | `step3d_scTenifoldKnk.R`, `step3d_EXOSC4_ko.R`, `step3i_GSE167363_scTenifoldKnk.R`, `step3i_extra_knockout.R`, `step3d_CellOracle_local.py`, `step3j_GSE167363_CellOracle.py`, `step3n_celloracle_stratified.py`, `step3n_ko_stratified.R`, `step3n_ko_stratified_patch.R`, `step3n_ko_stratified_compare.R`, `step22_KO_top10_enrichment.R` |
| S7 | Mendelian randomization | `extract_MR_4980.py`, `match_ards_4980.py`, `step1d_MR_4980.R`, `step48a_extract_reverse_MR.py`, `step48b_reverse_MR.R`, `step51_mr_power.R`, `prepare_ldsc.py` (LDSC), `convert_gwas_ma.py`, `align_coloc.py`, `step1d_coloc.R`, `step1d_heidi.R`. SMR+HEIDI proper: `step117_run_smr_heidi.sh` (runbook for the external SMR binary; archived outputs in `results/04_MR_SMR/`) |
| S11 | Protein-level validation | `step_hpa_validation.R`, `step_hpa_immune_cell.R` |
| S12 | Drug repositioning and molecular docking | `step116_cmap_enrichr_libraries.py` (primary Enrichr query), `step67_cmap_negative_control.py`, `step76_cmap_deg_null.py` (null backgrounds), `step68_opentargets_druggability.py`, `step115_chembl_known_drugs.py`, `step4_docking_site.py` (+ `step4_dock_one.py`, `step57_upp1_redock_fixed_box.py`), `step56_redocking.py`, `step66_alphafold_plddt.py`, `step70_mmgb_sasa.py`, `step69_plip_vnn1_luteolin.py`, `step5_network_ppi.R` |
| S13 | Single-gene GSEA | `step4a_single_gene_gsea.R` (two cohorts), `step46_single_gene_gsea_4cohorts.R` (all four cohorts) |
| S14 | Prognostic modeling in GSE65682 | `step105_survival_prognosis_GSE65682.R`, `step105b_pergene_KM_GSE65682.R`, `step106_extended_modeling_GSE65682.R`, `step109_three_statistical_facts.R` |
| — | Reproducibility self-check | `verify_manuscript.py` |

Superseded in-place revisions (kept for provenance, not executed): the GSE167363
annotation was finalized by `step3e_fix2_annotate.R` / `step3e_fix3_annotate.R`
superseding `step3e_GSE167363_annotate.R`; PBMC CellChat visualization was
finalized by `step3k_pbmc_figs.R` superseding the initial `step3b_cellchat_compare.R`
renders; UPP1 docking was corrected by `step57_upp1_redock_fixed_box.py`
superseding the original UPP1 box in `step4_docking_site.py`
(archived as `docking_site_results_OLD_UPP1box.csv`).

## Data sources (public, not redistributed here)

Full inventory with accession numbers, platforms, versions, sample counts, exclusions,
and access dates: [`DATA_INVENTORY.md`](DATA_INVENTORY.md) (and `DATA_INVENTORY.csv`).
Retrieval commands: [`FETCH_PUBLIC_DATA.md`](FETCH_PUBLIC_DATA.md).
Summary:

| Dataset | Source | Version / access date |
|---|---|---|
| GSE76293, GSE69528, GSE100159, GSE243066, GSE151263, GSE167363 | Gene Expression Omnibus | series matrices as deposited |
| GSE65682 (prognosis) | GEO | series matrix + family SOFT (GPL13667) |
| ieu-b-4980 (sepsis GWAS) | IEU OpenGWAS / GWAS Catalog GCST90281174 | as downloaded by `extract_MR_4980.py` |
| ARDS outcome GWAS | FinnGen | R11 (J10_ARDS) |
| eQTLGen whole-blood cis-eQTL | eqtlgen.org | 2019-12-11 probe-level release (BESD for SMR: 31,684 individuals) |
| Blood protein levels (CD177, VNN2; pan-disease Olink PEA) | Human Protein Atlas, Blood resource (Disease Blood Atlas) | **HPA version 24, files downloaded 2025-11-06** (v25, 71 diseases, released 2025-11-11, postdates access and was not used) |
| Target tractability | Open Targets Platform | accessed September 2026 (`step68`) |
| Known-drugs cross-check | ChEMBL web services | **release 37, accessed 2026-09-17** (`step115_chembl_known_drugs.py` reproduces the archived CSV byte-identically) |
| Drug-repositioning query | Enrichr API (Old CMap up/down, DSigDB, Drug Perturbations from GEO) | archived 2026-09; `step116_cmap_enrichr_libraries.py` re-issues the identical query |
| MSigDB | Broad GSEA website | c5.all.v7.5.1 (GSVA), h.all.v2024.1.Hs.symbols (GSEA, metabolism) |
| Protein structures | PDB / AlphaFold | UPP1: 3EUF (co-crystallized BAU); VNN1: 7SLV (co-crystallized 9TC); AlphaFold for CDK5RAP2 / CD177 / GRB10 / EXOSC4 |

`download_data.sh` / `download_data2.sh` fetch the eQTLGen archive, the FinnGen
R11 ARDS summary statistics, and the GSE243066 RAW tar. Everything else (GEO
series matrices, HPA TSVs, MSigDB GMTs, protein structures) must be fetched
manually from the sources above; `PATHS.md` lists the exact expected file
locations under `00_rawdata/`.

## Environment

R 4.6.1 and Python 3.13.12. Full inventories, extracted mechanically from the scripts:

- `environment/R_packages_used.txt` — 71 R packages (key ones: limma, WGCNA,
  clusterProfiler, TwoSampleMR, survival, rms, glmnet, ranger, GSVA, CellChat,
  monocle3 1.3.1, slingshot, Seurat, scTenifoldKnk, DESeq2, coloc, fgsea)
- `environment/python_packages_used.txt` — third-party Python modules (scanpy,
  anndata, celloracle, sklearn, OpenMM 8.5.1, rdkit, meeko, …)
- `environment/tool_versions.txt` — interpreter and tool pins (SMR 1.03,
  AutoDock Vina 1.2.5, OpenMM 8.5.1, PLINK 1.9, monocle3 1.3.1) and resource
  versions with access dates

Only the versions stated above were explicitly recorded; the remaining packages were
installed at their then-current releases and were not captured in a package-lock
snapshot. External binaries required by specific stages: AutoDock Vina 1.2.5 (S12);
SMR v1.03 + PLINK + LDSC (S7; see `analysis/step117_run_smr_heidi.sh`).

## Reproducibility

Scripts are deterministic apart from procedures with explicitly set seeds
(random survival forest, LASSO cross-validation, bootstrap optimism correction,
negative-control sampling); seeds are fixed in the scripts.
`verify_manuscript.py` re-checks the manuscript's word counts, citation
ordering, and numeric claim anchors against this repository.

## Data availability and materials

This study is a **secondary analysis of public data**: no new participants were
enrolled, no identifiable personal data were accessed, and no raw public data are
redistributed here (their licences and size prohibit it). The repository therefore
ships the complete *materials* set:

| Layer | Contents |
|---|---|
| Code | 217 scripts (R / Python / bash) + `run_all.sh` (17 stages → Supplementary Methods S1-S14) |
| Derived results | 369 files under `results/`, one folder per analysis module |
| Manuscript materials | manuscript text, Supplementary Methods S1-S14, figure legends, and all 27 supplementary tables as submitted |
| Input provenance | `DATA_INVENTORY.md` / `.csv` (19 input entries) and `FETCH_PUBLIC_DATA.md` |
| Environment | `environment/tool_versions.txt`, `R_packages_used.txt`, `python_packages_used.txt` |
| Integrity | `MANIFEST.txt` (size + SHA-256 per file), `code_archive.zip` |

Deliberately excluded and re-downloadable: raw public data (`00_rawdata/`),
large single-cell objects (`*.rds`, `*.h5ad`, `*.h5seurat`), and figure TIFFs
(shipped with the manuscript). See `.gitignore`.

## Permanent identifiers

- Code repository: **https://github.com/wuds2006/sepsis-ards-shared-hub-genes**
- Archived release (all versions): **https://doi.org/10.5281/zenodo.23668532** (Zenodo concept DOI)

## Cite this

Wu D, Yang H, Shen Y, Chen W, Wu J, Wu Z, Wu Q. *sepsis-ards-shared-hub-genes*
(version 1.1.0). Zenodo. https://doi.org/10.5281/zenodo.23668532
(`CITATION.cff` carries the same metadata in machine-readable form.)

## Licence

MIT. See `LICENSE`. Third-party data are governed by their own terms of use;
see `DATA_INVENTORY.md`.
