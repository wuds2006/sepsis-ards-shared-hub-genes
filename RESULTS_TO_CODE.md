# RESULTS_TO_CODE.md — which module produces which result

Every number in the manuscript traces back to one of thirteen analysis modules.
Each module is a directory in `data/frozen/`; the scripts that fill it live in
`R/`, `python/` and `shell/`, and the figure scripts in `R/figures/` read only
from `data/frozen/`. Reported values in the manuscript are the acceptance
targets when a module is re-run.

Stage numbering, the Supplementary Methods section each stage implements, and
the script-by-script breakdown are in [`docs/METHODS_TO_CODE.md`](docs/METHODS_TO_CODE.md);
this file is the module-level index.

| Module (`data/frozen/`) | Scripts | Representative outputs | Where it appears |
|---|---|---|---|
| `01_DEG` | `R/10_deg.R` | `GSE76293_DEG.csv`, `GSE69528_DEG.csv`, `*_DEG_all.csv`, group assignments | Figure 1A–B, Table S6 |
| `02_WGCNA` | `R/11_wgcna.R`, `R/12_wgcna_preservation.R` | `*_gene_module.csv`, `*_disease_module_genes.txt`, `modulePreservation_*`, `hub_gene_kME.csv` | Figure 1C–D, Table S10, Table S12 |
| `03_intersection` | `R/13_intersection.R`, `R/14_enrichment_go_kegg.R` | `Shared_CGs.csv`, `Shared_CGs_GO_BP.csv`, `Shared_CGs_KEGG.csv` | Figure 1E, Figure 2, Tables S7–S8 |
| `04_MR_SMR` | `R/50_mr_forward.R`, `51_mr_reverse.R`, `52_mr_power.R`, `53_colocalization.R`, `54_smr_heidi.R`; `python/02_mr_extract_forward.py`, `03_convert_gwas_sumstats.py`, `04_match_ards_against_sepsis.py`, `05_align_coloc_inputs.py`, `06_prepare_ldsc_inputs.py`, `07_mr_extract_reverse.py`; `shell/03_prepare_mr_inputs.sh`, `04_run_smr_heidi.sh` | `MR_*_results.csv`, `MR_power_by_OR.csv`, `coloc_*`, `SMR_*`, `HEIDI_*` | Figure 9, Tables S13–S14 |
| `05_ML_SHAP` | `R/20_ml_ards.R`, `21_ml_sepsis.R`, `22_ml_stability.R` | `ML_AUC_results.csv`, `SHAP_importance.csv`, `ML_key_genes.txt`, stability and permutation tables | Figure 2A–B, Figure 3, Table S9 |
| `06_external_validation` | `R/23_validation_ards.R`–`R/30_cd177_ablation.R` | `GSE243066_*`, `GSE100159_*`, `AUC_CI.csv`, `metaDE_single_gene_validation.csv` | Figure 6, Table 2, Tables S11, S15 |
| `07_hub_immune` | `R/35_cibersort_ssgsea.R`, `36_gsva_immune.R`, `37_hub6_gsva.R`, `39_ssgsea_4cohorts.R`, `40_cohort_pca.R` | `*_immune_diff.csv`, `*_hub_immune_cor.csv`, `*_hub6_GSVA.csv`, `*_hub6_ssGSEA.csv` | Figure 7, Figures S7, S9–S10, Table S16 |
| `08_singlecell` | `R/60_annotate_gse167363.R`–`R/70_donor_level.R`; `python/10_scanpy_cluster.py`–`17_cross_dataset.py` | `GSE151263_cell_annot.csv`, `GSE167363/`, `cellchat/`, `monocle3/`, `aucell/`, `knockout/`, `etosis_aucell/`, `hdwgcna/` | Figure 8, Figures S11–S18, S26–S27, Tables S17–S18, S25, S27–S28 |
| `09_docking_drug` | `python/20_docking_screen.py`–`26_plip_contacts.py`, `30_cmap_random_null.py`–`34_opentargets_druggability.py`; `R/figures/26_figS25_2d_contact_map.py`–`29_figS25_3d_scene.py` | `*_docking_*.csv`, `MMGBSA_*`, `*_hub_prediction.csv`, `CMap_*` | Figure 10, Figure S25, Tables S19–S21 |
| `10_HPA_validation` | `R/75_hpa_plasma_protein.R`, `76_hpa_immune_cells.R` | `blood_ms_concentration.tsv`, `HPA_immune_cell_specificity.csv`, `HPA_CD177_VNN2_infection_vs_healthy.csv` | Figure 10A |
| `11_network_pharma` | `R/15_ppi_network.R` | `PPI_nodes.csv`, `PPI_edges.csv`, `GO_BP_enrichment.csv`, `KEGG_2021_Human_enrichment.csv`, `Reactome_2022_enrichment.csv` | Figure 2, Figure 3D, Figure S23 |
| `12_single_gene_gsea` | `R/38_single_gene_gsea.R` | `single_gene_gsea_all_results.csv`, `single_gene_gsea_top2_*.csv` | Figure 4E, Figure S6 |
| `13_survival_prognosis` | `R/45_prognosis_cox.R`–`R/49_prognosis_statistics.R` | `tables/` (Cox, calibration, DCA, time-ROC, NRI/IDI, cut-point sensitivity), `models/` | Figure S20, Figures S21, Table S22–S24 |

## Conventions worth knowing

- **Frozen outputs are flat within each module.** The live pipeline writes into
  subdirectories (`08_singlecell/knockout/` etc.); the packaging step flattened
  some of them. Both views contain the same files.
- **Figure scripts never read raw data.** `R/figures/*` reads only
  `data/frozen/` (plus `data/derived/hub6_expression_long.csv`), which is why
  every figure can be redrawn in minutes without re-running the upstream work.
- **Model objects.** `13_survival_prognosis/models/` carries the fitted objects
  (`step106_data.rds`, `cox6_model.rds`, `step106_models.rds`) so the
  downstream figures and tables do not need to re-fit.
