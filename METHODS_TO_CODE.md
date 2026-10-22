# METHODS_TO_CODE.md — Manuscript Methods-to-Code Map

Every analytical claim in the manuscript *"Shared Hub Genes and Immune Signatures of
Sepsis and ARDS: A Multi-Omics Integration"* traces to the scripts and derived files
below. Scripts live in `analysis/` (run in numeric order within each stage); derived
outputs live in `results/`. Raw public data are not redistributed (see
`README.md → Data sources`); large single-cell objects (*.rds / *.h5ad) are regenerated
by the scripts and excluded via `.gitignore`.

## 1. Cohort curation and differential expression

| Manuscript step | Scripts | Key outputs |
|---|---|---|
| GEO download and series-matrix parsing | `download_data.sh`, `download_data2.sh`, `parse_geo_matrix.py` | `00_rawdata/` (excluded from archive) |
| Off-target sample removal and cohort curation (GSE76293 in-vitro neutrophils, GSE69528 uninfected T2D controls) | `step1a_DEG.R`, `step59_cohort_pca.R` | `results/01_DEG/`, cohort PCA diagnostics |
| limma differential expression (|log2 FC| > 1, BH-adjusted P < 0.05): 557 DEGs (ARDS), 1,458 (sepsis) | `step1a_DEG.R` | `results/01_DEG/`, Supplementary Table S6 |

## 2. WGCNA and the 179 shared candidates

| Manuscript step | Scripts | Key outputs |
|---|---|---|
| WGCNA module detection (13 ARDS / 11 sepsis disease-associated modules) | `step1b_WGCNA.R` | `results/02_WGCNA/` |
| Cross-cohort module preservation (Zsummary 19.5 / 21.4) | `step55_wgcna_preservation.R` | `results/02_WGCNA/` |
| DEG ∩ module intersection → 179 shared candidates | `step1c_intersection.R` | `results/03_intersection/`, Supplementary Table S7 |
| GO / KEGG enrichment (clusterProfiler) | `step2h_go_kegg_complete.R`, `step2i–2k_go_bp_cnetplot*.R` | `results/03_intersection/`, Figure 2 |

## 3. Machine learning, hub-gene selection, and validation

| Manuscript step | Scripts | Key outputs |
|---|---|---|
| 11 classifiers per discovery cohort, SHAP interpretation, LASSO/RF/SHAP concordance rule | `step2a_ML.R` (ARDS), `step2a2_ML_sepsis.R` (sepsis) | `results/05_ML_SHAP/`, Supplementary Table S9, Figure 3 |
| Label-permutation and small-cohort overfitting diagnostics | `step50_ml_stability_overfit.R` | `results/05_ML_SHAP/`, Supplementary Figure S5 |
| Cross-disease testing within discovery data (AUC 0.844–0.979; DYRK2/TCTN1 reversal) | `step2b2_GSE100159_candidates.R`, `step2b3_GSE243066_sepsis_candidates.R` | `results/06_external_validation/` |
| External validation, single-gene AUCs with CIs, meta-analysis (Hedges g, I²) | `step2b_validation_ARDS.R`, `step2b_GSE100159_validate.R`, `step2b4_GSE100159_ROC.R`, `step2b_AUC_CI.R` | `results/06_external_validation/`, Table 2, Supplementary Table S11 |
| MUC1 refutation sensitivity | `step75_muc1_sensitivity.R` | `results/06_external_validation/` |
| CD177 ablation (AUC 0.961 vs 0.984) | `step49_cd177_ablation.R` | `results/07_hub_immune/cd177_ablation_transfer_AUC.csv` |
| Ridge-regularized logistic nomograms | `step2e_nomogram.R`, `step63_fig4_nomogram_simplified.R` | `results/06_external_validation/`, Figure 5 |
| Exploratory ComBat batch-sensitivity analysis (transfer AUC 0.977) | `step52_validation_diagnostics.R` | `results/06_external_validation/`, Supplementary Figure S7 |

## 4. Prognostic modelling (GSE65682, 479 septic-shock patients)

| Manuscript step | Scripts | Key outputs |
|---|---|---|
| Per-gene and joint six-gene Cox models, risk score, KM, LRT vs clinical model, bootstrap optimism-corrected C index | `step105_survival_prognosis_GSE65682.R`, `step105b_pergene_KM_GSE65682.R` | `results/13_survival_prognosis/`, Supplementary Figure S20 |
| LASSO-Cox over 171 platform-represented candidates (1-SE vs minimum-deviance rules), random survival forest | `step106_extended_modeling_GSE65682.R` | `results/13_survival_prognosis/`, Supplementary Figure S21 |
| Prognostic figure assembly | `step107_sixgene_prognosis_figure.R` | figures |
| Three statistical facts: Cox PH test (global P=0.031), KM cut-point sensitivity (P 0.00026–0.0019), ssGSEA BH coverage | `step109_three_statistical_facts.R` | `results/13_survival_prognosis/tables/step109_*.csv`, Supplementary Table S24 (assembled by `step113`) |

## 5. Mendelian randomization and genetic causality

| Manuscript step | Scripts | Key outputs |
|---|---|---|
| GWAS summary-statistic extraction and harmonization (ieu-b-4980, FinnGen R11) | `extract_MR_4980.py`, `extract_MR_data.sh`, `extract_MR_py.py`, `convert_gwas_ma.py`, `match_ards_4980.py` | `results/04_MR_SMR/` |
| Forward IVW MR (β=−0.336, P=0.249) | `step1d_MR.R`, `step1d_MR2.R`, `step1d_MR_4980.R` | `results/04_MR_SMR/`, Figure 7A–B |
| Reverse MR (P=0.329) | `step48a_extract_reverse_MR.py`, `step48b_reverse_MR.R` | `results/04_MR_SMR/` |
| Post hoc power analysis (MDOR 2.26 forward / 1.04 reverse) | `step51_mr_power.R` | `results/04_MR_SMR/`, Figure 7C, Supplementary Table S13 |
| LD score regression (rg=−0.05; h²≤0.0022) | `prepare_ldsc.py` | `results/04_MR_SMR/` |
| SMR-HEIDI (five hub genes with cis-eQTL) | `step1d_heidi.R` | `results/04_MR_SMR/` |
| Colocalization | `step1d_coloc.R`, `align_coloc.py` | `results/04_MR_SMR/`, Figure 7D |

## 6. Immune microenvironment (bulk)

| Manuscript step | Scripts | Key outputs |
|---|---|---|
| CIBERSORT deconvolution + deconvolution-independent ssGSEA (18 lineage marker sets), BH correction | `step41_cibersort_ssgsea.R`, `step47_ssgsea_immune_4cohorts.R`, `step2c_hub_immune.R` | `results/07_hub_immune/`, Supplementary Table S2 |
| GSVA immune-pathway dysregulation, top-20/30/50 threshold sensitivity (complement) | `step2c2_GSVA_immune.R`, `step2c3_GSVA_immune_top20.R`, `step45_hub6_gsva.R` | `results/07_hub_immune/`, Figure 8G–H, Supplementary Table S15 |

## 7. Single-cell analyses (GSE151263 PBMC; GSE167363 CD14+ monocytes)

| Manuscript step | Scripts | Key outputs |
|---|---|---|
| QC, annotation (Seurat / Scanpy) | `step3a_singlecell.R`, `step3a_annotate.R`, `step3a_scanpy.py`, `step3e*_GSE167363_annotate.R` | `08_singlecell/` (objects excluded) |
| CellChat communication (MHC-I / MIF / MHC-II pathway strengths) | `step3b_cellchat.R`, `step3b_prepare_cellchat.py`, `step3b_cellchat_compare.R`, `step3f_GSE167363_cellchat.R`, `step3m_cellchat_supp.R` | `results/08_singlecell/` |
| netP strength/weight recomputation (pathway-level quantification reported in the text) | `step110_cellchat_netP_recalc.R` (GSE151263), `step111_cellchat_GSE167363_netP_recalc.R` (GSE167363) | `results/08_singlecell/`, Supplementary Table S25 |
| Pseudotime (slingshot, monocle3); sepsis monocytes at later positions | `step3c_slingshot.R`, `step3c_prepare_slingshot.py`, `step3h_GSE167363_slingshot.R`, `step3o_monocle3.R` | `results/08_singlecell/`, Figure 9C/G |
| Dual virtual knockout: CellOracle (CEPBP→UPP1) and scTenifoldKnk (rho 0.63–0.66 vs 0.44–0.50) | `step3d_CellOracle*.py`, `step3d_EXOSC4_ko.R`, `step3d_scTenifoldKnk.R`, `step3i_GSE167363_scTenifoldKnk.R`, `step3i_extra_knockout.R`, `step3j_GSE167363_CellOracle.py`, `step3n_ko_stratified*.R`, `step22_KO_top10_enrichment.R` | `results/08_singlecell/`, Figure 9D/H, Supplementary Figures S16–S18 |
| Per-donor verification (all five donors) | `step53_celloracle_per_donor.py` | `results/08_singlecell/` |
| Donor-level pseudobulk DESeq2 (CD177 log2FC 6.68, P=0.002) | `step54_donor_level.R`, `step60_donor_level_fig.R` | `results/08_singlecell/` |
| Metabolic program scoring (inverted-U), AUCell | `step3n_metabolism.py`, `step3g_aucell.R` | `results/08_singlecell/` |

## 8. Single-gene GSEA (four cohorts)

| Manuscript step | Scripts | Key outputs |
|---|---|---|
| Per-hub-gene GSEA in discovery and validation cohorts | `step4a_single_gene_gsea.R`, `step46_single_gene_gsea_4cohorts.R` | `results/12_single_gene_gsea/`, Figure 4E, Supplementary Figure S6 |

## 9. Protein-level evidence — Human Protein Atlas

| Manuscript step | Scripts | Key outputs |
|---|---|---|
| Blood Atlas pan-disease proteomics (CD177 9/10 infections up, ≤13.6-fold; VNN2 10/10) | `step_hpa_validation.R` | `results/10_HPA_validation/` |
| Immune-cell RNA specificity (CD177 → neutrophils; UPP1/VNN1 → classical monocytes) | `step_hpa_immune_cell.R` | `results/10_HPA_validation/` |

Data vintage: **HPA version 24 (Blood resource / Disease Blood Atlas, 59 diseases), files
downloaded 2025-11-06.** Raw TSVs (> 5 MB) are excluded from this archive; re-download from
<https://www.proteinatlas.org/about/download>. HPA v25 (released 2025-11-11, 71 diseases)
postdates the download and was not used.

## 10. In silico drug repurposing and docking

| Manuscript step | Scripts | Key outputs |
|---|---|---|
| Enrichr Old_CMAP_up/down query of the six-gene signature (dexamethasone rank 1, P=3.7×10⁻⁴) | `step67_cmap_negative_control.py` (live query + negative controls) | `results/09_docking_drug/cmap_hub6_realtime_query.csv` |
| Negative controls: 100 random non-hub sets, 100 transcriptome-wide sets | `step67_cmap_negative_control.py` | `results/09_docking_drug/cmap_negative_control_*.csv`, Supplementary Figure S3 |
| Stringency-matched background (genes up-regulated in both discovery cohorts; 2/100) | `step76_cmap_deg_null.py` | `results/09_docking_drug/cmap_deg_matched_null.csv` |
| AutoDock Vina docking of 10 natural-product ligands (luteolin/VNN1 −10.01 kcal/mol) | `step4_docking.py`, `step4_dock_one.py`, `step4_docking_site.py` | `results/09_docking_drug/`, Table 3, Supplementary Table S19 |
| Redocking validation (incl. UPP1 fixed box) | `step56_redocking.py`, `step57_upp1_redock_fixed_box.py` | `results/09_docking_drug/` |
| AlphaFold pLDDT confidence assessment | `step66_alphafold_plddt.py` | `results/09_docking_drug/alphafold_plddt_assessment.csv` |
| PLIP interaction analysis (VNN1–luteolin) | `step69_plip_vnn1_luteolin.py` | `results/09_docking_drug/`, Figure 10D |
| MM/GBSA rescoring (OpenMM; Gasteiger/UFF fallback disclosed in Methods) | `step70_mmgb_sasa.py` | `results/09_docking_drug/`, Supplementary Table S20 |
| Open Targets tractability profiling | `step68_opentargets_druggability.py` | `results/09_docking_drug/opentargets_*.csv` |
| ChEMBL known-drugs check (no approved/clinical small molecule for any of the six) | `step115_chembl_known_drugs.py` | `results/09_docking_drug/chembl_known_drugs.csv` |

Data vintage: **Open Targets Platform, accessed September 2026** (scripted, step 68);
**ChEMBL release 37, accessed 2026-09-17**. The ChEMBL result was originally produced by an
ad-hoc query; `step115` (written 2026-10-21, while release 37 was still live) reproduces the
archived CSV byte-identically. The five `*_hub_prediction.csv` files (DSigDB, DisGeNET,
Jensen DISEASES, Drug Perturbations from GEO 2014, Old CMAP) are Enrichr web-interface
exports from 2026-08-24; all quantitative CMap claims in the manuscript use the scripted
re-queries above (steps 67/76).

## 11. Network pharmacology (PPI)

| Manuscript step | Scripts | Key outputs |
|---|---|---|
| STRING PPI network and pathway annotation of the 179 candidates | `step5_network_ppi.R` | `results/11_network_pharma/` |

## 12. Figures, tables, and manuscript integrity

| Item | Scripts |
|---|---|
| Main-figure panel generation and composite assembly (Figures 1–10) | `step30–37`, `step42–44`, `step58`, `step61`, `step63–65`, `step71`, `step79–89` |
| Supplementary figures (S1–S21) | `step21–24`, `step27–29`, `step38`, `step42b`, `step46b`, `step47b`, `step74`, `step77–78`, `step90–93`, `step98–99`, `step102` |
| Supplementary tables (S1–S26; S24–S26 assembled 2026-10-19) | `step26_build_tables.py`, `step62_build_tables_S13_S17.py`, `step72_build_tableS21.py`, `step108b_supp_tables.py`, `step113_assemble_supplementary_tables_S24_S25_S26.py` |
| Figure-format export (PDF + PNG 300 dpi + TIFF 1000 dpi LZW) | `step15_export_figures.py`, `step112_export_missing_figure_formats.py` |
| Manuscript compliance check (word counts, citation order, style rules, paired-claim anchors) | `verify_manuscript.py` |
| Manuscript .docx rendering | `step95_make_revised_docx.py` |
| Code-archive build (this repository) | `make_code_archive.sh` |

## Notes

- R 4.6.1; Python 3. Key R packages: limma, WGCNA, clusterProfiler, TwoSampleMR, survival,
  rms, glmnet, ranger, CIBERSORT, GSVA, CellChat, monocle3, Seurat, DESeq2. Full pinned
  list: `environment/R_packages.csv`.
- Random seeds are fixed inside the scripts (random survival forest, LASSO CV, bootstrap,
  negative-control sampling: `random.seed(20260916)` in step 67).
- All statistical tests are two-sided with Benjamini–Hochberg adjustment where stated.
