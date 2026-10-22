# Supplementary Materials

This file accompanies the manuscript "Shared Hub Genes and Immune Signatures of Sepsis and ARDS: A Multi-Omics Integration." With the Path-A restructure of 2026-10-16, the complete Methods moved from the main text into the Online Supplement (Supplementary Methods S1–S14); the main-text METHODS is a summary pointing here. Tables and figures are numbered S1-S26 and S1-S21 by first-citation order in the revised main text (renumbering of 2026-09-25; see 05_Renumbering_Maps.md; Tables S24–S26 added 2026-10-19). Source paths refer to the analysis archive and retain their original directory names; the renumbered files accompanying this package are in figures/supplementary/ and figures/supplementary_tables/. Citations within this supplement use S-numbered references (see Supplementary References at the end of Supplementary_Methods.md).

## Contents

- **Supplementary Methods** (see `Supplementary_Methods.md`)
- **Supplementary Tables** (below)
- **Supplementary Figure Legends** (below)

---

## Supplementary Tables

### Table S1. Discovery-internal cross-disease screening of candidate hub genes

Single-gene AUCs of every machine-learning-nominated candidate tested within the discovery cohort of the opposite disease (ARDS candidates in GSE69528; sepsis candidates in GSE76293), with exact Mann-Whitney P values (floating-point floor 2.2 × 10⁻¹⁶) and the pass/fail status of the same-direction AUC>0.7 screen. Seven of nine ARDS-nominated candidates (the six final genes plus MUC1) and four of six sepsis-nominated candidates (METTL9, MTHFD2, BATF, PLP2) passed; DYRK2, TCTN1, and ZSCAN18 showed direction-reversed discrimination and were excluded. This screen uses no external validation data; external validation subsequently confirmed six of the seven ARDS-nominated candidates and refuted MUC1 (AUC=0.575 in GSE243066).
*Source: `05_ML_SHAP/discovery_crossdisease_AUC.csv`; script `scripts/step27_reviewer_fix_R.R`*

### Table S2. Immune-cell infiltration results (CIBERSORT)

Per-sample CIBERSORT-estimated relative abundances of 22 immune-cell subtypes, deconvolution quality metrics (permutation P value, Pearson correlation, and RMSE; 1,000 permutations), and two-sided Wilcoxon rank-sum group comparisons for all 22 subtypes in the four bulk cohorts (GSE76293, GSE69528, GSE243066, GSE100159). Note: estimates for the purified PMN cohort (GSE76293) are provided for transparency only; deconvolution of non-granulocyte lineages from PMN fractions is not interpretable and was not used to support any conclusion.
*Source: `07_hub_immune/cibersort/*_CIBERSORT_results.csv`, `*_celltype_wilcox.csv`*

### Table S3. Immune-cell infiltration results (ssGSEA)

Per-sample ssGSEA scores for 18 immune-cell types (marker sets assembled from published lineage markers; full list in the `marker_genes` sheet), two-sided Wilcoxon rank-sum group comparisons (with Benjamini–Hochberg-adjusted P values), and Spearman rank correlations between ssGSEA cell-type scores and the six hub genes, for all four bulk cohorts (GSE76293, GSE69528, GSE243066, GSE100159; 13 sheets including `marker_genes`). Note: consistent with the CIBERSORT caveat, ssGSEA estimates for the purified PMN cohort (GSE76293) are provided for transparency only and were not used to support any conclusion.
*Source: `07_hub_immune/ssgsea_immune/*_ssgsea_scores.csv`, `*_ssgsea_wilcox.csv`, `*_ssgsea_hub_spearman.csv`; script `scripts/step47_ssgsea_immune_4cohorts.R`*

### Table S4. CMap glucocorticoid-specificity analysis

Top ten compounds per Old CMAP library (up and down) for the six-hub-gene query signature with nominal P values and combined scores; random gene-set null results under three backgrounds (100 random six-gene sets from non-hub shared candidate genes, 100 from the transcriptome-wide measured universe, and 100 DEG-matched sets from the 123 genes up-regulated in both discovery cohorts with the hub genes excluded, all with an identical Enrichr query protocol), with per-set glucocorticoid top-10 indicator and best glucocorticoid-term P value. The DEG-matched background matches both disease relevance and coordinated up-regulation of the real signature: 2 of 100 such nulls matched or exceeded the real glucocorticoid P value (3.7 × 10⁻⁴) and 18 of 100 placed a glucocorticoid in the top ten.
*Source: `09_docking_drug/Old_CMAP_*_hub_prediction.csv`, `09_docking_drug/cmap_negative_control_*.csv`, `09_docking_drug/cmap_deg_matched_null.csv`; scripts `scripts/step67_cmap_negative_control.py`, `scripts/step76_cmap_deg_null.py`*

### Table S5. AlphaFold model-confidence assessment (pLDDT)

Per-structure pLDDT statistics (read from the B-factor column of each AlphaFold model) for the four targets without ligand-complexed experimental structures: global mean/median and fraction of residues at pLDDT ≥ 70 and ≥ 90, plus the same statistics restricted to the docking-box residues.
*Source: `09_docking_drug/alphafold_plddt_assessment.csv`; script `scripts/step66_alphafold_plddt.py`*

### Table S6. Complete differential expression results

Full DEG lists for GSE76293 (ARDS vs. control) and GSE69528 (sepsis vs. control), including log2 fold change, P value, and adjusted P value for all genes.
*Source: `01_DEG/GSE76293_DEG_all.csv`, `01_DEG/GSE69528_DEG_all.csv`*

### Table S7. The 179 shared candidate genes

The intersection of DEGs and disease-associated WGCNA module genes in both cohorts.
*Source: `03_intersection/`*

### Table S8. Machine-learning overfitting diagnostics

Four sections. (1) Label-permutation test of the 11 classifiers in the ARDS training cohort (GSE76293; 500 label shuffles per classifier; P=(1 + number of null AUCs at or above the observed AUC)/(1 + 500)): every classifier that could be refit under permutation exceeded its null distribution (permutation P=0.002 to 0.014; GBM could not be refit), showing that the cross-validated AUCs are not label-noise artifacts despite the small cohort. (2) Apparent versus cross-validated AUC for all 11 classifiers in both training cohorts, with the apparent-to-CV gap. (3) Exact one-sided 97.5% and 95% lower confidence bounds for the observed AUC=1.000 results (ARDS training 12 vs 12: lower bound 0.975; in-cohort refits in both validation cohorts: 0.993 and 0.991).
*Source: `05_ML_SHAP/stability/permutation_500x_results.csv`, `apparent_vs_cv_AUC_*.csv`, `auc_exact_lower_bounds.csv`; script `scripts/step50_ml_stability_overfit.R`*

### Table S9. Feature-selection stability of the machine-learning prioritization

Six sections. LASSO bootstrap selection frequencies (200 resamples) for every candidate gene in both training cohorts (final hub genes selected in 37% to 60% of ARDS resamples; the sepsis finalists were stable across 41% to 94%), random-forest seed stability (frequency of top-30 and final-candidate inclusion across 50 random seeds; the six final genes appeared in 94% to 100% of ARDS seeds), and Jaccard similarity of the candidate sets under plus-or-minus 10% perturbation of the two ranking thresholds (0.89 to 1.00 in both cohorts). *VNN1* was retained by random forest in every seed but fell below the top-20 mean-AUC criterion and entered the final list through the LASSO arm of the prespecified selection rule.
*Source: `05_ML_SHAP/stability/lasso_bootstrap_stability_*.csv`, `rf_seed_stability_*.csv`, `threshold_perturbation_jaccard_*.csv`; script `scripts/step50_ml_stability_overfit.R`*

### Table S10. Module membership and co-expression degree of the six hub genes

WGCNA module assignment, module membership (kME) within the top disease-associated module of each discovery cohort (turquoise in GSE76293, module-trait r=0.925, P=1.1 × 10⁻¹⁰; blue in GSE69528, r=0.752, P=1.8 × 10⁻²¹), and co-expression degree centrality (number of the 178 other shared candidate genes connected at |r| > 0.5 and |r| > 0.7) for each of the six hub genes in both discovery cohorts. All six genes have kME of 0.70 to 0.93 and 96 to 142 co-expression partners at |r| > 0.7, quantitatively supporting their intramodular hub status.
*Source: `02_WGCNA/hub_gene_kME.csv`; script `scripts/step27_reviewer_fix_R.R`*

### Table S11. Validation-cohort diagnostic performance and cross-platform sensitivity analyses

Five sections. (1) Single-gene AUCs with DeLong 95% confidence intervals in both validation cohorts. (2) Within-cohort z-standardized distribution overlap coefficients (OVL) between training and validation cohorts, standardized-mean-difference shifts, and within-cohort disease-versus-control effect sizes (Cohen d) in training and validation, with effect attenuation: within-cohort effects remained large in both validation cohorts (for example, *UPP1* Cohen d 2.72 in training versus 0.89 in validation), and standardized distributions overlapped substantially (OVL 0.66 to 0.88), explaining the coefficient-transfer shortfall. (3) Meta-analysis of per-gene standardized effects across the two validation cohorts (Hedges g under fixed- and random-effects DerSimonian-Laird models with I-squared): pooled g 1.14 to 2.32 (all P<0.00001), confirming the genes separate disease from control within each platform. (4) ComBat cross-platform sensitivity analysis of coefficient transfer (exploratory; pools discovery and validation data and therefore breaks the validation firewall): sepsis transfer AUC 0.708 to 0.977, ARDS 0.984 to 0.978, indicating a platform-scaling artifact as the dominant contributor to the sepsis transfer shortfall. (5) Seven-gene (including *MUC1*) versus six-gene model sensitivity analysis, bounding the optimistic bias introduced by using the validation cohorts to screen the seven ARDS-nominated candidates down to the final six: including *MUC1* left the ARDS validation AUC unchanged (0.984; 95% CI, 0.953 to 1.000) and slightly reduced the sepsis transfer AUC (0.708 to 0.699), so the reported validation performance was not inflated by the screening step. (6) *CD177* ablation sensitivity analysis: refitting the ridge-regularized ARDS model without *CD177* (five genes) transferred to GSE243066 at AUC=0.961 versus 0.984 for the full six-gene set (training AUC 1.000 for both), showing that transferred discrimination does not depend on any single gene.
*Source: `06_external_validation/reviewer_w5/single_gene_AUC_CI_validation.csv`, `distribution_overlap_train_vs_validation.csv`, `metaDE_single_gene_validation.csv`, `combat_sensitivity_model_transfer.csv`, `06_external_validation/muc1_sensitivity.csv`, `07_hub_immune/cd177_ablation_transfer_AUC.csv`; scripts `scripts/step52_validation_diagnostics.R`, `scripts/step75_muc1_sensitivity.R`*

### Table S12. Disease-versus-disease ROC and donor-level sensitivity analyses

Two sections. (1) Patient-level pseudobulk ROC (ARDS versus sepsis) in the single-cell PBMC cohort (GSE151263; three ARDS versus four sepsis donors): per-gene AUCs (0.33 to 0.75) and the combined six-gene z-score signature (AUC=0.50), with exact Mann-Whitney P values (all P>0.4), demonstrating that the signature does not discriminate between the two disease states. (2) Donor-level sensitivity analyses for the key single-cell contrasts: per-donor median pseudotime (Monocle3), per-donor median AUCell UPP1-program recovery scores in monocytes, and per-donor mean metabolic pathway scores (glycolysis, oxidative phosphorylation, fatty acid metabolism, xenobiotic metabolism) in both single-cell cohorts. For GSE167363, the two samples from the repeatedly sampled sepsis patient were collapsed to a single donor (four patients: two sepsis, two healthy controls). All donor-level statistics are directionally consistent with the cell-level contrasts reported in the main text; none reaches significance at the donor level because of the small number of donors (two to four per group).
*Source: `08_singlecell/disease_vs_disease_ROC_GSE151263.csv`, `08_singlecell/donor_sensitivity/donor_level_summary.csv`; script `scripts/step27_reviewer_fix_py.py`*

### Table S13. Mendelian randomization power analysis

Analytic power calculations (S11; mRnd-equivalent) and 5,000-draw simulations for both directions of the two-sample MR, using the IVW standard error under the multiplicative random-effects model as reported by TwoSampleMR. Forward analysis (sepsis to ARDS; 18 instruments, mean F=22.7): 28% power at a causal odds ratio of 1.5, and a minimum detectable odds ratio of 2.26 at 80% power, so the forward null is nearly uninformative for small effects. Reverse analysis (ARDS to sepsis; 7 instruments, mean F=21.8): minimum detectable odds ratio 1.04, so the reverse null does exclude small effects. All instrument-level inputs (effect sizes, standard errors, harmonization exclusions) are provided.
*Source: `04_MR_SMR/MR_power_summary.csv`, `MR_power_by_OR.csv`, `power_forward_sepsis_to_ARDS_IVs.csv`, `power_reverse_ARDS_to_sepsis_IVs.csv`; script `scripts/step51_mr_power.R`*

### Table S14. Mendelian randomization and SMR/colocalization results

Full MR estimates (IVW, MR-Egger, weighted median, weighted mode), heterogeneity, pleiotropy, and Steiger directionality, plus SMR and coloc results for the hub genes. SMR-HEIDI was applied to the five hub genes with an available blood cis-eQTL instrument in eQTLGen: *CDK5RAP2*, *CD177*, *GRB10*, *UPP1*, and *VNN1*; *EXOSC4* lacked a blood cis-eQTL and was not testable.
*Source: `04_MR_SMR/`*

### Table S15. Full GSVA results for 537 immune pathways

GSVA enrichment scores and limma differential results for all immune-related GO biological-process gene sets. Threshold sensitivity of the top-ranked shared pathways: complement activation is the only pathway shared by the top-20 dysregulated processes (ranked by adjusted P) of both discovery cohorts; relaxing the cutoff to the top 30 yields five shared pathways (adding ACUTE_INFLAMMATORY_RESPONSE, DEFENSE_RESPONSE_TO_FUNGUS, MACROPHAGE_ACTIVATION, and REGULATION_OF_PHAGOCYTOSIS), the top 50 yields seven (further adding KILLING_BY_HOST_OF_SYMBIONT_CELLS and POSITIVE_REGULATION_OF_DEFENSE_RESPONSE_TO_BACTERIUM), and the top 100 yields twenty-four; all shared pathways were concordant in direction (upregulated in disease).
*Source: `07_hub_immune/gsva_top20/*_GSVA_limma_all.csv`*

### Table S16. AUCell gene-set activity statistics for the six hub genes in both single-cell cohorts

Per-cell AUCell recovery scores for each hub-gene co-expression program (top-50 positively correlated genes plus the hub gene) and the composite HubSignature, with mean scores, two-sided Wilcoxon rank-sum statistics, and Benjamini–Hochberg-adjusted P values for each gene-set-by-cell-type comparison. GSE167363 (CD14+ monocytes) is compared as sepsis versus healthy controls (42 comparisons, 32 significant at FDR < 0.05); GSE151263 (PBMCs) is compared as ARDS versus sepsis (49 comparisons, 28 significant). The complete gene-set definitions are also provided.
*Source: `08_singlecell/GSE167363/aucell/GSE167363_AUCell_stats.csv`, `08_singlecell/aucell/GSE151263_AUCell_stats.csv`, `*_AUCell_genesets.csv`*

### Table S17. Donor-level analyses of the single-cell cohorts

Four sections. (1) Per-donor CellOracle analysis of the five GSE167363 donors, each with an independently built gene regulatory network (identical parameters): the CEBPB-knockout effect on *UPP1* was negative in all five donors (delta −0.33 to −2.08), with the largest magnitudes in the later and more severe sepsis samples. (2) Per-donor UPP1-pseudotime Spearman correlations and per-donor pseudobulk CEBPB and UPP1 expression (log2 counts per million) in GSE167363: four of five donors showed the same negative direction as the pooled analysis. (3) Pseudobulk DESeq2 (raw counts aggregated per donor and cell type) for GSE167363 (sepsis versus control, 2 versus 2 donors): *CD177* in classical monocytes was the only contrast significant at the donor level (log2 fold change 6.68, adjusted P=0.002). (4) Pseudobulk DESeq2 for GSE151263 (ARDS versus sepsis, 3 versus 4 donors).
*Source: `08_singlecell/celloracle_per_donor/CellOracle_CEBPB_KO_per_donor.csv`, `08_singlecell/donor_level/*.csv`; scripts `scripts/step53_celloracle_per_donor.py`, `scripts/step54_donor_level.R`*

---

### Table S18. Virtual-knockout results in CD14+ monocytes (GSE167363)

scTenifoldKnk differential-regulation results for the four knockable hub genes (UPP1, CDK5RAP2, EXOSC4, VNN1), including knockout feasibility (expressed-cell fraction and quality-control threshold), hub-gene cross-effect distances, and the CellOracle transcription-factor-by-gene knockout-effect matrix (delta of simulated versus observed expression, with the CEBPB–UPP1 effect of −0.912 as the strongest link).
*Source: `08_singlecell/GSE167363/knockout/`, `08_singlecell/GSE167363/celloracle/CellOracle_TF_KO_hub_effect.csv`*

### Table S19. Complete molecular-docking results (60 combinations)

Binding energies for ten natural-product ligands against six hub-gene protein targets.
*Source: `09_docking_drug/`*

### Table S20. MM/GBSA free-energy rescoring of docked poses

Single-point MM/GBSA binding free energies (kcal/mol) for the 50 target–ligand combinations of the five high-confidence structures (CDK5RAP2 excluded owing to moderate box-region pLDDT): complex, protein, and ligand single-point energies, ΔSASA, and the resulting ΔG, alongside the AutoDock Vina binding energies and their Spearman rank correlation. Protein: ff14SB with GBn2 implicit solvent; ligand: Gasteiger charges and UFF Lennard-Jones parameters; restrained local minimization (see Methods and Limitations).
*Source: `09_docking_drug/mmgb_*.csv`; script `scripts/step70_mmgb_sasa.py`*

### Table S21. Open Targets tractability annotations for the six hub genes

Open Targets Platform tractability buckets (modality, label) per gene, number of tractable buckets, subcellular-location evidence, and known-drug queries from ChEMBL (no approved or clinical-stage small-molecule drug recorded for any target).
*Source: `09_docking_drug/opentargets_druggability.csv`, `09_docking_drug/chembl_known_drugs.csv`; script `scripts/step68_opentargets_druggability.py`*

### Table S22. Prognostic model comparison in GSE65682

Discrimination and incremental-value metrics for the prognostic models fitted in GSE65682 (479 patients with septic shock; 28-day mortality 23.8%): apparent and bootstrap-optimism-corrected C index (300 resamples), time-dependent AUC at 7, 14, and 28 days, likelihood-ratio tests against the clinical model (age, sex, and Mars endotype), and net reclassification index and integrated discrimination improvement at 28 days. Models: Clinical (age, sex, endotype); Hub6 (six hub genes, z-standardized); Clinical + Hub6; LASSO-179 (penalized Cox signature over the 171 shared candidate genes represented on the platform; penalty at minimum 10-fold cross-validated deviance, 15 genes retained); Clinical + LASSO-179; and a random survival forest over the same predictors (out-of-bag and 10-fold out-of-fold estimates).
*Source: `13_survival_prognosis/tables/ext_model_performance.csv`, `ext_bootstrap_optimism.csv`, `ext_lrt.csv`, `ext_nri_idi.csv`; scripts `scripts/step105_survival_prognosis_GSE65682.R`, `scripts/step106_extended_modeling_GSE65682.R`*

### Table S23. LASSO-179 prognostic signature coefficients

Coefficients of the 15 genes retained by the LASSO Cox model at the penalty minimizing 10-fold cross-validated deviance in GSE65682, with hazard ratios per unit of log2 expression; none of the six hub genes was retained, and the more conservative one-standard-error rule retained no gene at all.
*Source: `13_survival_prognosis/tables/ext_lasso_coefficients.csv`; script `scripts/step106_extended_modeling_GSE65682.R`*

### Table S24. Sensitivity analyses of the six-gene prognostic model in GSE65682

Three sections. (1) Proportional-hazards diagnostics (cox.zph): the joint six-gene Cox model violated the proportional-hazards assumption globally (chi-square P=0.031) although no single gene did (minimum per-gene P=0.057, *EXOSC4*); the full model including age, sex, and Mars endotype passed (global P=0.131), so joint hub-gene estimates are read as time-averaged. (2) Kaplan-Meier cut-point sensitivity of the six-gene risk score: median split log-rank P=0.0019, tertile-extreme contrast P=0.00026, quartile-extreme contrast P=0.0015, continuous-score test P=0.00065, and three-group trend P=0.00049. (3) Multiplicity coverage of the immune-cell comparisons: CIBERSORT (22 subtypes, four cohorts) already carried Benjamini–Hochberg-adjusted P values (two cell types with all-zero fractions were not testable); ssGSEA (18 subtypes) originally lacked adjustment and was recomputed with Benjamini–Hochberg correction, with 7 of 18 (GSE76293) and 15 of 18 (GSE69528) subtype comparisons remaining significant.
*Source: `13_survival_prognosis/tables/step109_coxzph_hub6.csv`, `step109_coxzph_full.csv`, `step109_cutpoint_sensitivity.csv`, `07_hub_immune/GSE76293_immune_diff_BH.csv`, `GSE69528_immune_diff_BH.csv`; script `scripts/step109_three_statistical_facts.R`*

### Table S25. CellChat pathway-level communication strength comparison (recalculated 2026-10-19)

Pathway-level communication strength (netP$strength; total communication probability summed over cell-type pairs) and weight (number of significant ligand–receptor pairs, P≤0.05) recalculated for both single-cell cohorts, after confirming that netP$strength and netP$weight were absent from every archived object. Two comparisons are provided. (1) GSE151263 peripheral-blood mononuclear cells, sepsis versus ARDS (`CellChat_ARDS.rds`, `CellChat_Sepsis.rds`): MHC-I and MIF rank first and second in both disease states; aggregate MIF strength is 1.24 × 10⁻³ in sepsis versus 6.43 × 10⁻⁴ in ARDS (1.93-fold higher in sepsis). (2) GSE167363 CD14+ monocytes, sepsis versus healthy control (`CellChat_Sepsis.rds`, `CellChat_HC.rds`): MIF ranks first in sepsis (3.16 × 10⁻³ versus 2.85 × 10⁻³ in controls), whereas MHC-II leads the healthy network (1.34 × 10⁻² versus 1.87 × 10⁻³). For each cohort the full per-pathway differential network is given (both strengths, difference, ratio, dominant side, and interaction counts).
*Sources: `08_singlecell/cellchat/netP_strength_weight_ARDS.csv`, `08_singlecell/cellchat/netP_strength_weight_Sepsis.csv`, `08_singlecell/cellchat/pathway_diff_network_recalc_2026-10-19.csv`, `08_singlecell/GSE167363/cellchat/netP_strength_weight_Sepsis.csv`, `08_singlecell/GSE167363/cellchat/netP_strength_weight_HC.csv`, `08_singlecell/GSE167363/cellchat/pathway_diff_network_GSE167363_recalc_2026-10-19.csv`; scripts `scripts/step110_cellchat_netP_recalc.R`, `scripts/step111_cellchat_GSE167363_netP_recalc.R`*

### Table S26. Cohort metadata summary

Per-cohort sample composition, tissue, platform, sampling time, and ancestry. GSE151263: 13,996 cells from seven samples (ARDS 6,000 cells from three donors, 2,000 each; sepsis 7,996 cells from four donors, 1,998 to 2,000 each); cell-type composition: monocytes 4,050, CD4 T 2,979, CD8 T 1,821, platelets 1,700, NK 1,448, B 1,323, dendritic cells 675. GSE167363: 12,359 CD14+ monocytes from five samples (two healthy controls, two sepsis patients, one sampled twice). GSE65682: 479 patients with septic shock, PAXgene blood RNA at ICU admission; no SOFA or APACHE severity fields are available in the series matrix. GSE243066: blood collected within 24 hours of admission. Sampling time and ancestry for GSE76293, GSE69528, GSE100159, and GSE243066: [AUTHOR_INPUT_NEEDED: confirm from GEO records]; ancestry of the MR data sources (ieu-b-4980, FinnGen R11) and of the eQTLGen reference: European.
*Source: `08_singlecell/GSE151263_cell_annot.csv`, `00_rawdata/validation/GSE65682/GSE65682_series_matrix.txt`*

---

## Supplementary Figure Legends

### Figure S1. WGCNA soft-threshold selection and module dendrograms

Scale-free topology fit index and mean connectivity versus soft-thresholding power for GSE76293 and GSE69528. (A–C) GSE76293 (ARDS): scale-free topology fit R² versus soft-thresholding power (A; red dashed line, R²=0.8; selected power=10), mean connectivity versus power (B), and gene clustering dendrogram with detected module colors (C). (D–F) The same analyses for GSE69528 (sepsis; selected power=6). The top-3,000 genes by median absolute deviation were used with an unsigned network (minimum module size, 30; merge cut height, 0.25).
*Source: `05_supplementary/Supplementary_Figures_Tables/Figures/FigureS1_WGCNA_soft_threshold_dendrograms.{pdf,png,tiff}`; scripts `scripts/step23_figS1_wgcna.R`, `scripts/step25_assemble_supplementary.py`*

### Figure S2. Donor-level sensitivity analyses for key single-cell contrasts

Each point is one donor (sample); horizontal bars mark medians. Cell-level group comparisons in the main text treat cells as independent units and therefore overstate statistical precision; this figure shows the corresponding donor-level aggregations, which are directionally consistent with every key cell-level contrast but underpowered at two to four donors per group. (A) PBMC pseudotime (Monocle3), GSE151263, three ARDS versus four sepsis donors (per-donor medians 7.85 versus 7.16; cell-level P=0.055). (B) CD14+ monocyte pseudotime, GSE167363, two healthy-control versus two sepsis donors, with the two samples from the repeatedly sampled sepsis patient collapsed to a single donor (per-donor medians 2.43 versus 3.44; cell-level P=9.6 × 10⁻⁸¹). (C) AUCell UPP1-program recovery scores in monocytes, GSE151263 (per-donor medians 0.672 versus 0.715). (D) Per-donor mean glycolysis pathway scores in both single-cell cohorts (GSE167363: 0.0355 healthy controls versus 0.0554 sepsis; GSE151263: 0.0607 sepsis versus 0.0485 ARDS). All donor-level values with group assignments are provided in Supplementary Table S12.
*Source: `05_supplementary/Supplementary_Figures_Tables/Figures/FigureS2_donor_sensitivity.{pdf,png,tiff}`; script `scripts/step27_reviewer_fix_py.py`*

### Figure S3. CMap glucocorticoid enrichment against random gene-set nulls

Specificity of the glucocorticoid enrichment of the six-hub-gene signature. (A) Null distributions of the best glucocorticoid-term enrichment P value across random six-gene sets drawn from the 173 non-hub shared candidate genes (blue) and from the 123 genes up-regulated in both discovery cohorts with the hub genes excluded (DEG-matched background, orange; identical Enrichr query protocol); the red line marks the real signature (best glucocorticoid P=3.7 × 10⁻⁴), matched or exceeded by 13 of 100 candidate-background nulls (two strictly exceeded it; 11 tied the same discrete enrichment P) and by 2 of 100 DEG-matched nulls. (B) Frequency with which a glucocorticoid term appears in the top ten hits: 100% for the real signature versus 17% of candidate-background random sets, 18% of DEG-matched random sets, and 0% of transcriptome-wide random sets.
*Source: `05_supplementary/Supplementary_Figures_Tables/Figures/FigureS18_CMap_null_distribution.{pdf,png,tiff}` (`09_docking_drug/cmap_negative_control_*`, `09_docking_drug/cmap_deg_matched_null.csv`); scripts `scripts/step67_cmap_negative_control.py`, `scripts/step76_cmap_deg_null.py`, `scripts/step77_figS18_three_backgrounds.py`*

### Figure S4. Redocking validation of the docking protocol

Validation of the AutoDock Vina protocol by redocking the co-crystallized ligands into their cognate chain-A receptor structures (22 Å ligand-centroid box, exhaustiveness 16, five output modes). (A) UPP1 (PDB 3EUF, co-crystallized ligand BAU): heavy-atom root-mean-square deviation (RMSD) of each output mode from the crystallographic pose; the top-ranked mode reproduced the crystallographic pose at 1.54 Å. (B) VNN1 (PDB 7SLV, co-crystallized ligand 9TC): the top-ranked mode reproduced the crystallographic pose at 0.69 Å. Dashed lines mark the 2.0 Å acceptance criterion; Vina affinity values are annotated below each mode. RMSDs are in-place (no alignment) and symmetry-corrected.
*Source: `05_supplementary/Supplementary_Figures_Tables/Figures/FigureS19_Redocking_validation.{pdf,png,tiff}` (`09_docking_drug/redocking/redocking_results.csv`); script `scripts/step74_redocking_figS19.py`*

### Figure S5. Individual machine-learning classifier performance

ROC curves and confusion-matrix summaries for each of the 11 classifiers. (A, B) Receiver-operating-characteristic curves of the 11 classifiers (LASSO, ridge, elastic net, random forest, support vector machine, k-nearest neighbors, naive Bayes, neural network, logistic regression, XGBoost, and decision tree) built on the 179 shared candidate genes, plotted from pooled out-of-fold predictions of five-fold stratified cross-validation in GSE76293 (ARDS; A) and GSE69528 (sepsis; B). (C, D) Confusion matrices (classification threshold, 0.5) for each classifier in the two cohorts.
*Source: `05_supplementary/Supplementary_Figures_Tables/Figures/FigureS3_ML_ROC_confusion.{pdf,png,tiff}`, out-of-fold predictions `05_ML_SHAP/*_OOF_predictions.csv`; script `scripts/step24_figS2_ML_roc_cm.R`*

### Figure S6. Single-gene GSEA of the six hub genes across four cohorts

Running-enrichment-score plots (fgsea, MSigDB Hallmark gene sets v2024.1) for each of the six hub genes in all four bulk cohorts. For each gene, the top two pathways ranked by absolute normalized enrichment score (NES) are shown as running-enrichment-score curves (top) over the gene's co-expression rank, with the rank metric (Pearson correlation between the hub gene and every other gene) plotted underneath; genes belonging to either pathway are highlighted. (A–F) GSE76293 (ARDS discovery): all six genes are negatively enriched for MYC target genes and DNA repair, consistent with proliferative arrest of activated neutrophils. (G–L) GSE69528 (sepsis discovery): five hub genes (CD177, GRB10, UPP1, EXOSC4, VNN1) are positively enriched for TNF signaling via NF-κB and IL6/JAK/STAT3 signaling, whereas CDK5RAP2 is most strongly anti-correlated with MYC target genes (still positively enriched for IL6/JAK/STAT3). (M–R) GSE243066 (ARDS validation): all six genes are positively enriched for interferon-alpha and interferon-gamma responses. (S–X) GSE100159 (sepsis validation): the interferon axis recurs together with TNF and IL6/JAK/STAT3 signaling, and the MYC-target anti-correlation persists for GRB10.
*Source: `12_single_gene_gsea/FigureS4_single_gene_gsea.{pdf,png,tiff}` (original two-cohort render, superseded), rebuilt as `AJRCCM_submission/03_figures/Figures_ggplot2_NPG/composites/FigureS4_composite.pdf` from `single_gene_gsea_all_results_4cohorts.csv` and `single_gene_gsea_top2_4cohorts.csv`; script `scripts/step46_single_gene_gsea_4cohorts.R` (discovery cohorts reuse the archived results of `scripts/step4a_single_gene_gsea.R`)*

### Figure S7. Per-cohort principal-component analysis (batch-effect assessment)

Principal-component analysis of the four bulk cohorts, performed independently within each cohort on the 2,000 most variable genes (median absolute deviation), with samples colored by disease versus control status: (A) GSE76293 (ARDS discovery; Affymetrix HG-U133 Plus 2.0; 12 vs 12), (B) GSE69528 (sepsis discovery; Affymetrix HG-U133A), (C) GSE243066 (ARDS validation; Illumina NovaSeq 6000; 34 vs 15), and (D) GSE100159 (sepsis validation; Illumina HumanHT-12 v4; 33 vs 12). Within each cohort, samples separate primarily along the disease-control axis on the first two components, with no evidence of hidden batch structure; because the four cohorts were analyzed independently and combined only at the level of per-cohort statistics (differential expression, effect sizes, meta-analysis), cross-platform normalization was not required, and the ComBat analysis of Table S11 is reported as an exploratory sensitivity analysis only. Full sample coordinates are provided in the source CSV.
*Source: `03_figures/Figures_final/FigureS5_composite.{pdf,png,tiff}`, coordinates `06_external_validation/reviewer_w5/cohort_PCA_coords.csv`; script `scripts/step59_cohort_pca.R`*

### Figure S8. Mendelian randomization shows no independent causal effect

(A) Scatter plot of single-nucleotide polymorphism (SNP) effect sizes for sepsis (UK Biobank, ieu-b-4980; 11,643 cases) versus ARDS (FinnGen R11 J10_ARDS; 431 cases) for the 18 harmonized instrumental variables (instrument threshold P<5 × 10⁻⁶; mean F=22.8), with fitted slopes for inverse-variance-weighted (IVW), MR-Egger, weighted-median, and weighted-mode estimators. (B) Forest plot of per-SNP Wald ratios and the overall IVW and MR-Egger estimates. (C) Sensitivity analyses: all four estimators returned null causal estimates (IVW β=−0.336; 95% CI, −0.908 to 0.236; P=0.249), with no heterogeneity (Cochran Q=13.79, df = 17, P=0.682), no directional pleiotropy (MR-Egger intercept = 0.0035, P=0.954), and correct causal direction by Steiger filtering (P=2.65 × 10⁻³⁴). (D) SMR-HEIDI tests for the five hub genes with a cis-eQTL (eQTLGen blood): all ten gene-by-outcome tests were non-significant for both the SMR causal estimate (P_SMR > 0.05) and the HEIDI heterogeneity test (P_HEIDI > 0.05), providing no detectable causal signal for these gene-outcome pairs and arguing against pleiotropic confounding; *EXOSC4* lacked a blood cis-eQTL and was not testable. (E) Linkage-disequilibrium score regression (LDSC): observed-scale SNP heritability was low for both sepsis (h²=0.0022; SE, 0.0007) and ARDS (h²=0.0019; SE, 0.0006), and the genetic correlation was near zero and non-significant (rg=−0.05; SE, 0.25; P=0.83), indicating weak germline genetic architecture for both traits; under such low power, the null MR result is best read as absence of detectable causal signal rather than proof of no causal effect.
*Source: `04_MR_SMR/FigureS6_MR_validation.{pdf,png,tiff}` (panels A–B: `MR_sepsis4980_scatter.*`, `MR_sepsis4980_forest.*`; C–E generated from `MR_sepsis4980_*` result tables, `smr_*_heidi.smr`, and `ldsc_input/*.log`); script `scripts/step21_figS3_MR.py`*

### Figure S9. Hub-gene expression versus CIBERSORT immune-cell fractions

Spearman rank-correlation heatmaps between transcript levels of the six hub genes and CIBERSORT-estimated relative abundances of the 22 immune-cell subtypes in the four bulk cohorts: GSE76293 (A; ARDS discovery, purified PMN fractions), GSE69528 (B; sepsis discovery, whole blood), GSE243066 (C; ARDS validation), and GSE100159 (D; sepsis validation). Across the three whole-blood cohorts, all six hub genes correlated positively with M0 macrophage fractions (rho=0.15 to 0.72) and negatively with CD8 T-cell fractions (rho=−0.33 to −0.84) in all 18 gene-by-cohort comparisons. Grey cells denote undefined correlations for zero-variance cell types. Estimates from the PMN cohort (A) are shown for transparency only.
*Source: `07_hub_immune/cibersort/*_hub_ciber_cor.csv`; script `scripts/step42b_figS4_hub_ciber_cor.R`*

### Figure S10. ssGSEA immune-cell infiltration across the four bulk cohorts

Single-sample gene set enrichment analysis (ssGSEA) with marker gene sets of 18 immune-cell types assembled from published lineage markers, as a deconvolution-independent complement to CIBERSORT. (A–D) Per-sample ssGSEA scores (box plots with individual data points) in disease versus control groups for GSE76293 (A; ARDS, PMN), GSE69528 (B; sepsis, whole blood), GSE243066 (C; ARDS, PBMC), and GSE100159 (D; sepsis, whole blood); asterisks denote two-sided Wilcoxon rank-sum P values (*P<0.05, **P<0.01, ***P<0.001, ****P<0.0001). Th17 was not scored in GSE69528 because fewer than two marker genes were present on that platform. Monocyte and neutrophil scores were significantly higher in disease samples of the three whole-blood and PBMC cohorts, whereas CD8 T-cell, CD4 T-cell, natural-killer-cell, and B-cell scores were consistently lower. (E–H) Spearman rank correlations between ssGSEA cell-type scores and hub-gene expression in the same four cohorts (asterisks as above); hub genes correlated positively with monocyte and neutrophil scores and negatively with CD4 and CD8 T-cell scores, mirroring the CIBERSORT-based pattern (Supplementary Figure S9). Consistent with the CIBERSORT caveat, ssGSEA estimates from the purified PMN cohort (GSE76293) are shown for transparency only. ssGSEA scores quantify marker-program activity; elevated neutrophil-program scores in the PBMC cohort reflect myeloid-program activation, not neutrophil counts. Full statistics are provided in Supplementary Table S3.
*Source: `05_supplementary/Supplementary_Figures_Tables/Figures/FigureS8_ssgsea_immune_infiltration.{pdf,png,tiff}`; scripts `scripts/step47_ssgsea_immune_4cohorts.R`, `step47b_figS15_assembly.py`*

### Figure S11. Single-cell atlas of the GSE167363 CD14+ monocyte cohort

(A) Uniform manifold approximation and projection (UMAP) of 12,359 CD14+ monocytes from two healthy controls and two patients with sepsis (one sampled at two time points), after Harmony integration, annotated into seven subpopulations: classical (n=3,238), inflammatory (n=6,595), intermediate (n=341), interferon-stimulated gene-high (n=28), and non-classical monocytes (n=134), monocyte-derived dendritic cells (n=131), and plasmacytoid dendritic cells (n=1,892). (B) Subpopulation proportions by condition (healthy control versus sepsis). (C) Dot plot of lineage and subpopulation marker genes used for cluster-level annotation, including the inflammatory (*IL1B*, *CCL3*, *SAA1*, *TNF*) and interferon-stimulated gene (*ISG15*, *IFIT1*, *MX1*) module scores. (D–F) Expression of the six hub genes shown as dot plot (D), feature plots on the UMAP (E), and grouped violin plots (healthy control versus sepsis; F). *UPP1* was expressed in 42.7% of CD14+ monocytes, whereas *CDK5RAP2* (5.8%), *EXOSC4* (3.0%), *VNN1* (1.1%), *GRB10* (0.6%), and *CD177* (0.4%) showed progressively lower detection rates, which determined knockout feasibility.
*Source: `05_supplementary/Supplementary_Figures_Tables/Figures/FigureS9_single_cell_atlas_GSE167363.{pdf,png,tiff}` (panels: `08_singlecell/GSE167363/GSE167363_UMAP.*`, `GSE167363_celltype_prop.*`, `GSE167363_marker_dotplot.*`, `GSE167363_hub_dotplot.*`, `GSE167363_hub_featureplot.*`, `GSE167363_hub_violin_group.*`), `GSE167363_hub_expression_pct.csv`; script `scripts/step25_assemble_supplementary.py`*

### Figure S12. CellChat pathway-level comparison

Ranked signaling-pathway comparison between ARDS and sepsis in PBMCs (GSE151263), and between sepsis and healthy controls in CD14+ monocytes (GSE167363; 638 vs. 454 ligand–receptor interactions, with PECAM2, CD80, GRN, NECTIN, PARs, and interferon-II increased in sepsis and collagen, TGF-β, IL-1, and MHC-II dominant in healthy controls), plus complete differential ligand–receptor interaction lists. (A, B) Relative information flow of signaling pathways in GSE151263 (ARDS versus sepsis; A) and GSE167363 (healthy control versus sepsis; B). (C, D) Differential ligand–receptor interaction bubble plots in GSE151263 (C) and GSE167363 (D).
*Source: `05_supplementary/Supplementary_Figures_Tables/Figures/FigureS10_CellChat_pathway_comparison.{pdf,png,tiff}` (panels: `08_singlecell/cellchat/rankNet_pathways.png`, `diff_LR_bubble.png`, `08_singlecell/GSE167363/cellchat/rankNet_pathways.*`, `diff_LR_bubble.*`), interaction lists `Supplementary_Figures_Tables/Tables/CellChat_LR_*.csv`; script `scripts/step25_assemble_supplementary.py`*

### Figure S13. Pseudotime, AUCell, and regulatory analysis of the GSE167363 CD14+ monocyte cohort

(A) slingshot pseudotime trajectories of the monocyte compartment (four lineages rooted in classical monocytes), shown on UMAP and principal-component space, with pseudotime-inferred gene-expression trends of the six hub genes. *UPP1* declined markedly along differentiation (Spearman rho=−0.519; P<1 × 10⁻³⁰⁰); the other five hub genes showed weak but significant associations (|rho| ≤ 0.10, all false-discovery rate < 1 × 10⁻⁵). (B) AUCell activity of the six hub-gene co-expression programs (each defined as the hub gene plus its top-50 positively correlated genes within the single-cell matrix) and the composite HubSignature, shown as violin plots per subpopulation, heatmap of mean scores, and UMAP embeddings; 32 of 42 gene-set-by-cell-type comparisons were significant (Benjamini–Hochberg-adjusted P<0.05, sepsis versus healthy controls), and the *UPP1* program was most active in classical monocytes. (C) CellOracle transcription-factor knockout effects on the four knockable hub genes: heatmap of simulated expression changes (delta of simulated versus observed expression) across the top-ranked upstream transcription factors, with per-gene effect plots; CEBPB knockout produced the strongest effect on *UPP1* (delta=−0.912, reducing simulated expression from 1.26 to 0.35), STAT1 on *CDK5RAP2*, ETS2 on *EXOSC4*, and ENO1 on *VNN1*.
*Source: `05_supplementary/Supplementary_Figures_Tables/Figures/FigureS11_pseudotime_AUCell_CellOracle.{pdf,png,tiff}` (panels: `08_singlecell/GSE167363/slingshot/`, `08_singlecell/GSE167363/aucell/`, `08_singlecell/GSE167363/celloracle/`); script `scripts/step25_assemble_supplementary.py`*

### Figure S14. Monocle3 trajectory analysis of both single-cell cohorts

(A–C) GSE151263 (PBMCs; ARDS versus sepsis): (A) Monocle3 trajectories colored by cell-type annotation (one root node per graph partition, rooted in monocytes); (B) pseudotime by disease group (ARDS median 8.08 versus sepsis 7.19; P=0.055, two-sided Wilcoxon rank-sum test), showing a marginal further shift of ARDS over sepsis; (C) Moran's I of the six hub genes on the transcriptomic manifold (*UPP1* highest at 0.352; other hub genes 0.011 to 0.038, all q < 0.001; *CD177* not expressed). (D–F) GSE167363 (CD14+ monocytes; sepsis versus healthy controls): (D) trajectories colored by subpopulation (rooted in classical monocytes); (E) pseudotime by group (sepsis median 3.07 versus healthy controls 2.43; P=9.6 × 10⁻⁸¹); (F) Moran's I of the six hub genes (*UPP1* 0.495, the highest in either cohort; *CD177* 0.320). Pseudotime scales are cohort-specific and cannot be compared quantitatively across cohorts.
*Source: `05_supplementary/Supplementary_Figures_Tables/Figures/FigureS12_monocle3_trajectories.{pdf,png,tiff}` (panels: `08_singlecell/monocle3/GSE151263_monocle3_{trajectory_celltype,pseudotime_group,hub_moransI_bar}.*` and the GSE167363 counterparts); scripts `scripts/step3o_monocle3.R`, `scripts/step26_assemble_supplementary2.py`*

### Figure S15. Metabolic pathway activity in both single-cell cohorts

(A–C) GSE151263 (PBMCs; ARDS versus sepsis): (A) mean Hallmark metabolic module scores by cell type and group; (B) group-level pathway differences (two-sided Wilcoxon rank-sum tests; ARDS scored lower than sepsis on glycolysis, P=9.9 × 10⁻⁷²; oxidative phosphorylation, P=1.6 × 10⁻⁵²; fatty acid metabolism, P=1.5 × 10⁻⁴⁶; and xenobiotic metabolism, P=9.6 × 10⁻²⁸); (C) Spearman correlations between hub-gene expression and pathway scores. (D–F) GSE167363 (CD14+ monocytes; sepsis versus healthy controls): (D) pathway scores by subpopulation and group; (E) group-level differences (sepsis upregulated glycolysis, P=1.1 × 10⁻⁸⁵; oxidative phosphorylation, P=1.9 × 10⁻²⁰; fatty acid metabolism, P=6.4 × 10⁻³¹; and xenobiotic metabolism, P=2.4 × 10⁻⁵⁰); (F) hub-gene-pathway correlations, in which *UPP1* was the only positive correlate of metabolic activity (glycolysis, rho=0.168; xenobiotic metabolism, rho=0.343) while the other hub genes correlated negatively with oxidative phosphorylation (rho=−0.51 to −0.67). Together with (A–C), the two cohorts place metabolic activity on an inverted-U trajectory along the normal-sepsis-ARDS spectrum.
*Source: `05_supplementary/Supplementary_Figures_Tables/Figures/FigureS13_metabolic_pathways.{pdf,png,tiff}` (panels: `08_singlecell/metabolism/GSE151263_metab_{pathway_heatmap,group_difference,hub_correlation}.*` and the GSE167363 counterparts); scripts `scripts/step3n_metabolism.py`, `scripts/step26_assemble_supplementary2.py`*

### Figure S16. Donor-level sensitivity analyses of the CellOracle and pseudobulk results

(A) Per-donor CellOracle analysis in GSE167363: an independent gene regulatory network was built for each of the five CD14+ donors, and the effect of CEBPB in silico knockout on simulated UPP1 expression (delta of simulated minus observed mean per-cell expression) was recomputed per donor; the effect was negative in all five donors (−0.33 to −2.08), with the largest magnitudes in the later and more severe sepsis samples. (B) Per-donor Spearman correlation between UPP1 expression and monocyte pseudotime (four of five donors negative, rho −0.08 to −0.38; one early-sepsis donor positive). (C) Donor-level pseudobulk DESeq2 (raw counts aggregated per donor and cell type; sepsis versus healthy controls, 2 versus 2 donors) shown as log2 fold-change heatmap of the six hub genes across cell types; asterisks mark adjusted P<0.05; the induction of CD177 in classical monocytes (log2 fold change 6.68, adjusted P=0.002) is the only contrast significant at the donor level. Full values in Supplementary Table S17.
*Source: `03_figures/Figures_final/FigureS14_composite.{pdf,png,tiff}`; script `scripts/step60_donor_level_fig.R`*

---

### Figure S17. scTenifoldKnk virtual-knockout results

(A–F) Top-10 genes affected by virtual knockout of each hub gene, ranked by scTenifoldKnk regulatory distance and shown as log2 fold change of network edge weight; asterisks mark genes significant at false-discovery rate < 0.05. Knockouts in the CD14+ monocyte cohort (GSE167363): *UPP1* (A), *CDK5RAP2* (B), *EXOSC4* (C), and *VNN1* (D; the largest affected gene set, n=32 at FDR < 0.05). Knockouts in the PBMC cohort (GSE151263): *UPP1* (E) and *EXOSC4* (F). Ribosomal protein genes (RPL and RPS families) were the most consistently perturbed genes across all six knockouts. (G–H) Over-representation analysis (clusterProfiler, Benjamini–Hochberg-adjusted) of each knockout's top-10 affected genes: GO biological-process (G) and KEGG pathway (H) dot plots, with dot size proportional to gene count and color encoding −log10 adjusted P value. Cytoplasmic translation (GO; adjusted P down to 1.5 × 10⁻¹⁸) and the ribosome pathway (KEGG; adjusted P down to 9.1 × 10⁻¹⁷) were significantly enriched for every knockout, confirming that hub-gene knockout predominantly rewires translational machinery; *VNN1* knockout additionally affected actin-cytoskeleton and oxidative-stress–related genes (TMSB4X, GSTO1, PFN1, PRDX5).
*Source: `08_singlecell/GSE167363/knockout/FigureS15_KO_top10_enrichment.{pdf,png,tiff}`, `KO_top10_logFC_summary.csv`, `KO_top10_GO_BP.csv`, `KO_top10_KEGG.csv`, individual panels `KO_{cohort}_{gene}_top10_logFC.*`; scripts `scripts/step22_KO_top10_enrichment.R`, `scripts/step22b_patch_composite.R`*

### Figure S18. Group-stratified virtual knockout and transcription-factor perturbation

(A–C) GSE151263 (PBMCs; ARDS versus sepsis): (A) Gene Ontology biological-process enrichment of the top-150 genes perturbed by *UPP1* knockout, per disease group; (B) top-10 perturbed genes per group (direction and distance of differential regulation); (C) CellOracle group-stratified transcription-factor knockout effects on hub-gene expression (CEBPB knockout effect on *UPP1*: Δ = −0.63 in sepsis versus Δ = −0.31 in ARDS). (D–F) GSE167363 (CD14+ monocytes; sepsis versus healthy controls): (D) GO biological-process enrichment of *CDK5RAP2*-knockout-perturbed genes per group; (E) top-10 perturbed genes per group; (F) CellOracle group-stratified effects (CEBPB knockout effect on *UPP1*: Δ = −1.16 in sepsis versus Δ = −0.36 in healthy controls). Knockout fold-change profiles correlated at rho=0.66 (*UPP1*) and 0.63 (*EXOSC4*) between ARDS and sepsis, but only 0.50 (*UPP1*), 0.48 (*EXOSC4*), and 0.44 (*CDK5RAP2*) between sepsis and healthy controls (all P<1 × 10⁻¹⁵⁰), indicating that the two disease states share a more similar perturbed network architecture than either shares with health. *VNN1* knockout had the largest number of significantly affected genes (n=32 at false-discovery rate < 0.05); the *VNN1* program in healthy controls was not amenable to knockout (expressed in 0.3% of cells).
*Source: `05_supplementary/Supplementary_Figures_Tables/Figures/FigureS16_stratified_knockout.{pdf,png,tiff}` (panels: `08_singlecell/knockout_stratified/GSE151263_KO_stratified_{GO_BP,top10}_UPP1.*`, `GSE167363_KO_stratified_{GO_BP,top10}_CDK5RAP2.*`, `08_singlecell/celloracle_stratified/GSE{151263,167363}_CellOracle_stratified_hubdelta.*`); scripts `scripts/step3n_ko_stratified.R`, `step3n_ko_stratified_patch.R`, `step3n_ko_stratified_compare.R`, `step3n_celloracle_stratified.py`, `step26_assemble_supplementary2.py`*

### Figure S19. HPA protein-level validation of VNN2

VNN2 differential expression across ten infectious diseases versus healthy controls (Human Protein Atlas Blood Atlas, Olink proximity-extension assay), with the hub gene CD177 shown for comparison.
*Source: `05_supplementary/Supplementary_Figures_Tables/Figures/FigureS17_HPA_validation.{pdf,png,tiff}` (`10_HPA_validation/HPA_CD177_VNN2_infection_logFC.*`)*

### Figure S20. Six-hub-gene prognostic performance in the independent septic shock cohort GSE65682

479 patients with septic shock (MARS cohort; 114 deaths within 28 days). (A) Kaplan-Meier curves for patients grouped at the median of the six-hub-gene risk score (28-day survival 0.83 low-risk versus 0.70 high-risk; log-rank P=0.002; shaded bands, 95% confidence bands). (B) Kaplan-Meier curves for each hub gene individually (median split; log-rank P shown per panel); none of the six genes stratified mortality significantly. (C) Time-dependent receiver-operating-characteristic curves of the risk score at 7, 14, and 28 days (inverse-probability-of-censoring weighting; the 28-day evaluation was placed at 27.99 days to avoid the administrative censoring boundary). (D) 28-day calibration of the six-gene Cox model (rms calibrate, hare method, 400 bootstrap resamples); the bias-corrected curve follows the diagonal. (E) Multivariable Cox model adjusted for age and sex: risk-score hazard ratio 2.77 per unit (95% CI 1.57 to 4.90; P=0.0004), age hazard ratio 1.02 per year (P=0.004), and sex not significant (P=0.327).
*Source: `13_survival_prognosis/figures/FigureS20_sixgene_prognosis_GSE65682.{pdf,png,tiff}`; scripts `scripts/step105_survival_prognosis_GSE65682.R`, `scripts/step107_sixgene_prognosis_figure.R`*

### Figure S21. Extended prognostic modeling over the 179 shared candidate genes in GSE65682

(A) Ten-fold cross-validated partial-likelihood deviance of the LASSO Cox model over the 171 shared candidate genes represented on the platform; the penalty at minimum deviance (lambda.min, 15 genes retained) and the one-standard-error penalty (lambda.1se, no gene retained) are marked. (B) Coefficients of the 15 genes retained at lambda.min (none of the six hub genes was retained). (C) Time-dependent receiver-operating-characteristic curves at 28 days for the clinical model (age, sex, Mars endotype), the six-hub-gene model, their combination, the LASSO signature alone and with clinical covariates, and the random survival forest (out-of-fold predictions). (D) Five-group calibration of predicted versus Kaplan-Meier-observed 28-day mortality for the same models. (E) Decision-curve analysis (net benefit versus threshold probability, 0.05 to 0.60). (F) Incremental value over the clinical model: change in bootstrap-optimism-corrected C index, net reclassification index, and integrated discrimination improvement (x100) at 28 days. (G) Corrected C-index comparison including the random survival forest (out-of-bag). (H) Random-survival-forest variable importance (impurity) for the 12 most important predictors; age dominates, and no hub gene appears.
*Source: `13_survival_prognosis/figures/FigureS21_extended_prognostic_modeling.{pdf,png,tiff}`; script `scripts/step106_extended_modeling_GSE65682.R`*

## File Manifest

| File | Description |
|------|-------------|
| `Supplementary_Methods.md` | Complete methods relocated from the main text by the Path-A restructure of 2026-10-16 (S1-S14; S1-S7 follow first-citation order; citations converted to S-numbered references on 2026-10-19) |
| `Supplementary_Materials.md` | This index (Tables S1-S26 and Figure legends S1-S21) |
| `figures/supplementary/FigureS1-S21.*` | Supplementary figures (renumbered), PDF/PNG/TIFF triple export (PDF+PNG+TIFF verified complete on 2026-10-19, including S2, S4, S5, and S17) |
| `figures/supplementary_tables/TableS1-S26.*` | Supplementary tables (renumbered; S22-S26 as CSV; S24-S26 assembled 2026-10-19) |
| `figures/supplementary_tables/CellChat_LR_*.csv` | Complete differential ligand-receptor interaction lists |
| `figures/supplementary_tables/TableS24_prognostic_statistical_checks.csv` | Cox proportional-hazards diagnostics, Kaplan-Meier cut-point sensitivity, and multiplicity coverage of the immune-cell comparisons |
| `figures/supplementary_tables/TableS25_CellChat_pathway_strength_recalc.csv` | Recalculated pathway-level CellChat communication strength for both single-cell cohorts (GSE151263 sepsis versus ARDS; GSE167363 sepsis versus healthy control) |
| `figures/supplementary_tables/TableS26_cohort_metadata.csv` | Cohort metadata; ancestry and sampling time point for the four bulk cohorts flagged `[AUTHOR_INPUT_NEEDED]` |
