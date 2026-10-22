#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 103 —— 构建可公共归档的分析代码仓库（GitHub + Zenodo DOI）
按 02_Manuscript_Revised.md METHODS 的 14 个小节组织最终分析脚本；
剔除图片拼装/文稿工具类脚本（step15-102 中的 figure/manuscript tooling）；
被取代的早期版本在 README 中注明，不拷贝。
"""
import os, shutil

R    = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
SRC  = f"{R}/scripts"
REPO = f"{R}/AJRCCM_submission/10_analysis_code"

# METHODS 小节 -> (目标子文件夹, [脚本], 一句话说明)
MAP = [
    ("00_setup", ["install_R_packages.R", "common.R", "ajrccm_theme.R"],
     "R 包安装与共享主题/工具函数"),
    ("01_data_download", ["download_data.sh", "download_data2.sh", "parse_geo_matrix.py"],
     "GEO/GWAS/eQTLGen 原始数据下载与 series_matrix 探针级解析（含脱靶样本识别）"),
    ("02_preprocessing_DEG", ["step1a_DEG.R"],
     "METHODS-Data Preprocessing and Differential Expression Analysis：背景校正、分位数标准化、探针合并、limma 差异表达"),
    ("03_WGCNA", ["step1b_WGCNA.R", "step55_wgcna_preservation.R"],
     "METHODS-WGCNA：签名网络构建、软阈值选择、模块-性状相关；跨队列 modulePreservation Zsummary"),
    ("04_shared_genes_enrichment", ["step1c_intersection.R", "step2h_go_kegg_complete.R",
                                    "step2g_go_sankey_chord.R", "step2f_sankey_chord.R",
                                    "step5_network_ppi.R"],
     "METHODS-Identification and Enrichment Analysis：DEG∩WGCNA 交集、clusterProfiler GO(BP/MF/CC)+KEGG、桑基气泡图与弦图、共享基因网络注释"),
    ("05_machine_learning", ["step2a_ML.R", "step2a2_ML_sepsis.R",
                             "step50_ml_stability_overfit.R", "step27_reviewer_fix_R.R"],
     "METHODS-Machine-Learning Prioritization：11 分类器+SHAP（ARDS/Sepsis 分别训练）、选择稳定性（LASSO bootstrap/RF 50 seeds/阈值扰动 Jaccard/500 次标签置换）、发现集内跨疾病候选检验 + kME/共表达度中心性"),
    ("06_expression_GSEA", ["step2d_expression_ROC.R", "step2d2_boxplot_ROC.R",
                            "step46_single_gene_gsea_4cohorts.R"],
     "METHODS-Expression Patterns and GSEA：6 hub 基因四数据集表达比较 + 单基因 GSEA（fgsea, MSigDB Hallmark v2024.1, Pearson 相关排序）"),
    ("07_prediction_model", ["step2e_nomogram.R", "step49_cd177_ablation.R",
                             "step75_muc1_sensitivity.R"],
     "METHODS-Disease Prediction Model：rms/lrm（ridge penalty=1）列线图 + ROC + 校准 + DCA；CD177 消融五基因迁移；MUC1 七基因敏感性"),
    ("08_external_validation", ["step2b_GSE100159_validate.R", "step2b_validation_ARDS.R",
                                "step2b2_GSE100159_candidates.R", "step2b3_GSE243066_sepsis_candidates.R",
                                "step2b4_GSE100159_ROC.R", "step2b_AUC_CI.R",
                                "step52_validation_diagnostics.R", "step59_cohort_pca.R",
                                "step27_reviewer_fix_py.py"],
     "METHODS-External Validation and Diagnostic Performance：两独立队列单基因 ROC/AUC + DeLong CI + 精确单侧下限、OVL/Cohen d/Hedges g 合并、ComBat 探索性重评估、队列 PCA 批量评估、GSE151263 患者级 pseudobulk 疾病-疾病 ROC + 供体分层敏感性"),
    ("09_Mendelian_randomization", ["step1d_MR_4980.R", "step48a_extract_reverse_MR.py",
                                    "step48b_reverse_MR.R", "step51_mr_power.R",
                                    "step1d_heidi.R", "step1d_coloc.R", "align_coloc.py",
                                    "prepare_ldsc.py", "convert_gwas_ma.py",
                                    "extract_MR_4980.py", "match_ards_4980.py"],
     "METHODS-Mendelian Randomization：双向 TwoSampleMR（ieu-b-4980 ↔ FinnGen J10_ARDS, P<5e-6, 1Mb clump, IVW+Egger/WM/模式+Steiger）、Brion/mRnd 功效、SMR-HEIDI（eQTLGen）、coloc 共定位、LDSC 数据准备"),
    ("10_immune_infiltration", ["step41_cibersort_ssgsea.R", "step47_ssgsea_immune_4cohorts.R",
                                "step45_hub6_gsva.R", "step2c2_GSVA_immune.R",
                                "step2c3_GSVA_immune_top20.R", "step42b_figS4_hub_ciber_cor.R",
                                "step2c_hub_immune.R"],
     "METHODS-Immune Infiltration：CIBERSORT（LM22, nu-SVR, 1000 置换, 四队列）、18 细胞 ssGSEA（四队列 + hub 基因相关）、六基因 GSVA 签名评分、免疫 GO-BP GSVA top20（limma 差异）"),
    ("11_single_cell/01_processing", ["step3a_scanpy.py", "step3n_seurat_integration.R",
                                      "step3e_GSE167363_annotate.R", "step3e_fix_annotate.R",
                                      "step3e_fix2_annotate.R", "step3e_fix3_annotate.R"],
     "METHODS-Single-Cell：GSE151263 Scanpy+Harmony+Leiden；GSE167363 Seurat+Harmony 与标志物模块注释（含三轮修正链）、RPCA 整合对照"),
    ("11_single_cell/02_cellchat", ["step3b_prepare_cellchat.py", "step3b_cellchat.R",
                                    "step3b_cellchat_compare.R", "step3f_GSE167363_cellchat.R"],
     "METHODS-Single-Cell：两队列 CellChat 细胞通讯与信息流排序"),
    ("11_single_cell/03_trajectory", ["step3c_prepare_slingshot.py", "step3c_slingshot.R",
                                      "step3h_GSE167363_slingshot.R", "step3h_finish_slingshot.R",
                                      "step3o_monocle3.R"],
     "METHODS-Single-Cell：slingshot + Monocle3 双拟时序、Moran's I 空间自相关"),
    ("11_single_cell/04_aucell_metabolism", ["step3g_aucell.R", "step3n_metabolism.py",
                                             "step3l_cross_dataset_comparison.py"],
     "METHODS-Single-Cell：AUCell 共表达程序活性（基因+top50 正相关）、7 条 Hallmark 代谢通路评分、跨数据集六维度比较"),
    ("11_single_cell/05_virtual_knockout", ["step3d_scTenifoldKnk.R", "step3d_EXOSC4_ko.R",
                                            "step3i_GSE167363_scTenifoldKnk.R", "step3i_extra_knockout.R",
                                            "step3d_CellOracle.py", "step3d_CellOracle_local.py",
                                            "step3j_GSE167363_CellOracle.py",
                                            "step3n_ko_stratified.R", "step3n_ko_stratified_patch.R",
                                            "step3n_ko_stratified_compare.R", "step3n_celloracle_stratified.py",
                                            "step22_KO_top10_enrichment.R"],
     "METHODS-Single-Cell：scTenifoldKnk + CellOracle 双虚拟敲除（含分组分层、低表达补敲、KO top10 富集）"),
    ("11_single_cell/06_donor_level", ["step53_celloracle_per_donor.py", "step54_donor_level.R"],
     "METHODS-Single-Cell：供体级 CellOracle（每供体独立 GRN）、供体级 pseudobulk DESeq2 与轨迹相关"),
    ("12_protein_validation", ["step_hpa_validation.R", "step_hpa_immune_cell.R"],
     "METHODS-Protein-Level Validation：HPA Blood Atlas（Olink PEA）CD177/VNN2 泛感染性疾病蛋白差异、HPA 免疫细胞 RNA 特异性"),
    ("13_drug_repurposing_docking", ["step67_cmap_negative_control.py", "step76_cmap_deg_null.py",
                                     "step4_docking_site.py", "step4_dock_one.py",
                                     "step56_redocking.py", "step57_upp1_redock_fixed_box.py",
                                     "step66_alphafold_plddt.py", "step69_plip_vnn1_luteolin.py",
                                     "step70_mmgb_sasa.py", "step68_opentargets_druggability.py"],
     "METHODS-In Silico Drug Repurposing and Docking：Enrichr Old CMAP/DSigDB/GEO 扰动库查询 + 三种背景随机零分布、AutoDock Vina 定点对接（10 配体）、共结晶回对接验证、UPP1 box 修正、AlphaFold pLDDT、PLIP 相互作用、OpenMM MM/GBSA、Open Targets/ChEMBL 可成药性"),
    ("90_supplementary_tables", ["step26_build_tables.py", "step62_build_tables_S13_S17.py",
                                 "step72_build_tableS21.py"],
     "补充表 S1-S21 的数据装配（派生数据，非分析本体）"),
]

n_copied = 0
missing = []
for sub, files, desc in MAP:
    dst_dir = os.path.join(REPO, sub)
    os.makedirs(dst_dir, exist_ok=True)
    for f in files:
        s = os.path.join(SRC, f)
        if not os.path.exists(s):
            missing.append(f); continue
        shutil.copy2(s, os.path.join(dst_dir, f))
        n_copied += 1
    # 每文件夹一个简短 README
    with open(os.path.join(dst_dir, "README.md"), "w", encoding="utf-8") as fh:
        fh.write(f"# {sub}\n\n{desc}\n\n## Scripts\n\n")
        for f in files:
            fh.write(f"- `{f}`\n")

print(f"copied {n_copied} scripts into {REPO}")
if missing:
    print("MISSING:", missing)
