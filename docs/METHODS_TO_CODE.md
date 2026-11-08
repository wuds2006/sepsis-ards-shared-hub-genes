# Supplementary Methods → 代码对照

S1–S14 的编号与 `Supplementary Material.docx` 一致。表中“输出”指脚本落盘的
中间结果（多数在 `.gitignore` 排除的模块目录里，成品表另存于 `data/frozen/`），
“图”指该步结果进入的图号（v3 号系）。

工作目录为仓库根 `Code/`；每行脚本按编号顺序执行。

## S1 数据获取与样本清洗

| 脚本 | 内容 | 输出 |
|---|---|---|
| `shell/01_fetch_geo_series.sh` | 下载 7 个 GEO 系列的 series matrix / RAW.tar | `00_rawdata/` |
| `shell/02_fetch_annotation_and_gwas.sh` | 平台注释（GPL570/GPL10558/GPL6884）、OpenGWAS 与 FinnGen 汇总统计 | `00_rawdata/annotation/`、`04_MR_SMR/` |
| `python/01_parse_geo_matrix.py` | 系列矩阵解析、样本剔除（GSE76293 去 70 例体外培养、GSE69528 去 27 例糖尿病对照） | 各队列表达矩阵与分组表 |

剔除规则写在脚本顶部的注释里，与 Table 1 的样本数一致（12v12、83v28、33v12、34v15）。

## S2 机器学习与 SHAP

| 脚本 | 内容 | 图 |
|---|---|---|
| `R/20_ml_ards.R` | 11 种分类器（LR/LASSO/Ridge/ENet/RF/SVM/XGBoost/GBM/DT/KNN/NB/NN）在 GSE76293 训练，5 折 CV + SHAP | F2、F3、S5、S22 |
| `R/21_ml_sepsis.R` | 同套流程在 GSE69528 | F3、S5、S22 |
| `R/22_ml_stability.R` | 重复抽样稳定性、标签置换与过拟合诊断 | S5 |

一致性规则（LASSO ∩ RF ∩ SHAP，允许不完全一致）在 `R/20`、`R/21` 的
`consensus` 段落；`R/22` 给出置换零分布。

## S3 免疫细胞解卷积

| 脚本 | 内容 | 图 |
|---|---|---|
| `R/35_cibersort_ssgsea.R` | CIBERSORT（e1071 nu-SVR + preprocessCore 分位标准化）22 类分数、组间 Wilcoxon、22×22 Spearman；六基因 ssGSEA | F8、S9、S24 |
| `R/39_ssgsea_4cohorts.R` | 18 套 lineage marker 的 ssGSEA，四队列 | F8、S10 |
| `R/36_gsva_immune.R` | GSVA 免疫通路，top-20 差异通路与共享通路 | F8、S15 |

M0 分数的解读（髓系富集而非巨噬细胞绝对量）写在 `R/35` 的输出注释与图注中。

## S4 单细胞处理

| 脚本 | 内容 | 图 |
|---|---|---|
| `python/10_scanpy_cluster.py` | GSE151263 PBMC：QC、harmony 整合、leiden 聚类、UMAP | F9、S11、S13 |
| `R/60_annotate_gse167363.R` | GSE167363 CD14+ 单核：聚类注释与标志基因打分 | F9、S11 |
| `python/11_cellchat_prepare.py` | 导出 CellChat 输入（按组、按亚群） | — |
| `R/61_cellchat_pbmc.R`、`R/62_cellchat_netp_pbmc.R` | PBMC 通讯网络与 netP 重算 | F9、S12 |
| `R/63_cellchat_monocytes.R`、`R/64_cellchat_netp_monocytes.R` | 单核细胞通讯与 netP 重算 | F9、S12 |
| `R/65_monocle3.R` | monocle3 拟时序、graph_test、Moran's I（含 R 4.6 下 `my.moran.test` 的热修复） | F9、S13、S14 |
| `R/66_aucell.R` | 六基因 AUCell 打分 | F9、S13 |
| `R/70_donor_level.R` | donor 水平 pseudobulk（DESeq2）与一致性核查 | S2、S16 |
| `python/16_metabolism.py` | 代谢通路打分（scMetabolism 风格） | S15 |
| `python/17_cross_dataset.py` | 两套单细胞数据的跨数据集一致性 | S16 |
| `R/figures/46_supp_S26.R` | 单核细胞 ETosis（19 基因）与 NOX（17 基因）程序 AUCell 打分：亚群间、donor 配对与 ARDS/sepsis 边界的汇总图（读冻结表） | S26、Table S28 |
| `R/figures/47_supp_S27.R` | 单核细胞 hdWGCNA 共表达模块：soft power 4（R² = 0.84）、2 639 个 metacell、模块-性状与模块-ETosis 相关 | S27 |

## S5 虚拟敲除

| 脚本 | 内容 | 图 |
|---|---|---|
| `R/67_sctenifold_six_hub.R` | scTenifoldKnk，六个 hub 基因 | S17 |
| `R/68_sctenifold_monocytes.R` | GSE167363 内逐基因敲除 | S17 |
| `R/69_sctenifold_stratified.R` | HC/sepsis/ARDS 分层敲除与谱相关性 | S18 |
| `python/12_celloracle_pbmc.py`、`python/13_celloracle_monocytes.py` | CellOracle 基础 GRN 与 in silico KO | F9、S13 |
| `python/14_celloracle_stratified.py`、`python/15_celloracle_per_donor.py` | 分层与逐 donor 的 CEBPB→UPP1 依赖 | F9、S13 |
| `data/frozen/08_singlecell/hdwgcna/KO_vs_module_fisher.csv`、`KO_vs_etosis_fisher.csv` | 虚拟敲除扰动谱与 hdWGCNA 模块、ETosis 基因集的 Fisher 交叠（四个可敲除 hub 富集于核糖体生成模块、避开自身模块） | Table S27 |

## S6 外部验证

| 脚本 | 内容 | 图 |
|---|---|---|
| `R/23_validation_ards.R` | GSE243066 迁移模型 | F5、F6 |
| `R/24_validation_sepsis_candidates.R`、`R/25_validation_ards_candidates.R` | 候选基因逐一在验证队列的表现 | F3、F6 |
| `R/26_validation_sepsis_roc.R`、`R/27_validation_auc_ci.R` | 单基因 ROC 与 DeLong CI | F6、Table 2 |
| `R/28_validation_diagnostics.R` | 批次/平台诊断、ComBat 敏感性 | S7 |
| `R/29_muc1_sensitivity.R` | MUC1 被否决的敏感性分析 | F3 |
| `R/30_cd177_ablation.R` | 去 CD177 后的迁移表现 | Table S11 |

Hedges g 合并由 `R/27` 输出（`metaDE_single_gene_validation.csv`）：
FE 合并 g 1.14–2.32，CD177/UPP1/VNN1 的 I² 0.60–0.91。

## S7 孟德尔随机化

| 脚本 | 内容 | 图 |
|---|---|---|
| `shell/03_prepare_mr_inputs.sh` | 暴露/结局汇总统计的获取与格式统一 | `04_MR_SMR/` |
| `python/02_mr_extract_forward.py`、`python/07_mr_extract_reverse.py` | 正向 / 反向工具变量抽取 | 同上 |
| `python/03_convert_gwas_sumstats.py`、`python/04_match_ards_against_sepsis.py` | 等位基因对齐、SNP 匹配 | 同上 |
| `R/50_mr_forward.R` | IVW + MR-Egger + 加权中位数等敏感性分析 | F7 |
| `R/51_mr_reverse.R` | 反向 MR | F7 |
| `R/52_mr_power.R` | 事后功效与最小可检测 OR | F7、Table S13 |
| `python/05_align_coloc_inputs.py`、`R/53_colocalization.R` | eQTL×GWAS 共定位 | F7 |
| `R/54_smr_heidi.R`、`shell/04_run_smr_heidi.sh` | SMR + HEIDI（外部二进制运行手册） | F7 |
| `python/06_prepare_ldsc_inputs.py` | LDSC 输入（求和统计、LD score 路径） | Table S14 |

正向 β=−0.336（95% CI −0.908–0.236，P=0.249），反向 P=0.329。

## S8 芯片预处理与探针映射

| 脚本 | 内容 |
|---|---|
| `python/01_parse_geo_matrix.py` | 系列矩阵读取、`///` 多基因拆分、背景校正与 log2 |
| `R/10_deg.R` | 探针→基因折叠（同名取平均表达最高探针）、分位标准化、limma |
| `R/40_cohort_pca.R` | 四队列 PCA 与批次可视化（S7） |

## S9 差异表达与 WGCNA

| 脚本 | 内容 | 图 |
|---|---|---|
| `R/10_deg.R` | limma，\|log2FC\|>1 且 adj.P<0.05，得到 557（ARDS）/1458（sepsis） | F1、S6 |
| `R/11_wgcna.R` | soft-threshold、模块识别、模块-性状相关，13/11 个疾病相关模块 | F1、S1 |
| `R/12_wgcna_preservation.R` | 跨队列模块保守性（Zsummary 19.5 / 21.4） | S1 |
| `R/13_intersection.R` | DEG ∩ 模块基因 → 179 个共享候选 | F1、S1 |
| `R/14_enrichment_go_kegg.R` | GO BP / KEGG 富集 | F2、S23 |
| `R/15_ppi_network.R` | STRING PPI 与网络模块 | F3 |

## S10 GSVA 免疫通路

| 脚本 | 内容 | 图 |
|---|---|---|
| `R/36_gsva_immune.R` | GSVA，top-20 通路，共享通路随阈值（20/30/50）变化 | F8、S15 |
| `R/37_hub6_gsva.R` | 六基因与 GSVA 通路的关联 | F8 |
| `R/39_ssgsea_4cohorts.R` | ssGSEA 作为解卷积无关的交叉验证 | F8、S10 |

## S11 蛋白水平验证

| 脚本 | 内容 | 图 |
|---|---|---|
| `R/75_hpa_plasma_protein.R` | HPA Blood Atlas：CD177、VNN2 在 10 种感染性疾病中的血浆蛋白变化 | F10、S19 |
| `R/76_hpa_immune_cells.R` | 免疫细胞 RNA 表达特异性 | S19 |

## S12 药物重定位与分子对接

| 脚本 | 内容 | 图 |
|---|---|---|
| `python/30_cmap_random_null.py` | CMap 上调查询 vs 随机基因集零分布（13/100、0/100、2/100） | F10、S3 |
| `python/31_cmap_deg_matched_null.py` | 严格匹配背景（两队列共同上调基因） | S3 |
| `python/32_cmap_enrichr.py` | 多库交叉核对 | S3 |
| `python/21_docking_box_definition.py` | 结合口袋定义（文献位点 + fpocket） | — |
| `python/20_docking_screen.py` | AutoDock Vina 1.2.5 初筛：6 蛋白 × 10 天然产物，exhaustiveness=32 | F10、S4 |
| `python/22_redock_upp1.py` | UPP1 固定 box 重对接 | S4 |
| `python/23_redocking_validation.py` | 共晶配体重对接 RMSD 校验 | S4 |
| `python/24_mmgbsa_rescoring.py` | OpenMM 单点 MM/GBSA（ff14SB + GBn2/OBC-II） | F10、Table S19 |
| `python/25_alphafold_plddt.py` | AlphaFold pLDDT（整模型 vs 对接 box 区域） | F10 |
| `python/26_plip_contacts.py` | PLIP 接触残基与距离 | F10 |
| `python/33_chembl_known_drugs.py`、`python/34_opentargets_druggability.py` | ChEMBL 37 与 Open Targets 可成药性 | Table 3 |

Vina 最优：VNN1–luteolin −10.01 kcal/mol；两种打分在靶点层面一致、配体排序不一致，
因此不作配体层面优先级。

Figure S25 的接触图另由 `R/figures/26_figS25_2d_contact_map.py`、
`27_figS25_anchor.py`、`28_figS25_2d_coords.py`、`29_figS25_3d_scene.py`
从对接姿态重建。这些脚本需要 pose（`VNN1_luteolin_out.pdbqt`）、配体
（`luteolin.sdf`）与受体（`7SLV.pdb`），三者均不随仓库分发，运行时用
`DOCK_ROOT` 指向本地的 `09_docking_drug/` 目录（默认即该相对路径）。

## S13 单基因 GSEA

| 脚本 | 内容 | 图 |
|---|---|---|
| `R/38_single_gene_gsea.R` | 六基因在四队列分别做单基因 GSEA（Hallmark），NES 与 FDR | F4、S6 |

## S14 GSE65682 预后建模

| 脚本 | 内容 | 图 |
|---|---|---|
| `R/45_prognosis_cox.R` | 逐基因 Cox、联合六基因 Cox、bootstrap 乐观校正 | S20 |
| `R/46_prognosis_pergene_km.R` | 逐基因 KM | S20 |
| `R/47_prognosis_extended_models.R` | LASSO-Cox（1SE 与 min-deviance）、RSF、临床模型、LR 检验、NRI/IDI、DCA | S20、S21 |
| `R/48_prognosis_figure.R` | 预后图形汇总 | S20、S21 |
| `R/49_prognosis_statistics.R` | 三个统计口径的复核（C-index、AUC、LR P） | S20、S21 |

UPP1 是联合模型中唯一显著项（HR 0.64，95% CI 0.44–0.95，P=0.026）；
六基因评分 log-rank P=0.002，但相对临床模型无增量（LR P=0.16）。

## 绘图（不属于 S1–S14）

| 脚本 | 图 |
|---|---|
| `R/00_setup.R` | 版面、字体、配色、四格式导出、`out_dir()` |
| `R/figures/01_prep_tables.R` | 六基因表达长表（`data/derived/`） |
| `R/figures/01_fig1.R` … `10_fig10.R` | 主图 1–10 |
| `R/figures/20_supp_S1.R` … `44_supp_S24.R` | 补充图 S1–S24 |
| `R/figures/99_check_figures.R` | 四格式与字体字号、物理尺寸校验 |

校验项：四种格式齐备；SVG 只用 Helvetica；字号只出现 7 与 9 pt；SVG 内
`<text>` 节点 ≥10 个（保证文字未被转曲）；PNG 宽 180 mm、300 dpi；
PDF MediaBox 高度 200–220 mm。
