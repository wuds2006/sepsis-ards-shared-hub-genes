# 数据来源与版本

所有输入均为公共资源。原始矩阵不随仓库分发，可用 `shell/01_…` 与
`shell/02_…` 重新获取；下表列出访问号、用途与版本。

## 表达谱与单细胞

| 数据集 | 平台 / 类型 | 纳入样本 | 用途 |
|---|---|---|---|
| GSE76293 | GPL570 芯片（中性粒细胞） | 12 ARDS / 12 对照（剔除 70 例体外培养） | ARDS 发现队列 |
| GSE69528 | GPL10558 芯片（全血） | 83 脓毒症 / 28 对照（剔除 27 例 2 型糖尿病未感染） | 脓毒症发现队列 |
| GSE100159 | GPL6884 芯片（全血） | 33 脓毒症 / 12 对照 | 脓毒症外部验证 |
| GSE243066 | RNA-seq（全血） | 34 ARDS / 15 对照 | ARDS 外部验证 |
| GSE151263 | scRNA-seq PBMC | 13,996 个细胞（3 ARDS + 4 脓毒症供者） | 单细胞通讯、拟时序 |
| GSE167363 | scRNA-seq CD14+ 单核 | 12,359 个细胞 | 单细胞定位、虚拟敲除 |
| GSE65682 | GPL13667 芯片（感染性休克） | 479 例，28 天死亡 114 例 | 预后建模 |

## 遗传学

| 资源 | 版本 / 访问号 | 用途 |
|---|---|---|
| IEU OpenGWAS | `ieu-b-4980`（GWAS Catalog GCST90281174），脓毒症 | 正向 MR 暴露 |
| FinnGen | R11，`J10_ARDS` | MR 结局 |
| eQTLGen | 2019-12-11 发布，cis-eQTL（FDR<0.05） | SMR、HEIDI、共定位 |
| 1000 Genomes | Phase 3，EUR | LD 参考（PLINK `--make-bld`） |

## 蛋白、药物与结构

| 资源 | 版本 / 访问日期 | 用途 |
|---|---|---|
| Human Protein Atlas Blood Atlas | v24，2025-11-06 | CD177 / VNN2 血浆蛋白 |
| CMap（CLUE） | 2026-09 快照 | 药物重定位与零分布 |
| ChEMBL | release 37，2026-09-17 | 已知化合物与临床阶段 |
| Open Targets | 2026-09 | 可成药性 |
| PDB | 3EUF（UPP1）、7SLV（VNN1） | 对接受体 |
| AlphaFold DB | 与 UniProt 对应的预测模型 | pLDDT 与缺失残基核查 |

## 仓库内的 `data/frozen/`

13 个模块的落盘结果，共约 400 个文件、90 MB。绘图脚本只读这一层。

| 目录 | 内容 |
|---|---|
| `01_DEG` | 差异表达全表、折叠后的基因表达矩阵、分组表 |
| `02_WGCNA` | 模块基因、模块-性状相关、跨队列保守性 |
| `03_intersection` | 179 共享候选、GO/KEGG 富集表 |
| `04_MR_SMR` | 工具变量、IVW 与敏感性分析、SMR/HEIDI、共定位、LDSC |
| `05_ML_SHAP` | 11 分类器 AUC、SHAP 重要性、种子稳定性、置换零分布 |
| `06_external_validation` | 单基因 AUC 与 CI、迁移模型、Hedges g 合并（`metaDE_single_gene_validation.csv`） |
| `07_hub_immune` | CIBERSORT 22 类分数、组间检验、ssGSEA、GSVA |
| `08_singlecell` | 聚类注释、CellChat 通路强度、拟时序、Moran's I、虚拟敲除、代谢、CellOracle |
| `09_docking_drug` | Vina 打分、重对接 RMSD、MM/GBSA、PLIP 接触 |
| `10_HPA_validation` | 血浆蛋白与免疫细胞 RNA 表达 |
| `11_network_pharma` | STRING PPI、节点/边表、功能富集 |
| `12_single_gene_gsea` | 四队列单基因 GSEA 全结果 |
| `13_survival_prognosis` | Cox、KM、校准、DCA、NRI/IDI、LASSO/RSF、时间依赖 AUC |

## 重新获取

```bash
bash shell/01_fetch_geo_series.sh      # GEO 系列矩阵
bash shell/02_fetch_annotation_and_gwas.sh   # 平台注释、OpenGWAS、FinnGen
bash shell/03_prepare_mr_inputs.sh     # MR 输入整理
```

下载脚本会把文件放到 `00_rawdata/`；GEO 偶有限流，失败时按脚本里的
URL 手动下载到同一路径即可。
