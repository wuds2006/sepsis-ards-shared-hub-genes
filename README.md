# Shared hub genes and immune signatures of sepsis and ARDS — analysis code

Code and derived tables behind the manuscript *"Shared Hub Genes and Immune
Signatures of Sepsis and ARDS: A Multi-Omics Integration"* (submitted to
*AJRCCM*).

本仓库收录两段内容：一是从公共数据得到全部结果的分析脚本（R / Python / shell，
按 Supplementary Methods S1–S14 分阶段），二是把结果画成主图 1–10 与补充图
S1–S27 的绘图脚本（R，四格式导出：300 dpi PNG、矢量 PDF、LZW TIFF、
保留原生 `<text>` 的 SVG）。

---

## 1. 目录

```
Code/
├── run_all.sh                 全流程驱动（可按阶段单独跑）
├── R/
│   ├── 00_common.R            分析阶段共享：路径、配色、主题、工具函数
│   ├── 00_setup.R             绘图阶段共享：180 mm 版面、Helvetica 7/9 pt、四格式导出
│   ├── 10_…30_                S8/S9/S2/S6 预处理、差异表达、WGCNA、交集、机器学习、外部验证
│   ├── 35_…40_                S3/S10/S13 免疫浸润、GSVA、单基因 GSEA、队列 PCA
│   ├── 45_…49_                S14 GSE65682 预后建模
│   ├── 50_…54_                S7 孟德尔随机化、coloc、SMR-HEIDI
│   ├── 60_…70_                S4/S5 单细胞：注释、CellChat、monocle3、AUCell、虚拟敲除、donor 水平
│   ├── 75_…76_                S11 蛋白水平验证（HPA Blood Atlas）
│   └── figures/               主图 1–10、补充图 S1–S27、四格式校验
├── python/                    GEO 矩阵解析、GWAS 整理、scanpy/CellOracle/对接/CMap/ChEMBL
├── shell/                     公共数据下载、MR 输入准备、SMR-HEIDI 运行手册
├── data/
│   ├── frozen/                13 个模块的结果表（脚本读这些表出图）
│   └── derived/               由图件层从原始矩阵整理出的长表
├── environment/               R / Python 版本、包清单、安装脚本
└── docs/                      METHODS_TO_CODE.md、DATA_SOURCES.md
```

补充图 S26、S27 与补充表 S27、S28 对应的三项单细胞分析：

| 分析 | 绘图脚本 | 冻结结果 | 产出 |
|---|---|---|---|
| 单核细胞 ETosis 程序活性（AUCell；ETosis 19 基因集、NOX 17 基因集） | `R/figures/46_supp_S26.R` | `data/frozen/08_singlecell/etosis_aucell/` | Figure S26；Table S28 取自 `*_ET_genesets.csv` |
| 单核细胞共表达模块（hdWGCNA；2 639 个 metacell，soft power 4，R² = 0.84） | `R/figures/47_supp_S27.R` | `data/frozen/08_singlecell/hdwgcna/` | Figure S27 |
| 虚拟敲除扰动谱与模块的交叠（scTenifoldKnk 扰动谱 × Fisher 精确检验） | — | `data/frozen/08_singlecell/hdwgcna/KO_vs_*.csv` | Table S27 取自 `KO_vs_module_fisher.csv`、`KO_vs_etosis_fisher.csv` |

第三项是对已算出的扰动谱做集合交叠统计，没有独立的绘图脚本；两张 Fisher
表本身就是补充表的数值来源。

## 2. 怎么跑

工作目录固定为 `Code/`。

```bash
# 只重画图（不需要原始表达矩阵，data/frozen 已足够）
Rscript R/figures/06_fig6.R
Rscript R/figures/99_check_figures.R Figure6

# 全部主图 + 补充图 + 四格式校验
bash run_all.sh --stage 11

# 全流程（需先下载原始数据，见 docs/DATA_SOURCES.md）
bash shell/01_fetch_geo_series.sh
bash run_all.sh

# 看阶段清单
bash run_all.sh --list
```

图件输出到 `../Figures/{main,supplementary}`；用 `FIG_ROOT` 可以改到别处：

```bash
FIG_ROOT=/path/to/figures Rscript R/figures/01_fig1.R
```

单张图约 10–40 s，27 张补充图合计约 9 min（S1 的 WGCNA soft-threshold 最慢）。

## 3. 数据从哪来

原始数据全部为公共资源，未纳入本仓库（体积大、可重新下载）：GEO
GSE76293、GSE69528、GSE100159、GSE243066、GSE151263、GSE167363、GSE65682；
IEU OpenGWAS `ieu-b-4980`、FinnGen R11、eQTLGen、1000 Genomes EUR；
Human Protein Atlas Blood Atlas v24、ChEMBL 37、Open Targets、CMap。
下载与版本见 `docs/DATA_SOURCES.md`。

`data/frozen/` 是 13 个分析模块落盘的结果表（约 96 MB、421 个文件），
绘图脚本只读这里，因此**不重跑上游分析也能复现全部图件**。
`data/derived/hub6_expression_long.csv` 是从原始矩阵抽出的六基因表达长表，
已随仓库提供；只有在重建时才需要把 `ARDS_ROOT` 指到含原始矩阵的项目目录。

## 4. 复现边界

诚实说明，避免误读：

| 范围 | 状态 |
|---|---|
| 主图 1–10、补充图 S1–S27 | 由 `data/frozen` 一键重画，四格式校验通过 |
| 差异表达、WGCNA、交集、ML/SHAP、外部验证、免疫浸润、GSVA/GSEA、预后建模 | R 脚本齐全，输入为 GEO 原始矩阵，可重跑 |
| 单细胞（注释、CellChat、monocle3、AUCell、虚拟敲除、代谢、CellOracle） | R/Python 脚本齐全，需 Seurat/scanpy/CellOracle 环境与数小时至数天算力 |
| 孟德尔随机化 | 脚本齐全；SMR 1.03 与 LDSC 为外部二进制，用法见 `shell/04_run_smr_heidi.sh` |
| 分子对接 | AutoDock Vina 1.2.5 + Meeko + OpenMM；受体/配体结构从 PDB、AlphaFold 获取，本仓库不存构象 |
| CMap / ChEMBL / Open Targets | 脚本齐全，需联网查询；结果为 2026 年 9 月快照 |

上游脚本保留了当初的目录名作为输出位置（`01_DEG/`、`05_ML_SHAP/` 等），
这些目录已在 `.gitignore` 中排除，重跑不会污染仓库。

## 5. 环境

R 4.6.1、Python 3.13.12，包清单见 `environment/`。
R 包可用 `Rscript environment/install_R_packages.R` 补齐（只装缺的，不动既有环境）。

```bash
Rscript -e 'source("environment/install_R_packages.R")'
```

字体统一 Helvetica；PDF 用基础 `pdf()` 设备而非 `cairo_pdf`，文字保持可编辑
（本机 Arial 的 CID 查找会失败，已避开）。

## 6. 术语

- 六基因签名：*CDK5RAP2, CD177, GRB10, UPP1, EXOSC4, VNN1*
- ARDS：acute respiratory distress syndrome；sepsis：脓毒症
- S1–S14：Supplementary Methods 的编号，与 `docs/METHODS_TO_CODE.md` 一致

## 7. 许可与引用

MIT，见 `LICENSE`。引用格式见 `CITATION.cff`。
代码按“能跑通、能对上数”整理，未做通用化封装；换数据请自行核对样本分组与平台注释。
