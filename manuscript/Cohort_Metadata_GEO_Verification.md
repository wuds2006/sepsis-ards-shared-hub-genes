# 63. 队列元数据核查记录：采样时点与供体血统（2026-11-04）

**目的**：清零 `Supplementary_Methods.md` S1／S11／S12、`Supplementary_Materials.md` Table S26 说明、`TableS26_cohort_metadata.csv` 中的 `[AUTHOR_INPUT_NEEDED]` 占位符，并确认数据库 release 版本与基金编号。

**检索方式**：NCBI GEO SOFT 接口逐系列拉取 `!Series_*`（`targ=self&form=text&view=brief`）与全部样本 `!Sample_characteristics_ch1`（`targ=gsm&form=text&view=brief`）；采样时点在 GEO 未标注时，回查源文献全文（Europe PMC / PMC 全文，或出版社摘要）。

---

## 1. 七个转录组队列

| 数据集 | 提交机构 / 国家（GEO `contact_institute` / `contact_country`） | 样本构成（GEO） | 采样时点 | 供体血统 |
|---|---|---|---|---|
| GSE76293 | GSK / United Kingdom | 94 样本：12 ARDS + 12 HVT + 70 体外培养 PMN（已排除） | **诊断后 48 小时内**（源文献：AJRCCM 2016, PMID 27064380，"Patients underwent venipuncture, fiberoptic bronchoscopy (FOB), and bronchoalveolar lavage (BAL) within 48 hours of diagnosis"）；GEO 无该字段 | GEO 无字段；源文献未报告族群（仅年龄、性别、ARDS 严重度分层） |
| GSE69528 | Benaroya Research Institute / USA（样本取自泰国） | 138 样本：29 melioidosis + 54 其他病原脓毒症 + 28 健康 + 27 T2DM（T2DM 已排除） | GEO 无该字段；源文献仅记"All blood samples were obtained at the Khon Kaen Regional Hospital, Khon Kaen, Thailand"，未注明相对入院的时点；招募期 2007–2008 年 | **全部为亚洲人**（源文献：Genome Biol 2009, PMID 19903332，"all subjects were Asian"）；泰国东北部人群 |
| GSE100159 | Baylor Research Institute / USA | 47 样本：18 melioidosis + 15 其他感染 + 12 对照 + 2 恢复期（分析用 45） | **GEO 无该字段**（仅 tissue / dataset / sample class 三个 characteristics） | **GEO 无字段**，源文献未报告；记录中只有提交机构，无采集地点 |
| GSE243066 | Osaka University（NGS core facility, RIMD）/ Japan | 49 样本：34 ARDS + 15 健康 | **入院 24 小时内**（GEO `overall_design` 明示） | GEO 无字段；队列招募于日本 |
| GSE151263 | University of Pittsburgh / USA | 7 PBMC 样本（3 ARDS + 4 脓毒症），13,996 细胞 | **插管/机械通气启动后 24 小时内**（源文献：JCI Insight 2020, PMID 32554932，"Peripheral blood samples were collected ... within 24 hours of initiation of mechanical ventilation"） | GEO 无字段；源文献未报告族群 |
| GSE167363 | University of California, Riverside / USA | 12 样本 / 6 供体：2 健康 + 2 脓毒症（1 例采两次，T0 与 T6），12,359 CD14+ 单核细胞 | **脓毒症诊断后 6 小时内**（源文献：J Leukoc Biol 2021, PMID 34558746，"PBMC collected during the critical period between sepsis diagnosis and 6 h"；GEO 样本名 S1 T0 / S1 T6） | GEO 无字段（仅有 Sex、age range）；源文献未报告族群 |
| GSE65682 | University of Malta / Malta（提交者 Scicluna、van der Poll；MARS 队列本身为荷兰两家 ICU） | 479 例脓毒症休克 + 对照，PAXgene 全血 | **ICU 入院时，并在 ICU 住院期间**（GEO `overall_design` 明示） | GEO 无字段；已发表 MARS 队列描述中 White 占 86–88%（n=522 子集 White 449/522；汇总表 n=2407 White 2106） |

### 结论（写入稿件时的口径）

- **七个 GEO 系列的样本记录中均无 ancestry / ethnicity / race 字段**。因此稿件只能写 "ancestry is not annotated in the GEO records"，不得写 "all cohorts were of European ancestry"。
- 唯一有族群陈述的是 GSE69528 的源文献（all Asian）。GSE65682 的族群构成来自 MARS 队列的已发表描述（White 86–88%），不是 GEO 字段。
- GSE100159 的采样时点与采集地点在 GEO 与源文献中均无记录，只能写"未记录"。
- 转录组队列的国别可按提交机构或源文献写：英国、泰国（东北部）、美国、日本、荷兰。

---

## 2. 遗传学与外部数据库版本

| 项目 | 原占位符 | 核实结果 | 依据 |
|---|---|---|---|
| 脓毒症 GWAS | —（正文 Data Availability 已有，S7 未标 GWAS Catalog 号） | IEU OpenGWAS `ieu-b-4980` = GWAS Catalog **GCST90281174**（11,643 例 / 474,841 对照） | 用户 2026-11-04 确认；`44_Fact_Citation_Verification` 亦记为一致 |
| ARDS 结局 GWAS | — | FinnGen **R11**，endpoint `J10_ARDS`（431 例 / 493,301 对照） | 用户确认 |
| HPA Blood Atlas | `[AUTHOR_INPUT_NEEDED: Blood Atlas version and date]` | **version 24，accessed November 6, 2025** | 用户确认；与正文 METHODS 一致 |
| Open Targets Platform | `[AUTHOR_INPUT_NEEDED: confirm release number and exact date]` | **release 26.09，发布日期 24 September 2026**（dataVersion 26.09 / apiVersion 26.9.0） | Open Targets 官方 release notes 与 26.09 发布博客 |
| ChEMBL | `[AUTHOR_INPUT_NEEDED: release and date]` | **release 37（2026 年 5 月发布），accessed 17 September 2026** | EMBL-EBI ChEMBL 下载页：CHEMBL37 = May 2026，DOI 10.6019/CHEMBL.database.37 |
| 蛋白结构 | —（S12 已写明） | PDB **3EUF**（UPP1，共晶配体 BAU）、**7SLV**（VNN1，共晶配体 9TC）；无实验结构的四个靶点用 AlphaFold 模型 | 用户确认；`44_Fact_Citation_Verification` 记为一致 |

---

## 3. 基金编号

| 项目 | 原值 | 确认值 | 依据 |
|---|---|---|---|
| 福建省财政厅项目 | `{2025}881`（正文 FUNDING、中文译文、51_ 内部文档） | **`{2025}876`** | 用户 2026-11-04 确认；花括号 `{2025}` 为正式编号写法 |

注：`Manuscript_Submission.docx` 标题页 Funding 段此前已是 876（两 run 拆分：`…{2025}8` + `76`），正文 FUNDING 段仍为 881，本次已统一。

---

## 4. 落盘清单

| 文件 | 改动 |
|---|---|
| `Supplementary_Methods.md` S1 | 采样时点按 GEO／源文献逐队列写明；血统改为"GEO 未标注 + 仅 GSE69528 报告 all Asian + 各队列招募国" |
| `Supplementary_Methods.md` S7 | 补 GWAS Catalog `GCST90281174` |
| `Supplementary_Methods.md` S11 | HPA version 24，accessed November 6, 2025 |
| `Supplementary_Methods.md` S12 | Open Targets release 26.09（24 Sep 2026）；ChEMBL release 37（May 2026），accessed 17 Sep 2026 |
| `Supplementary_Materials.md` Table S26 说明 | 整段重写为"每队列一次"的统一格式（平台／组织／样本构成／采样时点／国别／血统），删除占位符 |
| `Supplementary_Materials.md` 来源行 | 描述改为"七个转录组队列 + 遗传参考数据集" |
| `TableS26_cohort_metadata.csv`（`figures/` 与 `code_archive/` 两份） | 7 行 ancestry 与 sampling_time_point 回填；GSE65682 由未核实的 "European (as reported…)" 改为 "Predominantly White (86%…)" |
| `02_Manuscript_Revised.md`、`03_Manuscript_Translation_CN.md` | FUNDING `{2025}881` → `{2025}876` |
| `Manuscript_Submission.docx` | 正文 FUNDING 段 `{2025}881` → `{2025}876`（仅该段变化，385 段／斜体 130／分页符 44 不变） |
| `Supplement Material/Supplementary Material.docx` | 4 段（S1、S26 说明、S11、S12）按 run 级重写；251 段／TOC 域 300／斜体 134／分页符 27 不变 |
| `51_Title_Page_CoverLetter_AI_Disclosure_2026-10-19.md` | 基金编号与 P1 待办置为已解决 |

**备份**：`/tmp/ards_backup_20261104/`（`Supplementary_Methods.md.bak2`、`Supplementary_Materials.md.bak`、`02_Manuscript_Revised.md.bak`、`TableS26_cohort_metadata.csv.bak`、`Supplementary Material.docx.bak`、`Manuscript_Submission.docx.bak2`）。

**GEO 原始记录留档**：`/tmp/geo/{GSE*}.self.txt`、`/tmp/geo/{GSE*}.gsm.txt`（会话临时目录，如需长期保留应移入 `00_rawdata/`）。
