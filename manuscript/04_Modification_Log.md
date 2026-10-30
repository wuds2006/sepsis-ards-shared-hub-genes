# Part 3 — Modification Log（修改日志）

对应：01_Review_Report.md 的审查结论 → 02_Manuscript_Revised.md 的具体落实。  
原则：less is more。原文已经历多轮去 AI 化处理，整体自然度较高，本轮仅做审查确认的必改项与机械规范化，未做全篇重写。

## A. 内容性修改（对应审稿报告必改/应改项）

| #   | 位置                    | 修改内容                                                                                                                                       | 依据                                    |
| --- | --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------- |
| A1  | Abstract, M\&M        | MR 句改为区分正向（效能不足，MDOR=2.26）与反向（有信息量，MDOR=1.04）                                                                                              | 审查 1.6；原文"null but uninformative"以偏概全 |
| A2  | Abstract, M\&M        | "complement activation was the only shared pathway" 加上"top-ranked dysregulated immune processes"限定                                         | 审查 2.2；结果是 top-20 限定的秩约束陈述            |
| A3  | Introduction, 第3段     | 对文献 6-8（原14-16）的笼统指控加 "to our knowledge" 限定                                                                                                | 审查 2.1；未逐篇核实不得全称断言                    |
| A4  | Introduction, 末段      | 删除结果预告（6 基因、验证 AUC、-10.01kcal/mol 等），压缩为研究前提与设界策略两句                                                                                        | 审查 1.2；AJRCCM 引言以目的收尾                 |
| A5  | Methods, Data Sources | 增补 GSE100159 剔除 2 例恢复期样本                                                                                                                   | 审查 3.2；使 Table 1 的 33+12+2=47 可由正文推导  |
| A6  | Methods, 机器学习段        | 删除泄漏进 Methods 的结果（9/6 候选、MUC1 AUC=0.575、七基因敏感性数值），改为指向 Results 的规则性陈述                                                                      | 审查 1.1                                |
| A7  | Methods, 标题           | "Immune infiltrations" 改为 "Immune Infiltration"                                                                                            | 审查 1.4；标题大小写与单复数一致性                   |
| A8  | Results, 外部确认段        | 精确单侧 95% 置信下限 0.992 更正为 0.993（GSE100159；源值 0.9925）                                                                                         | 审查 3.1；与 Table 2 一致                   |
| A9  | Discussion, 第1段       | "ARDS monocytes sit beyond the sepsis state" 改为 "disease-group monocytes occupy later pseudotime positions than their respective controls" | 审查 2.4；跨队列无共同对照，不能确立次序                |
| A10 | Discussion, MR 段      | 删除与 Results 重复的效能数字（28%、OR 2.26），改为引用 Results                                                                                              | 审查 1.3                                |
| A11 | Discussion, CMap 段    | 删除对 Results 零分布数字的逐一枚举（13/100、17/100、87/100、2/100 重复叙述），压缩为解读性三句                                                                           | 审查 1.3                                |
| A12 | Discussion            | 两个加粗段中小标题去加粗，保留为小标题行                                                                                                                       | 约束(4)：不用粗体表强调；结构小标题保留                 |

## B. 机械规范化（约束 1-17 的落实）

| #  | 规则             | 落实                                                                                                                 |
| -- | -------------- | ------------------------------------------------------------------------------------------------------------------ |
| B1 | 统计格式（约束 3、6）   | 全文 P=/AUC=/n=/β=/rg=/rho=/NES=/Zsummary= 等统计符号与数值间空格删除（如 "P = 0.249"→"P=0.249"）；" × 10"→"×10"；"95% CI, "→"95% CI " |
| B2 | 破折号（约束 2、7）    | 全文本无 em dash/——，保持为零；范围号保留 en dash                                                                                 |
| B3 | 半角符号（约束 8）     | 全文符号已为半角，无需改动                                                                                                      |
| B4 | AI 高频连接词（约束 9） | 无 "Therefore/Moreover/Furthermore/not only...but also/In conclusion" 句首模板；"In addition" 仅 1 处且自然，保留                |
| B5 | 粗体/斜体（约束 4）    | 删除唯二的两处加粗（讨论小标题）；基因名斜体为 AMA 命名规范而非强调，予以保留                                                                          |
| B6 | 词汇规范化（约束 1）    | 全文无 leverage/delve into/underscore/remarkably/dramatically/strikingly，无需替换                                         |
| B7 | 长句（约束 11）      | 审查确认的超长句已通过 A4/A6/A11 的压缩顺带拆解；其余方法学长句为枚举性技术描述，拆分会损失精确性，按"less is more"保留                                           |
| B8 | 未证实内容（约束 13）   | 未新增任何无出处内容；A3 对原有全称断言作了保守化                                                                                         |

## C. 编号系统重排（约束 8、9）

| 系列                          | 处理                                                          | 结果                                          |
| --------------------------- | ----------------------------------------------------------- | ------------------------------------------- |
| 主图 Figure 1-10              | 首次引用顺序本就有序（图2→图10），未变动                                      | 无需重排                                        |
| 参考文献                        | 原引用序在（5）后跳至（14-16），再跳回（6-9），（34-36）先于（17-33）出现；已按首次出现顺序整体重排 | 正文引用序号与文献列表同步更新（映射见 05_Renumbering_Maps.md） |
| Supplementary Figure S1-S19 | 原 S18、S19 在 Methods 中先于 S3 被引用；已按首次出现顺序重排                   | S18→S3、S19→S4，其余顺移                          |
| Supplementary Table S1-S21  | 原 S18、S19 提前引用，S20/S21 尾部逆序；已重排                             | S18→S4、S19→S5、S21→S20、S20→S21，其余顺移          |
| Supplementary Methods S1-S9 | 原 S4 先于 S2/S3 引用；已重排                                        | S4→S2、S7→S3、S8→S4、S9→S5、S5→S6、S6→S7         |

## D. 未处理项（需作者决定）

1. Data Availability 仍为占位符 [AUTHOR_INPUT_NEEDED: URL and DOI]，等待 GitHub/Zenodo 信息。
2. Figure 4 面板字母序（4F 先于 4A-D 被引）未重排：重拼 composite 会改动已定的图文件与图注，建议投稿前一并处理（审查 1.5）。
3. Table 2 表头合并单元格重复问题属于 docx 排版层，已记录在案，未在 markdown 层改动。

## 2026-09-25 (Plan B): main-figure citations rewired to RESULTS_NPG set

- All main-figure citations in the revised manuscript (EN + CN) were rewired  
  from the old Figures_final panel letters to the 03_figures_RESULTS_NPG set  
  (step96, 22 replacements each + Figure 3C addition); docx regenerated.
- New main-figure citations added where the old text cited none: Figure 4E  
  (single-gene GSEA NES heatmap), Figure 7A-D (MR forests, power, SMR-HEIDI),  
  Figure 9 whole-figure, Figure 10A (HPA Blood Atlas).
- Figure-number first-appearance order verified monotonic 1 to 10.
- Main figure legends for this set installed at figures/Figure_legends.md.

---

## 2026-10-07：Introduction 按 7 步框架重写（去防御性表达）

按用户指定的 7 步框架（研究背景 → 研究现状 → research gap → 研究问题 → commitment/论证思路合并段 → 贡献）将 Introduction 由 7 段重写为 5 段；应用户要求将"研究问题"与"commitment 与论证思路"合并为第 4 段。引用 1–19 全部保留且首引顺序连续，无需重排；无新增事实。

防御性表达处理：

- 删除 "the contribution lies not in nominating these genes de novo, but in extending..." 自贬结构，改为正向陈述（末段 "the analysis extends ... and adds ... depth"）。
- 删除 "framed, from the outset, as in silico drug repurposing prediction rather than drug repositioning in the validated sense" 预防性划界，改为正向描述所做工作及对照。
- "we bounded every downstream claim ..." 总结式对冲融入四条正交证据线的具体描述。
- "Part of the reason may be simple" 弱断言改为直接陈述（两病共享炎症机制却分开研究）。
- "to our knowledge without systematic off-target sample removal" 不确定措辞改为肯定句。
- "has rarely been asked" 明确为 "has rarely been examined at the level of coordinately regulated genes"。

同步更新：02_Manuscript_Revised.md、03_Manuscript_Translation_CN.md、manuscript_revised.docx（step95 重生成）。脚本：scripts/step100_update_introduction.py。

---

## 2026-10-07（二）：DISCUSSION 去防御性表达（14 处，用户确认后执行）

按已确认清单修改 Discussion：①重复堆叠免责各保留一处（混合解读 caveat、GO 分组非实验证据、donor 数、MR 正向零结果、CMap 循环性）；②6 处 "should be read as" 指导式框架保留 limitations 引导句 1 处、其余改直接陈述；③删预防性术语声明（"hub gene" 操作性含义整句）、弱断言（presumably）、防御性收尾（P4/P8 两遍 "rather than" 改正向）、泛化免责声明（"everything here is computational"）。所有具体 limitation（效能、LDSC、VNN2 选择性、MM/GBSA、n=24、非感染对照未测、nomogram 0.708、无临床协变量）全部保留。不新增事实、不改数据、不动引用编号。Discussion 字数变化见脚本输出。同步更新 02_Manuscript_Revised.md、03_Manuscript_Translation_CN.md、manuscript_revised.docx（step95 重生成）。脚本：scripts/step101_discussion_deai.py。

---

## 2026-10-15：Discussion 补入"共享基因名单分析情境依赖"段（CX3CR1/PID1/PTGDS 核对）

背景：用户要求核对 Jiang et al.（Front Immunol 2024;15:1480542，本文参考文献 6）提名的 *CX3CR1*、*PID1*、*PTGDS* 是否在本文 179 个共享候选基因集中，并将核对结果补入 Discussion。

核对数据（01_DEG/GSE69528_DEG_all.csv、GSE76293_DEG_all.csv；候选集 03_intersection/Shared_CGs.txt）：

| 基因       | GSE69528（脓毒症）                        | GSE76293（ARDS）                              | 是否入 179 候选集 |
| -------- | ------------------------------------ | ------------------------------------------- | ----------- |
| *CX3CR1* | log2FC=−1.13，校正P=0.000195（Down DEG）  | log2FC=−0.25，校正P=0.691（NS）                  | 否           |
| *PID1*   | log2FC=−1.10，校正P=2.97×10⁻⁸（Down DEG） | log2FC=−0.79，校正P=0.006（显著但未达 |log2FC|>1，NS） | 否           |
| *PTGDS*  | log2FC=−1.26，校正P=1.48×10⁻⁵（Down DEG） | log2FC=−0.27，校正P=0.309（NS）                  | 否           |

三基因均死于"两病交集"的 ARDS 侧 DEG 门（阈值 |log2FC|>1 且校正P<0.05）；其中 *CX3CR1*、*PTGDS* 在 GSE76293 连名义显著都不满足。FIMMU 文（文献 6）的检出依赖其不同的 ARDS 数据集（GSE32707、GSE10474）与宽松阈值（|fold change|>0.3、P<0.05）。

修改：在 Discussion 第一段（六基因总述）与第二段（生物学关联）之间新增一段，以文献 6 为外部对照，说明共享基因名单对队列构成与阈值选择的依赖，并衔接本文 curation-first 设计（预先固定阈值、检验前剔除脱靶样本、独立验证）的合理性。引用编号不变（复用文献 6，无新增、无重排）。

同步更新：02_Manuscript_Revised.md、03_Manuscript_Translation_CN.md、manuscript_revised.docx（step95 重生成）。脚本：scripts/step104_add_context_dependence_para.py。

## 2026-10-15：GSE65682 生存/预后分析与扩展建模（整合）

### 触发

在四个数据集上完善生存分析、预后模型构建与评价。核查后仅 GSE65682（MARS 脓毒性休克队列，479 例，28 天死亡 114 例）含生存终点；GSE76293（体外 PMN 实验）、GSE100159（类鼻疽）、GSE243066（ARDS vs 健康）均无生存字段。

### 分析（13_survival_prognosis/）

| 脚本                                    | 内容                                                                                |
| ------------------------------------- | --------------------------------------------------------------------------------- |
| step105_survival_prognosis_GSE65682.R | 6 hub 基因 Cox、风险评分 KM、timeROC、C 指数、校准、nomogram、调整后森林图                              |
| step105b_pergene_KM_GSE65682.R        | 单基因 KM                                                                            |
| step106_extended_modeling_GSE65682.R  | 171 基因 LASSO-Cox（10 折 CV）、弹性网/岭回归对照、随机生存森林、bootstrap 乐观校正（B=300）、NRI/IDI、DCA、校准比较 |
| step107_sixgene_prognosis_figure.R    | Figure S20 组装                                                                     |

### 主要结果

- 六基因：单基因均不显著；多因素 Cox 中 UPP1 显著（HR 0.64, P=0.026）；风险评分 KM P=0.002（28 天生存 0.70 vs 0.83）；校正 C 指数 0.566；在年龄+性别+endotype 模型上无增量（LRT P=0.161）。
- 扩展建模：LASSO 1-SE 规则未选基因（信号弱），lambda.min 保留 15 基因（无 hub 基因）；校正 C 指数 0.622（联合临床 0.638）；28 天 AUC 0.698；NRI 0.132（P=0.026）、IDI 0.094；RSF 无增益（OOB C 0.610）。
- 解读：诊断标志物与预后信息可分离；hub 基因为疾病状态标志物而非预后标志物。

### 手稿修改（英文+中文同步）

- Abstract：Methods 句加预后评价；M\&M 加预后结果句；Conclusions 改为"疾病状态标志物而非预后标志物"。
- Methods：Data Sources 加 GSE65682；新增小节"Prognostic Modeling in an Independent Septic Shock Cohort"。
- Results：新增小节"Prognostic Performance in an Independent Septic Shock Cohort"（两段）。
- Discussion：新增段落"Prognostic information is separable from the diagnostic signal"；limitations 段加预后局限句。
- DATA AVAILABILITY：加 GSE65682。
- REFERENCES：新增 39-43（Davenport/Blanche/Friedman/Wright/Pencina）。
- Supplementary_Materials.md：新增 Table S22、S23 条目与 Figure S20、S21 图注；File Manifest 更新。
- Supplementary_Methods.md：新增 S14（GSE65682 预后建模）。
- 新增补充图：FigureS20_sixgene_prognosis_GSE65682、FigureS21_extended_prognostic_modeling（PDF/PNG/TIFF 300dpi）。
- 新增补充表：TableS22_prognostic_model_comparison.csv、TableS23_lasso_prognostic_signature.csv。

## 2026-10-16：路径 A 重构（AJRCCM 合规改造：Methods 移入 Online Supplement + Results/Discussion 大幅压缩）

### 触发与依据

AJRCCM Original Research 限值（C 级多源一致记录，官方 OUP 页被反爬拦截；投稿前经投稿系统 checklist 复核）：正文 ≤3,500 词、摘要 ≤250 词、标题 ≤100 字符、关键词 ≤5 个、参考文献 ≤50 条。现稿实测正文 10,939 词、摘要 373 词、标题 155 字符、关键词 7 个。经用户确认执行方案 A（计划表 `13_Restructure_Plan_PathA_2026-10-16.md`）：标题选候选 1（84 字符）、关键词保留 machine learning、其余按计划表执行。

### 修改（02_Manuscript_Revised.md，三轮压缩）

- 标题：换为 "Shared Hub Genes and Immune Signatures of Sepsis and ARDS: A Multi-Omics Integration"（84 字符）；Running head 不变。
- 摘要：373 → 250 词（五段式结构不变，保留成对数字与反向 MR 句）。
- 关键词：7 → 5（sepsis; ARDS; hub genes; machine learning; Mendelian randomization）。
- INTRODUCTION：719 → 455 词（保留三段式限制陈述与四条正交验证线）。
- METHODS：3,412 → 300 词摘要版；全文移入 Online Supplement，结尾指针 "Complete methods... are provided in the Online Supplement (Supplementary Methods S1–S14)"。
- RESULTS：3,955 → 1,578 词；11 小节保留（5 个小节标题精简：GO/KEGG、Machine Learning Identifies Six Hub Genes、Expression Patterns and Enrichment、Prognostic Performance in Septic Shock、MR Under Limited Power、Single-Cell Localization）。
- DISCUSSION：2,853 → 1,161 词；CMap 解构段并入 druggability 段，局限段保留 8 句骨干。
- 参考文献：43 → 22 条（1–19 序号不变；39→20 Davenport、37→21 Schalkwijk、38→22 Liu C；工具类 17–33 与预后方法类 40–43 随 Methods 移入补充材料 author-year 制）。映射详见 05_Renumbering_Maps.md 新增小节。
- 压缩均为编辑性删减：被删数字在补充图表/表中有据可查；以下红线项全数保留——唯一 T5 句原文不动；成对引用 7 对（KM P=0.002 与 LRT P=0.16、MDOR 2.26/1.04、UPP1 HR 0.64 与单变量 P=0.640、真无因果+表型误分类、0.753–0.998+病对病 0.50 等）；阴性留存（E23 病对病、E24 单基因全阴、E29 无增量、E31 RSF、E10 MUC1 否证、E42 Treg 不解读、M0 降级读法等）。

### 核验（scripts/verify_manuscript.py，新增）

正文 3,494 ≤3,500 ✓；摘要 250 ✓；标题 84 ✓；关键词 5 ✓；22 条引用全在且首现严格递增 ✓（范围引用 (6–8)/(9–12)/(13–15) 展开逻辑修复）；主图 1–10、主表 1–3 全在 ✓；em dash 0、P= 无空格、robust/coherent 等 0 ✓；T5 句唯一 ✓；四组成对引用 True ✓。

### 同步更新

- Supplementary_Methods.md：承接正文 Methods 全文，合并去重补入 S2（选择稳定性 200 bootstrap/50 种子/标签置换、rms lrm 岭正则预测模型+nomogram+DCA）、S3（ssGSEA 18 套谱系标志集）、S6（Hedges g meta 分析+DeLong+ComBat 探索性）、S7（Brion 型功效分析 MDOR 2.26/1.04、LDSC、pQTL 可及性）、S9（clusterProfiler GO/KEGG、Zsummary 模块保留、kME）、S12（CMap Enrichr 三重空值背景、Vina 1.2.5 回对接 RMSD 1.54/0.69 Å、OpenMM 8.5.1 MM/GBSA 全参数、Open Targets/ChEMBL）；S14 内数字引用 (40)/(42)/(43) 转为 author-year（Blanche/Wright/Pencina）。头部加 Methods 移驻说明。
- Supplementary_Materials.md：标题同步为新标题；File Manifest 的 Supplementary_Methods 行更新为 "Complete methods relocated from the main text (Path-A, 2026-10-16)"。
- 05_Renumbering_Maps.md：新增 43→22 映射小节。
- 03_Manuscript_Translation_CN.md：按新版英文稿全文重译（见当日翻译同步记录）。
- manuscript_revised.docx：step95 重生成，SUBHEADS 更新为压缩后小节标题。


## 2026-10-16（第二批）：Discussion 修订计划表 ①–⑥ + 段7 削减包（#10/#11 出批，作者批准后执行）

### 修改内容（02_Manuscript_Revised.md Discussion）

| #   | 位置 | 修改                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | 依据                                 |
| --- | -- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------- |
| ①   | 段1 | "corroborated at the protein level for CD177" → "corroborated for CD177 protein across nine infectious contexts"                                                                                                                                                                                                                                                                                                                                                                                                                                    | E53 可用范围限定（#10 修订①）                |
| ②   | 段4 | "CEBPB remained the dominant upstream regulator" → "CellOracle again placed CEBPB as the dominant upstream regulator"                                                                                                                                                                                                                                                                                                                                                                                                                               | CellOracle 计算推断 placed 语气（#10 修订②） |
| ③   | 段8 | "every key finding" → "nearly every key finding"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | E45 供体 4/5 同向 1 例外（#10 修订③）        |
| ④   | 段9 | "six coordinately regulated hub genes" → "six co-expressed hub genes"                                                                                                                                                                                                                                                                                                                                                                                                                                                                               | E11 共表达证据不越 D9/D10（#10 修订④）        |
| ⑤   | 段6 | 插入 "A random survival forest added no gain, so the advantage reflects variable selection, not nonlinearity."                                                                                                                                                                                                                                                                                                                                                                                                                                        | E31 RSF 阴性留存补回（#11 修订⑤）            |
| ⑥   | 段2 | 尾句从句置换为 "and the prespecified standard proved real: *MUC1* was rejected at external validation, and *DYRK2* and *TCTN1* reversed direction across cohorts."                                                                                                                                                                                                                                                                                                                                                                                         | E10/E09 阴性留存补回（#11 修订⑥）            |
| 削减包 | 段7 | 195→174 词（−21）：任一→某个；删 direct/conceptually/known/still；"co-crystallized ligand"→"co-crystallized"；"is inhibited with nanomolar potency by"→"by the nanomolar, orally active"；"the protective effect of ... knockout in a mouse model of"→"... knockout protection in mouse"；"disease-relevant signatures of the same size"→"size-matched disease signatures"；"whose reproducibility has been questioned"→"of questioned reproducibility"。保护项不动：S1 druggability 句、S4 D18 成对（无批准药+VNN2 未测）、2/100 零分布锚点、(18,19)(21)(22) 引用、S6 hypothesis-generating 收口 | #09 段7 弹性 + #11 三句合并               |

### 核验（scripts/verify_manuscript.py）

正文 3,498 ≤3,500 ✓（Discussion 1,165）；摘要 250 ✓；22 条引用全在且首现严格递增 ✓；主图 1–10、主表 1–3 ✓；em dash 0、P= 无空格、禁用词 0 ✓；T5 句唯一 ✓；五组成对引用全 True ✓；MUC1/DYRK2/TCTN1/RSF 四项阴性发现在 Discussion 补回确认 ✓。

### 同步更新

- 03_Manuscript_Translation_CN.md：讨论九段对应七处同步（九种感染背景／MUC1 否证+方向反转／CellOracle 再次置于／随机生存森林句／段7 压缩／几乎每一／共表达），中英关键 token 对齐核验通过、旧句全清。
- manuscript_revised.docx：step95 重生成；基因斜体含 MUC1/DYRK2/TCTN1 各 2 处斜体 run；em dash 0、P 空格 0；六项修订字符串全在。

### 未执行项（沿用待办）

段4 "metabolic decoupling in established ARDS" 文献 (10) 归属核实（待查①，未批准删除，维持原文）；E59 字面点名（可选 +9 词，未批准）；Open Targets/HPA 查询日期、GitHub URL + Zenodo DOI、2026 年文献 Zotero 核对。

## 2026-10-16（第三批）：结构连贯性修补「组合 A」（#13 出批，作者批准后执行）

依据：`19_Structure_Coherence_Audit_and_Final_Skeleton_2026-10-16.md` 第三节词数账本组合 A（⑦甲＋⑧＋⑨＋⑪甲＋⑫甲＋⑬＋⑮甲＋⑱，净 −7）。

| #  | 位置    | 修改                                                                                                                                                                                                                                                                                        | 依据（审计编号）                               |
| -- | ----- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------- |
| ⑦甲 | 段2    | 段首加答句 "The difference from earlier screens is threshold- and cohort-dependent rather than a refutation:"；删段尾重复半句 "The divergence is not a refutation; it quantifies how strongly shared-gene lists depend on cohort composition and threshold choice,"（保留 #11 修订⑥ 的 MUC1/DYRK2/TCTN1 阴性句） | C2／#12 遗留项（−6）                         |
| ⑨  | 段3    | "rather than a random list" → ", not an arbitrary list"                                                                                                                                                                                                                                   | B1 概念强度与段9 "co-expressed" 对齐（−1）       |
| ⑬  | 段4    | "These data fit ARDS" → "These blood data fit ARDS"                                                                                                                                                                                                                                       | B2 数据层级标记（+1）                          |
| ⑱  | 段4／段8 | E23 定界句重定位：段8 删 "The signature does not discriminate ARDS from sepsis at the patient level."（−12）；段4 末句前插入 "The six-gene score did not, however, discriminate the two diseases at the patient level."（+14）                                                                                | E1 位置错置（净 +2）                          |
| ⑪甲 | 段5    | "That reading fits the shared genes" → "That reading, together with null SMR-HEIDI tests for the five genes with cis-eQTL, fits the shared genes"                                                                                                                                         | **A1 唯一 🔴 项**：疾病层 MR null 缺基因层桥接（+11） |
| ⑧  | 段8    | 主题句 "These findings should be interpreted against several constraints." → "The constraints that follow concern scope: comparison group, cohort, and clinical context."                                                                                                                    | C1 空泛框架句＋7→8 过渡（实测 +4，非 #12 所记 +5）     |
| ⑮甲 | 段9    | 删半循环修辞句 "For a comorbidity this common and this lethal, markers of disease state are likely to come from that intersection."                                                                                                                                                              | A4（−19）                                |
| ⑫甲 | 段9    | "and computational screening nominates VNN1 as the lead candidate." → "and tractability profiling nominates VNN1 and UPP1 for first testing."                                                                                                                                             | A3 与段7「按可成药性非结合能排序」及 E59 排序不一致冲突（+1）   |

### 核验（scripts/verify_manuscript.py）

正文 **3,491** ≤3,500 ✓（Discussion 1,158；Intro 455／Methods 300／Results 1,578）；摘要 250 ✓；标题 84 字符 ✓；关键词 5 ✓；22 条引用全在且首现严格递增 ✓；主图 1–10、主表 1–3 ✓；em dash 0、P= 无空格、禁用词 0 ✓；T5 句唯一 ✓；五组成对引用全 True ✓（含新增 AUC 0.753–0.998 与病对病 0.50 成对）。

### 同步更新

- 03_Manuscript_Translation_CN.md：讨论段3/4/5/6/9/10 六处对应改写（阈值与队列依赖答句／非任意名单／血液数据／E23 定界句／SMR-HEIDI 桥接／适用范围主题句／可成药性提名），中英七组关键 token 对齐核验通过、六组旧句全清、CN em dash 0。
- manuscript_revised.docx：step95 重生成；9 个基因斜体 run 齐备（UPP1 15、CD177 13、VNN1 11 等）；em dash 0、P 空格 0；七项新增字符串全在、三组旧串全清。

### 未执行项（沿用待办）

⑭（段1/段9 补体激活补 top-20 截断限定，+5）【待作者核实】截断稳健性后决定；⑯（段6 单队列重复句 −21/−14）、⑰（段8 对接重复句 −13）为可选腾挪，本轮未动；待查① 文献 (10) 归属；Open Targets/ChEMBL 检索时点；E59 字面点名（可选 +9）；HPA 查询日期、GitHub URL + Zenodo DOI、2026 年文献 Zotero 核对。

## 2026-10-16（第四批）：AJRCCM 文献对话插入 AJ-1～AJ-4（作者批准后执行）+ ⑰/⑯ 腾挪启用

依据：`21_Verified_References_and_Pending_Citations_2026-10-16.md` 第四节 AJ-1～AJ-4（均已 PubMed 核验）；词数支付按 #14 既定规则「若新增文献插入，须先从段7/段6 可选腾挪（⑰ −13、⑯ −21）中支付」。

### 插入内容（02_Manuscript_Revised.md Discussion）

| 标签        | 位置          | 插入句                                                                                                                    | 证据编号 / 用途                                             |
| --------- | ----------- | ---------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------- |
| AJ-3→(21) | 段1 末句追加     | "…it characterizes ARDS blood as well, where leukocyte expression separates reactive from uninflamed phenotypes (21)." | E38–E41 的外部锚点；回答 LIT-1：ARDS 血液表型可由白细胞表达分离             |
| AJ-1→(22) | 段4 CEBPB 句前 | "Most of the sepsis transcriptome is shared across infection sources (22)."                                            | E51/E46 的方法学先例；支撑「共享成分占主导」而非重复前人                      |
| AJ-4→(23) | 段5 末句追加     | "…rather than established causal targets, whereas a genetically supported intermediate exists for this axis (23)."     | E32–E37 的对照定位（既有 sRAGE 阳性中介 vs 本文六基因 null）；精确回答 LIT-7 |
| AJ-2→(26) | 段7 收口句追加    | "…strictly hypothesis-generating predictions, as befits endotype-dependent corticosteroid responses (26)."             | E54–E56 的降档依据（VANISH 事后分析：SRS 内型调节氢化可的松反应）            |

### 启用腾挪（既有的两处重复句）

- ⑰：段9（局限）删 "Docking and MM/GBSA provide predicted poses and single-point estimates interpretable for ranking only."（−13）——同一局限已由段7 收口句「interpretable for ranking only / hypothesis-generating」覆盖。
- ⑯：段7（预后）删 "The estimates are internal to one cohort, so the penalized signature requires an independent septic shock cohort before any prognostic use;"（−21），"if pursued, it should start from" → "Any prognostic development should start from"（0）——单队列局限已由段9 "The prognostic analysis rests on a single cohort with internal validation." 覆盖。
- 追加微压缩各 −1～−2：段4 "link the signature to"→"link it to"；段7 "*GRB10* would require targeting protein-protein interfaces"→"*GRB10* requires protein-protein interface targeting"；段9 末句 "No clinical covariates were available for adjustment…in these cohorts"→"No clinical covariates were available…unknown here"；两段插入句各精简 2 词。

### 文献重编号（22 → 26）

新增 21 Bos 2019、22 Burnham 2017、23 Jones 2020；旧 21→**24** Schalkwijk、旧 22→**25** Liu C；末条新增 26 Antcliffe 2019。映射已记入 `05_Renumbering_Maps.md`（Entries 1–20 不变）。

### 核验（scripts/verify_manuscript.py，n=26）

正文 **3,497** ≤3,500 ✓（Intro 455／Methods 300／Results 1,578／Discussion 1,164）；摘要 250 ✓；标题 84 ✓；关键词 5 ✓；**26 条引用全在且首现严格递增 1→26 ✓**；主图 1–10、主表 1–3 ✓；em dash 0、P= 无空格、禁用词 0 ✓；T5 句唯一 ✓；五组成对引用全 True ✓。

### 同步更新

- 03_Manuscript_Translation_CN.md：讨论四处插入句同步（反应型/非炎症型表型、脓毒症转录组跨感染源共享、同一疾病轴的遗传学中介、内型依赖皮质激素反应）＋(21)(22)→(24)(25) 重编号＋⑰/⑯ 对应删除与三处微压缩；中英七组关键 token 对齐核验通过、六组旧串全清、CN em dash 0。
- manuscript_revised.docx：step95 重生成；9 基因斜体 run 齐备；26 条参考文献全在；六项新字符串全在；旧句 "would require targeting" 已清；em dash 0、P 空格 0。

### 未执行项（沿用待办）

R-4（段7 (19) 归属改写三案待批）、R-5（段4 代谢脱耦短语删与否）、R-6（(14) CEBPB 全文核实）、R-7（(22→25) Vnn1 knockout 全文核实）、R-1～R-3（(13)/(15)/(16) 元数据修正）、R-8/R-9（补充方法 Bulik-Sullivan 歧义与文末书目）；Open Targets 检索时点、HPA 查询日期、GitHub URL + Zenodo DOI、2026 年文献 Zotero 核对。

## 2026-10-16（第五批）：M-0 重复句删除 + U5 判别对象限定句插入（作者批准后执行）

依据：`23_Conflict_Interpretation_Matrix_2026-10-16.md` 第五节 M-0（L109 整句重复缺陷）与第四节写法库 U5；作者批准「修 M-0 后只取一条 U5」。

### 修改内容（02_Manuscript_Revised.md Discussion）

| 项   | 位置                  | 改动                                                                                                                                                                                                                                                                                                                                                        | 词数  |
| --- | ------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --- |
| M-0 | 段6（L109）            | 删除 "Diagnostic and prognostic information are therefore largely separable in this gene set." 的**第二次出现**（同一行连续重复，编辑残留）；保留第一次与后句 "Any prognostic development should start from the broader candidate set rather than the hub genes."                                                                                                                        | −12 |
| U5  | 段4（L105）病对病阴性句并入前半句 | 原 "The six-gene score did not, however, discriminate the two diseases at the patient level." → "The external-validation estimates concern discrimination from health, not between the two diseases; the six-gene score did not, however, discriminate the two diseases at the patient level."（把 I-5 的「判别对象不同」写足：0.753–0.998 是对健康对照的判别，病对病 AUC=0.50 不与之矛盾） | +12 |

**词数账**：3,497 −12 +12 = **3,497** ≤3,500 ✓（Discussion 1,164 → 1,164）。与批准方案的偏差 +1（U5 从 11 词变为 12 词）：批准原句以 "These estimates" 起句，在段4 无明确先行词，改为 "The external-validation estimates" 以闭合指代；仍在限值内。

### 核验（scripts/verify_manuscript.py，n=26）

正文 **3,497** ≤3,500 ✓；摘要 250 ✓；标题 84 ✓；关键词 5 ✓；26 条引用全在且首现严格递增 ✓；主图 1–10、主表 1–3 ✓；em dash 0、P= 无空格、禁用词 0 ✓；T5 唯一 ✓；五组成对引用全 True ✓（含 AUC 0.753–0.998 与病对病 0.50 成对）。

### 同步更新

- 03_Manuscript_Translation_CN.md：段4 病对病句前补 "外部验证估计针对的是疾病与健康对照的区分，而非两种疾病之间的区分；"（M-0 无需改动——中文稿 L115 该句本就只有 1 次）。
- manuscript_revised.docx：step95 重生成；重复句计数 = 1（已去重）、U5 句在案。

### 本轮未动项

R-1～R-9 全部沿用悬置（R-4/R-5/R-6/R-7 优先）；脚本默认路径问题 M-1 未执行（须作者批准后改脚本）；写法库 U1–U4、U6–U14 未采纳（本轮只取 U5）。

## 2026-10-17（第六批）：#35 结构精修包 + #36 术语证据一致性包 + R-4 甲案 + 31_ N1（作者批准后整批执行）

依据：`41_Structural_Refinement_Discussion_2026-10-16.md`（W1v2/W2/W3/D-1a）、`42_Terminology_Evidence_Consistency_2026-10-17.md`（TE-A1～A7、TE-R1～R7、TE-D1～D10、TE-S1）、`31_Unexpected_Null_Results_Paragraphs_2026-10-16.md`（N1）。作者裁决：推荐路径整批执行；R-4 按甲案；34 词余量启用 N1。

### 摘要（02_Manuscript_Revised.md Abstract；净 0 词，250→250）

| ID    | 改动                                                                                                                                                                              | 词数 |
| ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -: |
| TE-A1 | "and immune microenvironment" → "and blood immune microenvironment"                                                                                                             | +1 |
| TE-A2 | "discriminated disease from control" → "discriminated disease from healthy controls"                                                                                            | +1 |
| TE-A3 | MR 句 → "Reverse Mendelian randomization constrained causal effects of ARDS on sepsis risk (minimum detectable odds ratio 1.04); the forward analysis was uninformative (2.26)." | −2 |
| TE-A4 | "the six-gene score stratified mortality ... yet added nothing to" → "the six-gene risk score stratified mortality ... yet added no increment beyond"                           | +2 |
| TE-A5 | "a penalized signature from the 179 shared candidates" → "a penalized 179-candidate signature"                                                                                  | −4 |
| TE-A6 | "the preferred flavonoid target" → "the preferred predicted binding partner"                                                                                                    | +1 |
| TE-A7 | "define ... lymphocyte-depleted ... rather than prognostic markers" → "mark ... cytotoxic-lymphocyte-depleted ... rather than independent prognostic markers"                   | +1 |


### Results（净 −3 词，1,578→1,575）

| ID       | 位置                       | 改动                                                                                                                                                                                    | 词数 |
| -------- | ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -: |
| TE-R1+R2 | External Confirmation 末句 | "(combined AUC=0.50 ... shared disease-response markers" → "(combined six-gene signature AUC=0.50 ... shared disease-state markers"                                                   | +2 |
| TE-R3    | Prognostic 第3句           | "added nothing to a model of" → "added no increment beyond a model of"                                                                                                                | +1 |
| TE-R4    | Prognostic 第5句           | "outperformed the six-gene score" → "outperformed the six-gene risk score"                                                                                                            | +1 |
| TE-R5    | MR 第3句                   | "meaningfully excludes causal effects" → "meaningfully constrains causal effects"                                                                                                     |  0 |
| TE-R6    | Single-Cell 第4句          | "knockout identified CEBPB" → "knockout placed CEBPB"                                                                                                                                 |  0 |
| TE-R7    | Therapeutic 第4句          | "agreed that the panel bound VNN1 favorably ... beyond identifying VNN1 as the preferred target" → "favored VNN1 at the target level ... no ligand-level prioritization is supported" | −7 |

### Discussion（1,164→1,149）

| 项                | 段落     | 改动                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| ---------------- | ------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| W1v2+TE-D1       | 段1 第3句 | "The shared immune phenotype was algorithm-independent: ... reproduced ... lymphoid depletion" → "Two distinct readouts reproduced this phenotype: ... returned ... cytotoxic-lymphocyte depletion"                                                                                                                                                                                                                                                                                                                                   |
| TE-D2            | 段2     | "the prespecified standard proved real" → "was enforced"                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| W2+TE-D3         | 段3 整段  | 重写为 "The six genes have annotation-level links rather than arbitrary membership: ... thicken the intersection without evidencing co-regulation, which requires experimental testing."（实测 −12 含 TE-D3；41_ 报价 W2 −8、42_ 复核 −9，以脚本口径为准）                                                                                                                                                                                                                                                                                                  |
| D-1a（R-5）        | 段4     | 删 "metabolic decoupling in established ARDS and the"（−6）                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| TE-D4            | 段4     | "the six-gene score did not" → "the combined six-gene signature did not"（+1）                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| N1+TE-D6+TE-D7   | 段6 首句  | 补 E24 并与 TE-D6 合并为 "Individually, none of the six was associated with 28-day mortality in univariable models; jointly the six-gene risk score stratified survival (log-rank P=0.002) but added no increment beyond ..."；"prognostic content" → "prognostic information"（合并句相对旧首句实测 +15，其中 TE-D6 +2、N1 边际 +13；31_ 对 N1 报价 +17 为未与 TE-D6 合并的口径，合并后避免重复 "stratified/28-day mortality"，少 4）                                                                                                                                               |
| W3+TE-D8（R-4 甲案） | 段7 整段  | 重写：druggability 句、"uridine analogs inhibit UPP1 ... RR6 inhibits VNN1, abolishing ... (24) and reproducing ... (25)"、删 "with narrow therapeutic windows"、CMap 句改为 "matched only 2 of 100 stringency-matched disease signatures and returned no glucocorticoid in the down-regulated query, while its direction reflects ... compendia of questioned reproducibility (18) and context-agnostic matching (19)"、"Both it and the docking energies are hypothesis-generating ... (26)"（实测 −13；41_ 报价 W3 −13、42_ 复核 −16，差异源于 R-4 甲案 +3 并入） |
| TE-D9            | 段8     | "No clinical covariates" → "No severity covariates"（0）                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| TE-D10           | 段9     | "operating against complement activation" → "set against complement activation"（0）                                                                                                                                                                                                                                                                                                                                                                                                                                                    |

### 补充方法与脚本

- TE-S1：`09_review_deAI_2026-09-25/Supplementary_Methods.md` 与 `05_supplementary/Supplementary_Methods.md` 两份同步——"played no role in candidate selection" → "played no role in discovery-stage nomination"；"served strictly as confirmation" → "served as final-set confirmation"，并在 MUC1 句末追加 "; because this step reduced seven ARDS-nominated candidates to six, downstream validation estimates carry a mild optimistic bias bounded by the seven-gene sensitivity analysis."（与正文 Methods 的七筛六声明一致；不计正文词数）。
- `scripts/verify_manuscript.py`：成对检查第 5 组正则由 `combined AUC=0\.50` 更新为 `combined six-gene signature AUC=0\.50`，以追踪 TE-R2 批准措辞；仅为核验规则跟随，不动任何数据口径。

### 核验（scripts/verify_manuscript.py，n=26，执行后实测）

正文 **3,479** ≤3,500 ✓（Intro 455／Methods 300／Results 1,575／Discussion 1,149；余量 21 词）；摘要 **250** ✓；标题 84 ✓；关键词 5 ✓；26 条引用全在且首现严格递增 ✓；主图 1–10、主表 1–3 ✓；em dash 0、P= 无空格、禁用词 0 ✓；T5 句唯一 ✓；五组成对引用全 True ✓。  
词数对账：Results −3 与 42_ 报价一致；Discussion 实测 1,149 = 42_ 推荐稿 1,136 + N1 合并实测 +13（42_ 第九节预估 N1 +17 为未合并口径，实际并入 TE-D6 句后净 +13，正文总计由预估 3,483 降为 3,479，余量 21）。

### 同步更新

- 03_Manuscript_Translation_CN.md：摘要 3 处（血液免疫微环境／健康对照／反向 MR 约束因果效应／6基因风险评分／179候选惩罚标签／优先预测结合伙伴；结论句"标记／细胞毒性淋巴细胞减少／独立预后标志物"）、结果 5 处（合并6基因标签AUC、疾病状态标志物、6基因风险评分、约束因果效应、置于上游调节因子、不支持配体水平优先级）、讨论 8 处（两种读出重现、注释水平关联段重写、删代谢解耦短语、合并6基因标签、E24 单变量阴性并入段6首句、段7 整段重译含 (18)/(19) 归属、无严重度协变量、段9"处于……背景之中／是否有任何一个构成药物靶点"）。
- manuscript_revised.docx：step95 重生成。

### 本轮结清与悬置

- **结清**：R-4（甲案已入稿）、R-5（D-1a 已删）；42_ 一致性包与 #35 精修包全部落地。
- **悬置**：R-6（(14) 是否含 CEBPB【待作者核实】）、R-7（(25) 是否含 *Vnn1* knockout 保护实验与 RR6 给药臂【待作者核实】）、三项统计事实（Cox PH／KM 切点敏感性／免疫细胞比较 BH 覆盖面【待作者核实】）、R-1～R-3、R-8/R-9；Open Targets 检索时点、HPA 查询日期、GitHub URL + Zenodo DOI。
- 余量 21 词未再分配（33_ X1 +26 超限、34_ C2c +18 可用但未批准）。

---

## 第七批（2026-10-17）：#37 语言精修 + #38 事实校正 + #39 期刊适配（三批合并执行）

**批准依据**：作者 2026-10-17 裁决——①三批合并执行（#37 → #38 → #39）；②摘要 AUC 与孤立 (2.26) 两项均按零词方案修；③Introduction **展开 WGCNA、不压词**（与 JA-02 Intro −6 互斥，取 JA-11）。  
**输入**：`43_English_Language_Refinement_2026-10-17.md`、`44_Fact_Citation_Verification_2026-10-17.md`、`45_Journal_Adaptation_Verification_2026-10-17.md`、`/Users/mac/Desktop/Author guidelines.docx`（AJRCCM 官方作者指南）。

### 一、#39 期刊适配（结构，本批核心）

| 项         | 位置                             | 改动                                                                                                                                                                  |                    Δ词 |
| --------- | ------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------: |
| **JA-01** | DISCUSSION 段9 → 新增一级标题         | 段9 整体移出 Discussion，前插 `CONCLUSION` 一级标题（段内文字零改动）                                                                                                                    | 0（Disc −47／Concl +47） |
| **JA-11** | INTRODUCTION 段4                | "intersected with WGCNA disease-associated modules" → "intersected with disease-associated WGCNA modules (weighted gene co-expression network analysis)"（首次出现处展开缩写） |                **+5** |
| JA-03a    | 摘要                             | "(AUC 0.753 to 0.998)" → "(areas under the curve 0.753 to 0.998)"                                                                                                   |                    +3 |
| JA-03b    | 摘要                             | "genes intersected with co-expression modules" → "genes intersecting co-expression modules"                                                                         |                    −1 |
| JA-03c    | 摘要                             | "but did not distinguish sepsis from ARDS" → "but not sepsis from ARDS"                                                                                             |                    −2 |
| JA-04a    | 摘要                             | "was uninformative (2.26)" → "was uninformative (minimum detectable odds ratio 2.26)"（消除被读作引文的风险；期刊明令摘要不列参考文献）                                                      |                    +4 |
| JA-04b    | 摘要                             | "from healthy controls in external validation" → "from controls in external validation"                                                                             |                    −1 |
| JA-04c    | 摘要                             | "Both diseases showed monocyte expansion" → "Both showed monocyte expansion"                                                                                        |                    −1 |
| JA-04d    | 摘要                             | "In 479 patients with septic shock," → "In 479 septic shock patients,"                                                                                              |                    −1 |
| JA-04e    | 摘要                             | "cytotoxic-lymphocyte-depleted blood compartment" → "cytotoxic-lymphocyte-depleted compartment"                                                                     |                    −1 |
| **JA-05** | REFERENCES 6                   | "Song X, Xu JF." → "Song X, et al."（期刊：最多 6 位作者，超出用 et al.）                                                                                                         |                     0 |
| JA-00     | `scripts/verify_manuscript.py` | 分节正则加入 `CONCLUSION`；正文求和加入 CONCLUSION；打印列表加入 CONCLUSION；成对检查 body 加入 CONCLUSION。**否则新增 Conclusion 的 47 词不计入正文，会造成假通过**                                              |                     0 |

### 二、#38 事实校正（沿用 44_ A 类）

| 项      | 位置                       | 改动                                                                                                                                                                                                                             | Δ词 |
| ------ | ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | -: |
| FA-01  | 摘要                       | "a penalized 179-candidate signature" → "171-candidate"（原始日志：171 基因 × 479 样本；179 为脚本内命名 LASSO-179 的误写）                                                                                                                         |  0 |
| FA-05  | Results ML 段             | "was stable but not unanimous (Supplementary Table S1)" → "(Supplementary Table S9)"（映射：旧 S7 → 新 S9）                                                                                                                           |  0 |
| FA-06丙 | Introduction 段3          | "inherits the reproducibility limits of cancer-cell-line reference compendia (18, 19)" → "inherits the reproducibility limits (18) and context-agnostic matching (19) of cancer-cell-line reference compendia"（与段7 R-4 甲案归属一致） | +3 |
| FA-02  | REFERENCES 13            | "*EBioMedicine*. 2025;123:106100" → "2026;123:106100"                                                                                                                                                                          |  0 |
| FA-03  | REFERENCES 15            | "*Adv Sci (Weinh)*. 2026:e76236" → "2026;13(52):e76236"                                                                                                                                                                        |  0 |
| FA-04  | REFERENCES 16            | "*Cell Death Discov*. 2026. doi:…" → "2026;12(1):385. doi:…"                                                                                                                                                                   |  0 |
| —      | FA-07（13-fold → 14-fold） | **未采纳**（可选，换算口径待作者裁定）                                                                                                                                                                                                          |  0 |

### 三、#37 语言精修（17 条推荐项，LE-09 取甲案）

- A 类空泛词/自评修饰/模板过渡：删 "Groundwork exists:"（−2）、删 "stringently"（−1）、删 "that follow"（−2，Discussion 段8）、删 "in short,"（−2，现 CONCLUSION 段）、"decisive" → "clear"（0）。
- B 类时态：asks→asked、result is→was、extends→extended、yields→yielded（均 0）。
- C 类搭配与断句："than control" → "than in controls"（+1）；"told a different story" → "gave a different result"（0）；段2 "…are not differentially expressed, and the prespecified standard was enforced:" → 断为两句（−1）；"read qualitatively and consistent with" → "read qualitatively that is consistent with"（+1）；"fits the shared genes as" → "is compatible with the shared genes acting as"（+3）。
- D 类连接："Two caveats bound" → "limit"（0）；"as befits" → "consistent with"（0）；段7 首句加 "nevertheless"（+1）。
- 未采纳：LE-17（摘要 (2.26) 消歧 +5，已由 JA-04a 零词方案实现）、LE-18（删 "key"，可选）。

### 四、核验（scripts/verify_manuscript.py，n=26，执行后实测）

| 指标                                          | 结果                                                               |
| ------------------------------------------- | ---------------------------------------------------------------- |
| 摘要                                          | **250 / 250** ✓                                                  |
| Introduction                                | **460**（软限 <450 未达，因作者选择展开 WGCNA；该软限为 "with rare exceptions" 条款） |
| Methods / Results / Discussion / Conclusion | 300 / 1,576 / 1,102 / 47                                         |
| **正文合计**                                    | **3,485 / 3,500**（余量 **15**）✓                                    |
| 引用                                          | 1–26 全在、首现严格递增 ✓                                                 |
| T5 句 / em dash / P= 空格 / 禁用词                | 1 ✓ / 0 ✓ / 0 ✓ / 0 ✓                                            |
| 五组成对约束                                      | 全 True ✓                                                         |
| 主图 1–10、主表 1–3、标题 84 字符、关键词 5               | ✓                                                                |

**诚实更正**

1. **JA-01 的"零词数"仅在 JA-00 同步执行后成立**——未打补丁前脚本会把 Conclusion 的 47 词漏计（显示 3,438 而非 3,485），属假通过；已按 JA-00 修补。
2. **JA-04d（"top-ranked" → "top"）实测为 0 词**（连字符不影响 whitespace 分词），与初判 −1 不符；改用 JA-04e（删 "blood"）补上，摘要逐步实测回到 250。
3. **作者取舍后果**：选择"展开 WGCNA（+5）"而放弃"Intro 压词（−6）"，Introduction 由 455 变为 **460**，正文余量由 27 降为 **15**。二者互斥已在 45_ 中明示。
4. FA-06 丙案实测 +3（44_ 报价 +3，一致）。

### 五、同步更新

- `03_Manuscript_Translation_CN.md`：#37 译法 12 处（删"前人已打下相当的基础："、"经严格整编"→"经整编"、"以下约束"→"这些约束"、删"简言之"、"决定性"→"明确"、"另一番图景"→"不同结果"、"高于对照"→"高于对照组"、段2 断句、"与……相一致"、"契合将"→"与……作为……相一致"、"相称"→"相一致"、"若仍要将"）；#38 三处（171 个候选基因、补充表 S9、(18)/(19) 归属）；#39 八处（**新增"## 结论"章节**、WGCNA 全称、曲线下面积、最小可检测比值比 2.26、疾病与对照、两者均、删"血液区室"的"血液"）。
- `manuscript_revised.docx`：step95 重生成；新串 7/7 就位、旧串 0/9、斜体 run 99、参考文献 26 条、em dash 0、三条元数据修正（2026;123／13(52):e76236／12(1):385）齐备。

### 六、本轮结清与悬置

- **结清**：#37 全部 17 条推荐项、#38 A 类 6 项（FA-07 除外）、#39 结构项 JA-01/03/04/05/11 + 工具项 JA-00。
- **转为作者待办（期刊硬性，见 45_ 第五节 A/B 清单）**：A1 作者贡献、A2 利益冲突段、A3 Funding、A4 Acknowledgments、A5 伦理声明、A6 正文字数（**3,485**）、A7 影响陈述、**A8 AI 使用披露**（投稿系统＋标题页＋稿件内三处）、A9 GitHub URL + Zenodo DOI（占位符仍在）、A10 统计方法参考文献（补入将触发 6 处编号重排）、A11 图件 Alt text（`Figure_legends.md` 现为 0）、A12 补充材料 S 编号参考文献（现为 author-year）、A13 摘要页底标注词数 250；B1–B10 投稿信/系统项。
- **沿用悬置**：R-6（(14) CEBPB）、R-7（(25) *Vnn1* 敲除与 RR6 给药臂）、三项统计事实（Cox PH／KM 切点／BH 覆盖面）、(15) 作者数、HPA/Open Targets 查询日期、SMR-HEIDI 5 基因名单、GSE151263 供体明细。
- 余量 15 词未再分配。

---

## 第八批（2026-10-17）：#40 后续——E59 补回 + 启用 JA-02（推翻此前"不压 Intro"裁决）

**触发**：作者质询「重写后单细胞+双虚拟敲除／网络药理／分子对接是否被删除」。核查结论：三类内容**未删除**，为「搬家（Results 段20/段22 + 补充材料）＋ 压缩（Discussion 只留判定锚点）＋ 删技术细节（进 Supplementary）」，依据为 Path-A 重构方案 13_ 第五节逐段删留表。核查同时发现 3 项真丢失：MIF 通讯（E44）、五基因转移 AUC 0.961、E59 在 Discussion 侧缺失。  
**作者裁决**：三项≥全补 → 但 MIF 与 0.961 因事实未核实**暂缓**，本批仅执行 E59 + 支付杠杆 L2（JA-02）。


### 一、逐条改动

| #             | 级别 | 位置                                 | 原文 → 改文                                                                                     | 依据                                      |      Δ词 |
| ------------- | -- | ---------------------------------- | ------------------------------------------------------------------------------------------- | --------------------------------------- | ------: |
| **X-01**      | 🔴 | DISCUSSION 段7（VNN2 旁系句之后、CMap 句之前） | （无）→ 插入 "Vina and MM/GBSA disagreed on ligand ranking, so no ligand priority is supported." | #06 3.2 成对约束第④对（E57⇄E59 不得只引一侧）；A6 固定句式 | **+13** |
| **JA-02a**    | 🟡 | INTRODUCTION                       | "(ARDS) (2, 3), **a syndrome of** acute hypoxemic respiratory failure" → 删 3 词              | 45_ JA-02（作者本批推翻原裁决后启用）                 |      −3 |
| **JA-02b**    | 🟡 | INTRODUCTION                       | "is well characterized **in isolation**, from" → 删 2 词                                      | 同上                                      |      −2 |
| **JA-02c**    | 🟡 | INTRODUCTION                       | "The premise that **genuinely** shared genes" → 删 1 词                                       | 同上                                      |      −1 |
| 中文 X-01       | —  | 讨论段7                               | 插入"Vina与MM/GBSA在配体排序上不一致，因此不支持任何配体优先级。"                                                     | 英中同步                                    |       — |
| 中文 JA-02a/b/c | —  | 引言                                 | ARDS 定义句精简／删"在孤立状态下"／删"真正"                                                                  | 英中同步                                    |       — |

全部断言"唯一命中"；英文 4 处、中文 4 处均 1 命中后写入。

### 二、执行后实测（`scripts/verify_manuscript.py` n=26）

| 项                                           | 结果                                                                    |
| ------------------------------------------- | --------------------------------------------------------------------- |
| 摘要                                          | **250 / 250** ✓                                                       |
| Introduction                                | **454**（软限 <450 仍差 **4 词**——JA-02 设计时针对 455，叠加 JA-11 的 +5 后为 460→454） |
| Methods / Results / Discussion / Conclusion | 300 / 1,576 / **1,115** / 47                                          |
| **正文合计**                                    | **3,492 / 3,500**（余量 **8**）✓                                          |
| 引用                                          | 1–26 全在、首现严格递增 ✓                                                      |
| T5 句 / em dash / 禁用词                        | 1 ✓ / 0 ✓ / 0 ✓                                                       |
| 五组成对约束                                      | 全 True ✓                                                              |

### 三、诚实更正

1. **启用 L2 后 Introduction 仍为 454 词，未回到软限内**——JA-02 原设计把 455 压到 449，但 #39 的 JA-11（展开 WGCNA，+5）先发生，故现值 460−6=454。若要 ≤450，需再从 Introduction 删 4 词（另提方案，本批未做）。
2. **MIF 与 0.961 未入稿，且各有未决事实**：`08_singlecell/GSE151263_vs_GSE167363_comparison.md` 记「脓毒症状态下 MIF/CCL/TNF 主导」，与证据卡 E44 的「CellChat 以 MIF 主导（ARDS vs sepsis）」冲突 → MIF 归属【待作者核实】；0.961 在 05\_/06\_/08_ 文本输出中未检索到 → 出处【待作者核实】。二者核实前不补，避免把方向写反或把无出处数字写入稿件。
3. 本批未采用 L1/L3/L4–L7 去词（因 MIF、0.961 暂缓后不再需要），段8 首句与段7 首句维持原样。

### 四、同步更新

- `03_Manuscript_Translation_CN.md`：4 处（E59 中文句插入讨论段7；引言 ARDS 定义句精简、删"在孤立状态下"、删"真正"）。
- `manuscript_revised.docx`：step95 重生成；新串 5/5（Vina 句、ARDS 句、is well characterized, from、The premise that shared genes、CONCLUSION）、旧串 0/3、斜体 run 99、参考文献 26 条。

### 五、本轮结清与悬置

- **结清**：#06 成对约束第④对（E57⇄E59）由未过转为通过；#40 X-01 执行完毕。
- **暂缓（作者核实后再补）**：MIF 通讯句（甲案 +9／乙案 +13，取决于归属）、五基因转移 AUC 0.961 句（+12）。届时正文 3,492 + 9~~13 + 12 = 3,513~~3,517，**均超 3,500** → 补回前必须先落实去词（可用 L1 −6／L3 −3／L4–L7 各 −1，共 −13，仍不敷）。
- **沿用悬置**：R-6、R-7、三项统计事实、(15) 作者数、HPA/Open Targets 日期、SMR-HEIDI 5 基因名单、GSE151263 供体明细、FA-07、A1–A13 与 B1–B10 期刊项。
- 余量 **8** 词未再分配。

---

## 第十批（2026-10-19）：27 条修改意见全量执行 + netP 重分析 + 差异网络重算

**触发**：作者五项指令——(1) 阈值依赖性写成结果本身（+23 词方案）；(2) 核查 08_singlecell 的 CellChat netP$strength/weight，缺失则重分析并写入；(3) 用桌面 ARDS/Sepsis interactions CSV 重算差异网络；(4) 27 条全进正文，允许 +170~200 词；(5) Q-7/Q-8/Q-10 按杂志标准修改，不明确处用占位符。

### 一、分析层（先于写稿）

1. **netP 核查与重分析**：三份 rds（CellChat_ARDS/Sepsis/merged）的 netP$strength 与 netP$weight 均为 NULL（netP$prob 为 7x7x14 数组存在）。脚本 `scripts/step110_cellchat_netP_recalc.R` 从对象自身补算两槽位并写回 rds；导出 `08_singlecell/cellchat/netP_strength_weight_ARDS.csv`、`netP_strength_weight_Sepsis.csv`。
2. **差异网络重算**：桌面 CSV 与本地副本逐字节一致。通路级强度对比导出 `pathway_diff_network_recalc_2026-10-19.csv`。三种算法交叉一致（数组法 = LR 求和法 = 既有 pathway_strength_comparison.csv）。
3. **MIF 归属裁定（结清 Q-1）**：MIF 通路强度 ARDS 6.43e-4 vs 脓毒症 1.24e-3，脓毒症侧高 1.93 倍；MHC-I 与 MIF 在两侧均列前两名。原图 9B 图例 "enrichment of MIF signaling in ARDS" 无数值支撑（9B 实为细胞对级差分互作计数热图），已修正。
4. **文献外部核验（结清 Q-3）**：Lim N, Pavlidis P. Sci Rep. 2021;11(1):17624（PMID 34475469）；Kamimoto K, et al. Nature. 2023;614(7949):742-751（CellOracle）；Osorio D, et al. Patterns (N Y). 2022;3(3):100434（scTenifoldKnk）。

### 二、稿件逐条改动（02_Manuscript_Revised.md）

| 意见     | 位置                               | 改动                                                                                                        | Δ词（实测净）   |
| ------ | -------------------------------- | --------------------------------------------------------------------------------------------------------- | --------- |
| 11     | Results 免疫段末句                    | 阈值依赖写成结果本身（top-30 五条、top-50 七条、方向一致）                                                                      | +23       |
| 12+13  | Results 预后段                      | 补 HR 0.92-1.16、CI 全跨 1、名义 P 与 BH 校正 P=0.45-0.64                                                           | +19       |
| 6      | Results 外部确认段                    | 删 CD177 后五基因转移 AUC 0.961 对 0.984（S11 第 6 节）                                                               | +14       |
| 3(1-3) | Discussion 约束段                   | 替换「未正式评估」句：PH 全局违反 P=0.031（单项最小 EXOSC4 0.057）、KM 切点 P=0.00026-0.0019、ssGSEA 已补 BH                         | +30       |
| 5      | Results 单细胞段                     | 新增 CellChat 句：MHC-I/MIF 两侧前二；MIF 强度脓毒症 1.24e-3 高于 ARDS 6.43e-4（图 9B；图 S12）                                | +30       |
| 4      | Discussion 段1 / 结论               | 意义提升句（内型归属未经检验）；贡献陈述句（交集边界清晰、可测量、开放于实验检验）                                                                 | +10 +10   |
| 10     | Results 蛋白段                      | 13-fold → 13.6-fold                                                                                       | 0         |
| 14     | Methods 段1 + Results MR 段        | GSE151263 13,996 细胞（3+4 供体）、GSE167363 12,359；SMR 五基因名单显化（EXOSC4 无 cis-eQTL 除外）                            | +15 +5    |
| 16     | Discussion MR 段                  | 祖源声明：MR 数据源与 eQTLGen 参照均为欧洲祖源                                                                             | +12       |
| 18     | Methods 段1                       | GSE65682 入 ICU 时采血                                                                                        | +4        |
| 19/20  | Methods 段2 + Data Availability   | HPA 检索日期占位符；Open Targets 记 2026 年 9 月；ChEMBL 日期占位符                                                        | +10       |
| 23+7   | Methods/Results/Discussion + 文献表 | 21 条统计与方法学文献入稿（21-41 首现于 Methods/Results；47 Lim & Pavlidis 于 Discussion；原 21-26 改为 42-46、48）；文献总数 26 → 48 | +21 token |
| 去词杠杆   | Discussion                       | L1（约束段首句 -7）、L3（删 nevertheless -1）、L6（删 that is -2）、L7（删 therefore -1）                                    | -11       |

### 三、执行后实测（verify_manuscript.py，n=48）

- 摘要 250/250；Intro 454；Methods 347；Results 1,672；Discussion 1,157；Conclusion 57；**正文合计 3,687**（净 +195，在 +170~200 预算内）
- 引文 1-48 全引、无缺、首现严格升序；em dash 0；"P=" 空格 0；五组成对约束全 True
- 唯一 "robust" 命中来自 Newman CIBERSORT 论文标题，属文献条目不可改

### 四、同步更新

- `03_Manuscript_Translation_CN.md`：全部改动逐段同步（CellChat MIF 句、阈值依赖句、PH/切点句、祖源句、47/48 编号）。
- `manuscript_revised.docx`：step95 重生成（本批将 CONCLUSION 补入 HEADS 标题集，此前为普通段落）。
- `figures/Figure_legends.md`：图 9B 图例修正（删无支撑的 MIF-in-ARDS 措辞，写入重算强度并指向 S12/S25）；图 8G/H 加阈值依赖注；GSE151263 细胞数 14,000 → 13,996；追加 Figure 1-10 + Graphical Abstract 共 11 段 Alt text（意见 24）。
- `Supplementary_Materials.md`：S11 增第 6 节（CD177 消融）；S14 增 SMR 五基因名单；S15 增阈值稳健性明细；新增 S24（预后三项统计事实）、S25（CellChat netP 重算差异网络）、S26（队列元数据）；manifest 更新至 S1-S26。
- `Supplementary_Methods.md`：author-year 全部转为 S 编号书目（意见 25；新增 Supplementary References S1-S11，S4 为 LDSC 方法学篇 PMID 25642630 并注明与 atlas 篇的区分，意见 8）；S4 细胞数 13,996；S12 补 "approved or clinical-stage" 操作定义与 ChEMBL 日期占位；S11 补 HPA 版本/日期占位；S1 补采样时点说明与占位。
- 新文件 `51_Title_Page_CoverLetter_AI_Disclosure_2026-10-19.md`：标题页各段、ATS 影响陈述、AI 三重披露文本、投稿信草稿、系统勾选项核对表——占位符待作者填充（意见 21/22/27）。

### 五、结清与悬置

- **结清**：意见 3/4/5/6/7/8/9/10/11/12/13/14（细胞数部分）/16/17/18（GSE65682 部分）/23/24/25/26；Q-1（MIF，按重算证据写脓毒症侧）、Q-2（词数预算）、Q-3、Q-11、Q-12。
- **仍待作者**：Q-4 SMR 五基因名单签核（已按排除逻辑写入，需确认）；Q-5 四个 bulk 队列祖源/采样时点（S26 占位）；Q-7 HPA 版本与日期、ChEMBL 日期（稿件内占位符）；Q-8 标题页六段内容；Q-9 AI 披露措辞确认；Q-13 投稿信细节。
- **诚实说明**：意见 11 预期（top-30/50 下补体仍唯一）实测不成立，已按作者裁定把阈值依赖性写成结果本身；意见 5 用户初判（ARDS 侧）与三重独立重算证据相反，稿件按证据写为脓毒症侧，若另有差异度量输出可再议。

---

## 第十一批（2026-10-19）：节点⑦意义提升 / 节点⑧贡献陈述 深度升级（D1→D3）

**触发**：作者质询"意见4是否完成"，核查结论为第十批仅植入 D1 种子句（余量 5 词展不开）。产出计划册 `52_Node7-8_Depth_Upgrade_Plan_2026-10-19.md`，作者批准：A1 三钩子全要、B1–B7 删词全部放行、Conclusion 基因名补斜体。

### 一、新增（+70 词）

1. **A1 · Discussion 收束段（+64 词）**：约束段之后、CONCLUSION 之前新起一段，三个"So what"全部落地：
   - 分层变量：病对病 AUC=0.50 → 共享宿主反应状态标志，适合同时入组两病的试验分层，不是鉴别诊断工具；
   - 反对直接靶向：MR/SMR 阴性 → 直接治疗性靶向缺乏因果遗传学证据支持；
   - 平台泛化教训：转移 0.708 经 ComBat 恢复 0.977 → 限制此类血液标签的是平台泛化而非发现统计。
2. **A2 · Conclusion 末两句改写（+6 词）**：枚举交付物——预设定（prespecified）、外部与蛋白水平确认、与因果/预后主张切割；保留"是否成药靶点需实验验证"对冲句；*VNN1*、*UPP1* 补斜体。

### 二、对冲删词（B1–B7，共 −68 词）

| #  | 位置               | 删除内容                                                                                                                   | 依据                              |
| -- | ---------------- | ---------------------------------------------------------------------------------------------------------------------- | ------------------------------- |
| B1 | Results · ML 段末  | "so validation estimates carry a mild optimistic bias bounded by the seven-gene sensitivity analysis"（−14）             | 与 Methods 末段逐字重复，Methods 保留     |
| B2 | Discussion 段1    | "with M0 fractions read as myeloid enrichment and regulatory T-cell direction method-dependent and uninterpreted"（−14） | 两事 Results 免疫段已首次声明             |
| B3 | Results · 预后段末   | "These estimates are internal to a single cohort and require independent confirmation."（−12）                           | 与约束段重复                          |
| B4 | Results · 外部确认段末 | "so the six genes behave as shared disease-state markers rather than disease-specific discriminators"（−13）             | 解释性结论移至新 A1 段（Discussion 层级更恰当） |
| B5 | Discussion 段2    | 措辞压缩（would have entered...where / pre-set → none enters...: /）（−6）                                                     | 事实不动                            |
| B6 | Results · 单细胞段   | "given its cross-cohort origin"（−4）                                                                                    | 限定语 Discussion 段5 已有            |
| B7 | Introduction 段4  | "(weighted gene co-expression network analysis)" 展开与 "eleven machine-learning classifiers" 细节（−6）                      | 与 Methods 重复                    |

### 三、执行后实测（verify_manuscript.py，n=48）

- Intro 448 / Methods 347 / Results 1,629 / Discussion 1,201 / Conclusion 64；**正文合计 3,689**（净 +2，与计划预估一致，预算上限 3,692 内）
- 引文 1–48 升序无缺；em dash 0；T5 句计数 1（未破坏）；五组配对论断全 True（含病对病 0.50：Results 数值 + A1 段双保险）
- 分节核对：Results −43（B1+B3+B4+B6），Discussion +44（−14−6+64），Conclusion +7（A2），Intro −6（B7）——合计 +2 ✓

### 四、同步更新

- `03_Manuscript_Translation_CN.md`：九处逐段同步（含收束段与结论直译），头部说明更新为第十一批/3,689 词。
- `manuscript_revised.docx`：step95 重生成。
- `52_Node7-8_Depth_Upgrade_Plan_2026-10-19.md`：本批计划册（已执行完毕，可归档）。

### 五、遗留提示（非阻塞）

1. **WGCNA 首现展开**随 B7 从 Introduction 删除，全文现无展开式。若作者希望符合 AMA 首现展开惯例，可在 Methods "WGCNA (22)" 处补 "weighted gene co-expression network analysis (WGCNA)"（+5 词，正文将达 3,694，超预算上限 2 词；或以 A1 段 "shows, finally, that"→"shows that"（−1）等再对冲）。待作者定夺。
2. 正文 3,689 仍超 AJRCCM 刊定 3,500 词 189 词——第十批"+170~200"裁决的既定后果，本批净 +2 未扩大。

---

## 第十二批（2026-10-19）：C 组文献查证 → B 组裁决落稿 → D 组资产补齐 → 定稿

**触发**：作者 11 项指令（C 组文献查证优先；B5 不压缩按现状投；执行 D 组；整理建仓与 DOI 清单；伦理措辞与批件判断；关键词 MeSH 化；刷新 12_ 约束卡；D3 不构建；E 组两侧 CellChat 定量；B1–B4/B6–B8 逐项裁决；关键词再确认）。

### 一、C 组文献查证（5 条，PDF 全文 + 外部核验）

产出 `54_Reference_Verification_Report_2026-10-19.md`。

| 项  | 编号                | 结论                                                                                                                                         | 稿件处置                                       |
| -- | ----------------- | ------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------ |
| C1 | (14) 含 CEBPB      | ✅ 证实；六基因 Model-score 特征含 CEBPB，与严重度相关，n=1,265 多中心，Neu_CD177 亚群 Model score 与坏死性凋亡评分最高，qPCR CEBPB p<0.01                                    | 不改                                         |
| C2 | (46) 含 *Vnn1* 敲除臂 | ✅ 证实；VNN1⁻/⁻ 小鼠 + RR6 双臂，knockout 显著减轻肺损伤（PMID 41575600）                                                                                   | 不改                                         |
| C3 | (18)(19) 可否支撑归因   | ✅ 均成立；(18) 明写"recent studies have questioned the reliability"且为 cancer cell lines；(19) 明写 CMap 方法"agnostic to … multiple disease contexts" | 不改                                         |
| C4 | (15) 作者人数         | ✅ 41 位作者 → "Jung M, Beer A, Ghosh S, et al." 符合 AMA                                                                                        | 不改                                         |
| C5 | (16) 模型细节         | ✅ 有 shRNA 敲低与药理回补，但模型为肺腺癌（LUAD/TAM），且明确无条件敲除                                                                                               | **+5 词补 "in a lung adenocarcinoma model"** |

**关键说明**：原预案"C1/C2 若为'不含'则释放约 18 词以支付 B 组"未发生——两条承重引用均被原文坐实，无需删改。另如实记录：初次检索因 PDF 自动换行断词（"cell-  
lines"、"reliabil-  
ity"）而产生漏判，归一化后重检修正，以归一化结果为准。

### 二、E 组 · 两侧 CellChat 定量支撑（指令 9）

`08_singlecell/GSE167363/cellchat/` 下存在独立对象（CellChat_Sepsis / CellChat_HC / CellChat_merged），  
其 `netP$strength/weight` 同样为 NULL。新增脚本 `scripts/step111_cellchat_GSE167363_netP_recalc.R`  
从 `netP$prob`（7×7×35 数组）补算并写回 rds，导出 `netP_strength_weight_*.csv` 与  
`pathway_diff_network_GSE167363_recalc_2026-10-19.csv`。

结果（CD14+ 单核细胞，脓毒症 vs 健康）：

| 通路     | 脓毒症强度     | 健康强度      | 比值   | 主导侧 |
| ------ | --------- | --------- | ---- | --- |
| MIF    | 3.16×10⁻³ | 2.85×10⁻³ | 1.11 | 脓毒症 |
| MHC-II | 1.87×10⁻³ | 1.34×10⁻² | 0.14 | 健康  |
| CypA   | 9.76×10⁻⁴ | 5.88×10⁻⁵ | 1.66 | 脓毒症 |
| MHC-I  | 3.95×10⁻⁴ | 9.67×10⁻⁶ | 40.9 | 脓毒症 |


Results 单细胞段新增一句（+30 词），锚定 (Figure 9F; Supplementary Table S25)。  
至此两侧 CellChat（GSE151263 脓毒症 vs ARDS、GSE167363 脓毒症 vs 健康）均有定量支撑。

### 三、B 组裁决落稿（指令 10）

| 项  | 位置            | 改动                                                                                                    | Δ词  |
| -- | ------------- | ----------------------------------------------------------------------------------------------------- | --- |
| B1 | Results MR 段  | "Power analysis" → "**Post hoc** power analysis"                                                      | +2  |
| B2 | Discussion 段4 | 病对病句末补 "(AUC=0.50; n=7 patients)"                                                                     | +3  |
| B3 | Discussion 段6 | "transcriptomic endotype" 后补 (20) 锚点                                                                  | +1  |
| B4 | Discussion 段2 | 补筛选偏倚句 "Validation cohorts narrowed seven candidates to six, so those estimates carry mild optimism." | +13 |
| B6 | Introduction  | WGCNA 首现恢复展开 "(weighted gene co-expression network analysis)"                                         | +5  |
| B7 | 全稿            | ⑦⑧ 维持 D3；①–⑥ 深度维持                                                                                     | 0   |
| B8 | 全稿            | 三项口径（RQ4 混合判定／H2 后补定位／默认项）确认                                                                          | 0   |
| C5 | Discussion 段3 | 补 "in a lung adenocarcinoma model"                                                                    | +5  |

**B5（最大未决）已裁决**：正文不压缩、按现状投；字数标注刷新为 **3,747**；定稿后由作者统一刷新。

### 四、D 组资产（指令 3）

| 项   | 内容                                                    | 状态                                                          |
| --- | ----------------------------------------------------- | ----------------------------------------------------------- |
| D1a | S24 实体化：Cox PH 诊断 + KM 切点敏感性 + ssGSEA BH 覆盖面          | ✅ `TableS24_prognostic_statistical_checks.csv`（14 行，三区块）    |
| D1b | S25 实体化：两侧 CellChat 通路强度重算                            | ✅ `TableS25_CellChat_pathway_strength_recalc.csv`（57 行，两队列） |
| D1c | S26 占位：队列元数据，祖源/采样时点标 `[AUTHOR_INPUT_NEEDED]`         | ✅ `TableS26_cohort_metadata.csv`（10 行）                      |
| D2  | 图件格式补齐：主图 1/2/10、补充图 S2/S4/S5/S7…                     | ✅ PNG@300 dpi + TIFF@1000 dpi(LZW)；现**全部三件套齐全**             |
| D4  | `05_Renumbering_Maps.md` 追加附录：26→48 文献映射 + S24–S26 登记 | ✅                                                           |
| D5  | `Supplementary_Materials.md` manifest 更新至与实体一致        | ✅                                                           |
| D3  | Graphical Abstract                                    | ⛔ **按指令不构建**（Alt text 已备，待素材）                               |

汇编/导出脚本：`scripts/step112_export_missing_figure_formats.py`、  
`scripts/step113_assemble_supplementary_tables_S24_S25_S26.py`。

### 五、期刊标准项

| 项                       | 处置                                                                                                                                                                                   |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| 关键词（指令 6/11）            | 改为 MeSH 化且不与标题重复：**Transcriptome; Machine Learning; Mendelian Randomization Analysis; Single-Cell Analysis; Drug Repositioning**；备选 Respiratory Distress Syndrome (D012128) 已写入 51\_ |
| 伦理声明（指令 5）              | 新增独立小节 `ETHICS APPROVAL AND CONSENT TO PARTICIPATE`："Not applicable." + 二次分析／去标识／无新入组说明；并给出"不需要伦理批件、无需批件号"的判断与"waived"取舍理由                                                           |
| Data availability（指令 4） | 新建 `55_Code_Archive_and_DOI_Plan_2026-10-19.md`（英文）：归档目录、脚本清单、8 步建仓与 Zenodo 取 DOI 流程；归档已构建 `code_archive/` 53 MB、600 文件、`code_archive.zip` 18 MB、`MANIFEST.txt` 含逐文件 SHA-256         |
| 12_ 约束卡（指令 7）           | 刷新：标题 155→**84**、摘要 373→**249**、关键词 7→**5**、图表"未查得"→**官方 no limit**                                                                                                                  |
| D3（指令 8）                | 不构建 Graphical Abstract                                                                                                                                                               |

### 六、执行后实测（verify_manuscript.py，n=48）

- 摘要 **249**（≤250）；Intro **453** / Methods **347** / Results **1,660** / Discussion **1,223** / Conclusion **64**
- **正文合计 3,747**（Intro..Conclusion）；另有 ETHICS 63 词、DATA AVAILABILITY 109 词不计入
- 引文 1–48 全引、无缺、首现严格升序；em dash 0；"P=" 后空格 0；T5 句计数 1；五组配对论断全 True
- 唯一 "robust" 命中来自 Newman CIBERSORT 论文标题（文献条目，不可改）

### 七、同步更新

- `02_Manuscript_Revised.md`：C5 + B1 + B2 + B3 + B4 + B6 + E 组 + 伦理小节
- `03_Manuscript_Translation_CN.md`：逐段同步（含关键词、WGCNA 展开、GSE167363 CellChat 句、事后功效、段2/段4/段6、伦理小节、头部字数刷新）
- `manuscript_revised.docx`：step95 重生成
- `figures/Figure_legends.md`：Alt text 11 段已在第十批备齐
- `Supplementary_Materials.md` / `Supplementary_Methods.md`：S24–S26 登记、manifest 与实体一致
- 新建 `54_`（文献查证报告）、`55_`（建仓与 DOI 英文手册）
- `scripts/step109–113`：三项统计事实、两侧 netP 重算、图件格式补齐、S24–S26 汇编

### 八、仍待作者提供（提交前）

1. GitHub 仓库 URL 与 Zenodo **concept DOI**（按 `55_` 的 8 步执行；勿用 version DOI）
2. HPA Blood Atlas 版本与检索日期；ChEMBL release 与检索日期（稿件内仍为占位符）
3. 标题页六段：作者贡献 / 利益冲突 / Funding / Acknowledgments / 伦理备案确认 / 影响陈述
4. AI 披露措辞确认（草稿在 `51_`）
5. 投稿信细节：建议与回避审稿人、ORCID、预印本声明
6. S26 的 4 个 bulk 队列祖源与采样时点
7. SMR 五基因名单签核（已按"EXOSC4 无 cis-eQTL 故排除"写入，需确认）

### 九、口径声明

- C 组查证未引入任何新的 DOI/PMID 进正文；(15) 的 PMID 42325214 仅备查。
- 本批未新增任何数据、文献或论证强度；S24/S25 为既有产物的汇编，未重算科学结论  
  （唯一例外：`step111` 为 GSE167363 补算此前缺失的 `netP$strength/weight` 槽位，属对象修复而非新分析）。
- 正文 3,747 词超刊定 3,500 词 247 词，系作者 2026-10-19 "不压缩、按现状投"裁决的既定后果。

---

## 第十三批 · 定稿后复核（2026-10-20）

**触发**：第十二批 11 项定稿后，作者指令「定稿后我再统一刷新」之前的完整性复核。

### 1. 终检脚本口径修正（工具层，不影响稿件）

- **问题**：`scripts/verify_manuscript.py` 的风格检查（em dash / `P=` 后空格 / AI 风格词）作用于全文，误报 `robust` ×1。
- **定位**：命中位于**参考文献 34 的论文原题**（Newman AM, et al. *Robust enumeration of cell subsets from tissue expression profiles.* Nat Methods. 2015;12(5):453-457）。文献标题须逐字转录，改词即构成杜撰文献，触 #06 红线。
- **处置**：引入 `_style_text = text[:text.find("\nREFERENCES")]`，风格检查改为**仅扫正文（标题..DATA AVAILABILITY/ETHICS），排除 REFERENCES 段**，并在输出标题标注 `== style (body only; reference list excluded) ==`。
- **结果**：重跑后 `em dash 0` / `P= 后空格 0` / `banned AI words 全 0`。

### 2. 约束卡 12_ 过时值补刷（指令 7 的未竟部分）

`12_Journal_Constraints_Card_2026-10-16.md` 此前仅刷新了第二节，第三、五、七节仍留着重构前旧值，与第二节自相矛盾。已补完：

- **第三节**：旧值「正文 10,939 词」「摘要 373 词」标注 ~~作废~~ 并保留可追溯；补记实际处置（重构至 3,747 词，作者 2026-10-19 裁决不压缩按现状投；摘要已合规 249 词）。
- **新增 3.4 节**「现行合规状态一览（2026-10-19 终检）」：11 行对照表，含正文 3,747 ⚠️、摘要 249 ✅、参考文献 48 ✅、标题 84 ✅、关键词 5 ✅、图/表 ✅、em dash 0 ✅、AI 词 0 ✅（含扫描范围说明）、T5 唯一句 1 ✅、五组配对论断 5/5 ✅。
- **第五节**：8 项待确认改为表格，标注 1–5 已关闭、6–8 仍待。
- **第七节**：改为「唯一阻塞投稿的是第 5 组素材」8 项表（5.1 GitHub URL、5.2 Zenodo concept DOI、5.3 HPA 版本日期、5.4 ChEMBL release、5.5 S26 祖源时点、5.6 SMR 五基因签核、5.7 标题页六段、5.8 投稿信细节），每项附现稿占位与执行指引。
- **末尾**：加版本记录 v1/v2/v3。

### 3. 资产完整性复核（D 组）

- **图件三件套**：主图 Figure1/2/10_composite 与补充图 S2/S4/S5/S17 的 PDF + PNG@300 dpi + TIFF@1000 dpi(LZW) **全部 =1，齐备**。
- **补充表 S24–S26**：实体位于 `figures/supplementary_tables/`
  - `TableS24_prognostic_statistical_checks.csv`（15 行）
  - `TableS25_CellChat_pathway_strength_recalc.csv`（58 行）
  - `TableS26_cohort_metadata.csv`（11 行）
- **代码归档**：`code_archive/` 600 条目（analysis 212 / results 365 / manuscript 19 / environment 2），`code_archive.zip` 18 MB，MIT LICENSE 就位。

### 4. 四联时序

02\_（23:18）→ 03\_（00:01）→ 04\_（00:03）→ docx（00:04），同步无倒置。**本批仅改工具脚本与内部约束卡，稿件四联未动，docx 无需重生成。**

### 5. 本批未做的事（避免误报）

- 未改动 `02_Manuscript_Revised.md` 任何内容（robust 位于参考文献原题，属须逐字保留项）。
- 未执行建仓与 Zenodo 归档（按约定由作者执行 `55_` 的 8 步）。

---

## 第十四批 · 引言改写与润色（2026-10-21）

**范围**：仅 INTRODUCTION。Methods / Results / Discussion / Conclusion / 摘要 / 参考文献列表一律未动。  
**调用技能**：academic-writing-assistant、academic-humanization-low-ai、anti-defensive-writing（均按路径加载）。

### 六节点决策记录

**节点 1 · 文献综述切入（矛盾线）**  
候选线：①疾病分开研究 vs 共享机制；②既往筛选方法学缺陷；③单分子研究碎片化；④计算重定位无验证。  
选定：**①为主线，②③作支线**——"两病共享炎症机制却长期被分开研究"是本领域最根本且最易被审稿人认同的矛盾，②③ 具体化该矛盾，④ 收束于方法空白。原稿已按此组织，予以保留并压缩表述。

**节点 2 · 研究空白呈现（由矛盾反推）**  
从三条支线分别反推：共享机制已知 → 但筛选建立在未整编矩阵（脱靶样本）；共享基因已报道 → 但同数据内提名加验证；单分子已报道 → 但因果与假重复未处理、重定位继承癌细胞系局限。三段式保留，删去"First/Second/Third"的机械列举感，改为"Three gaps follow"后直接陈述，减少层级词。

**节点 3 · 研究问题聚焦**  
原问句含三层（再整编、外部验证、正交因果与单细胞检验），颗粒度合适，保留。删去原句尾部冗余的"and remain supported after"并列，压缩为"and survive re-curation, external validation, and orthogonal causal and single-cell testing"。

**节点 4 · 研究贡献界定**  
原贡献段"curated molecular intersection"四点列举保留，压缩为同位语结构。删除"well-supported"等自评形容词（anti-defensive），以"established sepsis markers"客观表述替代。

**节点 5 · 研究设计衔接**  
"Four orthogonal lines then tested their consistency..."保留四线名称（与 Methods/Results 严格对应），仅删"whether they behave"等虚词，未改动任何方法学指称。

**节点 6 · 完成**  
引言 454 → **349 词**（−105）。引用 1–48 全引、无缺号、首现严格升序；em dash 0；AI 风格词 0；T5 唯一句 1；五组配对论断 5/5。

### 逐段改动（英文）

| 段         | 原词  | 新词 | 主要处理                                                                                                                           |
| --------- | --- | -- | ------------------------------------------------------------------------------------------------------------------------------ |
| 1 临床与共享机制 | 88  | 58 | 删 ARDS 定义从句（方法与结果已充分界定）、"much of their inflammatory machinery"→直接列出效应器                                                         |
| 2 既往研究    | 105 | 84 | 拆为两句；删"clinical""biphasic septic response of""myeloid-derived suppressor-cell expansion""Independent work has converged on"等冗余 |
| 3 三个空白    | 84  | 79 | "Three limitations persist. First…Second…Third…"→"Three gaps follow."直接陈述；删"reference compendia"的 reference                    |
| 4 研究问题与设计 | 119 | 96 | 删"with cross-disease testing"展开、"pathway-level"的 -level、"quantified"→"against nulls"                                           |
| 5 贡献      | 56  | 32 | 删"surviving re-curation and""immune and regulatory""For CD177 and UPP1, the analysis extended…"→ 直接陈述                          |

### 人类化与反防御性处理

- 句长人为制造差异：段 2 首句缩至 5 词（"Shared-gene screens exist (6–8)."），后接长句。
- 删除"the analysis extended two well-supported sepsis markers"中的自评式修饰。
- 未使用"综上所述/值得注意的是/需要澄清的是"类模板句；未新增任何 hedge。
- 保留的必要限定：WGCNA 首现展开（B6 决策，AMA 惯例）、"confined to discovery data"（方法学范围，影响结论解读）。

### 字数结果（如实记录）

| 项                       | 改写前   | 改写后       |
| ----------------------- | ----- | --------- |
| Introduction            | 454   | **349**   |
| 正文合计（Intro..Conclusion） | 3,747 | **3,643** |
| 距 AJRCCM 上限 3,500       | 超 247 | **超 143** |

**结论：仅改引言无法达标。** 引言现有 349 词，要令正文落至 3,500 需再减 143 词，即引言须压到 206 词；而引言承载引用 1–19 的全文首现锚点（其中 1、2、4、5、17 为**唯一**锚点，正文他处不再出现），压到 206 词必然删掉大量文献述评，且违反该刊"Introduction 应给出背景与清晰研究理由"的要求。已试算极简电报式版本（240 词），对应正文仍为 3,533 词，依旧超限 33 词且质量不可用。**故本次停在质量底线 349 词，差额 143 词须由作者裁决处置**（见交付说明的三个方案）。

### 四联同步

02_ 英文稿 ✅｜03_ 中文稿 ✅（引言逐段同步 + 头部字数刷新为 3,643）｜manuscript_revised.docx ✅（step95 重生成）｜04_ 日志 ✅（本条目）。


---

## 第十六批 · Results/Discussion 精简至 3,500 以内（2026-10-21）

**触发**：作者指令「仍想压到 3,500，动 Results/Discussion 的 225 词，不能改动文章的内容与意义」。
**调用技能**：anti-defensive-writing、academic-humanization-low-ai、academic-writing-assistant。
**范围**：仅 RESULTS 与 DISCUSSION。Introduction（431 词）、Methods、Conclusion、摘要、参考文献列表一律未动。

### 削减原则（三选一，全部落在「不删信息」区间）

1. **删重复**：Results 与 Discussion 表述同一事实处，只保留一处；
2. **删防御性套语**：只用于防备假想质疑、不提供证据/范围/逻辑的句子（anti-defensive 检测清单命中项）；
3. **压虚词**：去冗余修饰、合并从句，不动任何数值、引用锚点、图/表锚点、阴性结果与结论。

### 逐段削减明细（英文）

| 段落 | 削减 | 主要处理 |
|---|---|---|
| R 共享候选基因 | −8 | "179 shared candidate genes"→"179 shared candidates"；删冗余修饰 |
| R GO/KEGG | −15 | 删末句 "This was a focused search…not a comprehensive survey."（纯免责，无证据增量） |
| R 机器学习 | −7 | 合并选择器列举；"label-permutation testing"→"label permutation" |
| R 表达模式 | −6 | 删与"all four cohorts"重复的 "consistently" |
| R 预测模型 | −8 | 删 "The sepsis transfer decrement was partly explained:" 引导语 |
| R 外部确认 | −10 | 压缩 meta 分析与病对病对比表述 |
| R 预后 | −23 | 删末句随机生存森林（数值 0.610 移入 Discussion D6，避免两处重复） |
| R MR | −8 | "weak germline genetic signal" 精简；删冗余引导 |
| R 免疫微环境 | −19 | 删 "estimates…reported for transparency only"（防御性套语）；删 "and was not interpreted" |
| R 单细胞 | −13 | 删末句 "all single-cell findings are mechanistic signals requiring confirmation"（D8 已述局限） |
| R 蛋白/治疗 | −28 | 删末句 "All findings in this section are computational predictions…"（D7 已述 hypothesis-generating）；压缩对接句 |
| D 段1 主要发现 | −24 | 删 ssGSEA 复现句（Results 已述）；删 "Whether this intersection refines endotype assignment is untested here." |
| D 段2 与既往筛选 | −7 | 删 "The prespecified standard was enforced:" 引导语（事实已由 *MUC1*/*DYRK2*/*TCTN1* 承载） |
| D 段3 注释关联 | −13 | 删末句（与 Results "annotation-based, not demonstrated co-regulated modules" 重复） |
| D 段4 两病关系 | −28 | 删 CellOracle 定位复述（Results 已述），改为直接接外部文献 (14) |
| D 段5 MR 解读 | −3 | 删 "together with"→"with"；删 "independent" |
| D 段6 预后 | −1 | 承接 Results 移入的 0.610 |
| D 段7 成药性 | −16 | 删 "no approved or clinical-stage small molecule exists…"（Results 已述）；删 "complete small-molecule tractability profile" |
| D 段8 局限 | −19 | 删段首套语 "The constraints concern scope."；压缩供体水平句 |
| D 段9 临床收束 | −4 | 压缩分层表述 |

**合计 −235 词**（Results 1,660→1,535，Discussion 1,223→1,113）。

### 校验结果（全部通过）

- 正文合计 **3,490** 词（Intro 431 / Methods 347 / Results 1,535 / Discussion 1,113 / Conclusion 64），**≤3,500 ✅，余量 10 词**
- 引用 1–48 全引、无缺号、首现严格升序 True
- em dash 0；`P=` 后空格 0；AI 风格词全 0；T5 唯一句 1；五组配对论断 5/5 True
- 主图 1–10、主表 1–3 锚点完整

### 过程中自检发现并修复的两个失误（如实记录）

1. 首轮把 D5 的 "**does** constrain causal effects of ARDS on sepsis risk" 误删为 "constrains"，致 T5 唯一句计数归零 → 已恢复。
2. 首轮把 R12 的 "combined **six-gene signature** AUC=0.50" 误缩为 "combined signature AUC=0.50"，致配对论断第 5 项转 False → 已恢复。
两处恢复共 +3 词，另从 D1/D2/D4/D5/R22 再压 8 词，最终 3,490。

### 未删项（守住的红线）

- 全部数值与统计量（557/1,458、Zsummary、MDOR 2.26/1.04、AUC 0.753–0.998、0.961/0.984、PH P=0.031、切点 0.00026–0.0019、BH 0.45–0.64、HR 0.64、C 指数、rho、13.6-fold、log2FC 6.68 等）
- 全部阴性结果（MR 双向零结果、单基因预后不显著、病对病 AUC=0.50、下调 CMap 无糖皮质激素、一倍标准误未保留基因）
- 全部引用锚点（尤其 (38)–(48) 等唯一锚点）与图/表锚点
- 结论与三个 "So what"（分层变量 / 反对直接靶向 / 平台泛化）

### ⚠️ 提交前风险提示

Methods 中仍有 `[AUTHOR_INPUT_NEEDED]` 占位（HPA 版本与检索日期、ChEMBL release 与检索日期），回填后正文将增加约 15–25 词，**可能突破 3,500**。当前余量仅 10 词。建议回填后重跑 `scripts/verify_manuscript.py` 复核；若超限，从 Methods 或本批保留的冗余处再压。

### 四联同步

02_ 英文稿 ✅｜03_ 中文稿 ✅（结果讨论逐段同步 + 头部字数刷新为 3,490）｜manuscript_revised.docx ✅（step95 重生成）｜04_ 日志 ✅（本条目）。


---

## 第十七批 · HPA/ChEMBL 版本与检索日期回填（2026-10-21）

**触发**：作者确认按修改表回填两项数据源版本与检索日期。

### 查证依据（已另行记入工作日志）

- **HPA**：本地 8 个 TSV 下载于 2025-11-06（mtime 与 zip 内嵌时间戳互证）；数据含 59 疾病（10 感染类），对应 **HPA v24**（2024-10-22 发布的 Disease Blood Atlas，59 疾病 Olink PEA）；v25 于 2025-11-11 发布（71 疾病），晚于下载日且疾病数不符，排除。
- **ChEMBL**：chembl_known_drugs.csv 创建于 2026-09-17 00:28，当日工作日志记录「ChEMBL 六靶点均无临床期药物」；2026-09-17 时线上版本为 **ChEMBL 37**（2026-05 发布，38 未发布）。

### 回填内容（三处）

1. Methods：Human Protein Atlas Blood Atlas **(accessed November 6, 2025)**
2. Methods：ChEMBL **(release 37, accessed September 17, 2026)**
3. DATA AVAILABILITY：Human Protein Atlas Blood Atlas **(version 24, accessed November 6, 2025)**

### 执行说明（如实记录）

作者于 23:32 先行手动完成了 02_ 英文稿三处回填（与本批修改表一致）；本批完成 03_ 中文稿三处同步、docx 重生成、头部字数标注刷新。

### 校验结果

- 正文 **3,494** 词（+4，≤3,500，余量 6）；摘要 249 ✅
- 引用 1–48 全引、无缺号、首现严格升序 True
- em dash 0；AI 风格词全 0；T5 唯一句 1；五组配对论断 5/5 True
- 02_ 与 docx 中 `[AUTHOR_INPUT_NEEDED]` 仅剩 GitHub URL/DOI 一处（待建仓回填）

### 四联同步

02_ 英文稿（作者回填）✅｜03_ 中文稿 ✅｜manuscript_revised.docx ✅｜04_ 日志 ✅（本条目）。

### 遗留

- `chembl_known_drugs.csv` 无归档生成脚本（ad-hoc 查询的溯源缺口，结果 CSV 已入 code_archive）；如需补齐可写 step115 固化查询。

---

## 第十八批 · GitHub 建仓 + Zenodo DOI 正式回填（2026-10-24 / 10-26 落稿）

**触发**：作者完成 GitHub 账号认证与仓库创建，取得 Zenodo concept DOI `10.5281/zenodo.23668532`。这是 E-J4（Data Availability 须含可访问代码库与持久 DOI）的最后一项阻塞素材。

### 建仓与认证过程（如实记录）

1. **账号确认**：GitHub 连接器已授权，`get_me` 返回登录名 **wuds2006**（id 107284368），因此无需浏览器输入密码或设备码。
2. **MCP 建库失败（403）**：`create_repository` 返回 `Resource not accessible by integration` —— 连接器 token 无建库/写权限，故不能依赖 MCP 建库或推文件，改用「网页建库 + SSH push」。
3. **SSH 密钥**：生成 `~/.ssh/id_ed25519_github`（ed25519、无密码、指纹 `SHA256:8opYZMWhOGBchkj2SWDga9Yg/zq9+PkZFjHWntupGTI`），新建 `~/.ssh/config` 把 `Host github.com` 绑定到该密钥。
4. **授权踩坑**：首次公钥被添加为 **Signing Key**（仅用于提交签名，不能用于 push 认证），`ssh -T` 报 `Permission denied (publickey)`；改为 **Authentication Key** 后通过，返回 `Hi wuds2006!`。
5. **push 踩坑**：后台命令运行在沙箱内，`~/.ssh` 被拦截读取，导致 git push 卡在 host key 确认提示；改为前台执行 + 显式 `GIT_SSH_COMMAND`（指定 `-i` 私钥与 `UserKnownHostsFile`）后一次成功。

**结果**：仓库 <https://github.com/wuds2006/sepsis-ards-shared-hub-genes> 已上线，分支 `main` 位于 `f910022`，613 个文件，annotated tag `v1.0.0` 已推送（GitHub API 验证 `refs/tags/v1.0.0` 在线）。

### Zenodo 网络故障排查（2026-10-26）

作者报「Zenodo 连不上」。诊断结论：**不是 Zenodo 故障，也不是账号问题，而是本地 DNS 污染**。

| 检查项 | 结果 |
|---|---|
| `curl https://zenodo.org` | HTTP 000（0.15 s 立即失败） |
| `nslookup zenodo.org`（路由器 192.168.0.1） | **Address: 0.0.0.0** ← 根因 |
| 外部 DNS（223.5.5.5 / 119.29.29.29 / 114 / 1.1.1.1）及 Google、Cloudflare DoH | 全部超时或同样被污染 |
| `curl --noproxy '*' --resolve "zenodo.org:443:188.185.42.113"` | **HTTP 200** ← IP 层从未被封 |

真实 IP 经阿里 DoH 的 JSON 端点取得：
`curl -s "https://dns.alidns.com/resolve?name=zenodo.org&type=A"`（注意是 `/resolve`，不是 `/dns-query`）。

**修复**：经作者同意，用 `osascript ... with administrator privileges` 弹系统授权框，向 `/etc/hosts` 追加三行：

```
# Zenodo direct resolution - local router DNS returns 0.0.0.0 for zenodo.org
188.185.42.113 zenodo.org
188.185.45.63 zenodo.org
137.138.121.78 zenodo.org
```

并 `dscacheutil -flushcache && killall -HUP mDNSResponder`。原文件备份于 `/tmp/hosts.bak.<时间戳>`，随时可撤。
复查：模拟浏览器（无代理变量）访问 `https://zenodo.org/` 返回 HTTP 200（5.8 s）。

### 稿件回填内容（DATA AVAILABILITY，一处）

原文：

> All analysis scripts are available in a public GitHub repository and archived with a persistent DOI on Zenodo **[AUTHOR_INPUT_NEEDED: URL and DOI]**.

改为：

> All analysis scripts are available at https://github.com/wuds2006/sepsis-ards-shared-hub-genes and archived at Zenodo (https://doi.org/10.5281/zenodo.23668532).

该句比原句**少 10 词**，且 DATA AVAILABILITY 按字数规则本就不计入正文，故正文仍为 **3,494** 词，余量 6。

### 配套资产

- `CITATION.cff` / `.zenodo.json`：作者 **Wu, Dansen**，ORCID **0000-0002-8429-168X**，`repository-code` 与 `url` 均为上述 GitHub 地址。
- 桌面「建仓材料」：`QUICKSTART.md` 改写为「已登录/已推送版」（11 项状态板、两次浏览器点击、push/tag/release、Zenodo 五步、三条 fallback），并新增 **5b 节**记录本次 DNS 污染的完整诊断与修复方法。
- 代码主管线 `analysis/run_all.sh`（17 stages = S1–S14 + 校验）与 `RESULTS_TO_CODE.md` 逐数字映射表均已入档。

### 校验结果（step95 后重跑 verify_manuscript.py）

- 正文 **3,494** 词（≤3,500）✅；摘要 249 ✅；标题 84 字符 ✅；关键词 5 ✅
- 引用首现严格升序 True；正文引用 22 条无缺号
- em dash 0；AI 风格词全 0；T5 唯一句 1；五组配对论断 5/5 True
- **稿件内 `[AUTHOR_INPUT_NEEDED]` 已清零**（GitHub URL + Zenodo DOI 为最后一处）

### 四联同步

02_ 英文稿 ✅｜03_ 中文稿 ✅｜manuscript_revised.docx ✅｜04_ 日志 ✅（本条目）。

### 遗留（不影响投稿，但建议补齐）

- 标题页六段内容（作者姓名/单位/ORCID、通讯作者、贡献、利益冲突、基金、致谢）
- 投稿信细节（建议/回避审稿人、ORCID、预印本声明）
- S26 四个 bulk 队列的祖源与采样时点（现为占位）
- SMR 五基因名单签核
- Zenodo 侧：若 concept DOI 页面元数据有误，可在 Zenodo 网页直接编辑（DOI 本身不可改）

---

## 第十九批 · 参考文献时效性与替换（2026-10-29）

**触发**：作者复核 docx 参考文献，列出 23 条发表较早的文献，要求判断是否为「经典 + 高影响 + 密切相关」；
若否，检索 PubMed 近 5 年真实高质量文献并替换，明确禁止虚假替换。

### 方法

- PubMed E-utilities（esearch/esummary/efetch），日期过滤 2021–2026；
- 每条候选再用 Crossref `api.crossref.org/works/<DOI>` 核验卷期页与作者。
- 发现 NCBI `efetch` 的 DOI 字段存在错配（clusterProfiler 4.0 被标为 `10.1109/TCBB.2019.2956708`、
  Sanderson 2022 被标为 `10.1038/s41467-019-12424-x`），故**一律以 Crossref 为准**。

### 替换结果（4 条，编号不变，总数仍 48）

| 编号 | 原文献 | 新文献 | PMID / DOI |
|---|---|---|---|
| 10 | Hotchkiss 2013 *Nat Rev Immunol* | van der Poll T, Shankar-Hari M, Wiersinga WJ. The immunology of sepsis. *Immunity*. 2021;54(11):2450–2464 | 34758337 / 10.1016/j.immuni.2021.10.012 |
| 11 | Landelle 2010 *Intensive Care Med* | Monneret G, et al. Monocyte HLA-DR expression in septic shock patients: insights from a 20-year real-world cohort of 1023 cases. *Intensive Care Med*. 2025;51(10):1820–1832 | 40986015 / 10.1007/s00134-025-08110-w |
| 23 | Yu 2012 *OMICS*（clusterProfiler 1.0） | Wu T, et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *Innovation (Camb)*. 2021;2(3):100141 | 34557778 / 10.1016/j.xinn.2021.100141 |
| 30 | Smith & Ebrahim 2003 *Int J Epidemiol* | Sanderson E, Glymour MM, Holmes MV, et al. Mendelian randomization. *Nat Rev Methods Primers*. 2022;2:6 | 37325194 / 10.1038/s43586-021-00092-5 |

编号不变的原因：4 条新文献的**首现位置与原文献完全相同**（10、11 在 Introduction 的 (9–12)；
23、30 在 METHODS 同一处），故正文引用编号无需重排。

### 判定为「不可替代」的 18 条及理由要点

- **方法学/软件/算法原始出处（15 条）**：21 limma、22 WGCNA、24 Random forests、26 glmnet、27 Cox、
  28 Cox-LASSO、29 RSF、31 LDSC、32 SMR/GSMR、33 COLOC、34 CIBERSORT、35 GSVA、37 BH、
  38 ComBat、41 DESeq2 —— 引用方法原始文献是出版规范，替换为应用/综述类文献会被判引用不当。
- **论断或数据集绑定（3 条）**：
  - 9 Boomer 2011 *JAMA*：免疫抑制的原创人体证据（淋巴细胞凋亡），换成综述会使论证从原始证据降级为转述；
  - 20 Davenport 2016 *Lancet Respir Med*：GSE65682 的 Mars1–Mars4 endotype 定义文献；
  - 45 Schalkwijk 2014：正文"RR6 消除大鼠血浆 vanin 活性"论断的唯一支撑，近 5 年 VNN1 抑制剂文献
    报道的是其他化合物，替换即构成虚假引用。

### 留档的未采用候选（6 条，均经双源核验）

Giamarellos-Bourboulis 2024 *Nat Immunol*（38168953）、Yao 2023 *Mil Med Res*（37337301）、
Skrivankova 2021 STROBE-MR *JAMA*（34698778）、Guan 2023 VNN1（37608823）、
Meng 2026 VNN1-ROS-NETs（42473476）、Lv 2025 VNN1 综述（41090263）。

### 产出与校验

- 新增评估报告：`56_Reference_Currency_Review_2026-10-29.md`（23 条逐条判定 + 双源核验记录 + 备选池）。
- 参考文献仍 **48 条**（≤50）；2021 年前 29 条中 23 条为方法学/定义性文献。
- `verify_manuscript.py` 重跑：正文 **3,494** 词 ✅；摘要 249 ✅；标题 84 字符 ✅；关键词 5 ✅；
  引用首现严格升序 True；em dash 0；AI 风格词 0；T5 唯一句 1；配对论断 5/5 True。

### 四联同步

02_ 英文稿 ✅｜03_ 中文稿（无参考文献列表，编号未变，无需改动）✅｜
manuscript_revised.docx（step95 重生成）✅｜04_ 日志 ✅（本条目）

---

## 第二十批 · AJRCCM 相关文献（12 篇 PDF）引用评估与纳入（2026-10-30）

**触发**：作者提供 12 篇 AJRCCM 及相关杂志文献，要求按 DISCUSSION 实际论断判断可否引用，
核验真实性/支持关系/论断强度后引用，总数控制在 50 以内。

### 核验方法

- 12 篇 PDF 文本提取（pypdf）后逐篇阅读；正刊论文再经 PubMed + Crossref 双源核验卷期页与 DOI。
- 三步判定：真实性 → 能否支撑 Discussion 某一句具体论断 → 证据等级决定措辞强度（一律不写成因果）。
- 起点为 48 条，**净增额度只有 2 条**。

### 已引用（2 条，均为 AJRCCM 正刊原创）

| 编号 | 文献 | 支撑的论断 | 核验 |
|---|---|---|---|
| **43** | Aguilar OA, et al. *MICB* genomic variant is associated with NKG2D-mediated acute lung injury and death. *AJRCCM* 2024;209(1):70–82 | DISCUSSION 首段末："血液 NK 减少不排除 NK 在肺内的效应作用"；原文结论为 "support a role for"，故措辞用 "implicates ... in lung injury" | PMID 37878820 / 10.1164/rccm.202303-0472OC，双源一致 |
| **44** | Sinha P, et al. Molecular phenotypes of ARDS in the ROSE trial. *AJRCCM* 2024;209(7):816–828 | DISCUSSION 第 5 段："ARDS 本身含转录组轨迹各异的分子表型"；因其表型由临床+蛋白 LCA 定义，措辞未写成"转录组定义表型" | PMID 38345571 / 10.1164/rccm.202308-1490OC，双源一致 |

编号重排：新 43、44；原 43–48 → **45–50**（Burnham、Jones、Schalkwijk、Liu、Lim、Antcliffe），英文稿与中文稿同步。

### 未引用（10 篇）及理由

- **正刊 2 篇（额度已满，列为首选备用）**：Kiernan E 2025;211(9):1681–1688（脓毒症亚表型对液体复苏反应不同，
  可支撑"试验分层"论断，建议替换参考文献 17 或 8 后纳入）；Matthay MA 2024;209(1):37–47（新全球 ARDS 定义，
  与本文 Berlin 入组口径冲突，仅在局限性段可用）。
- **会议摘要 8 篇**：证据等级不足（无同行评议全文、样本量小），不用于支撑论断；其中
  Jones A5372（脓毒症相关 ARDS 的 STAT3/IL-6 遗传通路，与本文 IL6/JAK-STAT3 富集直接契合）与
  Alipanah-Lechner A3131（多组学重现 ARDS 亚表型、线粒体功能障碍）**正式发表后应优先纳入**；
  Evrard A4576 的论点已由 Sinha 2024 覆盖；ILD/RIARDS/肌钙蛋白 3 篇与本文论断无接口。
- 完整逐篇判定见 `57_AJRCCM_PDF_Reference_Assessment_2026-10-30.md`。

### 字数处理

新增两句后正文 **3,517 词**；经 13 处不改变意义的精简（"in this gene set"、"alongside platform scaling"、
"is unknown here"、"rather than a transcriptionally distinct entity" 等冗余限定语）回到 **3,499 词**（≤3,500）。

### 校验

参考文献 **50 条**（达上限）；正文引用首现**严格升序 True**；50 条全部被引用、无孤儿条目；
摘要 249；标题 84 字符；关键词 5；em dash 0；AI 风格词 0；T5 唯一句 1；配对论断 5/5 True。

### 四联同步

02_ 英文稿 ✅｜03_ 中文稿（新增两句译文 + 编号映射 43→45 … 48→50）✅｜
manuscript_revised.docx ✅｜04_ 日志 ✅（本条目）

---

## 第二十一批 · 备选池替换与 Kiernan 2025 纳入（2026-10-30）

在参考文献已满 50 条的前提下，按 2026-10-29 评估留下的备选池执行替换，总数仍为 50。

### 五处处置

| 编号 | 处置 | 文献 |
|---|---|---|
| 8 | 替换旧 9（Boomer 2011 *JAMA*） | Giamarellos-Bourboulis EJ, et al. *Nat Immunol* 2024;25(1):19–28（PMID 38168953） |
| 12 | 新增 | Yao RQ, et al. *Mil Med Res* 2023;10(1):27（PMID 37337301，单细胞 HLA-DR^low S100A^high 单核细胞） |
| 17 | 替换旧 17（Chen H *Medicine*） | Kiernan E, et al. *AJRCCM* 2025;211(9):1681–1688（PMID 40668865，SP1/SP2 液体复苏反应不同） |
| 30 | 替换旧 30（Sanderson 2022） | Skrivankova VW, et al. STROBE-MR. *JAMA* 2021;326(16):1614–1621（PMID 34698778） |
| 48 | 替换旧 48（Liu C *Inflammation*） | Meng S, et al. *J Inflamm Res* 2026;19:591174（PMID 42473476，VNN1-ROS-NETs 脓毒症肺损伤） |

删除旧 8（Liu T *J Inflamm Res* 2025）以腾出额度；(6–8) 改为 (6, 7)。
编号映射：旧 9→8、10→9、11→10、12→11，13 及以后不变。

### 论断边界（避免虚假引用）

- Yao 2023 用 resolves，不作因果断言。
- STROBE-MR 为报告规范，故 METHODS 措辞改为 "following STROBE-MR (30)"，不冒充方法学出处。
- Kiernan 2025 的亚表型由**血浆蛋白**模型（Ang-1/2 + sTNFR-1）定义，故只写"反应不同"，未写成转录组定义亚表型；数据：SP2 41% vs 27%、SP1 9% vs 9%、交互 P=0.02。
- Meng 2026 未指明抑制剂种类，故不把 RR6 与该研究绑定；Schalkwijk 2014（新 47）保留。
- 已知代价：Boomer 2011 的原创人体证据被权威综述取代（用户指定的取舍）。

### 字数

3,499 → 替换与新增后 3,516 → 11 处精简 3,500 → 恢复 2 处配对论断限定语 → 再 5 处微调 **3,495**（目标区间 3,490–3,497）。
注意：`verify_manuscript.py` 的守卫要求"true absence of causality"与"established causal targets"原文保留，误删会使配对论断转 False，已恢复。

### 校验

正文 **3,495**；摘要 249；参考文献 **50 条**；首现严格升序 True、50 条全部引用无孤儿；
主图 1–10 / 表 1–3 全引用；em dash 0；AI 词 0；T5 唯一句 1；配对论断 5/5 True。

### 四联同步

02_ 英文稿 ✅｜03_ 中文稿（编号映射 + 三处新译文 + 删除"两种读出重现了该表型"）✅｜
manuscript_revised.docx（step95 重生成，新 5 篇在、旧 3 篇已清，条目 50）✅｜04_ 日志 ✅（本条目）

完整核验与理由见 `58_Reference_Swap_Round2_2026-10-30.md`。

---

## 第二十二批（2026-10-30）：第 7 条反向替换——恢复 Boomer 2011（淋巴细胞凋亡的原创人体证据）

### 触发

作者决定撤回第二十一批对 8 号位的取舍：Boomer 2011 *JAMA* 是"脓毒症淋巴细胞凋亡 + 单核细胞 HLA-DR 下调"的**原创人体证据**（脾与肺组织），
权威综述不可替代其出处地位；故以第 7 条（Wei 2025 *Front Mol Biosci*）为其腾位。

### 处置（1 处替换，编号不动）

| 编号 | 处置 | 文献 |
|---|---|---|
| 7 | 替换（Wei J, et al. *Front Mol Biosci* 2025;12:1665387 移出） | Boomer JS, To K, Chang KC, Takasu O, Osborne DF, Walton AH, et al. Immunosuppression in patients who die of sepsis and multiple organ failure. *JAMA*. 2011;306(23):2594–2605. doi:10.1001/jama.2011.1829（PMID 22187279） |

元数据双源核验：PubMed efetch（PMID 22187279）与 Crossref（10.1001/jama.2011.1829）一致；
注意 10.1001/jama.2011.**1828** 在 Crossref 上指向另一篇（Needlestick Injuries），**不可用**。

### 正文连带改动（论断与编号对齐）

- INTRODUCTION："Shared-gene screens spanning the two diseases have been reported (6, 7)." → "(6)"（共享基因筛选现仅 Jiang 2024 一条，与 DISCUSSION 中"（6）提名 CX3CR1/PID1/PTGDS"口径一致）。
- INTRODUCTION 与 DISCUSSION 两处免疫抑制论断："(8–11)" → "(7–11)"，Boomer 2011 居该组首位，即组内最早的原创人体证据。
- 首现顺序仍为 6 → 7（同句内 6 先于 7），严格升序 True。

### 论断边界

- Boomer 2011 为人体脾/肺组织的观察性证据，故仍限定为"以淋巴细胞凋亡、MDSC 扩增与单核细胞 HLA-DR 下调为特征"，
  未改写成因果或预后断言。
- 第 7 条不再承担"共享基因筛选"论断——该论断归 (6)；Boomer 未被用于任何它不支持的句子。

### 字数

3,495 → **3,494**（"(6, 7)"→"(6)" 减 1 词；参考文献列表不计入正文）。区间 3,490–3,497，余量 6 词。

### 校验

正文 **3,494**；摘要 249；参考文献 **50 条**；全部引用无孤儿、首现严格升序 True；
主图 1–10 / 表 1–3 全引用；em dash 0；AI 词 0；T5 唯一句 1；配对论断 5/5 True。

### 四联同步

02_ 英文稿 ✅｜03_ 中文稿（"（6, 7）"→"（6）"、两处"（8-11）"→"（7-11）"，并在引言补"最早在人体脾与肺组织中得到证明（7）"；头部数据说明更新为 3,494 词 / 50 条）✅｜
manuscript_revised.docx（step95 重生成；Boomer 在、Wei J 已清，条目 50）✅｜
投稿包 `~/Desktop/投稿材料`（step114 + safe-copy shim 重建：主图 10、补充图 21、补充表 27、内嵌图 10、正文表 3，包内 50 条）✅｜
桌面建仓仓库 `~/Desktop/建仓材料/repo/manuscript/02_Manuscript_Revised.md` 已同步（未提交，待作者决定是否 commit/push）✅

完整核验与理由见 `59_Reference_Swap_Round3_Boomer_2026-10-30.md`。

---

## 第二十三批（2026-10-30）：作者信息与声明回填 ＋ 防御性写作第二轮清理 ＋ B 方案加词 ＋ 审稿人建议

### 触发

作者提供完整的作者署名、单位、ORCID、邮箱、基金、CRediT、利益冲突、AI 披露、Impact statement 与通讯作者信息，
并要求：① 据此更新稿件；② 用 anti-defensive-writing 规则对全文做第二轮防御性写作审查；③ 采用 B 方案加词（+8，不再为字数精简）；
④ 腾位候选更新为第 6 条 Jiang 2024、Wei 2025 退回备选池；⑤ 依据研究内容与所引 *AJRCCM* 论文作者群推荐审稿人。

### 一、作者信息与声明回填（02_ 声明小节 ＋ 51_ 标题页）

| 项 | 处置 |
|---|---|
| AUTHOR CONTRIBUTIONS（CRediT） | 按作者原文逐句写入 02_ `AUTHOR CONTRIBUTIONS` 与 51_ 标题页（含串行 Oxford comma 统一） |
| FUNDING | 写入 02_ `FUNDING` 与 51_：三项基金（0060112505 → Xiuling Shang；2019HSJJ09、{2025}881 → Dansen Wu） |
| CONFLICT OF INTEREST | 写入 02_ 与 51_："The authors reported no relevant financial relationships with ineligible companies." |
| AI USAGE DISCLOSURE | 三处统一为同一句（投稿系统勾选 ↔ 标题页 ↔ 稿件内小节） |
| Impact statement | 写入 51_ 投稿系统字段（作者原文，未改） |
| 作者与单位 | 7 位作者（Huobao Yang / Yanghui Shen / Wenjun Chen / Jiali Wu / Zhimin Wu / Qihua Wu / Dansen Wu），通讯作者置末位；单位 1–4；ORCID 与邮箱全部填入 51_ |
| 投稿信 | `step114` 中 (5)(6)(7) 条与落款改为实名实单位（Dansen Wu、ORCID、邮箱），(7) 指向 `61_` |

**⚠️ P0 一致性问题（未擅自改）**：CRediT 写 "Y.S., W.C. and J.W. contributed to performing the experiments"，
与稿件 `ETHICS APPROVAL AND CONSENT TO PARTICIPATE` 中 "no human or animal experiments were performed by the authors" 直接矛盾。
建议改为 "contributed to data curation, formal analysis, and development of methodology"。待作者裁决。
另：作者提供的字数 3,747（Introduction 453 / Methods 347 / Results 1,660 / Discussion 1,223 / Conclusion 64）为压缩前口径，
与当前稿 3,486 不一致，51_ 已按当前稿改写并说明。

### 二、B 方案加词（已执行）

DISCUSSION 免疫抑制论断后加：**"first shown in human spleen and lung tissue (7)"**（+8 词），
点名 Boomer 2011 为淋巴细胞凋亡的原创人体证据。作者裁决"不再精简"，故未因超限回删。

### 三、防御性写作第二轮清理（A 类 6 处，−17 词）

| # | 位置 | 删改 | 依据 |
|---|---|---|---|
| A1 | DISCUSSION | "…is threshold- and cohort-dependent, **not a refutation**" → "…**reflects thresholds and cohorts**" | 预防"是否定前人"的假想质疑 |
| A2 | DISCUSSION | "Validation cohorts narrowed seven candidates to six, so those estimates carry mild optimism." 整句删除 | 与 METHODS 末句同一限定重复 |
| A3 | DISCUSSION | "**Two caveats limit inference:**" 删标签 | 先示警再举证的免责框架 |
| A4 | RESULTS | "screening did not inflate performance" → "the screening step left performance unchanged" | "not X"式防御 → 中性陈述 |
| A5 | DISCUSSION | "Direct therapeutic targeting is unsupported…" → "Causal-genetic evidence supports their use as markers, not as therapeutic targets." | 负面框架 → 正面说明证据支持什么 |
| A6 | RESULTS | 删冗余转折 "but" | 前后为并列证据 |

B 类 12 处（真实选择偏倚边界、过拟合标注、注释分组非模块、完全分离假象、零结果表述、血液数据边界等）**明确保留**，不得删除。

**字数**：3,503（B 方案加词后）→ **3,486**（−17）。B 方案加词完整保留，且终值回到上限内（余量 14）。

### 四、腾位候选更新

后续若需增引，腾位候选正式定为第 **6** 条 Jiang 2024 *Front Immunol*；Wei 2025 *Front Mol Biosci* 退回备选池（已写入 `59_` 第八节）。

### 五、审稿人建议

新文件 `61_Reviewer_Suggestions_2026-10-30.md`：主名单 8 名（Hotchkiss、Calfee、Bos、Davenport、van der Poll、Sinha、Shankar-Hari、Uhel）＋ 方法学备选 3 名（Scicluna、Burgess、Alipanah-Lechner），
给出三套组合（推荐"丙：Calfee、Bos、Davenport、Hotchkiss"）与回避名单。**邮箱一律不臆造**，仅 Calfee 来自实验室公开页，其余标 `[需核实]`。

### 校验（verify_manuscript.py，n_refs=50）

正文 **3,486**（引言 445 / Methods 353 / Results 1,529 / Discussion 1,095 / Conclusion 64）；摘要 **249**；
参考文献 **50 条**（全引用、无孤儿、首现严格升序 True）；主图 1–10 / 表 1–3 全引用；
em dash 0；AI 词 0；T5 唯一句 1；配对论断 **5/5 True**。

### 四联同步

02_ 英文稿 ✅｜03_ 中文稿 ✅（A1–A6 对应中文句同步改写，头部数据更新为 3,486 词 / 50 条）｜
`manuscript_revised.docx`（step95 重生成；Boomer 在、"spleen and lung" 在、"not a refutation"/"Two caveats" 已清、声明与基金在内、50 条）✅｜
投稿包 `~/Desktop/投稿材料`（step114 + safe-copy shim 重建：主图 10、补充图 21、补充表 27、内嵌图 10、正文表 3；Cover_Letter.docx 落款已实名化）✅｜
桌面建仓仓库 `~/Desktop/建仓材料/repo/manuscript/02_Manuscript_Revised.md` 已同步（**未 commit/push**）✅

完整依据见 `60_Defensive_Writing_Audit_2026-10-30.md` 与 `61_Reviewer_Suggestions_2026-10-30.md`。

---

## 第二十四批（2026-10-30）：CRediT 措辞裁决落地（P0 关闭）

### 触发

第二十三批标记的两处待裁决，作者已裁定。

### 处置

| # | 裁决 | 执行 |
|---|---|---|
| 1 | CRediT 改为 "contributed to data curation, formal analysis, and development of methodology" | 02_ `AUTHOR CONTRIBUTIONS` 与 51_ 标题页两处同步替换；51_ 第六节 P0 行标记为 ✅ 已解决；第一节省略提示同步改写 |
| 2 | 正文 **3,486** 词为准 | 51_ 字数口径维持 3,486（Introduction 445 / Methods 353 / Results 1,529 / Discussion 1,095 / Conclusion 64），摘要 249；作者提供的 3,747 旧口径作废 |

替换后与 `ETHICS APPROVAL AND CONSENT TO PARTICIPATE` 的 "no human or animal experiments were performed by the authors" 不再冲突，P0 关闭。

### 校验

正文 **3,486**（声明小节不计入，改动不影响字数）；摘要 249；参考文献 **50 条**、无孤儿、首现严格升序 True；
主图 1–10 / 表 1–3 全引用；em dash 0；AI 词 0；T5 唯一句 1；配对论断 **5/5 True**。

### 四联同步

02_ ✅｜03_ 中文稿（不含声明小节，无需改）✅｜`manuscript_revised.docx`（step95 重生成；"performing the experiments" 已清、"data curation, formal analysis" 已在）✅｜
投稿包 `~/Desktop/投稿材料`（step114 + safe-copy shim 重建：主图 10、补充图 21、补充表 27、内嵌图 10、正文表 3；包内 Manuscript_Submission.docx 同步）✅｜
51_ 已复制入 `05_内部备查_不投稿`；桌面建仓仓库 `02_`/`04_` 已同步（**未 commit/push**）✅

### 剩余待作者（非阻塞）

署名顺序确认、基金编号 `{2025}881` 花括号、通讯作者电话、Acknowledgments、预印本有无、审稿人最终 3–5 名与邮箱核实。

---

## 第二十四批（2026-10-30）：CRediT 措辞裁决落地（P0 关闭）

### 触发

第二十三批标记的两处待裁决，作者已裁定。

### 处置

| # | 裁决 | 执行 |
|---|---|---|
| 1 | CRediT 改为 "contributed to data curation, formal analysis, and development of methodology" | 02_ `AUTHOR CONTRIBUTIONS` 与 51_ 标题页两处同步替换；51_ 第六节 P0 行标记为 ✅ 已解决；第一节省略提示同步改写 |
| 2 | 正文 **3,486** 词为准 | 51_ 字数口径维持 3,486（Introduction 445 / Methods 353 / Results 1,529 / Discussion 1,095 / Conclusion 64），摘要 249；作者提供的 3,747 旧口径作废 |

替换后与 `ETHICS APPROVAL AND CONSENT TO PARTICIPATE` 的 "no human or animal experiments were performed by the authors" 不再冲突，P0 关闭。

### 校验

正文 **3,486**（声明小节不计入，改动不影响字数）；摘要 249；参考文献 **50 条**、无孤儿、首现严格升序 True；
主图 1–10 / 表 1–3 全引用；em dash 0；AI 词 0；T5 唯一句 1；配对论断 **5/5 True**。

### 四联同步

02_ ✅｜03_ 中文稿（不含声明小节，无需改）✅｜`manuscript_revised.docx`（step95 重生成；"performing the experiments" 已清、"data curation, formal analysis" 已在）✅｜
投稿包 `~/Desktop/投稿材料`（step114 + safe-copy shim 重建：主图 10、补充图 21、补充表 27、内嵌图 10、正文表 3；包内 Manuscript_Submission.docx 同步）✅｜
51_ 已复制入 `05_内部备查_不投稿`；桌面建仓仓库 `02_`/`04_` 已同步（**未 commit/push**）✅

### 剩余待作者（非阻塞）

署名顺序确认、基金编号 `{2025}881` 花括号、通讯作者电话、Acknowledgments、预印本有无、审稿人最终 3–5 名与邮箱核实。
