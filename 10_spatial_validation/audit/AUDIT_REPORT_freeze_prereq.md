# GSE208253 空间验证模块 — 结果冻结前数值一致性审计报告

**审计性质**：read-only / summary-only（8.docx 授权一次性执行）
**审计日期**：2026-10-02
**审计脚本**：`spatial_analysis/audit_numeric_consistency.py`（可复跑）
**审计结论**：**PASS** — 全部 30 项检查通过（修复后）

---

## 1. 审计范围与方法

| 项目 | 内容 |
|---|---|
| 真值来源 | step1–step7 原始 TSV（`score/` `moran/` `effects/` `robustness/` `quality/` `input/`） |
| 比对对象 | `report/GSE208253_spatial_validation_report.html`、`report/GSE208253_validation_master_tables.md`、图 caption |
| 核对字段 | 所有 meta p、depth-adjusted p、per-sample Moran's FDR、sample_3 spot counts、67-gene detection counts、sign test p、Cliff's δ、Spearman ρ |
| 禁止事项 | 未手工挑选数字；未删除基因；未重筛 gene set；未换 score；未改 region/comparator；未运行 DEG 或任何结果驱动探索 |
| 方法 | 脚本从 TSV 重新推导每个数字（含交叉样本汇总的原始统计定义），再与文档逐一比对 |

---

## 2. 发现的缺陷（3 项，全部为 reporting defect，不影响任何科学结论）

### 缺陷 1（严重）🔴 core vs nc 的 meta p 指数错读

| 项 | 内容 |
|---|---|
| 位置 | HTML §6.1 结论 callout；HTML §8 汇总表；MD §2.1 汇总行 |
| 文档原值 | `3.5 × 10⁻⁹`（HTML） / `3.5e-9`（MD） |
| TSV 真值 | **`3.467035e-18`** |
| 偏差 | **差 9 个数量级**（科学计数法指数 −18 被误写为 −9） |
| 性质 | 抄录/转换错误。数值方向（极度显著）未变，但 provenance 不可接受 |
| 处理 | 已统一更正为 **3.5 × 10⁻¹⁸** |

### 缺陷 2（严重）🔴 depth-adjusted core vs nc 的 p 用了错误的统计量

| 项 | 内容 |
|---|---|
| 位置 | MD §5（`0.351`）；WB5 转述（`0.351`） |
| 文档原值 | `0.351` |
| 正确值 | **`0.3709`** |
| 根本原因 | `0.351` 是把 t 统计量按**标准正态分布**算出的近似值（`math.erf`），而非该检验应有的 **Student-t 分布 df=11** |
| 检验定义（权威） | `step6_depth_check.R` 第 44 行 `t.test(sub$d_adj)` = 对 12 例 per-sample 深度校正后 Cliff's δ 的**单样本 t 检验**（H₀: mean δ = 0），df = 11 |
| 交叉验证 | R `t.test()` 独立复算 p = **0.370901**；Python `scipy.stats.t.sf(0.9329, 11)*2` = **0.370896** — 二者一致 |
| 两套定义明细 | core vs nc: t=0.9329 → t 分布 **0.3709** / 正态近似 0.3509（错）<br>edge vs nc: t=3.1196 → t 分布 **0.0098** / 正态近似 0.0018<br>transitory vs nc: t=3.5240 → t 分布 **0.0048** / 正态近似 0.0004 |
| 处理 | MD/HTML 统一为 **0.3709**；报告新增「两套定义并列」表格，明确标注检验定义，杜绝再次混淆。`0.351` 已废弃 |

### 缺陷 3（轻微）🟡 HTML 的 depth-adjusted p 数值正确但算法不透明

| 项 | 内容 |
|---|---|
| 位置 | HTML §10.2、§13.2 |
| 现象 | 原用 `math.erf` 正态近似渲染，恰好得出 0.371（四舍五入后与正确的 0.3709 接近），但算法错误 |
| 风险 | 数值对而算法错 → 一旦样本数或 δ 分布变化就会失配，属不可持续实现 |
| 处理 | 已改为精确 Student-t 尾概率（scipy，含 incomplete-beta 兜底），显式输出 t、df、p 三件套 |

### 附带修复
- 措辞按 8.docx 要求收紧：HTML §5「**真实的空间聚集**」→「**可重复的正空间自相关**」，并新增措辞边界说明段（空间自相关不能区分生物学组织与技术性空间结构）。
- §13.2 补入英文字样 *"The 67-gene score showed reproducible positive spatial autocorrelation across all 12 samples."*

---

## 3. 审计中发现的审计工具自身缺陷（已修）

| 缺陷 | 现象 | 修复 |
|---|---|---|
| markdown 裸 `>` 破坏 tag 剥离 | MD 表头 `\| δ>0 \| δ<0 \|` 被 `<[^>]+>` 正则跨单元格吞掉，导致其后 meta p 值"消失"，误报 MISS | ① 审计脚本改用**具名 HTML 标签白名单**正则；② MD 源已把表头 `>` 转义为 `&gt;` |
| 废弃值自引用假阳性 | 审计在"确认旧值已移除"时，匹配到**说明修复的注释本身**（注释合法引用了旧值） | 新增 `strip_deprecation_notes()`，先在剔除修复说明的前提下再检查旧值是否残留 |

> 这两条本身也是审计价值的一部分：一个会误报的审计脚本，与一个会漏报的审计脚本同样危险。

---

## 4. 逐项核对结果（修复后）

### 4.1 Moran's I
| 检查项 | TSV 真值 | 文档 | 状态 |
|---|---|---|---|
| 样本数 | 12 | 12/12 | ✅ |
| I > 0 | 12 | 12/12 | ✅ |
| FDR < 0.05 | 12 | 12/12 | ✅ |
| I 范围 | 0.0571 – 0.2304 | 0.057–0.230 | ✅ |
| I 中位 | 0.1260 | 0.126 | ✅ |
| s10 与 spdep 一致 | 0.128172079425（逐位相同） | 已声明 | ✅ |

### 4.2 region 效应量
| 检查项 | TSV 真值 | 文档 | 状态 |
|---|---|---|---|
| core vs nc δ<0 | 12/12 | 12/12 | ✅ |
| core meta δ | −0.2754 | −0.275 | ✅ |
| core meta p | **3.467e-18** | **3.5 × 10⁻¹⁸** | ✅（已修） |
| core meta I² | 77.0% | 77% | ✅ |
| core leave-s3-out meta p | 1.057e-16 | 1.1 × 10⁻¹⁶ | ✅ |
| edge vs nc δ<0 | 10/12 | 10/12 | ✅ |
| edge meta p | 1.121e-4 | 1.1 × 10⁻⁴ | ✅ |
| transitory vs nc δ<0 | 11/12 | 11/12 | ✅ |
| transitory meta p | 1.854e-6 | 1.9 × 10⁻⁶ | ✅ |

### 4.3 pathology
| 检查项 | TSV 真值 | 文档 | 状态 |
|---|---|---|---|
| SCC vs non-SCC δ<0 | 10/12 | 10/12 | ✅ |
| meta δ | −0.1327 | −0.133 | ✅ |
| meta p | 7.216e-6 | 7.2 × 10⁻⁶ | ✅ |
| leave-s3-out meta p | 3.702e-5 | 3.7 × 10⁻⁵ | ✅ |
| s9 SCC spots | 2494 | 2494 | ✅ |
| s10 SCC spots | 2132 | 2132 | ✅ |

### 4.4 robustness（meanZ）
| 检查项 | TSV 真值 | 文档 | 状态 |
|---|---|---|---|
| Spearman 负相关样本数 | 12/12 | 12/12 | ✅ |
| Spearman 均值 | −0.244 | −0.244 | ✅ |
| 方向一致性 | 5/48 = 10.4% | 5/48 (10.4%) | ✅ |
| Primary 未切换 | — | 声明 UCell | ✅ |

### 4.5 depth-adjusted（核心争议项）
| 检查项 | TSV 真值 | 文档 | 状态 |
|---|---|---|---|
| 原始 mean δ | −0.2869 | −0.287 | ✅ |
| 校正后 mean δ | +0.036504 | +0.037 | ✅ |
| 校正后 δ>0 例数 | 8/12 | 8/12 | ✅ |
| t 统计量 | 0.9329 | 0.9329 | ✅ |
| **p（t 分布 df=11，正确）** | **0.370901** | **0.3709** | ✅（已修） |
| p（正态近似，已废弃） | 0.350877 | 已标注为备查 | ✅ |

### 4.6 spot / detection 计数
| 检查项 | TSV 真值 | 文档 | 状态 |
|---|---|---|---|
| in-tissue spots | 26,371 | 26,371 | ✅ |
| annotated spots | 24,399 | 24,399 | ✅ |
| s3 in-tissue / annotated | 969 / 476 | 969 / 476 | ✅ |
| 67-gene 检出 | 67 | 67/67 | ✅ |
| 检出 <5% 基因数 | 8 | 8/67 | ✅ |
| 检出 <1% 基因数 | 2 | 2/67 | ✅ |
| CASP5 | 0.49% | 0.49% | ✅ |
| GLS2 | 0.60% | 0.60% | ✅ |

---

## 5. 三层科学结论（冻结版）

1. **空间结构性 — SUPPORTED**
   12/12 样本 Moran's I > 0 且 FDR < 0.05 → *reproducible positive spatial autocorrelation across all 12 samples*。
   ⚠️ 措辞边界：支持的是 non-random spatial organization，**不等于**已证明生物学空间聚集。

2. **肿瘤特异定位 — NOT CONFIRMED**
   预先规定 core vs nc（12/12 δ<0）与独立 SCC vs non-SCC（10/12 δ<0）均未验证预期的肿瘤高定位。两套独立区域定义指向相似现象，故不可简单称为噪声。

3. **原始 core 负向 — 不得解释为生物学耗竭**
   core/nc library depth 差 2.2–8.7 倍；深度校正后 core vs nc δ 由 −0.287 塌缩至 +0.037，单样本 t 检验 df=11，p = 0.3709（NS）。原始方向性在很大程度上是**测序深度混杂**。

**meanZ 身份**：仅预先规定 robustness，12/12 与 UCell 负相关、方向一致性 10.4% 如实保留；**不升级为救援证据**，**不切换 Primary**。

---

## 6. 未授权事项（冻结后禁止）

- ❌ 删除任何基因 / 重筛或"优化" gene set
- ❌ 更换 Primary score（UCell 恒为 Primary）
- ❌ 修改 region_4class 定义或 non-SCC comparator
- ❌ 运行 DEG / CellChat / trajectory / ML / drug prediction
- ❌ 任何以"修复"为名的结果驱动 GSE208253 探索（寻找最佳 score / 最佳区域阈值）

---

## 7. 交付与可复现

| 文件 | 说明 |
|---|---|
| `audit/numeric_consistency_audit.txt` | 人可读审计报告（本文件同源） |
| `audit/numeric_consistency_audit.json` | 机器可读审计结果 |
| `audit_numeric_consistency.py` | 可复跑审计脚本（read-only） |
| `report/GSE208253_spatial_validation_report.html` | 修正后主报告 |
| `report/GSE208253_validation_master_tables.md` | 修正后主表（含溯源头注） |

复跑方式：
```bash
cd D:/GSE208253_spatial/spatial_analysis
python audit_numeric_consistency.py
```
预期输出：`VERDICT: PASS`
