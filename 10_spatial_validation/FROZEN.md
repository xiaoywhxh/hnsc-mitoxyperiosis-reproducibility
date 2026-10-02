# GSE208253 空间验证模块 — FROZEN / CLOSED

**冻结日期**：2026-10-02
**正式结案**：2026-10-02（9.docx 裁决：接受审计 PASS + 接受 FROZEN.md，模块正式 CLOSED）
**冻结依据**：预先规定空间验证验收通过（8.docx）+ 结果冻结前数值一致性审计 PASS（`VERDICT: PASS`）
**审计报告**：`audit/AUDIT_REPORT_freeze_prereq.md`（`audit_numeric_consistency_audit.txt` 同源）
**下一阶段**：仅授权 manuscript / figure / supplement integration 与引用一致性检查；**本模块分析定义与结果自此不可变更**。

---

## 权威数值锁定（不可再变）

| 指标 | 权威值 |
|---|---|
| core-vs-nc REML meta δ | **−0.2754** |
| 95% CI | **−0.3375 … −0.2134** |
| meta p | **3.467 × 10⁻¹⁸** |
| depth-adjusted core-vs-nc mean δ | **+0.0365**（8/12 为正） |
| depth-adjusted one-sample t | **0.9329**，df = **11** |
| depth-adjusted p | **0.3709** |
| `0.351` | ⚠️ **已废弃** normal-approximation 历史值，仅存于 audit provenance，**不得再作为科学结果引用** |

---

## 冻结的三层科学结论

### ① 空间结构性 — SUPPORTED
67-gene UCell score 在 12/12 样本中呈正 Moran's I 且 FDR < 0.05（I 0.057–0.230，中位 0.126）。
表述：*The 67-gene score showed reproducible positive spatial autocorrelation across all 12 samples.*
**措辞边界**：支持 non-random spatial organization；不等同于已证明生物学空间聚集（空间自相关不区分生物学组织与技术性空间结构）。

### ② 肿瘤特异定位 — NOT CONFIRMED
- 预先规定 Primary：core vs nc，δ < 0 在 12/12（meta δ = **−0.2754**，95% CI **−0.3375 … −0.2134**，p = **3.467 × 10⁻¹⁸**，I² = 77%；sign test p = 4.88×10⁻⁴）
- 独立验证：SCC vs 预先定义 non-SCC，δ < 0 在 10/12（meta δ = −0.133，p = 7.2×10⁻⁶，I² = 91%）
- 两套独立区域定义指向相似现象；但方向与"肿瘤高表达"预期相反。

### ③ 原始 core 负向 — 不得解释为生物学耗竭
core spot 中位 library size 为 nc 的 2.2–8.7 倍；深度校正后 core vs nc mean δ 由 −0.2869 → **+0.0365**（8/12 为正），
单样本 t 检验 df = 11，t = **0.9329**，**p = 0.3709（NS）**。
→ 原始方向性在很大程度上是**测序深度混杂**，不得解释为"program 在肿瘤核心生物学耗竭"。

---

## 冻结的方法学身份

| 项 | 状态 |
|---|---|
| Primary score | **UCell**（rank-based，per-sample 独立计算）— 恒为唯一 Primary |
| 67 genes | **全部保留**，未删、未重筛、未优化 |
| meanZ | **仅预先规定 robustness**；12/12 与 UCell 负相关（均值 −0.244），方向一致性 5/48 (10.4%) — 如实保留，**不升级为救援证据** |
| 低检出基因 | CASP5 (0.49%)、GLS2 (0.60%) 等 8 个 < 5% 基因**保留**，作技术限制报告 |
| sample_3 | Moran's I 用全 969 spots；annotation 分析用 476 可靠 spots；未插补；leave-s3-out 已做 |
| 跨样本推断 | 每例独立效应量 + sign test + random-effects meta；**未** pooled spot-level 推断 |
| pathology comparator | 预先定义（打分前固定 include/exclude 清单） |

---

## 冻结后禁止事项（未授权）

- ❌ 删除任何基因 / 重筛或"优化" gene set
- ❌ 更换 Primary score（UCell 恒为唯一 Primary；meanZ 仅 robustness）
- ❌ 修改 region_4class 定义或 non-SCC comparator
- ❌ 运行 DEG / CellChat / trajectory / ML / drug prediction
- ❌ 任何以"修复"为名的结果驱动 GSE208253 探索
- ❌ 为获得阳性定位而尝试新的 score 或阈值
- ⚠️ 不再返回"冻结前复审"；下一道门 = 论文整合后的 **manuscript-level consistency review**

---

## 本模块产物清单（只读冻结）

```
spatial_analysis/
├── input/spot_table_all.tsv              26,371 × (sample,barcode,coords,region,pathology)
├── score/spot_scores.tsv                 26,371 × (ucell, meanz, library, nGene)
├── moran/morans_I_by_sample.tsv          12 × (I, E[I], z, p, FDR, neighbours)
├── effects/
│   ├── region_effect_sizes.tsv           36 rows (3 comparisons × 12 samples)
│   ├── region_cross_sample_summary.tsv
│   ├── region_leave_s3_out_summary.tsv
│   ├── pathology_effect_sizes.tsv        12 rows
│   ├── pathology_cross_sample_summary.tsv
│   └── pathology_leave_s3_out_summary.tsv
├── robustness/
│   ├── robustness_spearman.tsv
│   ├── robustness_direction_agreement.tsv
│   └── depth_adjusted_region_direction.tsv
├── quality/gene_detection_by_sample.tsv  67 × 12 detection matrix
├── plots/                                6 figures × (PNG + PDF)
├── report/                               HTML + master tables MD
└── audit/                                audit report + json + script
```

---

## 变更历史

| 日期 | 事件 |
|---|---|
| 2026-10-02 | 预先规定空间验证执行完成（step1–step7） |
| 2026-10-02 | 8.docx 验收；**数值一致性审计 PASS**；修正 3 项 reporting defect（meta p 指数错读、depth-adj p 统计量、措辞收紧）；**模块冻结** |
| 2026-10-02 | 9.docx 裁决：**接受审计 PASS 与 FROZEN.md，模块正式 CLOSED**；权威数值锁定；FROZEN.md 内 meta p 由 `3.5×10⁻¹⁸` 统一为 **`3.467×10⁻¹⁸`**；转段 manuscript integration |
