# Canonical Scoring Specification — v1.0

> **锁定日期**: 2026-09-27 ｜ **适用范围**: HNSC mitoxyperiosis 项目全部队列/site/subgroup/HPV/purity/stress/specificity/oral 敏感性脚本
> **性质**: 全项目统一评分口径的唯一权威定义；此后所有脚本必须引用本 specification，不得各自隐式实现。
> **背景**: 2026-09-27 审计发现旧 S12 链（73 基因）与旧 S20 链（TCGA 漏 log2 + 子集 z）分别偏离预规定义，已统一并锁定为以下四条。

## 四条锁定规则

1. **Gene set — 固定 67-gene common program**
   跨平台评分一律使用 `program_common_genes_3OS.csv` 中的 67 个共同基因；不使用 73-gene 候选集（73 仅用于生物学定义与 LASSO 候选池）。被排除的 6 基因 = PRKN, AKT2, MIEF2, MARCHF5, STING1, CGAS（GSE65858 平台不可测）。

2. **TCGA 表达变换 — DESeq2-normalized counts → log2(x+1)**
   `TCGA_HNSC_tumor_full_expr.Rds` 为 DESeq2 归一化原始计数（实测值域 0 ~ 31,485），评分前**必须** `log2(x + 1)`，随后基因级 z。任何 TCGA 脚本若直接对 raw matrix 做 `t(scale(t(expr)))` 均属违规。

3. **GEO 表达变换 — 不再 log**
   GSE41613 / GSE42743 / GSE65858 表达矩阵已为 log 尺度（实测 1.7~15.6），直接 `t(scale(t(expr)))`，不得再 log（防双重变换）。

4. **Subgroup Cox 的两层 SD 口径（关键）**
   - **gene-level z 参考 = 父队列**：对 anatomical-subsite 等 subgroup 分析，基因级 Z-score 标准化以**完整父队列**为参照总体（TCGA 501 / GSE65858 270 / GSE42743 74），**先**算 z **再**按 site 切子集。
   - **program-score per-1-SD = 分析子集内**：Cox 中的 `scale(score)` 在**所分析子集内**标准化到 SD=1，故 HR 表示「该子集内 score 每增加 1 SD」的风险比。
   - 这两层 SD 不可混淆（正是旧 S20 的 interpretation ambiguity 来源）。

## 亚位/子集 subsetting 强制写法

逻辑索引含 NA 会产生「幽灵行」（见 GSE65858 4 个 site 未映射样本：`pd[pd$site=="Oral",]` 得 87 行而 Cox complete-case 仅 83）。所有 subgroup 子集化一律用 `which()` 或显式 `!is.na(...)`：
```r
idx <- which(df$site_group == st & !is.na(df$os_event))
sub <- df[idx, ]
```

## 与旧实现的对应关系

| 旧脚本 | 偏差 | 处置 |
|---|---|---|
| `V5_P1_2_oral_external.R`（S20/Fig S13） | TCGA 漏 log2 + 子集 z | → `..._LEGACY_20260927.R`（存档禁运行），由 `..._CANONICAL_20260927.R` 取代 |
| `site_subgroup_evaluation.R`（S12） | 73 基因 + 潜在 ghost-row | 由 `regenerate_canonical_tables.R` 的 67g + which() 口径重算 |

## 锁定后的统一结果（父队列 z 口径）

- TCGA overall 67g：HR 1.1923（1.0342–1.3747），P=0.0154
- TCGA oral：HR 1.3407（1.1138–1.6139），P=0.0019
- GSE41613（100% oral）：HR 1.2039（0.9059–1.6001），P=0.201
- GSE42743 overall：HR 0.9606（0.6895–1.3382），P=0.812
- GSE42743 oral：HR 0.8893（0.6258–1.2637），P=0.513
- GSE65858 oral：HR 0.6800（0.4518–1.0235），P=0.0645
