# HNSC Mitoxyperiosis — Reproducibility Repository

> **Manuscript**: Cross-Cohort Characterization of Mitoxyperiosis-Related Transcriptional Programs in Head and Neck Squamous Cell Carcinoma
> **Target journal**: Scientific Reports
> **Status**: DATA LOCKED / SCIENTIFIC CONTENT FROZEN（2026-09-27）
> **Modules**: `01_cohort_qc` … `09_figures` cover the bulk/single-cell analyses; `10_spatial_validation` covers the orthogonal GSE208253 spatial validation (frozen 2026-10-02).

**This repository contains the canonical analysis code underlying the locked manuscript results.** The repository was assembled without recomputing or modifying the frozen analyses; users may rerun the scripts after obtaining the required public datasets and configuring local data paths.

---

## 1. Canonical Scoring Definition（权威评分口径，v1.0）

全项目统一评分口径的唯一权威定义，见 `04_external_validation/Canonical_scoring_specification_v1.0.md`。四条锁定规则：

| # | 规则 | 说明 |
|---|---|---|
| 1 | **固定 67-gene common program** | 跨平台评分一律用 `locked_gene_sets/program_common_genes_3OS.csv` 的 67 基因；73-gene 仅用于生物学定义与 LASSO 候选池 |
| 2 | **TCGA 变换 = DESeq2 counts → log2(x+1)** | `TCGA_HNSC_tumor_full_expr.Rds` 是 DESeq2 归一化计数（值域 0~31,485），评分前必须 log2(x+1) 再基因级 z |
| 3 | **GEO 变换 = 不再 log** | GSE41613/GSE42743/GSE65858 已为 log 尺度，直接 z，防双重变换 |
| 4 | **两层 SD 口径** | 基因级 z 以**父队列**为参照；Cox 的 per-1-SD 在**分析子集内**标准化 |

**一句话 canonical 流程**：fixed 67 genes → TCGA log2(x+1) → parent-cohort gene-wise z → subgroup per-1-SD Cox。

## 2. 关键锁定结果（Sentinel Values）

| 指标 | 值 |
|---|---|
| TCGA overall | N=501/218 events, HR **1.1923** (1.0342–1.3747), P=0.0154 |
| TCGA oral | N=309/143, HR **1.3407** (1.1138–1.6139), P=0.0019 |
| GSE41613 | HR 1.2039 (0.9059–1.6001), P=0.201 |
| GSE42743 overall / oral | 0.9606 / 0.8893 |
| GSE65858 oral | HR **0.6800** (0.4518–1.0235), P=0.0645 |
| external pooled | ≈1.00 (0.85–1.18) |
| CD274 β | 0.58（模型 = CD274 ~ program + IFNG + purity + site）|
| single-cell | 19 evaluable patients |

## 3. 目录结构

```
reproducibility_repository/
├── README.md                     # 本文件
├── environment.md                # R 版本 + 包清单 + 数据源
├── locked_gene_sets/             # 锁定基因集（唯一权威）
│   ├── program_common_genes_3OS.csv   # 67-gene common program（锁定）
│   └── program_external_common73.csv  # 73-gene 生物学候选集
├── functions/                    # 共享函数
│   ├── 00_shared_functions.R
│   └── 01_gene_sets.R
├── 01_cohort_qc/                 # TCGA 501 队列构建
├── 02_gene_program/              # 67-gene program 定义 + 评分 + 证据表
├── 03_model_development/         # 11-gene LASSO + bootstrap 内部验证
├── 04_external_validation/       # canonical 外部验证（含 oral provenance）
├── 05_site_HPV/                  # 亚位点 / HPV 敏感性
├── 06_single_cell/               # GSE103322 单细胞定位 + null 模型
├── 07_immune_purity/             # CIBERSORT / ESTIMATE purity / 药物敏感性
├── 08_specificity_null/          # 特异性 benchmark + matched null
└── 09_figures/                   # 主图/补充图生成
```

## 4. Legacy vs Canonical Provenance（oral 分析）

`V5_P1_2_oral_external.R` 存在两处偏离（TCGA 漏 log2 + 子集 z），2026-09-27 已修正：

| 文件 | 状态 |
|---|---|
| `04_external_validation/V5_P1_2_oral_external_CANONICAL_20260927.R` | ✅ **canonical**（唯一可运行版本） |
| `04_external_validation/V5_P1_2_oral_external_LEGACY_20260927.R` | ⚠️ legacy 存档，**禁止运行**（头部有 supersedes 说明） |

同理 `05_site_HPV/site_subgroup_evaluation.R`（S12 旧，73 基因 + ghost-row 隐患）已被 canonical 口径重算；其脚本保留仅作 provenance。

## 5. 数据源（本 repository 不含数据）

脚本保留锁定分析时使用的原始本地数据路径 `D:/HNSC_mitoxyperiosis_positron/data`（以及输出路径 `.../V2_Reanalysis/...`）作为 provenance。**执行前，用户需将各脚本头部的 `DATA_DIR` / `OUT_DIR` 替换为对应的本地路径**（例如 `DATA_DIR <- file.path("data")`）。数据均为公开来源：

| 数据 | 来源 | 关键文件 |
|---|---|---|
| TCGA-HNSC | GDC / TCGAbiolinks（501 tumors） | `data/TCGA_full/TCGA_HNSC_tumor_full_expr.Rds`, `data/TCGA_HNSC_clinical.Rds` |
| GSE41613 | GEO（GPL570, n=97 oral） | `data/GSE41613_expr.Rds`, `data/GSE41613_surv.Rds` |
| GSE42743 | GEO（GPL570, n=74） | `data/GSE42743_processed/` |
| GSE65858 | GEO（Illumina HT-12, n=270） | `data/GSE65858_expr.Rds`, `data/GSE65858_pdata.Rds` |
| GSE27020 | GEO（GPL96, n=109 laryngeal, DFS） | `data/GSE27020_series_matrix.txt.gz` |
| GSE103322 | GEO（scRNA-seq HNSCC atlas） | `data/GSE103322_HNSCC_all_data.txt.gz` |
| DepMap / CCLE | DepMap portal / CCLE | `data/DepMap/` |

> ⚠️ **隐私声明**：本 repository 只含代码与基因集，**不含患者级数据**。TCGA/GEO 数据由第三方下载后置于本地 `data/`（未纳入本 repository）。基因集为公开基因符号，无隐私风险。

## 6. 复现步骤（概览）

1. 准备数据到本地 `data/`（按上表）。
2. 将各脚本头部 `DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"` 改为本地路径。
3. 按模块顺序运行：01_cohort_qc → 02_gene_program → 03_model_development → 04_external_validation → 05_site_HPV / 06_single_cell / 07_immune_purity / 08_specificity_null → 09_figures。
4. 核对输出与第 2 节 sentinel values 一致（参考 `Canonical_scoring_specification_v1.0.md` 的统一结果）。

## 7. 环境

R 4.3.3 + 22 个包（版本见 `environment.md`）。无 renv.lock（原始项目未启用 renv）；`environment.md` 已记录完整包版本清单供复现。

## 8. 禁词（除非否定语境）

validated prognostic signature / robust biomarker / HPV-independent / site-specific prognostic effect / mitoxyperiosis activity / activated mitoxyperiosis / clinical utility

---

_生成：2026-09-27 ｜ 由 Final Reviewer Simulation v1（通过）+ 四层 XLSX QC（2594 cells 零 mismatch）支撑_
