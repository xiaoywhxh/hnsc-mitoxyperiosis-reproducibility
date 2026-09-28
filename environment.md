# Environment — HNSC Mitoxyperiosis Reproducibility

> 生成：2026-09-27（由 Rscript 实查版本，非手填）

## R 版本

- **R 4.3.3 (2024-02-29 ucrt)** — Windows x64
- 无 renv.lock（原始项目未启用 renv）；下表为实查的完整包版本清单

## R 包清单（实查版本）

| 包 | 版本 | 用途 |
|---|---|---|
| survival | 3.5.8 | Cox 比例风险（全项目核心） |
| ggplot2 | 4.0.3 | 绘图 |
| timeROC | 0.4 | 时间依赖 ROC |
| glmnet | 5.0 | LASSO 惩罚 Cox（11-gene 探索模型） |
| GSVA | 1.50.5 | 基因集变异分析 |
| singscore | 1.22.0 | 单样本评分（program score 备选） |
| data.table | 1.18.4 | 数据处理 |
| estimate | 1.0.13 | ESTIMATE 纯度（CD274 β 回归协变量） |
| pheatmap | 1.0.12 | 热图 |
| patchwork | 1.3.2 | 拼图（Figure S21） |
| survminer | 0.5.0 | 生存曲线 |
| org.Hs.eg.db | 3.18.0 | 基因 ID 注释 |
| hgu133plus2.db | 3.13.0 | GPL570 探针注释（GSE41613/GSE42743） |
| hgu133a.db | 3.13.0 | GPL96 探针注释（GSE27020） |
| reshape2 | 1.4.4 | 数据整形 |
| dplyr | 1.2.1 | 数据操作 |
| R.utils | 2.13.0 | 工具 |
| png | 0.1.9 | PNG 读写（Figure 组装） |
| grid | 4.3.3 | 图形底层 |
| ggpubr | 1.0.0 | 出版级图 |

> CIBERSORT 使用本地 R 端口（LM22 reference），非 CRAN 包；具体实现见 `07_immune_purity/` 相关脚本。

## 数据源（本 repository 不含数据，需第三方自行获取）

| 数据 | 来源 | 关键文件 |
|---|---|---|
| TCGA-HNSC | GDC / TCGAbiolinks（501 tumors） | `data/TCGA_full/TCGA_HNSC_tumor_full_expr.Rds`, `data/TCGA_HNSC_clinical.Rds` |
| GSE41613 | GEO（GPL570, n=97 oral cavity） | `data/GSE41613_expr.Rds`, `data/GSE41613_surv.Rds` |
| GSE42743 | GEO（GPL570, n=74） | `data/GSE42743_processed/` |
| GSE65858 | GEO（Illumina HT-12, n=270） | `data/GSE65858_expr.Rds`, `data/GSE65858_pdata.Rds` |
| GSE27020 | GEO（GPL96, n=109 laryngeal, DFS） | `data/GSE27020_series_matrix.txt.gz` |
| GSE103322 | GEO（scRNA-seq HNSCC atlas） | `data/GSE103322_HNSCC_all_data.txt.gz` |
| DepMap / CCLE | DepMap portal / CCLE | `data/DepMap/` |

## 路径说明

脚本内硬编码绝对路径 `D:/HNSC_mitoxyperiosis_positron/data`（数据）与 `D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/...`（输出）。**复现时需改为本地路径**。这是原始项目的既有约定，本 repository 保留脚本原样以忠实 provenance。
