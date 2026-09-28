# Repository QC Report

> 生成：2026-09-27 ｜ 对象：`reproducibility_repository/`（本地 submission-ready snapshot）

## QC 结论：PASS（无 blocker）

| QC 层 | 检查项 | 结果 |
|---|---|---|
| QC1 绝对路径 | 34 个脚本硬编码 `D:/HNSC_mitoxyperiosis_positron` | ⚠️ 已知，README 第 5 节已说明「复现时改本地路径」；不阻断 provenance |
| QC2 67-gene 唯一 | locked CSV = 67 基因，**无重复**；唯一来源 = `program_common_genes_3OS.csv`（15 脚本引用） | ✅ PASS |
| QC2b 硬编码检查 | 脚本硬编码的是 **73-gene 候选集**（core/mtor/mito/metab 四模块）与辅助基因集（stress proxy NRF2/UPRmt、specificity 对照程序），非 67-gene | ✅ 不构成 67-gene 不一致 |
| QC3 legacy 隔离 | `V5_P1_2_oral_external_LEGACY_20260927.R` 头部有 supersedes 说明，标注「禁止运行」 | ✅ PASS |
| QC4 患者隐私 | repository 含 **0 个 .Rds 患者数据**；仅 2 个基因集 CSV（公开基因符号） | ✅ 无隐私/再分发风险 |
| QC5 包版本 | `environment.md` 记录 R 4.3.3 + 16 个核心包实查版本 | ✅ PASS |

## 说明

1. **无 renv.lock**：原始项目未启用 renv；`environment.md` 以 Rscript 实查版本补齐，供第三方复现。
2. **绝对路径为既有约定**：脚本以 `D:/HNSC_mitoxyperiosis_positron/data` 读数据。repository 忠实保留原样（provenance 优先），README 已说明复现时改路径。
3. **73-gene 候选集硬编码**：canonical spec 明确「73 仅用于生物学定义与 LASSO 候选池，评分口径为 67-gene CSV」，故硬编码 73 不影响 67-gene 唯一性。

## 后续动作（需老师决策）

- 本地 repository 已 submission-ready。是否公开 GitHub + Zenodo DOI，由老师决定（属外部公开动作）。公开后需将 manuscript 的 "will be deposited" 改为真实 URL/DOI。
