## B1: TCGA-HNSC 全队列重建（502 tumor + 44 normal）
## 数据来源: 旧电脑迁移的 TCGA-HNSC.txt（548 样本 × 60,660 基因）
## 目标: 1) 解析完整表达矩阵 2) 与现有 226 例对比 3) 重建 tumor-only 队列 4) attrition flow

suppressMessages({
  library(data.table)
  library(survival)
})

DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/01_cohort_qc"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. 读取完整表达矩阵（60,660 基因 × 548 样本）
# ============================================================
cat(">>> 读取完整 TCGA-HNSC 表达矩阵...\n")
expr_full <- fread(file.path(DATA_DIR, "TCGA_full", "TCGA-HNSC_full.txt"),
                   header = TRUE, sep = "\t", check.names = FALSE)
gene_ids <- expr_full[[1]]
expr_full <- as.matrix(expr_full[, -1])
rownames(expr_full) <- gene_ids
cat("完整矩阵:", nrow(expr_full), "genes x", ncol(expr_full), "samples\n")

# 样本类型 (TCGA barcode: TCGA-XX-XXXX-01A-... 第4段前2位 = 01 tumor / 11 normal / 06 metastatic)
sample_ids <- colnames(expr_full)
parse_type <- function(s) {
  parts <- strsplit(s, "-")[[1]]
  if (length(parts) >= 4) substr(parts[4], 1, 2) else "NA"
}
sample_types <- sapply(sample_ids, parse_type)
cat("样本类型分布:\n")
print(table(sample_types))

tumor_idx <- which(sample_types == "01")
normal_idx <- which(sample_types == "11")
cat("\nTumor(01):", length(tumor_idx), " Normal(11):", length(normal_idx), "\n")

# ============================================================
# 2. 与现有 226 例对比
# ============================================================
cat("\n>>> 与现有 226 例对比...\n")
tpm_old <- readRDS(file.path(DATA_DIR, "TCGA_HNSC_TPM_tumor.Rds"))
old_samples <- colnames(tpm_old)
cat("现有 TPM 样本数:", length(old_samples), "\n")

# 现有样本是短格式 TCGA-XX-XXXX 还是全格式?
cat("现有样本格式示例:", head(old_samples, 3), "\n")
new_tumor_ids <- sample_ids[tumor_idx]
# 标准化比较: 取前 12 位 patient barcode
old_pat <- substr(old_samples, 1, 12)
new_pat <- substr(new_tumor_ids, 1, 12)
overlap <- intersect(old_pat, new_pat)
cat("现有 226 例中与完整数据匹配:", length(overlap), " 不匹配:", length(old_pat) - length(overlap), "\n")
new_only <- setdiff(new_pat, old_pat)
cat("完整数据中新增患者:", length(new_only), "\n")

# ============================================================
# 3. 保存 tumor-only 全队列 (Rds, 供后续分析)
# ============================================================
cat("\n>>> 保存 tumor-only 全队列...\n")
tumor_mat <- expr_full[, tumor_idx]
cat("Tumor 矩阵:", nrow(tumor_mat), "x", ncol(tumor_mat), "\n")

# 去重复样本 (同一 patient 多个 tumor sample 时保留第一个)
tumor_pat <- substr(colnames(tumor_mat), 1, 12)
dup_pat <- tumor_pat[duplicated(tumor_pat)]
if (length(dup_pat) > 0) {
  cat("重复患者样本数:", length(dup_pat), "个 (每患者保留第一个)\n")
  keep <- !duplicated(tumor_pat)
  tumor_mat <- tumor_mat[, keep]
  cat("去重后 Tumor:", ncol(tumor_mat), "\n")
}

saveRDS(tumor_mat, file.path(DATA_DIR, "TCGA_full", "TCGA_HNSC_tumor_full_expr.Rds"))
write.csv(data.frame(sample = colnames(tumor_mat), patient = substr(colnames(tumor_mat), 1, 12)),
          file.path(OUT_DIR, "TCGA_full_tumor_manifest.csv"), row.names = FALSE)

# 保存 normal 矩阵
if (length(normal_idx) > 0) {
  normal_mat <- expr_full[, normal_idx]
  saveRDS(normal_mat, file.path(DATA_DIR, "TCGA_full", "TCGA_HNSC_normal_expr.Rds"))
  cat("Normal 矩阵:", nrow(normal_mat), "x", ncol(normal_mat), "\n")
}

# ============================================================
# 4. PRKN 表达验证（完整数据）
# ============================================================
cat("\n>>> PRKN 表达验证...\n")
prkn_hit <- which(rownames(expr_full) == "PRKN")
if (length(prkn_hit) > 0) {
  prkn_tumor <- expr_full[prkn_hit, tumor_idx]
  prkn_normal <- expr_full[prkn_hit, normal_idx]
  cat(sprintf("PRKN Tumor: n=%d, median=%.3f, mean=%.3f\n", length(prkn_tumor), median(prkn_tumor), mean(prkn_tumor)))
  cat(sprintf("PRKN Normal: n=%d, median=%.3f, mean=%.3f\n", length(prkn_normal), median(prkn_normal), mean(prkn_normal)))
  wt <- tryCatch(wilcox.test(prkn_tumor, prkn_normal), error = function(e) NULL)
  if (!is.null(wt)) cat(sprintf("Mann-Whitney P = %.3g\n", wt$p.value))
} else {
  cat("PRKN 未在 rownames 中找到! 检查:", grep("PRKN|PARK2", rownames(expr_full), value = TRUE)[1:5], "\n")
}

# ============================================================
# 5. 5-gene 覆盖检查
# ============================================================
cat("\n>>> 5-gene 覆盖检查...\n")
coef5 <- c("PRKN", "VDAC1", "SLC25A5", "HK1", "MAOB")
for (g in coef5) {
  hit <- which(rownames(expr_full) == g)
  if (length(hit) > 0) {
    vals <- expr_full[hit, tumor_idx]
    cat(sprintf("  %s: 存在, tumor 非零率 %.1f%%, 范围 %.3f-%.3f\n",
                g, mean(vals > 0) * 100, min(vals), max(vals)))
  } else {
    cat(sprintf("  %s: 缺失\n", g))
  }
}

cat("\n=== B1 数据准备 DONE ===")
