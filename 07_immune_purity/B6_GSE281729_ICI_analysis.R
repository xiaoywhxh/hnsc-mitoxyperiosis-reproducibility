## B6: GSE281729 — 真实 ICI 治疗队列分析
## 42 例可切除 HNSCC，nivolumab ± IDO 抑制剂，含 baseline RNA-seq / HPV / 病理反应
## 分析: baseline samples 中 mitoxyperiosis-related score vs responder/non-responder

suppressMessages({
  library(ggplot2)
})

DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data/GSE281729"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/12_ici"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. 解析 processed 文件
# ============================================================
lines <- readLines(gzfile(file.path(DATA_DIR, "GSE281729_processed.txt.gz"), open = "rt"), warn = FALSE)
hdr <- strsplit(lines[1], "\t")[[1]]
meta_rows <- 14  # 行 0-13 是元数据（含 ENSG VLOOKUP 标题行）

# 元数据
meta <- lapply(1:8, function(i) strsplit(lines[i], "\t")[[1]])
names(meta) <- c("Subject", "FirstDrug", "SecondDrug", "Doses", "TimePoint", "HPV", "Smoking", "Response")

n_samples <- length(meta[[1]]) - 1
cat("样本数:", n_samples, "\n")

# 基因表达矩阵（从第 meta_rows+1 行开始；第 1 列=ENSG, 第 2 列=symbol, 第 3+ 列=表达）
gene_rows <- (meta_rows + 1):length(lines)
expr_mat <- matrix(NA, nrow = length(gene_rows), ncol = n_samples)
gene_names <- character(length(gene_rows))
for (j in seq_along(gene_rows)) {
  parts <- strsplit(lines[gene_rows[j]], "\t")[[1]]
  gene_names[j] <- parts[2]  # symbol
  expr_mat[j, ] <- as.numeric(parts[3:(n_samples + 2)])
}
rownames(expr_mat) <- gene_names
colnames(expr_mat) <- meta$Subject[-1]
cat("表达矩阵:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")

# 样本标签
sample_df <- data.frame(
  Sample = meta$Subject[-1],
  Drug = meta$FirstDrug[-1],
  SecondDrug = meta$SecondDrug[-1],
  TimePoint = trimws(meta$TimePoint[-1]),
  HPV = trimws(meta$HPV[-1]),
  Response = trimws(meta$Response[-1]),
  stringsAsFactors = FALSE
)
cat("\nTimePoint 分布:", paste(names(table(sample_df$TimePoint)), table(sample_df$TimePoint), sep = "=", collapse = ", "), "\n")
cat("Response 分布:", paste(names(table(sample_df$Response)), table(sample_df$Response), sep = "=", collapse = ", "), "\n")

# ============================================================
# 2. Risk score (PRKN 缺失 → 用重训练的 4-gene reduced model 系数)
# ============================================================
# 5-gene locked 系数
coef5 <- c(PRKN = 0.128131184842568, VDAC1 = 0.00177343236540608,
           SLC25A5 = 0.000299123798540387, HK1 = 0.00316219765097715,
           MAOB = 0.00354527605147633)
# 4-gene retrained reduced model 系数（TCGA 重训练，05_4gene_retrain.R）
coef4 <- c(VDAC1 = 0.00182238735085086, SLC25A5 = 0.000350720872393251,
           HK1 = 0.00351984025324572, MAOB = 0.00549044277512729)

gene5 <- names(coef5)
gene4 <- names(coef4)
found <- gene5[gene5 %in% rownames(expr_mat)]
missing <- setdiff(gene5, rownames(expr_mat))
cat("\n5-gene 在 GSE281729 覆盖:", paste(found, collapse = ", "), " 缺失:", ifelse(length(missing) == 0, "无", paste(missing, collapse = ", ")), "\n")

# PRKN 缺失时使用 4-gene reduced model
use_coef <- if ("PRKN" %in% missing) coef4 else coef5[found]
use_genes <- names(use_coef)
cat("使用模型:", if ("PRKN" %in% missing) "4-gene reduced (PRKN absent)" else "5-gene locked", "\n")

# 仅 baseline 样本
base_idx <- which(trimws(sample_df$TimePoint) %in% c("Pre", "Pre "))
cat("Baseline 样本数:", length(base_idx), "\n")

# 计算 risk score（baseline samples）
sub <- expr_mat[use_genes, base_idx, drop = FALSE]
# 剔除表达全 NA 的样本（如 51 H-I）
keep <- which(apply(sub, 2, function(x) !any(is.na(x))))
cat("剔除全 NA 样本:", sum(apply(sub, 2, function(x) any(is.na(x)))), "个, 保留:", length(keep), "\n")
base_idx <- base_idx[keep]
sub <- sub[, keep, drop = FALSE]

z <- t(scale(t(sub)))  # gene-level Z-score
risk <- colSums(use_coef * z)
names(risk) <- colnames(sub)

base_df <- sample_df[base_idx, ]
base_df$RiskScore <- risk[as.character(base_df$Sample)]

# ============================================================
# 3. Responder vs Non-responder 分析
# ============================================================
# Response 分级: CR/PR/R = responder; NR/min Resp = non-responder (病理反应)
base_df$RespBin <- ifelse(base_df$Response %in% c("CR", "PR", "R"), "Responder", "Non-responder")
cat("\nResponder 分布:", paste(names(table(base_df$RespBin)), table(base_df$RespBin), sep = "=", collapse = ", "), "\n")

# Wilcoxon rank-sum
res_comp <- tryCatch(wilcox.test(RiskScore ~ RespBin, data = base_df), error = function(e) NULL)
if (!is.null(res_comp)) {
  cat(sprintf("Responder vs Non-responder: P=%.4f\n", res_comp$p.value))
  meds <- tapply(base_df$RiskScore, base_df$RespBin, median)
  cat("Median risk: Responder =", round(meds["Responder"], 3),
      " Non-responder =", round(meds["Non-responder"], 3), "\n")
}

# logistic regression (adjust HPV) — 仅使用完整数据行
base_df$HPV_bin <- ifelse(base_df$HPV %in% c("Positive", "Positive "), 1, 0)
base_df$RespBin_num <- ifelse(base_df$RespBin == "Responder", 1, 0)
log_df <- base_df[complete.cases(base_df[, c("RiskScore", "HPV_bin", "RespBin_num")]), ]
fit_log <- tryCatch(glm(RespBin_num ~ RiskScore + HPV_bin, data = log_df, family = binomial),
                    error = function(e) NULL)
if (!is.null(fit_log)) {
  s <- summary(fit_log)
  cat("Logistic (risk + HPV): risk coef =", round(coef(fit_log)["RiskScore"], 3),
      " P =", round(s$coefficients["RiskScore", 4], 4), "\n")
} else {
  cat("Logistic: 拟合失败（样本量小或分离）\n")
}

# 图
p1 <- ggplot(base_df, aes(x = RespBin, y = RiskScore, fill = RespBin)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  scale_fill_manual(values = c("Responder" = "#E64B35", "Non-responder" = "#4DBBD5")) +
  labs(title = "GSE281729: mitoxyperiosis-related risk score by pathologic response",
       subtitle = "Neoadjuvant nivolumab ± IDO inhibitor, baseline samples",
       x = "", y = "5-gene risk score (Z)") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(OUT_DIR, "Fig_ICI_responder_boxplot.pdf"), p1, width = 5.5, height = 5)
ggsave(file.path(OUT_DIR, "Fig_ICI_responder_boxplot.png"), p1, width = 5.5, height = 5, dpi = 300)

# 保存
write.csv(base_df, file.path(OUT_DIR, "GSE281729_baseline_risk_scores.csv"), row.names = FALSE)

# ============================================================
# 4. 摘要输出
# ============================================================
summary_txt <- sprintf(
"GSE281729 ICI 队列分析摘要 (B6)
===============================
样本: %d (baseline %d)
模型: %s (PRKN %s)
Responder vs Non-responder:
  Wilcoxon P = %.4f
  Median risk: Responder %.3f vs Non-responder %.3f
Logistic (risk + HPV): risk coef %.3f, P = %.4f
",
n_samples, nrow(base_df),
ifelse("PRKN" %in% missing, "4-gene reduced", "5-gene locked"),
ifelse("PRKN" %in% missing, "缺失→4-gene", "覆盖"),
ifelse(is.null(res_comp), NA, res_comp$p.value),
ifelse(is.null(res_comp), NA, meds["Responder"]),
ifelse(is.null(res_comp), NA, meds["Non-responder"]),
ifelse(is.null(fit_log), NA, coef(fit_log)["RiskScore"]),
ifelse(is.null(fit_log), NA, summary(fit_log)$coefficients["RiskScore", 4]))
writeLines(summary_txt, file.path(OUT_DIR, "GSE281729_summary.txt"))
cat(summary_txt)

cat("\n=== B6 DONE. Outputs in:", OUT_DIR, "===\n")
