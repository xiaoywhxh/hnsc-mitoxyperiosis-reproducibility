## B1: 全队列 (502 tumor) 关键分析重跑 — 与 226 例对比
## 目标: 1) attrition flow 2) program Cox 3) 5-gene LASSO/性能 4) PRKN 相关

suppressMessages({
  library(data.table)
  library(survival)
  library(glmnet)
  library(timeROC)
})

DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/01_cohort_qc"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

set.seed(20260720)

# ============================================================
# 1. 加载全队列表达 + 临床
# ============================================================
tumor_mat <- readRDS(file.path(DATA_DIR, "TCGA_full", "TCGA_HNSC_tumor_full_expr.Rds"))  # 60660 x 502
clin <- readRDS(file.path(DATA_DIR, "TCGA_HNSC_clinical.Rds"))  # 528 patients
surv_old <- readRDS(file.path(DATA_DIR, "TCGA_HNSC_survival.Rds"))  # 225 (226 例队列)

# 用 GDC clinical 的生存数据 — 先检查 surv_old 能否扩展到全队列
# clinical 有 vital_status, days_to_death, days_to_last_followup
clin$OS_time <- ifelse(!is.na(clin$days_to_death), clin$days_to_death,
                       ifelse(!is.na(clin$days_to_last_followup), clin$days_to_last_followup, NA))
clin$OS_status <- ifelse(clin$vital_status == "Dead", 1, 0)

patient_ids <- substr(colnames(tumor_mat), 1, 12)
clin_match <- match(patient_ids, clin$bcr_patient_barcode)

# 完整队列生存数据
surv_full <- data.frame(
  sample = colnames(tumor_mat),
  patient = patient_ids,
  OS_time = clin$OS_time[clin_match],
  OS_status = clin$OS_status[clin_match],
  stringsAsFactors = FALSE
)
cat("502 tumor 中有生存数据:", sum(!is.na(surv_full$OS_time)), "\n")

# ============================================================
# 2. Attrition flow
# ============================================================
cat("\n>>> TCGA-HNSC Attrition Flow:\n")
n_total <- 528
n_tumor_expr <- 502
n_with_os <- sum(!is.na(surv_full$OS_time))
cat(sprintf("TCGA-HNSC 总病例: %d\n", n_total))
cat(sprintf("→ 有表达数据 (01 tumor): %d\n", n_tumor_expr))
cat(sprintf("→ 有完整 OS: %d\n", n_with_os))

flow_df <- data.frame(
  Step = c("All TCGA-HNSC cases", "Tumor samples with RNA-seq", "With complete OS data"),
  N = c(n_total, n_tumor_expr, n_with_os)
)
write.csv(flow_df, file.path(OUT_DIR, "TCGA_attrition_flow.csv"), row.names = FALSE)
print(flow_df)

# ============================================================
# 3. Program-level Cox（502 vs 226 对比）
# ============================================================
cat("\n>>> Program-level Cox (502 队列)...\n")
# 需要 ssGSEA 评分 — 用完整数据重算 (73 基因程序)
core_genes <- c("PRKN", "VDAC1", "VDAC2", "VDAC3", "BAX", "BAK1", "BID",
                "BBC3", "PMAIP1", "BCL2", "BCL2L1", "MCL1", "TSPO", "AIFM1")
mtor_genes <- c("MTOR", "RICTOR", "RPTOR", "MLST8", "MAPKAP1", "PRR5", "PRR5L",
                "DEPTOR", "TSC1", "TSC2", "RHEB", "AKT1", "AKT2", "AKT3",
                "PTEN", "PIK3CA", "PIK3CB", "PIK3CD", "RRAGA", "RRAGB", "RRAGC", "RRAGD")
mito_genes <- c("DNM1L", "FIS1", "MFF", "MIEF1", "MIEF2", "MFN1", "MFN2", "OPA1",
                "MARCHF5", "PINK1", "SLC25A3", "SLC25A5", "SLC25A6", "CYCS")
metab_genes <- c("HK1", "HK2", "PFKFB3", "PKM", "LDHA", "LDHB", "IDH1", "IDH2",
                 "MDH1", "MDH2", "CS", "GLS", "GLS2", "GOT1", "GOT2",
                 "NLRP3", "CASP1", "CASP4", "CASP5", "STING1", "CGAS", "MAOB", "ENDOG")
all73 <- c(core_genes, mtor_genes, mito_genes, metab_genes)

# 基因覆盖
found73 <- intersect(all73, rownames(tumor_mat))
cat("73 基因覆盖:", length(found73), "/73\n")

# 计算 overall program score（简单均值 Z-score，作为 ssGSEA 近似）
# 注: 完整 ssGSEA 需要 GSVA 包，这里用均值 Z 作为 program-level 代理
expr73 <- log2(tumor_mat[found73, ] + 1)
z73 <- t(scale(t(expr73)))  # gene-level Z
program_score <- colMeans(z73, na.rm = TRUE)

# Cox (502)
surv502 <- surv_full[!is.na(surv_full$OS_time), ]
score502 <- program_score[surv502$sample]
fit502 <- coxph(Surv(OS_time, OS_status) ~ scale(score502), data = surv502)
s502 <- summary(fit502)
cat(sprintf("502 队列 overall program: HR=%.3f (95%%CI %.3f-%.3f), P=%.4g\n",
            exp(coef(fit502)), exp(confint(fit502))[1], exp(confint(fit502))[2], s502$coefficients[5]))

# 对比 226 队列 (用现有 TPM 数据重算同样指标)
tpm_old <- readRDS(file.path(DATA_DIR, "TCGA_HNSC_TPM_tumor.Rds"))
surv_old2 <- surv_old
common_old <- intersect(colnames(tpm_old), surv_old2$sample)
expr73_old <- log2(tpm_old[intersect(found73, rownames(tpm_old)), common_old] + 1)
z73_old <- t(scale(t(expr73_old)))
score226 <- colMeans(z73_old, na.rm = TRUE)
fit226 <- coxph(Surv(OS_time, OS_status) ~ scale(score226), data = surv_old2[match(common_old, surv_old2$sample), ])
s226 <- summary(fit226)
cat(sprintf("226 队列 overall program: HR=%.3f (95%%CI %.3f-%.3f), P=%.4g\n",
            exp(coef(fit226)), exp(confint(fit226))[1], exp(confint(fit226))[2], s226$coefficients[5]))

# ============================================================
# 4. 5-gene LASSO (502 队列)
# ============================================================
cat("\n>>> 5-gene LASSO-Cox (502 队列)...\n")
X502 <- t(z73[, surv502$sample])  # samples x genes
y502 <- surv502
fit_lasso502 <- cv.glmnet(X502, Surv(y502$OS_time, y502$OS_status), family = "cox",
                          alpha = 1, nfolds = 10, foldid = rep(1:10, length.out = nrow(X502)))
coefs502 <- coef(fit_lasso502, s = "lambda.min")[, 1]
sel502 <- names(coefs502)[coefs502 != 0]
cat("502 队列 LASSO 选择:", length(sel502), ":", paste(sel502, collapse = ", "), "\n")
if (length(sel502) > 0) {
  risk502 <- as.numeric(X502[, sel502, drop = FALSE] %*% coefs502[sel502])
  cfit <- coxph(Surv(OS_time, OS_status) ~ risk502, data = y502)
  cat(sprintf("502 队列模型 C-index: %.3f, HR per SD: %.3f\n",
              summary(cfit)$concordance[1],
              exp(coef(summary(cfit))[1]) * sd(risk502)))
  # AUC
  roc502 <- timeROC(T = y502$OS_time, delta = y502$OS_status, marker = risk502, cause = 1,
                    times = c(365, 1095, 1825), iid = FALSE)
  cat("502 队列 AUC 1/3/5yr:", round(roc502$AUC, 3), "\n")
}

# ============================================================
# 5. 输出
# ============================================================
write.csv(data.frame(cohort = c("226", "502"), N = c(226, 502)),
          file.path(OUT_DIR, "cohort_comparison.csv"), row.names = FALSE)

cat("\n=== B1 全队列重跑 DONE ===")
